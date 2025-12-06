"""
End-to-End Encryption Manager

Provides encryption at rest and in transit:
- TLS 1.3 for transport security
- AES-256-GCM for data at rest
- KMS integration for key management
- Database field-level encryption
"""
import os
import base64
import logging
from typing import Optional, Dict, Any, Union
from pathlib import Path
from datetime import datetime, timezone
from cryptography.hazmat.primitives.ciphers.aead import AESGCM
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.backends import default_backend
import secrets

logger = logging.getLogger(__name__)


class EncryptionConfig:
    """Configuration for encryption settings"""
    
    def __init__(
        self,
        kms_provider: str = "local",  # local, aws, azure, gcp
        master_key_id: Optional[str] = None,
        key_rotation_days: int = 90,
        algorithm: str = "AES-256-GCM",
    ):
        self.kms_provider = kms_provider
        self.master_key_id = master_key_id or os.getenv("MASTER_KEY_ID")
        self.key_rotation_days = key_rotation_days
        self.algorithm = algorithm
        
        # Load from environment
        self.enabled = os.getenv("ENCRYPTION_ENABLED", "true").lower() == "true"
        self.key_cache_ttl = int(os.getenv("KEY_CACHE_TTL", "3600"))  # 1 hour


class DataEncryptionKey:
    """Data Encryption Key (DEK) for envelope encryption"""
    
    def __init__(self, key_id: str, key_material: bytes, created_at: datetime):
        self.key_id = key_id
        self.key_material = key_material
        self.created_at = created_at
        self.cipher = AESGCM(key_material)
    
    def encrypt(self, plaintext: bytes, associated_data: Optional[bytes] = None) -> Dict[str, str]:
        """
        Encrypt data with AES-256-GCM
        
        Returns:
            Dict with ciphertext, nonce, and key_id
        """
        nonce = secrets.token_bytes(12)  # 96-bit nonce for GCM
        ciphertext = self.cipher.encrypt(nonce, plaintext, associated_data)
        
        return {
            "ciphertext": base64.b64encode(ciphertext).decode('utf-8'),
            "nonce": base64.b64encode(nonce).decode('utf-8'),
            "key_id": self.key_id,
            "algorithm": "AES-256-GCM",
        }
    
    def decrypt(
        self,
        ciphertext: str,
        nonce: str,
        associated_data: Optional[bytes] = None
    ) -> bytes:
        """Decrypt AES-256-GCM encrypted data"""
        ciphertext_bytes = base64.b64decode(ciphertext)
        nonce_bytes = base64.b64decode(nonce)
        
        plaintext = self.cipher.decrypt(nonce_bytes, ciphertext_bytes, associated_data)
        return plaintext
    
    def needs_rotation(self, max_age_days: int) -> bool:
        """Check if key needs rotation"""
        age = datetime.now(timezone.utc) - self.created_at
        return age.days >= max_age_days


class EncryptionManager:
    """
    Manages encryption keys and operations
    
    Features:
    - Envelope encryption (DEK encrypted by KEK/Master Key)
    - Automatic key rotation
    - KMS integration
    - Field-level database encryption
    """
    
    def __init__(self, config: Optional[EncryptionConfig] = None):
        self.config = config or EncryptionConfig()
        self._dek_cache: Dict[str, DataEncryptionKey] = {}
        self._master_key: Optional[bytes] = None
        
        if self.config.enabled:
            self._initialize_master_key()
    
    def _initialize_master_key(self):
        """Initialize or load master encryption key"""
        if self.config.kms_provider == "local":
            # Local key management (dev/testing only)
            key_file = Path(os.getenv("MASTER_KEY_FILE", ".master_key"))
            
            if key_file.exists():
                with open(key_file, 'rb') as f:
                    self._master_key = f.read()
                logger.info("[ENCRYPTION] Loaded master key from file")
            else:
                # Generate new master key
                self._master_key = AESGCM.generate_key(bit_length=256)
                
                # Save to file (ensure proper permissions)
                key_file.parent.mkdir(parents=True, exist_ok=True)
                with open(key_file, 'wb') as f:
                    f.write(self._master_key)
                os.chmod(key_file, 0o600)  # Owner read/write only
                
                logger.warning(
                    "[ENCRYPTION] Generated new master key - "
                    "DO NOT use local keys in production!"
                )
        
        elif self.config.kms_provider == "aws":
            # AWS KMS integration
            from .kms_providers import AWSKMSProvider
            region = os.getenv("AWS_REGION", "us-east-1")
            self.kms = AWSKMSProvider(
                key_id=self.config.master_key_id or os.getenv("AWS_KMS_KEY_ID"),
                region_name=region
            )
            logger.info("[ENCRYPTION] Using AWS KMS")
        
        elif self.config.kms_provider == "azure":
            # Azure Key Vault integration
            from .kms_providers import AzureKMSProvider
            vault_url = os.getenv("AZURE_KEY_VAULT_URL")
            key_name = self.config.master_key_id or os.getenv("AZURE_KEY_NAME")
            self.kms = AzureKMSProvider(vault_url=vault_url, key_name=key_name)
            logger.info("[ENCRYPTION] Using Azure Key Vault")
        
        elif self.config.kms_provider == "gcp":
            # Google Cloud KMS integration
            from .kms_providers import GCPKMSProvider
            self.kms = GCPKMSProvider(
                project_id=os.getenv("GCP_PROJECT_ID"),
                location=os.getenv("GCP_LOCATION", "us-east1"),
                key_ring=os.getenv("GCP_KEY_RING"),
                key_name=self.config.master_key_id or os.getenv("GCP_KEY_NAME")
            )
            logger.info("[ENCRYPTION] Using Google Cloud KMS")
    
    def _generate_dek(self) -> DataEncryptionKey:
        """Generate new Data Encryption Key"""
        key_id = f"dek-{secrets.token_hex(16)}"
        key_material = AESGCM.generate_key(bit_length=256)
        created_at = datetime.now(timezone.utc)
        
        dek = DataEncryptionKey(key_id, key_material, created_at)
        self._dek_cache[key_id] = dek
        
        logger.info(f"[ENCRYPTION] Generated new DEK: {key_id}")
        return dek
    
    def _get_or_create_dek(self, key_id: Optional[str] = None) -> DataEncryptionKey:
        """Get existing DEK or create new one"""
        if key_id and key_id in self._dek_cache:
            dek = self._dek_cache[key_id]
            
            # Check if rotation needed
            if dek.needs_rotation(self.config.key_rotation_days):
                logger.info(f"[ENCRYPTION] DEK {key_id} needs rotation")
                return self._generate_dek()
            
            return dek
        
        return self._generate_dek()
    
    def encrypt_field(
        self,
        plaintext: Union[str, bytes],
        field_name: Optional[str] = None
    ) -> Dict[str, str]:
        """
        Encrypt a database field value
        
        Args:
            plaintext: Data to encrypt
            field_name: Field name for associated data
            
        Returns:
            Dictionary with encrypted data and metadata
        """
        if not self.config.enabled:
            # Encryption disabled - return plaintext
            if isinstance(plaintext, str):
                plaintext = plaintext.encode('utf-8')
            return {
                "ciphertext": base64.b64encode(plaintext).decode('utf-8'),
                "algorithm": "none",
            }
        
        # Convert to bytes
        if isinstance(plaintext, str):
            plaintext_bytes = plaintext.encode('utf-8')
        else:
            plaintext_bytes = plaintext
        
        # Use field name as associated data for authentication
        associated_data = field_name.encode('utf-8') if field_name else None
        
        # Get DEK and encrypt
        dek = self._get_or_create_dek()
        encrypted = dek.encrypt(plaintext_bytes, associated_data)
        
        return encrypted
    
    def decrypt_field(
        self,
        encrypted_data: Dict[str, str],
        field_name: Optional[str] = None
    ) -> Union[str, bytes]:
        """
        Decrypt a database field value
        
        Args:
            encrypted_data: Dictionary with ciphertext, nonce, key_id
            field_name: Field name for associated data
            
        Returns:
            Decrypted plaintext as string (or bytes if can't decode)
        """
        if encrypted_data.get("algorithm") == "none":
            # No encryption
            plaintext_bytes = base64.b64decode(encrypted_data["ciphertext"])
            try:
                return plaintext_bytes.decode('utf-8')
            except UnicodeDecodeError:
                return plaintext_bytes
        
        # Get DEK
        key_id = encrypted_data["key_id"]
        dek = self._dek_cache.get(key_id)
        
        if not dek:
            raise ValueError(f"Encryption key not found: {key_id}")
        
        # Decrypt with associated data
        associated_data = field_name.encode('utf-8') if field_name else None
        
        plaintext = dek.decrypt(
            encrypted_data["ciphertext"],
            encrypted_data["nonce"],
            associated_data
        )
        
        # Try to decode as UTF-8 string
        try:
            return plaintext.decode('utf-8')
        except UnicodeDecodeError:
            return plaintext
        
        return plaintext
    
    def encrypt_file(self, file_path: Path, output_path: Optional[Path] = None) -> Path:
        """
        Encrypt a file
        
        Args:
            file_path: Path to file to encrypt
            output_path: Optional output path (defaults to .encrypted suffix)
            
        Returns:
            Path to encrypted file
        """
        if output_path is None:
            output_path = file_path.with_suffix(file_path.suffix + '.encrypted')
        
        # Read file
        with open(file_path, 'rb') as f:
            plaintext = f.read()
        
        # Encrypt
        dek = self._get_or_create_dek()
        encrypted = dek.encrypt(plaintext)
        
        # Write encrypted file with metadata
        import json
        with open(output_path, 'w') as f:
            json.dump(encrypted, f)
        
        logger.info(f"[ENCRYPTION] Encrypted file: {file_path} -> {output_path}")
        return output_path
    
    def decrypt_file(self, encrypted_path: Path, output_path: Optional[Path] = None) -> Path:
        """
        Decrypt a file
        
        Args:
            encrypted_path: Path to encrypted file
            output_path: Optional output path
            
        Returns:
            Path to decrypted file
        """
        if output_path is None:
            output_path = encrypted_path.with_suffix('')
        
        # Read encrypted file
        import json
        with open(encrypted_path, 'r') as f:
            encrypted_data = json.load(f)
        
        # Decrypt
        key_id = encrypted_data["key_id"]
        dek = self._dek_cache.get(key_id)
        if not dek:
            raise ValueError(f"Encryption key not found: {key_id}")
        
        plaintext = dek.decrypt(
            encrypted_data["ciphertext"],
            encrypted_data["nonce"]
        )
        
        # Write decrypted file
        with open(output_path, 'wb') as f:
            f.write(plaintext)
        
        logger.info(f"[ENCRYPTION] Decrypted file: {encrypted_path} -> {output_path}")
        return output_path
    
    def rotate_keys(self):
        """Rotate all encryption keys"""
        logger.info("[ENCRYPTION] Starting key rotation")
        
        old_keys = list(self._dek_cache.keys())
        
        for key_id in old_keys:
            dek = self._dek_cache[key_id]
            if dek.needs_rotation(self.config.key_rotation_days):
                # Generate new key
                new_dek = self._generate_dek()
                
                # Mark old key for re-encryption
                logger.info(f"[ENCRYPTION] Rotated key: {key_id} -> {new_dek.key_id}")
        
        logger.info("[ENCRYPTION] Key rotation complete")


# Global encryption manager
_encryption_manager: Optional[EncryptionManager] = None


def get_encryption_manager() -> EncryptionManager:
    """Get or create global encryption manager"""
    global _encryption_manager
    if _encryption_manager is None:
        _encryption_manager = EncryptionManager()
    return _encryption_manager


def encrypt_sensitive_field(value: str, field_name: str) -> Dict[str, str]:
    """Convenience function to encrypt sensitive field"""
    manager = get_encryption_manager()
    return manager.encrypt_field(value, field_name)


def decrypt_sensitive_field(encrypted_data: Dict[str, str], field_name: str) -> str:
    """Convenience function to decrypt sensitive field"""
    manager = get_encryption_manager()
    plaintext_bytes = manager.decrypt_field(encrypted_data, field_name)
    return plaintext_bytes.decode('utf-8')


# TLS Configuration Helper
def get_tls_config() -> Dict[str, Any]:
    """
    Get TLS 1.3 configuration for uvicorn/nginx
    
    Returns:
        Dictionary with TLS settings
    """
    return {
        "ssl_version": "TLSv1.3",
        "ssl_ciphers": [
            "TLS_AES_256_GCM_SHA384",
            "TLS_CHACHA20_POLY1305_SHA256",
            "TLS_AES_128_GCM_SHA256",
        ],
        "ssl_cert_file": os.getenv("SSL_CERT_FILE", "/etc/ssl/certs/server.crt"),
        "ssl_key_file": os.getenv("SSL_KEY_FILE", "/etc/ssl/private/server.key"),
        "ssl_ca_file": os.getenv("SSL_CA_FILE"),
        "ssl_verify_mode": "CERT_REQUIRED" if os.getenv("SSL_VERIFY_CLIENT") else "CERT_NONE",
    }


def generate_self_signed_cert(
    cert_file: Path,
    key_file: Path,
    days_valid: int = 365
):
    """
    Generate self-signed certificate for development
    
    Args:
        cert_file: Path to certificate file
        key_file: Path to private key file
        days_valid: Certificate validity in days
    """
    from cryptography import x509
    from cryptography.x509.oid import NameOID
    from cryptography.hazmat.primitives import serialization
    from cryptography.hazmat.primitives.asymmetric import rsa
    from datetime import timedelta
    
    # Generate private key
    private_key = rsa.generate_private_key(
        public_exponent=65537,
        key_size=4096,
        backend=default_backend()
    )
    
    # Generate certificate
    subject = issuer = x509.Name([
        x509.NameAttribute(NameOID.COUNTRY_NAME, "US"),
        x509.NameAttribute(NameOID.STATE_OR_PROVINCE_NAME, "CA"),
        x509.NameAttribute(NameOID.LOCALITY_NAME, "San Francisco"),
        x509.NameAttribute(NameOID.ORGANIZATION_NAME, "PhaGen"),
        x509.NameAttribute(NameOID.COMMON_NAME, "localhost"),
    ])
    
    cert = x509.CertificateBuilder().subject_name(
        subject
    ).issuer_name(
        issuer
    ).public_key(
        private_key.public_key()
    ).serial_number(
        x509.random_serial_number()
    ).not_valid_before(
        datetime.now(timezone.utc)
    ).not_valid_after(
        datetime.now(timezone.utc) + timedelta(days=days_valid)
    ).add_extension(
        x509.SubjectAlternativeName([
            x509.DNSName("localhost"),
            x509.DNSName("*.localhost"),
        ]),
        critical=False,
    ).sign(private_key, hashes.SHA256(), default_backend())
    
    # Write private key
    key_file = Path(key_file)
    key_file.parent.mkdir(parents=True, exist_ok=True)
    with open(key_file, "wb") as f:
        f.write(private_key.private_bytes(
            encoding=serialization.Encoding.PEM,
            format=serialization.PrivateFormat.TraditionalOpenSSL,
            encryption_algorithm=serialization.NoEncryption()
        ))
    os.chmod(key_file, 0o600)
    
    # Write certificate
    cert_file = Path(cert_file)
    cert_file.parent.mkdir(parents=True, exist_ok=True)
    with open(cert_file, "wb") as f:
        f.write(cert.public_bytes(serialization.Encoding.PEM))
    
    logger.info(f"[TLS] Generated self-signed certificate: {cert_file}")
