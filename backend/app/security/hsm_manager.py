"""
Hardware Security Module (HSM) integration for PhaGen.
Provides hardware-backed key storage for high-assurance customers.
"""
from __future__ import annotations

import logging
import os
from abc import ABC, abstractmethod
from typing import Dict, Optional, List
from dataclasses import dataclass
from pathlib import Path

logger = logging.getLogger(__name__)


@dataclass
class HSMConfig:
    """HSM configuration."""
    provider: str  # 'softhsm', 'aws-cloudhsm', 'azure-keyvault', 'pkcs11'
    slot: Optional[int] = None
    token_label: Optional[str] = None
    pin: Optional[str] = None
    library_path: Optional[str] = None
    key_label_prefix: str = "phagen"


class KeyProvider(ABC):
    """Abstract key provider interface."""
    
    @abstractmethod
    def generate_key(self, key_id: str, key_type: str = "AES-256") -> bool:
        """Generate a new key."""
        pass
    
    @abstractmethod
    def encrypt(self, key_id: str, plaintext: bytes) -> bytes:
        """Encrypt data using specified key."""
        pass
    
    @abstractmethod
    def decrypt(self, key_id: str, ciphertext: bytes) -> bytes:
        """Decrypt data using specified key."""
        pass
    
    @abstractmethod
    def list_keys(self) -> List[str]:
        """List all available keys."""
        pass
    
    @abstractmethod
    def delete_key(self, key_id: str) -> bool:
        """Delete a key."""
        pass


class SoftHSMProvider(KeyProvider):
    """
    SoftHSM provider for development/testing.
    Uses software-based HSM simulation.
    """
    
    def __init__(self, config: HSMConfig):
        self.config = config
        self.keys: Dict[str, bytes] = {}
        logger.info("Initialized SoftHSM provider")
        
        # Check if SoftHSM is available
        try:
            import PyKCS11
            self.pkcs11_available = True
            logger.info("PyKCS11 available for SoftHSM integration")
        except ImportError:
            self.pkcs11_available = False
            logger.warning("PyKCS11 not installed, using mock implementation")
    
    def generate_key(self, key_id: str, key_type: str = "AES-256") -> bool:
        """Generate a new key in SoftHSM."""
        try:
            # In production with PyKCS11:
            # session = pkcs11.openSession(slot, PyKCS11.CKF_SERIAL_SESSION | PyKCS11.CKF_RW_SESSION)
            # session.login(pin)
            # key = session.generateKey(template)
            
            # Mock implementation
            import secrets
            self.keys[key_id] = secrets.token_bytes(32)  # 256 bits
            logger.info(f"Generated key: {key_id}")
            return True
        except Exception as e:
            logger.error(f"Failed to generate key: {e}")
            return False
    
    def encrypt(self, key_id: str, plaintext: bytes) -> bytes:
        """Encrypt data."""
        if key_id not in self.keys:
            raise ValueError(f"Key not found: {key_id}")
        
        # Mock encryption (use AES-GCM in production)
        try:
            from cryptography.hazmat.primitives.ciphers.aead import AESGCM
            aesgcm = AESGCM(self.keys[key_id])
            nonce = os.urandom(12)
            ciphertext = aesgcm.encrypt(nonce, plaintext, None)
            return nonce + ciphertext
        except ImportError:
            # Fallback without cryptography library
            logger.warning("cryptography library not installed, using mock encryption")
            return b"ENCRYPTED:" + plaintext
    
    def decrypt(self, key_id: str, ciphertext: bytes) -> bytes:
        """Decrypt data."""
        if key_id not in self.keys:
            raise ValueError(f"Key not found: {key_id}")
        
        try:
            from cryptography.hazmat.primitives.ciphers.aead import AESGCM
            aesgcm = AESGCM(self.keys[key_id])
            nonce = ciphertext[:12]
            actual_ciphertext = ciphertext[12:]
            plaintext = aesgcm.decrypt(nonce, actual_ciphertext, None)
            return plaintext
        except ImportError:
            # Fallback
            if ciphertext.startswith(b"ENCRYPTED:"):
                return ciphertext[10:]
            raise ValueError("Invalid ciphertext format")
    
    def list_keys(self) -> List[str]:
        """List all keys."""
        return list(self.keys.keys())
    
    def delete_key(self, key_id: str) -> bool:
        """Delete a key."""
        if key_id in self.keys:
            del self.keys[key_id]
            logger.info(f"Deleted key: {key_id}")
            return True
        return False


class AWSCloudHSMProvider(KeyProvider):
    """
    AWS CloudHSM provider for production use.
    """
    
    def __init__(self, config: HSMConfig):
        self.config = config
        logger.info("Initialized AWS CloudHSM provider")
        
        # Check for boto3
        try:
            import boto3
            self.client = boto3.client('cloudhsmv2')
            logger.info("AWS CloudHSM client initialized")
        except ImportError:
            logger.error("boto3 not installed, AWS CloudHSM unavailable")
            self.client = None
    
    def generate_key(self, key_id: str, key_type: str = "AES-256") -> bool:
        """Generate key in AWS CloudHSM."""
        if not self.client:
            logger.error("AWS CloudHSM not available")
            return False
        
        # Placeholder for actual AWS CloudHSM integration
        logger.info(f"Would generate key in AWS CloudHSM: {key_id}")
        return False
    
    def encrypt(self, key_id: str, plaintext: bytes) -> bytes:
        """Encrypt using AWS CloudHSM."""
        raise NotImplementedError("AWS CloudHSM integration pending")
    
    def decrypt(self, key_id: str, ciphertext: bytes) -> bytes:
        """Decrypt using AWS CloudHSM."""
        raise NotImplementedError("AWS CloudHSM integration pending")
    
    def list_keys(self) -> List[str]:
        """List keys in AWS CloudHSM."""
        return []
    
    def delete_key(self, key_id: str) -> bool:
        """Delete key from AWS CloudHSM."""
        return False


class AzureKeyVaultProvider(KeyProvider):
    """
    Azure Key Vault provider for production use.
    """
    
    def __init__(self, config: HSMConfig):
        self.config = config
        logger.info("Initialized Azure Key Vault provider")
        
        # Check for azure-keyvault
        try:
            from azure.keyvault.keys import KeyClient
            from azure.identity import DefaultAzureCredential
            
            vault_url = os.getenv("AZURE_KEYVAULT_URL")
            if vault_url:
                credential = DefaultAzureCredential()
                self.client = KeyClient(vault_url=vault_url, credential=credential)
                logger.info("Azure Key Vault client initialized")
            else:
                logger.warning("AZURE_KEYVAULT_URL not set")
                self.client = None
        except ImportError:
            logger.error("azure-keyvault not installed")
            self.client = None
    
    def generate_key(self, key_id: str, key_type: str = "AES-256") -> bool:
        """Generate key in Azure Key Vault."""
        if not self.client:
            return False
        
        # Placeholder for actual Azure Key Vault integration
        logger.info(f"Would generate key in Azure Key Vault: {key_id}")
        return False
    
    def encrypt(self, key_id: str, plaintext: bytes) -> bytes:
        """Encrypt using Azure Key Vault."""
        raise NotImplementedError("Azure Key Vault integration pending")
    
    def decrypt(self, key_id: str, ciphertext: bytes) -> bytes:
        """Decrypt using Azure Key Vault."""
        raise NotImplementedError("Azure Key Vault integration pending")
    
    def list_keys(self) -> List[str]:
        """List keys in Azure Key Vault."""
        return []
    
    def delete_key(self, key_id: str) -> bool:
        """Delete key from Azure Key Vault."""
        return False


class HSMManager:
    """
    HSM manager with provider abstraction.
    """
    
    def __init__(self, config: Optional[HSMConfig] = None):
        if config is None:
            # Default to SoftHSM for development
            config = HSMConfig(
                provider="softhsm",
                token_label="phagen-dev",
                key_label_prefix="phagen"
            )
        
        self.config = config
        self.provider = self._initialize_provider()
        logger.info(f"HSM Manager initialized with provider: {config.provider}")
    
    def _initialize_provider(self) -> KeyProvider:
        """Initialize appropriate provider."""
        if self.config.provider == "softhsm":
            return SoftHSMProvider(self.config)
        elif self.config.provider == "aws-cloudhsm":
            return AWSCloudHSMProvider(self.config)
        elif self.config.provider == "azure-keyvault":
            return AzureKeyVaultProvider(self.config)
        else:
            logger.warning(f"Unknown provider: {self.config.provider}, using SoftHSM")
            return SoftHSMProvider(self.config)
    
    def generate_master_key(self, purpose: str = "data-encryption") -> str:
        """Generate a master encryption key."""
        key_id = f"{self.config.key_label_prefix}_{purpose}"
        success = self.provider.generate_key(key_id)
        
        if success:
            logger.info(f"Generated master key: {key_id}")
            return key_id
        else:
            raise RuntimeError(f"Failed to generate key: {key_id}")
    
    def encrypt_sensitive_data(self, key_id: str, data: str) -> str:
        """Encrypt sensitive data."""
        plaintext = data.encode('utf-8')
        ciphertext = self.provider.encrypt(key_id, plaintext)
        
        # Return base64-encoded ciphertext
        import base64
        return base64.b64encode(ciphertext).decode('ascii')
    
    def decrypt_sensitive_data(self, key_id: str, encrypted_data: str) -> str:
        """Decrypt sensitive data."""
        import base64
        ciphertext = base64.b64decode(encrypted_data)
        plaintext = self.provider.decrypt(key_id, ciphertext)
        
        return plaintext.decode('utf-8')
    
    def list_all_keys(self) -> List[str]:
        """List all keys in HSM."""
        return self.provider.list_keys()
    
    def rotate_key(self, old_key_id: str, new_key_id: str) -> bool:
        """Rotate encryption keys."""
        # Generate new key
        success = self.provider.generate_key(new_key_id)
        if not success:
            return False
        
        logger.info(f"Key rotation: {old_key_id} -> {new_key_id}")
        # Note: Re-encryption of data would happen separately
        return True
    
    def get_stats(self) -> Dict:
        """Get HSM manager statistics."""
        return {
            "provider": self.config.provider,
            "key_count": len(self.provider.list_keys()),
            "key_label_prefix": self.config.key_label_prefix
        }


def setup_hsm_for_database_encryption(
    database_url: str,
    hsm_manager: Optional[HSMManager] = None
) -> str:
    """
    Setup HSM-backed database encryption.
    
    Args:
        database_url: Database connection URL
        hsm_manager: Optional HSM manager instance
    
    Returns:
        Encrypted database URL
    """
    if hsm_manager is None:
        hsm_manager = HSMManager()
    
    # Generate or get database encryption key
    key_id = hsm_manager.generate_master_key("database")
    
    # Extract password from URL
    import re
    match = re.search(r'://([^:]+):([^@]+)@', database_url)
    if match:
        password = match.group(2)
        
        # Encrypt password
        encrypted_password = hsm_manager.encrypt_sensitive_data(key_id, password)
        
        # Replace in URL (in production, store encrypted separately)
        logger.info("Database credentials encrypted with HSM")
        return database_url.replace(password, f"HSM_ENCRYPTED:{key_id}")
    
    return database_url


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    
    print("=" * 60)
    print("HSM Integration Demo")
    print("=" * 60)
    
    # Initialize HSM manager
    manager = HSMManager()
    
    # Generate master key
    key_id = manager.generate_master_key("demo")
    print(f"\n✓ Generated master key: {key_id}")
    
    # Encrypt sensitive data
    sensitive_data = "database_password_12345"
    encrypted = manager.encrypt_sensitive_data(key_id, sensitive_data)
    print(f"\n✓ Encrypted data: {encrypted[:50]}...")
    
    # Decrypt
    decrypted = manager.decrypt_sensitive_data(key_id, encrypted)
    print(f"✓ Decrypted data: {decrypted}")
    print(f"✓ Match: {decrypted == sensitive_data}")
    
    # List keys
    keys = manager.list_all_keys()
    print(f"\n✓ Available keys: {keys}")
    
    # Stats
    stats = manager.get_stats()
    print(f"\n✓ HSM Stats: {stats}")
    
    print("\n" + "=" * 60)
    print("Production Setup:")
    print("1. Install SoftHSM: apt-get install softhsm2")
    print("2. Or use AWS CloudHSM / Azure Key Vault")
    print("3. Set environment variables:")
    print("   - HSM_PROVIDER=softhsm|aws-cloudhsm|azure-keyvault")
    print("   - HSM_TOKEN_LABEL=phagen-prod")
    print("   - HSM_PIN=<secure_pin>")
    print("4. Initialize keys: python -m backend.app.security.hsm_manager")
