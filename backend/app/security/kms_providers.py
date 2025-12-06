"""
Key Management Service (KMS) providers for cloud platforms

Supports AWS KMS, Azure Key Vault, GCP Cloud KMS for master key management
"""
import os
import logging
from typing import Optional
from base64 import b64encode, b64decode

logger = logging.getLogger(__name__)


class KMSProvider:
    """Base class for KMS providers"""
    
    def encrypt(self, plaintext: bytes, key_id: str) -> bytes:
        """Encrypt data with master key"""
        raise NotImplementedError
    
    def decrypt(self, ciphertext: bytes, key_id: str) -> bytes:
        """Decrypt data with master key"""
        raise NotImplementedError
    
    def generate_data_key(self, key_id: str) -> tuple[bytes, bytes]:
        """
        Generate data encryption key
        
        Returns:
            (plaintext_key, encrypted_key)
        """
        raise NotImplementedError


class AWSKMSProvider(KMSProvider):
    """AWS Key Management Service provider"""
    
    def __init__(self, key_id: str, region_name: Optional[str] = None):
        self.key_id = key_id
        self.region_name = region_name or os.getenv("AWS_REGION", "us-east-1")
        
        try:
            import boto3
            self.client = boto3.client('kms', region_name=self.region_name)
            logger.info(f"[KMS] Connected to AWS KMS in {self.region_name}")
        except ImportError:
            raise ImportError("boto3 required for AWS KMS: pip install boto3")
    
    def encrypt(self, plaintext: bytes, key_id: Optional[str] = None) -> bytes:
        """Encrypt with AWS KMS"""
        response = self.client.encrypt(
            KeyId=key_id or self.key_id,
            Plaintext=plaintext
        )
        return response['CiphertextBlob']
    
    def decrypt(self, ciphertext: bytes, key_id: Optional[str] = None) -> bytes:
        """Decrypt with AWS KMS"""
        response = self.client.decrypt(
            KeyId=key_id or self.key_id,
            CiphertextBlob=ciphertext
        )
        return response['Plaintext']
    
    def generate_data_key(self, key_id: Optional[str] = None) -> tuple[bytes, bytes]:
        """Generate data key with AWS KMS"""
        response = self.client.generate_data_key(
            KeyId=key_id or self.key_id,
            KeySpec='AES_256'
        )
        return response['Plaintext'], response['CiphertextBlob']


class AzureKMSProvider(KMSProvider):
    """Azure Key Vault provider"""
    
    def __init__(self, vault_url: str, key_name: str):
        self.vault_url = vault_url
        self.key_name = key_name
        
        try:
            from azure.keyvault.keys.crypto import CryptographyClient
            from azure.keyvault.keys import KeyClient
            from azure.identity import DefaultAzureCredential
            
            credential = DefaultAzureCredential()
            key_client = KeyClient(vault_url=vault_url, credential=credential)
            key = key_client.get_key(key_name)
            
            self.crypto_client = CryptographyClient(key, credential=credential)
            logger.info(f"[KMS] Connected to Azure Key Vault: {vault_url}")
        except ImportError:
            raise ImportError(
                "azure-keyvault-keys required: pip install azure-keyvault-keys azure-identity"
            )
    
    def encrypt(self, plaintext: bytes, key_id: Optional[str] = None) -> bytes:
        """Encrypt with Azure Key Vault"""
        from azure.keyvault.keys.crypto import EncryptionAlgorithm
        
        result = self.crypto_client.encrypt(
            EncryptionAlgorithm.rsa_oaep_256,
            plaintext
        )
        return result.ciphertext
    
    def decrypt(self, ciphertext: bytes, key_id: Optional[str] = None) -> bytes:
        """Decrypt with Azure Key Vault"""
        from azure.keyvault.keys.crypto import EncryptionAlgorithm
        
        result = self.crypto_client.decrypt(
            EncryptionAlgorithm.rsa_oaep_256,
            ciphertext
        )
        return result.plaintext
    
    def generate_data_key(self, key_id: Optional[str] = None) -> tuple[bytes, bytes]:
        """Generate data key with Azure Key Vault"""
        import secrets
        
        # Generate random 256-bit key
        plaintext_key = secrets.token_bytes(32)
        
        # Encrypt with master key
        encrypted_key = self.encrypt(plaintext_key)
        
        return plaintext_key, encrypted_key


class GCPKMSProvider(KMSProvider):
    """Google Cloud KMS provider"""
    
    def __init__(self, project_id: str, location: str, key_ring: str, key_name: str):
        self.project_id = project_id
        self.location = location
        self.key_ring = key_ring
        self.key_name = key_name
        
        try:
            from google.cloud import kms
            self.client = kms.KeyManagementServiceClient()
            
            self.key_path = self.client.crypto_key_path(
                project_id, location, key_ring, key_name
            )
            logger.info(f"[KMS] Connected to GCP Cloud KMS: {self.key_path}")
        except ImportError:
            raise ImportError("google-cloud-kms required: pip install google-cloud-kms")
    
    def encrypt(self, plaintext: bytes, key_id: Optional[str] = None) -> bytes:
        """Encrypt with GCP Cloud KMS"""
        response = self.client.encrypt(
            request={
                'name': key_id or self.key_path,
                'plaintext': plaintext
            }
        )
        return response.ciphertext
    
    def decrypt(self, ciphertext: bytes, key_id: Optional[str] = None) -> bytes:
        """Decrypt with GCP Cloud KMS"""
        response = self.client.decrypt(
            request={
                'name': key_id or self.key_path,
                'ciphertext': ciphertext
            }
        )
        return response.plaintext
    
    def generate_data_key(self, key_id: Optional[str] = None) -> tuple[bytes, bytes]:
        """Generate data key with GCP Cloud KMS"""
        import secrets
        
        # Generate random 256-bit key
        plaintext_key = secrets.token_bytes(32)
        
        # Encrypt with master key
        encrypted_key = self.encrypt(plaintext_key)
        
        return plaintext_key, encrypted_key


def get_kms_provider(
    provider: str,
    key_id: Optional[str] = None,
    **kwargs
) -> KMSProvider:
    """
    Factory function to create KMS provider
    
    Args:
        provider: Provider type ('aws', 'azure', 'gcp')
        key_id: Master key identifier
        **kwargs: Provider-specific arguments
        
    Returns:
        KMSProvider instance
    """
    provider = provider.lower()
    
    if provider == "aws":
        return AWSKMSProvider(
            key_id=key_id or os.getenv("AWS_KMS_KEY_ID"),
            region_name=kwargs.get("region_name")
        )
    
    elif provider == "azure":
        return AzureKMSProvider(
            vault_url=kwargs.get("vault_url") or os.getenv("AZURE_KEY_VAULT_URL"),
            key_name=key_id or os.getenv("AZURE_KEY_NAME")
        )
    
    elif provider == "gcp":
        return GCPKMSProvider(
            project_id=kwargs.get("project_id") or os.getenv("GCP_PROJECT_ID"),
            location=kwargs.get("location") or os.getenv("GCP_LOCATION", "us-east1"),
            key_ring=kwargs.get("key_ring") or os.getenv("GCP_KEY_RING"),
            key_name=key_id or os.getenv("GCP_KEY_NAME")
        )
    
    else:
        raise ValueError(f"Unknown KMS provider: {provider}")
