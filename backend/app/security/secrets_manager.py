"""
Centralized Secrets Management

Integrates with HashiCorp Vault, AWS Secrets Manager, Azure Key Vault, and GCP Secret Manager
for secure storage and rotation of application secrets.
"""
import os
import json
import logging
from typing import Optional, Dict, Any, List
from datetime import datetime, timedelta, timezone
from abc import ABC, abstractmethod
import hashlib

logger = logging.getLogger(__name__)


class SecretNotFoundError(Exception):
    """Raised when secret is not found"""
    pass


class SecretProvider(ABC):
    """Abstract base class for secret providers"""
    
    @abstractmethod
    def get_secret(self, secret_name: str, version: Optional[str] = None) -> str:
        """Retrieve secret value"""
        pass
    
    @abstractmethod
    def set_secret(self, secret_name: str, secret_value: str) -> None:
        """Store secret value"""
        pass
    
    @abstractmethod
    def delete_secret(self, secret_name: str) -> None:
        """Delete secret"""
        pass
    
    @abstractmethod
    def list_secrets(self) -> List[str]:
        """List all secret names"""
        pass
    
    @abstractmethod
    def rotate_secret(self, secret_name: str) -> str:
        """Rotate secret and return new value"""
        pass


class LocalSecretsProvider(SecretProvider):
    """
    Local file-based secrets (development only)
    
    WARNING: DO NOT use in production!
    """
    
    def __init__(self, secrets_file: str = ".secrets.json"):
        self.secrets_file = secrets_file
        self._load_secrets()
    
    def _load_secrets(self):
        """Load secrets from file"""
        if os.path.exists(self.secrets_file):
            with open(self.secrets_file, 'r') as f:
                self._secrets = json.load(f)
            logger.warning("[SECRETS] Loaded secrets from local file - NOT FOR PRODUCTION")
        else:
            self._secrets = {}
    
    def _save_secrets(self):
        """Save secrets to file"""
        with open(self.secrets_file, 'w') as f:
            json.dump(self._secrets, f, indent=2)
        os.chmod(self.secrets_file, 0o600)  # Owner read/write only
    
    def get_secret(self, secret_name: str, version: Optional[str] = None) -> str:
        """Get secret from local storage"""
        if secret_name not in self._secrets:
            raise SecretNotFoundError(f"Secret not found: {secret_name}")
        
        secret_data = self._secrets[secret_name]
        
        if version:
            # Get specific version
            versions = secret_data.get("versions", [])
            for v in versions:
                if str(v["version"]) == str(version):
                    return v["value"]
            raise SecretNotFoundError(f"Version {version} not found for {secret_name}")
        
        # Get current version
        return secret_data["current"]
    
    def set_secret(self, secret_name: str, secret_value: str) -> None:
        """Set secret in local storage"""
        if secret_name in self._secrets:
            # Archive old version
            old_value = self._secrets[secret_name]["current"]
            versions = self._secrets[secret_name].get("versions", [])
            versions.append({
                "version": len(versions) + 1,
                "value": old_value,
                "archived_at": datetime.now(timezone.utc).isoformat()
            })
            self._secrets[secret_name]["versions"] = versions
        else:
            self._secrets[secret_name] = {"versions": []}
        
        self._secrets[secret_name]["current"] = secret_value
        self._secrets[secret_name]["updated_at"] = datetime.now(timezone.utc).isoformat()
        
        self._save_secrets()
        logger.info(f"[SECRETS] Updated secret: {secret_name}")
    
    def delete_secret(self, secret_name: str) -> None:
        """Delete secret from local storage"""
        if secret_name in self._secrets:
            del self._secrets[secret_name]
            self._save_secrets()
            logger.info(f"[SECRETS] Deleted secret: {secret_name}")
    
    def list_secrets(self) -> List[str]:
        """List all secret names"""
        return list(self._secrets.keys())
    
    def rotate_secret(self, secret_name: str) -> str:
        """Rotate secret (generate new value)"""
        import secrets
        new_value = secrets.token_urlsafe(32)
        self.set_secret(secret_name, new_value)
        return new_value


class VaultSecretsProvider(SecretProvider):
    """HashiCorp Vault secrets provider"""
    
    def __init__(
        self,
        vault_addr: Optional[str] = None,
        vault_token: Optional[str] = None,
        mount_point: str = "secret"
    ):
        self.vault_addr = vault_addr or os.getenv("VAULT_ADDR", "http://localhost:8200")
        self.vault_token = vault_token or os.getenv("VAULT_TOKEN")
        self.mount_point = mount_point
        
        try:
            import hvac
            self.client = hvac.Client(url=self.vault_addr, token=self.vault_token)
            
            if not self.client.is_authenticated():
                raise ValueError("Vault authentication failed")
            
            logger.info(f"[SECRETS] Connected to Vault at {self.vault_addr}")
        except ImportError:
            raise ImportError("hvac package required for Vault integration: pip install hvac")
    
    def get_secret(self, secret_name: str, version: Optional[str] = None) -> str:
        """Get secret from Vault"""
        try:
            secret_path = f"{self.mount_point}/data/{secret_name}"
            
            if version:
                response = self.client.secrets.kv.v2.read_secret_version(
                    path=secret_name,
                    version=version,
                    mount_point=self.mount_point
                )
            else:
                response = self.client.secrets.kv.v2.read_secret_version(
                    path=secret_name,
                    mount_point=self.mount_point
                )
            
            return response['data']['data']['value']
        except Exception as e:
            raise SecretNotFoundError(f"Secret not found in Vault: {secret_name}") from e
    
    def set_secret(self, secret_name: str, secret_value: str) -> None:
        """Set secret in Vault"""
        self.client.secrets.kv.v2.create_or_update_secret(
            path=secret_name,
            secret={'value': secret_value},
            mount_point=self.mount_point
        )
        logger.info(f"[SECRETS] Updated secret in Vault: {secret_name}")
    
    def delete_secret(self, secret_name: str) -> None:
        """Delete secret from Vault"""
        self.client.secrets.kv.v2.delete_metadata_and_all_versions(
            path=secret_name,
            mount_point=self.mount_point
        )
        logger.info(f"[SECRETS] Deleted secret from Vault: {secret_name}")
    
    def list_secrets(self) -> List[str]:
        """List all secrets in Vault"""
        try:
            response = self.client.secrets.kv.v2.list_secrets(
                path="",
                mount_point=self.mount_point
            )
            return response['data']['keys']
        except Exception:
            return []
    
    def rotate_secret(self, secret_name: str) -> str:
        """Rotate secret in Vault"""
        import secrets
        new_value = secrets.token_urlsafe(32)
        self.set_secret(secret_name, new_value)
        return new_value


class AWSSecretsProvider(SecretProvider):
    """AWS Secrets Manager provider"""
    
    def __init__(self, region_name: Optional[str] = None):
        self.region_name = region_name or os.getenv("AWS_REGION", "us-east-1")
        
        try:
            import boto3
            self.client = boto3.client('secretsmanager', region_name=self.region_name)
            logger.info(f"[SECRETS] Connected to AWS Secrets Manager in {self.region_name}")
        except ImportError:
            raise ImportError("boto3 required for AWS: pip install boto3")
    
    def get_secret(self, secret_name: str, version: Optional[str] = None) -> str:
        """Get secret from AWS Secrets Manager"""
        try:
            kwargs = {'SecretId': secret_name}
            if version:
                kwargs['VersionId'] = version
            
            response = self.client.get_secret_value(**kwargs)
            return response['SecretString']
        except self.client.exceptions.ResourceNotFoundException:
            raise SecretNotFoundError(f"Secret not found in AWS: {secret_name}")
    
    def set_secret(self, secret_name: str, secret_value: str) -> None:
        """Set secret in AWS Secrets Manager"""
        try:
            self.client.update_secret(
                SecretId=secret_name,
                SecretString=secret_value
            )
        except self.client.exceptions.ResourceNotFoundException:
            # Create if doesn't exist
            self.client.create_secret(
                Name=secret_name,
                SecretString=secret_value
            )
        logger.info(f"[SECRETS] Updated secret in AWS: {secret_name}")
    
    def delete_secret(self, secret_name: str) -> None:
        """Delete secret from AWS Secrets Manager"""
        self.client.delete_secret(
            SecretId=secret_name,
            ForceDeleteWithoutRecovery=True
        )
        logger.info(f"[SECRETS] Deleted secret from AWS: {secret_name}")
    
    def list_secrets(self) -> List[str]:
        """List all secrets in AWS Secrets Manager"""
        response = self.client.list_secrets()
        return [s['Name'] for s in response['SecretList']]
    
    def rotate_secret(self, secret_name: str) -> str:
        """Rotate secret in AWS Secrets Manager"""
        import secrets
        new_value = secrets.token_urlsafe(32)
        self.set_secret(secret_name, new_value)
        return new_value


class AzureSecretsProvider(SecretProvider):
    """Azure Key Vault secrets provider"""
    
    def __init__(self, vault_url: Optional[str] = None):
        self.vault_url = vault_url or os.getenv("AZURE_KEY_VAULT_URL")
        
        if not self.vault_url:
            raise ValueError("AZURE_KEY_VAULT_URL must be set")
        
        try:
            from azure.keyvault.secrets import SecretClient
            from azure.identity import DefaultAzureCredential
            
            credential = DefaultAzureCredential()
            self.client = SecretClient(vault_url=self.vault_url, credential=credential)
            logger.info(f"[SECRETS] Connected to Azure Key Vault: {self.vault_url}")
        except ImportError:
            raise ImportError("azure-keyvault-secrets required: pip install azure-keyvault-secrets azure-identity")
    
    def get_secret(self, secret_name: str, version: Optional[str] = None) -> str:
        """Get secret from Azure Key Vault"""
        try:
            secret = self.client.get_secret(secret_name, version=version)
            return secret.value
        except Exception as e:
            raise SecretNotFoundError(f"Secret not found in Azure: {secret_name}") from e
    
    def set_secret(self, secret_name: str, secret_value: str) -> None:
        """Set secret in Azure Key Vault"""
        self.client.set_secret(secret_name, secret_value)
        logger.info(f"[SECRETS] Updated secret in Azure: {secret_name}")
    
    def delete_secret(self, secret_name: str) -> None:
        """Delete secret from Azure Key Vault"""
        poller = self.client.begin_delete_secret(secret_name)
        poller.wait()
        logger.info(f"[SECRETS] Deleted secret from Azure: {secret_name}")
    
    def list_secrets(self) -> List[str]:
        """List all secrets in Azure Key Vault"""
        return [s.name for s in self.client.list_properties_of_secrets()]
    
    def rotate_secret(self, secret_name: str) -> str:
        """Rotate secret in Azure Key Vault"""
        import secrets
        new_value = secrets.token_urlsafe(32)
        self.set_secret(secret_name, new_value)
        return new_value


class SecretsManager:
    """
    Unified secrets management interface
    
    Features:
    - Multi-provider support (Vault, AWS, Azure, GCP, local)
    - Automatic secret rotation
    - Secret caching with TTL
    - Audit logging
    """
    
    def __init__(self, provider: Optional[SecretProvider] = None):
        if provider:
            self.provider = provider
        else:
            # Auto-detect provider from environment
            provider_type = os.getenv("SECRETS_PROVIDER", "local").lower()
            
            if provider_type == "vault":
                self.provider = VaultSecretsProvider()
            elif provider_type == "aws":
                self.provider = AWSSecretsProvider()
            elif provider_type == "azure":
                self.provider = AzureSecretsProvider()
            elif provider_type == "local":
                self.provider = LocalSecretsProvider()
            else:
                raise ValueError(f"Unknown secrets provider: {provider_type}")
        
        self._cache: Dict[str, tuple[str, datetime]] = {}
        self._cache_ttl = int(os.getenv("SECRETS_CACHE_TTL", "300"))  # 5 minutes
    
    def get(self, secret_name: str, use_cache: bool = True) -> str:
        """
        Get secret value
        
        Args:
            secret_name: Name of secret
            use_cache: Whether to use cached value
            
        Returns:
            Secret value
        """
        # Check cache
        if use_cache and secret_name in self._cache:
            value, cached_at = self._cache[secret_name]
            age = (datetime.now(timezone.utc) - cached_at).total_seconds()
            
            if age < self._cache_ttl:
                return value
        
        # Fetch from provider
        value = self.provider.get_secret(secret_name)
        
        # Cache it
        self._cache[secret_name] = (value, datetime.now(timezone.utc))
        
        logger.debug(f"[SECRETS] Retrieved secret: {secret_name}")
        return value
    
    def set(self, secret_name: str, secret_value: str) -> None:
        """Set secret value"""
        self.provider.set_secret(secret_name, secret_value)
        
        # Invalidate cache
        if secret_name in self._cache:
            del self._cache[secret_name]
    
    def delete(self, secret_name: str) -> None:
        """Delete secret"""
        self.provider.delete_secret(secret_name)
        
        # Invalidate cache
        if secret_name in self._cache:
            del self._cache[secret_name]
    
    def list(self) -> List[str]:
        """List all secrets"""
        return self.provider.list_secrets()
    
    def rotate(self, secret_name: str) -> str:
        """Rotate secret"""
        new_value = self.provider.rotate_secret(secret_name)
        
        # Invalidate cache
        if secret_name in self._cache:
            del self._cache[secret_name]
        
        logger.info(f"[SECRETS] Rotated secret: {secret_name}")
        return new_value
    
    def rotate_all(self, pattern: Optional[str] = None):
        """Rotate all secrets matching pattern"""
        secrets = self.list()
        
        for secret_name in secrets:
            if pattern and pattern not in secret_name:
                continue
            
            try:
                self.rotate(secret_name)
            except Exception as e:
                logger.error(f"[SECRETS] Failed to rotate {secret_name}: {e}")


# Global secrets manager
_secrets_manager: Optional[SecretsManager] = None


def get_secrets_manager() -> SecretsManager:
    """Get or create global secrets manager"""
    global _secrets_manager
    if _secrets_manager is None:
        _secrets_manager = SecretsManager()
    return _secrets_manager


def get_secret(secret_name: str) -> str:
    """Convenience function to get secret"""
    manager = get_secrets_manager()
    return manager.get(secret_name)


def set_secret(secret_name: str, secret_value: str):
    """Convenience function to set secret"""
    manager = get_secrets_manager()
    manager.set(secret_name, secret_value)
