"""
PII Redaction and Data Loss Prevention (DLP)

Protects sensitive data in pharmaceutical environments:
- Automatically detects and redacts PII/PHI
- Blocks exfiltration of sensitive identifiers
- Maintains data utility while ensuring compliance
"""
import re
import logging
import hashlib
from typing import Dict, Any, List, Set, Optional, Tuple
from enum import Enum

logger = logging.getLogger(__name__)


class PIIType(str, Enum):
    """Types of PII/PHI that require redaction"""
    EMAIL = "email"
    PHONE = "phone"
    SSN = "ssn"
    PATIENT_ID = "patient_id"
    MEDICAL_RECORD = "medical_record"
    NAME = "name"
    ADDRESS = "address"
    DATE_OF_BIRTH = "date_of_birth"
    IP_ADDRESS = "ip_address"
    CREDIT_CARD = "credit_card"


class RedactionPattern:
    """Pattern definitions for PII detection"""
    
    # Email addresses
    EMAIL = re.compile(
        r'\b[A-Za-z0-9._%+-]+@[A-Za-z0-9.-]+\.[A-Z|a-z]{2,}\b',
        re.IGNORECASE
    )
    
    # Phone numbers (various formats including Indian)
    # Matches: (555) 123-4567, 555.123.4567, 555-123-4567, +91-9876543210, 9876543210
    PHONE = re.compile(
        r'\b(?:\+?91[-.\s]?)?(?:\+?1[-.\s]?)?\(?\d{3}\)?[-.\s]?\d{3}[-.\s]?\d{4}\b|\b\+?91[-.\s]?\d{10}\b|\b[6-9]\d{9}\b'
    )
    
    # US Social Security Numbers
    SSN = re.compile(
        r'\b(?!000|666|9\d{2})\d{3}-(?!00)\d{2}-(?!0000)\d{4}\b'
    )
    
    # Medical Record Numbers (various formats)
    MEDICAL_RECORD = re.compile(
        r'\b(?:MR|MRN|MEDICAL[-_\s]?REC(?:ORD)?)[:\s#]*([A-Z0-9]{6,12})\b',
        re.IGNORECASE
    )
    
    # Patient IDs
    PATIENT_ID = re.compile(
        r'\b(?:PATIENT[-_\s]?ID|PID|PT[-_\s]?ID)[:\s#]*([A-Z0-9]{6,12})\b',
        re.IGNORECASE
    )
    
    # Dates of Birth (various formats)
    DATE_OF_BIRTH = re.compile(
        r'\b(?:DOB|DATE[-_\s]?OF[-_\s]?BIRTH)[:\s]*(\d{1,2}[/-]\d{1,2}[/-]\d{2,4})\b',
        re.IGNORECASE
    )
    
    # IPv4 addresses
    IP_ADDRESS = re.compile(
        r'\b(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\b'
    )
    
    # Credit card numbers (basic pattern)
    CREDIT_CARD = re.compile(
        r'\b(?:\d{4}[-\s]?){3}\d{4}\b'
    )
    
    # Common name patterns (surname, firstname format)
    NAME_PATTERN = re.compile(
        r'\b(?:Dr\.|Mr\.|Mrs\.|Ms\.)\s+[A-Z][a-z]+\s+[A-Z][a-z]+\b'
    )


class PIIRedactor:
    """
    PII/PHI redaction engine
    
    Features:
    - Pattern-based detection
    - Consistent redaction (same PII gets same hash)
    - Preserves data utility where possible
    - Audit trail of redactions
    """
    
    def __init__(
        self,
        redaction_salt: Optional[str] = None,
        preserve_domain: bool = True,
        enable_audit: bool = True,
    ):
        """
        Initialize PII redactor
        
        Args:
            redaction_salt: Salt for consistent hashing (should be secret)
            preserve_domain: Keep email domains visible (@domain.com)
            enable_audit: Log redaction statistics
        """
        self.redaction_salt = redaction_salt or "phagen_default_salt_change_in_prod"
        self.preserve_domain = preserve_domain
        self.enable_audit = enable_audit
        self.redaction_stats: Dict[PIIType, int] = {}
    
    def _hash_pii(self, value: str, pii_type: PIIType) -> str:
        """
        Create consistent hash of PII value
        
        Args:
            value: The PII value to hash
            pii_type: Type of PII
            
        Returns:
            Deterministic hash for this PII value
        """
        salted = f"{self.redaction_salt}:{pii_type.value}:{value}"
        hash_hex = hashlib.sha256(salted.encode()).hexdigest()[:8]
        return f"[{pii_type.value.upper()}-{hash_hex}]"
    
    def _redact_email(self, email: str) -> str:
        """Redact email address, optionally preserving domain"""
        if self.preserve_domain and '@' in email:
            local, domain = email.split('@', 1)
            hashed_local = self._hash_pii(local, PIIType.EMAIL)
            return f"{hashed_local}@{domain}"
        return self._hash_pii(email, PIIType.EMAIL)
    
    def redact_text(self, text: str, enabled_types: Optional[Set[PIIType]] = None) -> Tuple[str, Dict[PIIType, int]]:
        """
        Redact PII from text
        
        Args:
            text: Text to redact
            enabled_types: Set of PII types to redact (None = all)
            
        Returns:
            Tuple of (redacted_text, redaction_counts)
        """
        if not text:
            return text, {}
        
        enabled_types = enabled_types or set(PIIType)
        redacted = text
        counts: Dict[PIIType, int] = {}
        
        # Email redaction
        if PIIType.EMAIL in enabled_types:
            emails = RedactionPattern.EMAIL.findall(redacted)
            for email in emails:
                redacted = redacted.replace(email, self._redact_email(email))
                counts[PIIType.EMAIL] = counts.get(PIIType.EMAIL, 0) + 1
        
        # Phone number redaction
        if PIIType.PHONE in enabled_types:
            phones = RedactionPattern.PHONE.finditer(redacted)
            replacements = []
            for match in phones:
                phone = match.group(0)
                hashed = self._hash_pii(phone, PIIType.PHONE)
                replacements.append((phone, hashed))
            
            for phone, hashed in replacements:
                redacted = redacted.replace(phone, hashed, 1)
                counts[PIIType.PHONE] = counts.get(PIIType.PHONE, 0) + 1
        
        # SSN redaction
        if PIIType.SSN in enabled_types:
            ssns = RedactionPattern.SSN.findall(redacted)
            for ssn in ssns:
                hashed = self._hash_pii(ssn, PIIType.SSN)
                redacted = redacted.replace(ssn, hashed)
                counts[PIIType.SSN] = counts.get(PIIType.SSN, 0) + 1
        
        # Medical record number redaction
        if PIIType.MEDICAL_RECORD in enabled_types:
            mrns = RedactionPattern.MEDICAL_RECORD.findall(redacted)
            for mrn in mrns:
                hashed = self._hash_pii(mrn, PIIType.MEDICAL_RECORD)
                redacted = redacted.replace(mrn, hashed)
                counts[PIIType.MEDICAL_RECORD] = counts.get(PIIType.MEDICAL_RECORD, 0) + 1
        
        # Patient ID redaction
        if PIIType.PATIENT_ID in enabled_types:
            pids = RedactionPattern.PATIENT_ID.findall(redacted)
            for pid in pids:
                hashed = self._hash_pii(pid, PIIType.PATIENT_ID)
                redacted = redacted.replace(pid, hashed)
                counts[PIIType.PATIENT_ID] = counts.get(PIIType.PATIENT_ID, 0) + 1
        
        # DOB redaction
        if PIIType.DATE_OF_BIRTH in enabled_types:
            dobs = RedactionPattern.DATE_OF_BIRTH.findall(redacted)
            for dob in dobs:
                hashed = self._hash_pii(dob, PIIType.DATE_OF_BIRTH)
                redacted = redacted.replace(dob, hashed)
                counts[PIIType.DATE_OF_BIRTH] = counts.get(PIIType.DATE_OF_BIRTH, 0) + 1
        
        # IP address redaction
        if PIIType.IP_ADDRESS in enabled_types:
            ips = RedactionPattern.IP_ADDRESS.findall(redacted)
            for ip in ips:
                # Preserve first two octets for network analysis
                parts = ip.split('.')
                if len(parts) == 4:
                    redacted_ip = f"{parts[0]}.{parts[1]}.xxx.xxx"
                    redacted = redacted.replace(ip, redacted_ip)
                    counts[PIIType.IP_ADDRESS] = counts.get(PIIType.IP_ADDRESS, 0) + 1
        
        # Credit card redaction
        if PIIType.CREDIT_CARD in enabled_types:
            cards = RedactionPattern.CREDIT_CARD.findall(redacted)
            for card in cards:
                hashed = self._hash_pii(card, PIIType.CREDIT_CARD)
                redacted = redacted.replace(card, hashed)
                counts[PIIType.CREDIT_CARD] = counts.get(PIIType.CREDIT_CARD, 0) + 1
        
        # Update stats if audit enabled
        if self.enable_audit:
            for pii_type, count in counts.items():
                self.redaction_stats[pii_type] = self.redaction_stats.get(pii_type, 0) + count
        
        return redacted, counts
    
    def redact_dict(self, data: Dict[str, Any], enabled_types: Optional[Set[PIIType]] = None) -> Dict[str, Any]:
        """
        Recursively redact PII from dictionary values
        
        Args:
            data: Dictionary to redact
            enabled_types: Set of PII types to redact
            
        Returns:
            Dictionary with redacted values
        """
        redacted = {}
        
        for key, value in data.items():
            if isinstance(value, str):
                redacted[key], _ = self.redact_text(value, enabled_types)
            elif isinstance(value, dict):
                redacted[key] = self.redact_dict(value, enabled_types)
            elif isinstance(value, list):
                redacted[key] = [
                    self.redact_dict(item, enabled_types) if isinstance(item, dict)
                    else self.redact_text(item, enabled_types)[0] if isinstance(item, str)
                    else item
                    for item in value
                ]
            else:
                redacted[key] = value
        
        return redacted
    
    def get_redaction_stats(self) -> Dict[str, int]:
        """Get statistics on redactions performed"""
        return {pii_type.value: count for pii_type, count in self.redaction_stats.items()}
    
    def reset_stats(self) -> None:
        """Reset redaction statistics"""
        self.redaction_stats.clear()


class DLPPolicy:
    """
    Data Loss Prevention policy enforcement
    
    Prevents exfiltration of sensitive data through:
    - Blocking outputs containing unredacted PII
    - Monitoring for suspicious patterns
    - Alerting on policy violations
    """
    
    def __init__(
        self,
        block_on_pii: bool = True,
        max_pii_threshold: int = 5,
        alert_callback = None,
    ):
        """
        Initialize DLP policy
        
        Args:
            block_on_pii: Block operations if PII detected
            max_pii_threshold: Maximum allowed PII instances before blocking
            alert_callback: Function to call on policy violations
        """
        self.block_on_pii = block_on_pii
        self.max_pii_threshold = max_pii_threshold
        self.alert_callback = alert_callback
        self.redactor = PIIRedactor(enable_audit=True)
    
    def scan_for_pii(self, text: str) -> Tuple[bool, Dict[PIIType, int]]:
        """
        Scan text for PII without redacting
        
        Args:
            text: Text to scan
            
        Returns:
            Tuple of (has_pii, pii_counts)
        """
        _, counts = self.redactor.redact_text(text)
        total_pii = sum(counts.values())
        return total_pii > 0, counts
    
    def enforce_policy(self, text: str, operation: str = "unknown") -> Tuple[bool, str]:
        """
        Enforce DLP policy on text
        
        Args:
            text: Text to check
            operation: Description of operation (for logging)
            
        Returns:
            Tuple of (allowed, reason)
        """
        has_pii, counts = self.scan_for_pii(text)
        
        if not has_pii:
            return True, "No PII detected"
        
        total_pii = sum(counts.values())
        
        if total_pii > self.max_pii_threshold:
            reason = f"DLP policy violation: {total_pii} PII instances detected (threshold: {self.max_pii_threshold})"
            
            if self.alert_callback:
                self.alert_callback(operation, counts, reason)
            
            logger.warning(f"[DLP] {reason} in operation: {operation}")
            logger.warning(f"[DLP] PII breakdown: {counts}")
            
            if self.block_on_pii:
                return False, reason
        
        return True, f"PII detected but within threshold: {total_pii}/{self.max_pii_threshold}"
    
    def redact_and_enforce(self, text: str, operation: str = "unknown") -> Tuple[str, bool, str]:
        """
        Combined redaction and policy enforcement
        
        Args:
            text: Text to process
            operation: Description of operation
            
        Returns:
            Tuple of (redacted_text, allowed, reason)
        """
        # First redact
        redacted, counts = self.redactor.redact_text(text)
        
        # Then enforce policy on original text
        allowed, reason = self.enforce_policy(text, operation)
        
        return redacted, allowed, reason


# Global instances
_global_pii_redactor: Optional[PIIRedactor] = None
_global_dlp_policy: Optional[DLPPolicy] = None


def get_pii_redactor() -> PIIRedactor:
    """Get or create global PII redactor instance"""
    global _global_pii_redactor
    if _global_pii_redactor is None:
        _global_pii_redactor = PIIRedactor()
    return _global_pii_redactor


def get_dlp_policy() -> DLPPolicy:
    """Get or create global DLP policy instance"""
    global _global_dlp_policy
    if _global_dlp_policy is None:
        _global_dlp_policy = DLPPolicy()
    return _global_dlp_policy


# Convenience function for quick redaction
def redact_pii(text: str) -> str:
    """Quick PII redaction using global redactor"""
    redactor = get_pii_redactor()
    redacted, _ = redactor.redact_text(text)
    return redacted
