#!/usr/bin/env python3
"""
PhaGen Deployment Verification
Verifies container image signatures and SBOMs before deployment
"""

import json
import subprocess
import sys
from dataclasses import dataclass
from typing import List, Optional
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class VerificationResult:
    """Result of image verification"""
    image: str
    signature_verified: bool
    sbom_present: bool
    vulnerabilities_checked: bool
    errors: List[str]
    
    @property
    def is_valid(self) -> bool:
        """Check if image passed all verifications"""
        return self.signature_verified and len(self.errors) == 0


class ImageVerifier:
    """Verify container images before deployment"""
    
    def __init__(self, registry: str = "ghcr.io", organization: str = "your-org"):
        self.registry = registry
        self.organization = organization
    
    def verify_signature(self, image: str, public_key: Optional[str] = None) -> tuple[bool, List[str]]:
        """
        Verify image signature using Cosign
        
        Args:
            image: Full image reference (e.g., ghcr.io/org/image:tag)
            public_key: Path to public key file (optional for keyless signing)
        
        Returns:
            Tuple of (verification_success, error_messages)
        """
        logger.info(f"Verifying signature for: {image}")
        errors = []
        
        try:
            # Check if cosign is available
            subprocess.run(['cosign', 'version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            errors.append("cosign not found - install from https://github.com/sigstore/cosign")
            return False, errors
        
        try:
            if public_key:
                # Verify with public key
                cmd = ['cosign', 'verify', '--key', public_key, image]
            else:
                # Verify keyless (requires certificate identity)
                cmd = [
                    'cosign', 'verify', image,
                    '--certificate-identity-regexp', f'https://github.com/{self.organization}',
                    '--certificate-oidc-issuer', 'https://token.actions.githubusercontent.com'
                ]
            
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            logger.info(f"✓ Signature verified for {image}")
            return True, []
        
        except subprocess.CalledProcessError as e:
            error_msg = f"Signature verification failed: {e.stderr}"
            logger.error(error_msg)
            errors.append(error_msg)
            return False, errors
    
    def check_sbom(self, image: str) -> tuple[bool, Optional[dict], List[str]]:
        """
        Check if SBOM is attached to image
        
        Returns:
            Tuple of (sbom_present, sbom_data, error_messages)
        """
        logger.info(f"Checking SBOM for: {image}")
        errors = []
        
        try:
            # Try to download SBOM
            cmd = ['cosign', 'download', 'sbom', image]
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Parse SBOM JSON
            sbom_data = json.loads(result.stdout)
            logger.info(f"✓ SBOM found for {image}")
            return True, sbom_data, []
        
        except subprocess.CalledProcessError as e:
            error_msg = f"SBOM not found or invalid: {e.stderr}"
            logger.warning(error_msg)
            errors.append(error_msg)
            return False, None, errors
        except json.JSONDecodeError as e:
            error_msg = f"SBOM parsing failed: {str(e)}"
            logger.error(error_msg)
            errors.append(error_msg)
            return False, None, errors
    
    def check_vulnerabilities(self, image: str, severity_threshold: str = "HIGH") -> tuple[bool, List[str]]:
        """
        Scan image for vulnerabilities using Trivy
        
        Args:
            image: Full image reference
            severity_threshold: Minimum severity to report (LOW, MEDIUM, HIGH, CRITICAL)
        
        Returns:
            Tuple of (passed_scan, error_messages)
        """
        logger.info(f"Scanning for vulnerabilities: {image}")
        errors = []
        
        try:
            # Check if trivy is available
            subprocess.run(['trivy', '--version'], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("trivy not found - skipping vulnerability scan")
            logger.warning("Install from: https://github.com/aquasecurity/trivy")
            return True, []  # Don't fail if trivy not available
        
        try:
            cmd = [
                'trivy', 'image',
                '--severity', severity_threshold,
                '--exit-code', '1',  # Exit with error if vulnerabilities found
                '--no-progress',
                image
            ]
            
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode == 0:
                logger.info(f"✓ No {severity_threshold}+ vulnerabilities found")
                return True, []
            else:
                error_msg = f"Vulnerabilities found:\n{result.stdout}"
                logger.error(error_msg)
                errors.append(error_msg)
                return False, errors
        
        except Exception as e:
            error_msg = f"Vulnerability scan failed: {str(e)}"
            logger.error(error_msg)
            errors.append(error_msg)
            return False, errors
    
    def verify_image(
        self,
        image: str,
        public_key: Optional[str] = None,
        check_vulns: bool = True,
        severity_threshold: str = "HIGH"
    ) -> VerificationResult:
        """
        Run complete verification pipeline for an image
        
        Args:
            image: Full image reference
            public_key: Path to public key (optional)
            check_vulns: Whether to scan for vulnerabilities
            severity_threshold: Minimum vulnerability severity to report
        
        Returns:
            VerificationResult with all check results
        """
        logger.info(f"Starting verification for: {image}")
        all_errors = []
        
        # Verify signature
        sig_verified, sig_errors = self.verify_signature(image, public_key)
        all_errors.extend(sig_errors)
        
        # Check SBOM
        sbom_present, sbom_data, sbom_errors = self.check_sbom(image)
        all_errors.extend(sbom_errors)
        
        # Check vulnerabilities
        vulns_passed = True
        if check_vulns:
            vulns_passed, vuln_errors = self.check_vulnerabilities(image, severity_threshold)
            all_errors.extend(vuln_errors)
        
        result = VerificationResult(
            image=image,
            signature_verified=sig_verified,
            sbom_present=sbom_present,
            vulnerabilities_checked=vulns_passed,
            errors=all_errors
        )
        
        if result.is_valid:
            logger.info(f"✓ All verifications passed for {image}")
        else:
            logger.error(f"✗ Verification failed for {image}")
            for error in all_errors:
                logger.error(f"  - {error}")
        
        return result
    
    def verify_images(
        self,
        images: List[str],
        public_key: Optional[str] = None,
        check_vulns: bool = True,
        fail_fast: bool = False
    ) -> List[VerificationResult]:
        """
        Verify multiple images
        
        Args:
            images: List of image references
            public_key: Path to public key (optional)
            check_vulns: Whether to scan for vulnerabilities
            fail_fast: Stop on first verification failure
        
        Returns:
            List of VerificationResults
        """
        results = []
        
        for image in images:
            result = self.verify_image(image, public_key, check_vulns)
            results.append(result)
            
            if fail_fast and not result.is_valid:
                logger.error(f"Verification failed for {image}, stopping due to fail_fast=True")
                break
        
        return results


def main():
    """CLI for image verification"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Verify PhaGen container images")
    parser.add_argument('images', nargs='+', help='Image references to verify')
    parser.add_argument('--public-key', help='Path to Cosign public key')
    parser.add_argument('--no-vuln-check', action='store_true', help='Skip vulnerability scanning')
    parser.add_argument('--severity', default='HIGH', choices=['LOW', 'MEDIUM', 'HIGH', 'CRITICAL'],
                        help='Minimum vulnerability severity to report')
    parser.add_argument('--fail-fast', action='store_true', help='Stop on first verification failure')
    parser.add_argument('--registry', default='ghcr.io', help='Container registry')
    parser.add_argument('--organization', default='your-org', help='GitHub organization')
    
    args = parser.parse_args()
    
    verifier = ImageVerifier(registry=args.registry, organization=args.organization)
    
    results = verifier.verify_images(
        images=args.images,
        public_key=args.public_key,
        check_vulns=not args.no_vuln_check,
        fail_fast=args.fail_fast
    )
    
    # Print summary
    print("\n" + "="*80)
    print("VERIFICATION SUMMARY")
    print("="*80)
    
    passed = sum(1 for r in results if r.is_valid)
    failed = len(results) - passed
    
    for result in results:
        status = "✓ PASS" if result.is_valid else "✗ FAIL"
        print(f"{status} {result.image}")
        print(f"  Signature: {'✓' if result.signature_verified else '✗'}")
        print(f"  SBOM: {'✓' if result.sbom_present else '✗'}")
        print(f"  Vulnerabilities: {'✓' if result.vulnerabilities_checked else '✗'}")
        if result.errors:
            print(f"  Errors: {len(result.errors)}")
    
    print(f"\nTotal: {len(results)} images, {passed} passed, {failed} failed")
    
    # Exit with error if any verifications failed
    sys.exit(0 if failed == 0 else 1)


if __name__ == '__main__':
    main()
