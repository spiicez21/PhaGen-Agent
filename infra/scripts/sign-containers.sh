#!/usr/bin/env bash
# Container image signing and supply-chain verification pipeline
# Uses Cosign for signing and SBOM generation

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${GREEN}[INFO]${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[WARN]${NC} $1"
}

log_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

# Check prerequisites
check_prerequisites() {
    log_info "Checking prerequisites..."
    
    # Check for cosign
    if ! command -v cosign &> /dev/null; then
        log_error "cosign not found. Install from: https://github.com/sigstore/cosign"
        exit 1
    fi
    
    # Check for syft (SBOM generation)
    if ! command -v syft &> /dev/null; then
        log_warn "syft not found. Install from: https://github.com/anchore/syft"
        log_warn "SBOM generation will be skipped"
        SKIP_SBOM=true
    else
        SKIP_SBOM=false
    fi
    
    log_info "Prerequisites check complete"
}

# Generate SBOM for image
generate_sbom() {
    local image=$1
    local output_file=$2
    
    if [ "$SKIP_SBOM" = true ]; then
        log_warn "Skipping SBOM generation (syft not installed)"
        return
    fi
    
    log_info "Generating SBOM for $image..."
    syft "$image" -o spdx-json > "$output_file"
    log_info "SBOM saved to: $output_file"
}

# Sign container image with Cosign
sign_image() {
    local image=$1
    local key_file=${2:-"cosign.key"}
    
    log_info "Signing image: $image"
    
    # Check if key exists
    if [ ! -f "$key_file" ]; then
        log_warn "Key file not found: $key_file"
        log_info "Generating new keypair..."
        cosign generate-key-pair
    fi
    
    # Sign the image
    cosign sign --key "$key_file" "$image"
    
    log_info "Image signed successfully"
}

# Verify signed image
verify_image() {
    local image=$1
    local pub_key=${2:-"cosign.pub"}
    
    log_info "Verifying image signature: $image"
    
    if [ ! -f "$pub_key" ]; then
        log_error "Public key not found: $pub_key"
        return 1
    fi
    
    # Verify signature
    if cosign verify --key "$pub_key" "$image"; then
        log_info "✓ Image signature verified"
        return 0
    else
        log_error "✗ Image signature verification failed"
        return 1
    fi
}

# Attach SBOM to image
attach_sbom() {
    local image=$1
    local sbom_file=$2
    
    if [ ! -f "$sbom_file" ]; then
        log_error "SBOM file not found: $sbom_file"
        return 1
    fi
    
    log_info "Attaching SBOM to image..."
    cosign attach sbom --sbom "$sbom_file" "$image"
    log_info "SBOM attached successfully"
}

# Sign and verify pipeline
sign_and_verify_pipeline() {
    local image=$1
    local key_file=${2:-"cosign.key"}
    local pub_key=${3:-"cosign.pub"}
    
    log_info "Starting sign and verify pipeline for: $image"
    
    # Generate SBOM
    local sbom_file="${image//\//_}_sbom.json"
    sbom_file="${sbom_file//:/_}.json"
    generate_sbom "$image" "$sbom_file"
    
    # Sign image
    sign_image "$image" "$key_file"
    
    # Attach SBOM
    if [ -f "$sbom_file" ]; then
        attach_sbom "$image" "$sbom_file"
    fi
    
    # Verify signature
    verify_image "$image" "$pub_key"
    
    log_info "Pipeline complete for: $image"
}

# Generate provenance attestation
generate_provenance() {
    local image=$1
    local git_commit=$(git rev-parse HEAD 2>/dev/null || echo "unknown")
    local git_branch=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "unknown")
    local build_date=$(date -u +"%Y-%m-%dT%H:%M:%SZ")
    
    local provenance_file="${image//\//_}_provenance.json"
    provenance_file="${provenance_file//:/_}.json"
    
    log_info "Generating provenance attestation..."
    
    cat > "$provenance_file" <<EOF
{
  "image": "$image",
  "git_commit": "$git_commit",
  "git_branch": "$git_branch",
  "build_date": "$build_date",
  "builder": "$(whoami)@$(hostname)",
  "build_platform": "$(uname -s)/$(uname -m)"
}
EOF
    
    log_info "Provenance saved to: $provenance_file"
    
    # Attach provenance as attestation
    cosign attest --key cosign.key --predicate "$provenance_file" "$image" || log_warn "Attestation attach failed"
}

# Main function
main() {
    log_info "==================================================="
    log_info "PhaGen Container Signing Pipeline"
    log_info "==================================================="
    
    check_prerequisites
    
    # Default images to sign
    IMAGES=(
        "phagen-backend:latest"
        "phagen-frontend:latest"
        "phagen-rdkit-service:latest"
    )
    
    # Check if custom images provided
    if [ $# -gt 0 ]; then
        IMAGES=("$@")
    fi
    
    log_info "Images to sign: ${IMAGES[*]}"
    
    # Process each image
    for image in "${IMAGES[@]}"; do
        log_info ""
        log_info "Processing: $image"
        
        # Check if image exists locally
        if ! docker image inspect "$image" &> /dev/null; then
            log_warn "Image not found locally: $image"
            log_warn "Build the image first: docker build -t $image ."
            continue
        fi
        
        # Run pipeline
        sign_and_verify_pipeline "$image"
        
        # Generate provenance
        generate_provenance "$image"
    done
    
    log_info ""
    log_info "==================================================="
    log_info "Pipeline complete!"
    log_info "==================================================="
    log_info ""
    log_info "Next steps:"
    log_info "1. Push signed images: docker push <image>"
    log_info "2. Distribute public key (cosign.pub) to consumers"
    log_info "3. Verify on pull: cosign verify --key cosign.pub <image>"
}

# Script usage
usage() {
    cat <<EOF
Usage: $0 [image1] [image2] ...

Sign container images with Cosign and generate SBOM.

Examples:
  $0                                    # Sign default PhaGen images
  $0 phagen-backend:latest              # Sign specific image
  $0 phagen-backend:v1.0 phagen-frontend:v1.0  # Sign multiple images

Prerequisites:
  - cosign: https://github.com/sigstore/cosign
  - syft: https://github.com/anchore/syft (optional, for SBOM)
  - docker: Container images must be built locally

Environment:
  COSIGN_PASSWORD: Password for private key (optional)
  SKIP_SBOM: Set to 'true' to skip SBOM generation
EOF
}

# Handle help flag
if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    usage
    exit 0
fi

# Run main
main "$@"
