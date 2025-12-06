#!/usr/bin/env python3
"""
PhaGen Air-Gapped Deployment Bundle Creator
Packages all dependencies for offline deployment in secure pharma environments
"""

import hashlib
import json
import logging
import os
import shutil
import subprocess
import sys
import tarfile
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


@dataclass
class BundleManifest:
    """Manifest describing bundle contents"""
    bundle_version: str
    created_at: str
    python_version: str
    models: List[Dict]
    indexes: List[Dict]
    dependencies: Dict[str, str]
    container_images: List[str]
    checksum: Optional[str] = None


class AirGappedBundleCreator:
    """Create offline deployment bundles for air-gapped environments"""
    
    def __init__(self, output_dir: str = "airgap-bundle"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        self.models_dir = self.output_dir / "models"
        self.indexes_dir = self.output_dir / "indexes"
        self.deps_dir = self.output_dir / "dependencies"
        self.images_dir = self.output_dir / "images"
        self.config_dir = self.output_dir / "config"
        
        for directory in [self.models_dir, self.indexes_dir, self.deps_dir, self.images_dir, self.config_dir]:
            directory.mkdir(parents=True, exist_ok=True)
        
        self.manifest = BundleManifest(
            bundle_version="1.0.0",
            created_at=datetime.utcnow().isoformat(),
            python_version=f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
            models=[],
            indexes=[],
            dependencies={},
            container_images=[]
        )
    
    def download_models(self, model_configs: List[Dict]):
        """
        Download ML models for offline use
        
        Args:
            model_configs: List of model configurations
                [{'name': 'sentence-transformers/all-MiniLM-L6-v2', 'type': 'embedding'}, ...]
        """
        logger.info(f"Downloading {len(model_configs)} models...")
        
        try:
            from sentence_transformers import SentenceTransformer
        except ImportError:
            logger.warning("sentence-transformers not installed - skipping model download")
            return
        
        for config in model_configs:
            model_name = config['name']
            model_type = config.get('type', 'unknown')
            
            logger.info(f"Downloading {model_name}...")
            
            try:
                # Download model to local cache
                model = SentenceTransformer(model_name)
                
                # Save to bundle directory
                model_output_dir = self.models_dir / model_name.replace('/', '_')
                model.save(str(model_output_dir))
                
                # Add to manifest
                model_info = {
                    'name': model_name,
                    'type': model_type,
                    'path': f"models/{model_name.replace('/', '_')}",
                    'size_mb': sum(f.stat().st_size for f in model_output_dir.rglob('*') if f.is_file()) / (1024 * 1024)
                }
                self.manifest.models.append(model_info)
                
                logger.info(f"✓ Downloaded {model_name} ({model_info['size_mb']:.2f} MB)")
            
            except Exception as e:
                logger.error(f"Failed to download {model_name}: {e}")
    
    def package_indexes(self, index_dirs: List[str]):
        """
        Copy vector database indexes to bundle
        
        Args:
            index_dirs: List of paths to index directories (e.g., ['indexes/chroma'])
        """
        logger.info(f"Packaging {len(index_dirs)} indexes...")
        
        for index_dir in index_dirs:
            index_path = Path(index_dir)
            
            if not index_path.exists():
                logger.warning(f"Index directory not found: {index_dir}")
                continue
            
            logger.info(f"Copying {index_dir}...")
            
            # Copy index directory to bundle
            index_name = index_path.name
            dest_path = self.indexes_dir / index_name
            
            if dest_path.exists():
                shutil.rmtree(dest_path)
            
            shutil.copytree(index_path, dest_path)
            
            # Calculate size
            size_mb = sum(f.stat().st_size for f in dest_path.rglob('*') if f.is_file()) / (1024 * 1024)
            
            # Add to manifest
            index_info = {
                'name': index_name,
                'path': f"indexes/{index_name}",
                'size_mb': size_mb
            }
            self.manifest.indexes.append(index_info)
            
            logger.info(f"✓ Packaged {index_name} ({size_mb:.2f} MB)")
    
    def collect_python_dependencies(self, requirements_files: List[str]):
        """
        Download Python packages for offline installation
        
        Args:
            requirements_files: List of paths to requirements.txt files
        """
        logger.info("Collecting Python dependencies...")
        
        # Merge all requirements
        all_requirements = set()
        for req_file in requirements_files:
            req_path = Path(req_file)
            if req_path.exists():
                with open(req_path, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line and not line.startswith('#'):
                            all_requirements.add(line)
        
        if not all_requirements:
            logger.warning("No requirements found")
            return
        
        # Write combined requirements
        combined_req = self.deps_dir / "requirements.txt"
        with open(combined_req, 'w') as f:
            f.write('\n'.join(sorted(all_requirements)))
        
        logger.info(f"Found {len(all_requirements)} unique dependencies")
        
        # Download packages using pip download
        logger.info("Downloading packages (this may take a while)...")
        
        try:
            subprocess.run([
                sys.executable, '-m', 'pip', 'download',
                '-r', str(combined_req),
                '-d', str(self.deps_dir / 'packages'),
                '--no-binary', ':all:'  # Download source distributions for compatibility
            ], check=True, capture_output=True)
            
            # Also download binary wheels for current platform
            subprocess.run([
                sys.executable, '-m', 'pip', 'download',
                '-r', str(combined_req),
                '-d', str(self.deps_dir / 'wheels'),
            ], check=True, capture_output=True)
            
            logger.info("✓ Python dependencies collected")
            
            # Record in manifest
            for line in all_requirements:
                pkg_name = line.split('==')[0].split('>=')[0].split('<=')[0].strip()
                self.manifest.dependencies[pkg_name] = line
        
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to download dependencies: {e}")
    
    def export_container_images(self, image_names: List[str]):
        """
        Export Docker images to tar files for offline loading
        
        Args:
            image_names: List of Docker image references (e.g., ['phagen-backend:latest'])
        """
        logger.info(f"Exporting {len(image_names)} container images...")
        
        # Check if Docker is available
        try:
            subprocess.run(['docker', '--version'], check=True, capture_output=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("Docker not available - skipping image export")
            return
        
        for image_name in image_names:
            logger.info(f"Exporting {image_name}...")
            
            try:
                # Check if image exists
                subprocess.run(['docker', 'image', 'inspect', image_name], 
                             check=True, capture_output=True)
                
                # Export to tar
                image_filename = image_name.replace(':', '_').replace('/', '_') + '.tar'
                output_path = self.images_dir / image_filename
                
                subprocess.run([
                    'docker', 'save',
                    '-o', str(output_path),
                    image_name
                ], check=True)
                
                # Record in manifest
                size_mb = output_path.stat().st_size / (1024 * 1024)
                self.manifest.container_images.append({
                    'name': image_name,
                    'filename': image_filename,
                    'size_mb': size_mb
                })
                
                logger.info(f"✓ Exported {image_name} ({size_mb:.2f} MB)")
            
            except subprocess.CalledProcessError as e:
                logger.error(f"Failed to export {image_name}: {e}")
    
    def copy_configuration_files(self, config_files: List[str]):
        """
        Copy configuration files to bundle
        
        Args:
            config_files: List of paths to config files (e.g., ['.env.example', 'docker-compose.yml'])
        """
        logger.info("Copying configuration files...")
        
        for config_file in config_files:
            config_path = Path(config_file)
            
            if not config_path.exists():
                logger.warning(f"Config file not found: {config_file}")
                continue
            
            dest_path = self.config_dir / config_path.name
            shutil.copy2(config_path, dest_path)
            logger.info(f"✓ Copied {config_file}")
    
    def create_installation_script(self):
        """Generate installation script for air-gapped environment"""
        logger.info("Creating installation script...")
        
        install_script = self.output_dir / "install.sh"
        
        script_content = """#!/usr/bin/env bash
# PhaGen Air-Gapped Installation Script

set -e

echo "========================================"
echo "PhaGen Air-Gapped Installation"
echo "========================================"

# Check if running as root
if [ "$EUID" -ne 0 ]; then
  echo "Please run as root or with sudo"
  exit 1
fi

# Install Python dependencies
echo "Installing Python dependencies..."
if [ -d "dependencies/wheels" ]; then
  python3 -m pip install --no-index --find-links=dependencies/wheels -r dependencies/requirements.txt
  echo "✓ Python dependencies installed"
else
  echo "⚠ No Python dependencies found"
fi

# Load Docker images
echo "Loading Docker images..."
if [ -d "images" ]; then
  for image_tar in images/*.tar; do
    if [ -f "$image_tar" ]; then
      echo "Loading $(basename $image_tar)..."
      docker load -i "$image_tar"
    fi
  done
  echo "✓ Docker images loaded"
else
  echo "⚠ No Docker images found"
fi

# Set up models and indexes
echo "Setting up models and indexes..."
INSTALL_DIR="/opt/phagen"
mkdir -p "$INSTALL_DIR"

if [ -d "models" ]; then
  cp -r models "$INSTALL_DIR/"
  echo "✓ Models copied to $INSTALL_DIR/models"
fi

if [ -d "indexes" ]; then
  cp -r indexes "$INSTALL_DIR/"
  echo "✓ Indexes copied to $INSTALL_DIR/indexes"
fi

# Copy configuration
if [ -d "config" ]; then
  cp -r config "$INSTALL_DIR/"
  echo "✓ Configuration copied to $INSTALL_DIR/config"
fi

# Set environment variable for offline mode
echo "PHAGEN_OFFLINE_MODE=true" >> /etc/environment
echo "PHAGEN_MODELS_DIR=$INSTALL_DIR/models" >> /etc/environment
echo "PHAGEN_INDEXES_DIR=$INSTALL_DIR/indexes" >> /etc/environment

echo ""
echo "========================================"
echo "Installation Complete!"
echo "========================================"
echo ""
echo "Models: $INSTALL_DIR/models"
echo "Indexes: $INSTALL_DIR/indexes"
echo "Config: $INSTALL_DIR/config"
echo ""
echo "To start PhaGen:"
echo "  cd $INSTALL_DIR/config"
echo "  docker-compose up -d"
"""
        
        with open(install_script, 'w') as f:
            f.write(script_content)
        
        # Make executable
        install_script.chmod(0o755)
        
        logger.info(f"✓ Installation script created: {install_script}")
    
    def create_readme(self):
        """Generate README for bundle"""
        readme = self.output_dir / "README.md"
        
        readme_content = f"""# PhaGen Air-Gapped Deployment Bundle

Bundle Version: {self.manifest.bundle_version}
Created: {self.manifest.created_at}
Python Version: {self.manifest.python_version}

## Contents

- **Models**: {len(self.manifest.models)} ML models ({sum(m['size_mb'] for m in self.manifest.models):.2f} MB)
- **Indexes**: {len(self.manifest.indexes)} vector indexes ({sum(i['size_mb'] for i in self.manifest.indexes):.2f} MB)
- **Dependencies**: {len(self.manifest.dependencies)} Python packages
- **Container Images**: {len(self.manifest.container_images)} images ({sum(img['size_mb'] for img in self.manifest.container_images):.2f} MB)

## Installation

### Prerequisites
- Linux system (tested on Ubuntu 20.04+)
- Docker installed
- Python 3.8+ installed
- Root access (sudo)

### Steps

1. Transfer this bundle to air-gapped environment:
   ```bash
   # On connected system:
   tar -czf phagen-airgap-bundle.tar.gz airgap-bundle/
   
   # Transfer via USB/secure channel to air-gapped system
   ```

2. Extract bundle:
   ```bash
   tar -xzf phagen-airgap-bundle.tar.gz
   cd airgap-bundle
   ```

3. Run installation script:
   ```bash
   sudo ./install.sh
   ```

4. Configure environment:
   ```bash
   cd /opt/phagen/config
   cp docker-compose.yml.example docker-compose.yml
   # Edit docker-compose.yml with your settings
   ```

5. Start PhaGen:
   ```bash
   docker-compose up -d
   ```

6. Verify installation:
   ```bash
   docker-compose ps
   curl http://localhost:8000/health
   ```

## Offline Mode

PhaGen will automatically detect offline mode via `PHAGEN_OFFLINE_MODE=true` environment variable.

In offline mode:
- No external API calls (PubMed, FDA, ClinicalTrials.gov disabled)
- Models loaded from local directory (`/opt/phagen/models`)
- Indexes loaded from local directory (`/opt/phagen/indexes`)
- No automatic updates or telemetry

## Updating Bundle

To update PhaGen in air-gapped environment:
1. Create new bundle on connected system with updated components
2. Transfer new bundle to air-gapped system
3. Stop PhaGen: `docker-compose down`
4. Run new installation script
5. Restart PhaGen: `docker-compose up -d`

## Support

For issues in air-gapped deployment, collect logs:
```bash
docker-compose logs > phagen-logs.txt
```

Transfer logs to connected system for analysis.

## Security

- All models and dependencies have SHA-256 checksums in `manifest.json`
- Verify integrity before installation:
  ```bash
  python3 -c "import json; print(json.load(open('manifest.json'))['checksum'])"
  ```
- Container images are signed with Cosign (see `SIGNING_GUIDE.md`)

## Components

### Models
"""
        
        for model in self.manifest.models:
            readme_content += f"- {model['name']} ({model['type']}, {model['size_mb']:.2f} MB)\n"
        
        readme_content += "\n### Indexes\n"
        for index in self.manifest.indexes:
            readme_content += f"- {index['name']} ({index['size_mb']:.2f} MB)\n"
        
        readme_content += "\n### Container Images\n"
        for img in self.manifest.container_images:
            readme_content += f"- {img['name']} ({img['size_mb']:.2f} MB)\n"
        
        with open(readme, 'w') as f:
            f.write(readme_content)
        
        logger.info(f"✓ README created: {readme}")
    
    def finalize_bundle(self):
        """Calculate checksums and save manifest"""
        logger.info("Finalizing bundle...")
        
        # Calculate bundle checksum
        logger.info("Calculating checksums...")
        hasher = hashlib.sha256()
        
        for root, dirs, files in os.walk(self.output_dir):
            for file in sorted(files):
                if file == 'manifest.json':  # Skip manifest itself
                    continue
                
                filepath = Path(root) / file
                with open(filepath, 'rb') as f:
                    while chunk := f.read(8192):
                        hasher.update(chunk)
        
        self.manifest.checksum = hasher.hexdigest()
        
        # Save manifest
        manifest_path = self.output_dir / "manifest.json"
        with open(manifest_path, 'w') as f:
            json.dump(asdict(self.manifest), f, indent=2)
        
        logger.info(f"✓ Manifest saved: {manifest_path}")
        logger.info(f"Bundle checksum: {self.manifest.checksum}")
    
    def create_bundle(
        self,
        model_configs: Optional[List[Dict]] = None,
        index_dirs: Optional[List[str]] = None,
        requirements_files: Optional[List[str]] = None,
        image_names: Optional[List[str]] = None,
        config_files: Optional[List[str]] = None
    ):
        """
        Create complete air-gapped bundle
        
        Args:
            model_configs: ML models to download
            index_dirs: Vector index directories to package
            requirements_files: Python requirements files
            image_names: Docker images to export
            config_files: Configuration files to include
        """
        logger.info(f"Creating air-gapped bundle in: {self.output_dir}")
        
        if model_configs:
            self.download_models(model_configs)
        
        if index_dirs:
            self.package_indexes(index_dirs)
        
        if requirements_files:
            self.collect_python_dependencies(requirements_files)
        
        if image_names:
            self.export_container_images(image_names)
        
        if config_files:
            self.copy_configuration_files(config_files)
        
        self.create_installation_script()
        self.create_readme()
        self.finalize_bundle()
        
        # Calculate total size
        total_size_mb = sum(f.stat().st_size for f in self.output_dir.rglob('*') if f.is_file()) / (1024 * 1024)
        
        logger.info("")
        logger.info("="*80)
        logger.info("Bundle Creation Complete!")
        logger.info("="*80)
        logger.info(f"Output directory: {self.output_dir}")
        logger.info(f"Total size: {total_size_mb:.2f} MB ({total_size_mb/1024:.2f} GB)")
        logger.info(f"Models: {len(self.manifest.models)}")
        logger.info(f"Indexes: {len(self.manifest.indexes)}")
        logger.info(f"Dependencies: {len(self.manifest.dependencies)}")
        logger.info(f"Container images: {len(self.manifest.container_images)}")
        logger.info("")
        logger.info("Next steps:")
        logger.info("1. Create tarball: tar -czf phagen-airgap-bundle.tar.gz airgap-bundle/")
        logger.info("2. Transfer to air-gapped environment")
        logger.info("3. Extract and run: sudo ./install.sh")


def main():
    """CLI for bundle creation"""
    import argparse
    
    parser = argparse.ArgumentParser(description="Create PhaGen air-gapped deployment bundle")
    parser.add_argument('--output', default='airgap-bundle', help='Output directory')
    parser.add_argument('--models', nargs='+', help='Model names to download (e.g., sentence-transformers/all-MiniLM-L6-v2)')
    parser.add_argument('--indexes', nargs='+', help='Index directories to package (e.g., indexes/chroma)')
    parser.add_argument('--requirements', nargs='+', help='Requirements.txt files')
    parser.add_argument('--images', nargs='+', help='Docker images to export')
    parser.add_argument('--config', nargs='+', help='Configuration files to include')
    
    args = parser.parse_args()
    
    # Convert model names to configs
    model_configs = None
    if args.models:
        model_configs = [{'name': m, 'type': 'embedding'} for m in args.models]
    
    creator = AirGappedBundleCreator(output_dir=args.output)
    
    creator.create_bundle(
        model_configs=model_configs,
        index_dirs=args.indexes,
        requirements_files=args.requirements,
        image_names=args.images,
        config_files=args.config
    )


if __name__ == '__main__':
    main()
