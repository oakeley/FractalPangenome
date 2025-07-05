#!/usr/bin/env python3
"""
Startup script for Fractal Pangenome CLI
This script handles import issues and provides a clean entry point
"""

import sys
import os
from pathlib import Path

# Add current directory to Python path
current_dir = Path(__file__).parent.absolute()
sys.path.insert(0, str(current_dir))

# Now import and run the CLI
try:
    from genome_cli import cli
    
    if __name__ == "__main__":
        cli()
        
except ImportError as e:
    print(f"Import error: {e}")
    print("\nTrying to install missing dependencies...")
    
    # Check if required packages are installed
    missing_packages = []
    
    try:
        import click
    except ImportError:
        missing_packages.append("click")
    
    try:
        import rich
    except ImportError:
        missing_packages.append("rich")
    
    try:
        import yaml
    except ImportError:
        missing_packages.append("pyyaml")
    
    try:
        import docker
    except ImportError:
        missing_packages.append("docker")
    
    if missing_packages:
        print(f"Missing packages: {', '.join(missing_packages)}")
        print(f"Please install them with: pip install {' '.join(missing_packages)}")
        sys.exit(1)
    else:
        print("All dependencies appear to be installed.")
        print("The import error might be due to file structure.")
        print("Please check that all required files are in the same directory.")
        sys.exit(1)
