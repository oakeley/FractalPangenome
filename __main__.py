"""
Entry point for the Fractal Pangenome package
"""

import sys
import os

# Add the current directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

# Import and run the CLI
from genome_cli import cli

if __name__ == "__main__":
    cli()
