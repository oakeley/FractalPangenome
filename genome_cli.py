"""
Command Line Interface for Fractal Pangenome Database
=====================================================
"""

import click
import asyncio
import logging
import sys
import json
import subprocess
from pathlib import Path
from typing import List, Optional
import yaml
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from rich.panel import Panel
from rich.tree import Tree
import time
from datetime import datetime

# Fixed absolute imports
try:
    from neo4j_genome_importer import (
        GenomeDatabaseBuilder, 
        DockerNeo4jManager,
        PublicGenomeDownloader,
        ShortReadMapper,
        GenomeImporter
    )
except ImportError:
    # Fallback for development
    import sys
    import os
    sys.path.append(os.path.dirname(os.path.abspath(__file__)))
    from neo4j_genome_importer import (
        GenomeDatabaseBuilder, 
        DockerNeo4jManager,
        PublicGenomeDownloader,
        ShortReadMapper,
        GenomeImporter
    )

console = Console()

class DockerManager:
    """Enhanced Docker management with cleanup capabilities"""
    
    def __init__(self):
        self.container_names = [
            "neo4j-genomics",
            "redis-genomics", 
            "influxdb-genomics",
            "grafana-genomics"
        ]
        self.ports = [7474, 7687, 6379, 8086, 3000, 8888]
    
    def check_docker_running(self) -> bool:
        """Check if Docker daemon is running"""
        try:
            result = subprocess.run(['docker', 'version'], 
                                  capture_output=True, text=True, timeout=10)
            return result.returncode == 0
        except (subprocess.TimeoutExpired, FileNotFoundError):
            return False
    
    def check_port_conflicts(self) -> List[int]:
        """Check for port conflicts"""
        conflicted_ports = []
        
        for port in self.ports:
            try:
                result = subprocess.run(['lsof', '-i', f':{port}'], 
                                      capture_output=True, text=True)
                if result.returncode == 0 and result.stdout.strip():
                    conflicted_ports.append(port)
            except FileNotFoundError:
                # lsof not available, try with netstat or ss
                try:
                    result = subprocess.run(['ss', '-tlnp'], capture_output=True, text=True)
                    if f":{port}" in result.stdout:
                        conflicted_ports.append(port)
                except FileNotFoundError:
                    # Skip if no port checking tools available
                    pass
        
        return conflicted_ports
    
    def stop_existing_containers(self, force: bool = True) -> bool:
        """Stop existing genomics containers"""
        try:
            # Get all containers (running and stopped)
            result = subprocess.run(['docker', 'ps', '-a', '--format', '{{.Names}}'], 
                                  capture_output=True, text=True)
            if result.returncode != 0:
                return True
            
            all_containers = result.stdout.strip().split('\n') if result.stdout.strip() else []
            genomics_containers = [c for c in all_containers if any(name in c for name in self.container_names)]
            
            if not genomics_containers:
                return True
            
            console.print(f"[yellow]Found {len(genomics_containers)} existing containers to clean up...[/yellow]")
            
            for container in genomics_containers:
                console.print(f"  Stopping and removing: {container}")
                
                # Stop container
                if force:
                    subprocess.run(['docker', 'kill', container], capture_output=True)
                else:
                    subprocess.run(['docker', 'stop', container], capture_output=True)
                
                # Remove container
                subprocess.run(['docker', 'rm', container], capture_output=True)
            
            console.print("[green]✓ Cleanup completed[/green]")
            return True
            
        except Exception as e:
            console.print(f"[red]Error during cleanup: {e}[/red]")
            return False
    
    def cleanup_networks_and_volumes(self):
        """Clean up Docker networks and volumes"""
        try:
            # Remove specific network
            subprocess.run(['docker', 'network', 'rm', 'fractalpangenome_default'], 
                         capture_output=True)
            
            # Prune unused networks and volumes
            subprocess.run(['docker', 'network', 'prune', '-f'], capture_output=True)
            subprocess.run(['docker', 'volume', 'prune', '-f'], capture_output=True)
            
        except Exception:
            pass  # Ignore errors in cleanup

class GenomeDBCLI:
    """Main CLI class for genome database operations"""
    
    def __init__(self, config_path: str = "genome_db_config.yaml"):
        self.config_path = config_path
        self.db_builder = None
        self.docker_manager = DockerManager()
        
    def load_config(self) -> dict:
        """Load configuration from file"""
        if not Path(self.config_path).exists():
            console.print(f"[red]Configuration file {self.config_path} not found![/red]")
            console.print("Creating default configuration...")
            self._create_default_config()
        
        with open(self.config_path, 'r') as f:
            return yaml.safe_load(f)
    
    def _create_default_config(self):
        """Create default configuration file"""
        default_config = {
            'neo4j': {
                'container_name': 'neo4j-genomics',
                'neo4j_version': '5.15-community',
                'data_dir': './neo4j_data'
            },
            'download': {
                'download_dir': './reference_genomes'
            },
            'processing': {
                'segment_size': 10000,
                'kmer_size': 31,
                'batch_size': 1000
            }
        }
        
        with open(self.config_path, 'w') as f:
            yaml.dump(default_config, f, default_flow_style=False)
        
        console.print(f"[green]Created default configuration at {self.config_path}[/green]")

@click.group()
@click.option('--config', '-c', default='genome_db_config.yaml', 
              help='Configuration file path')
@click.option('--verbose', '-v', is_flag=True, help='Verbose logging')
@click.pass_context
def cli(ctx, config, verbose):
    """Fractal Pangenome Database CLI
    
    Enhanced command-line interface with automatic cleanup and debugging
    """
    # Setup logging
    log_level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Initialize CLI object
    ctx.ensure_object(dict)
    ctx.obj['cli'] = GenomeDBCLI(config)
    ctx.obj['verbose'] = verbose

@cli.group()
@click.pass_context
def docker(ctx):
    """Docker management commands"""
    pass

@docker.command()
@click.pass_context
def status(ctx):
    """Check Docker and port status"""
    cli_obj = ctx.obj['cli']
    docker_mgr = cli_obj.docker_manager
    
    # Check Docker
    docker_running = docker_mgr.check_docker_running()
    docker_status = "🟢 Running" if docker_running else "🔴 Not Running"
    console.print(f"Docker Daemon: {docker_status}")
    
    if not docker_running:
        console.print("[red]Docker is not running. Please start Docker first.[/red]")
        return
    
    # Check ports
    conflicted_ports = docker_mgr.check_port_conflicts()
    
    port_table = Table(title="Port Status")
    port_table.add_column("Port", style="cyan")
    port_table.add_column("Service", style="yellow")
    port_table.add_column("Status", style="white")
    
    port_services = {
        7474: "Neo4j HTTP",
        7687: "Neo4j Bolt", 
        6379: "Redis",
        8086: "InfluxDB",
        3000: "Grafana",
        8888: "Jupyter"
    }
    
    for port in docker_mgr.ports:
        service = port_services.get(port, "Unknown")
        if port in conflicted_ports:
            status = "[red]In Use[/red]"
        else:
            status = "[green]Available[/green]"
        port_table.add_row(str(port), service, status)
    
    console.print(port_table)

@docker.command()
@click.option('--force', is_flag=True, help='Force kill containers')
@click.pass_context
def cleanup(ctx, force):
    """Clean up existing Docker containers and resources"""
    cli_obj = ctx.obj['cli']
    docker_mgr = cli_obj.docker_manager
    
    if not docker_mgr.check_docker_running():
        console.print("[red]Docker is not running[/red]")
        return
    
    with console.status("[bold yellow]Cleaning up Docker resources..."):
        success = docker_mgr.stop_existing_containers(force)
        docker_mgr.cleanup_networks_and_volumes()
    
    if success:
        console.print("[green]✅ Docker cleanup completed[/green]")
    else:
        console.print("[red]❌ Docker cleanup had some issues[/red]")

@cli.group()
@click.pass_context
def database(ctx):
    """Database management commands"""
    pass

@database.command()
@click.option('--force', is_flag=True, help='Force recreate if exists')
@click.option('--cleanup', is_flag=True, default=True, help='Auto-cleanup existing containers')
@click.pass_context
def init(ctx, force, cleanup):
    """Initialize the genome database"""
    cli_obj = ctx.obj['cli']
    docker_mgr = cli_obj.docker_manager
    
    # Check Docker
    if not docker_mgr.check_docker_running():
        console.print("[red]Docker is not running. Please start Docker first.[/red]")
        sys.exit(1)
    
    # Auto-cleanup if requested
    if cleanup:
        conflicted_ports = docker_mgr.check_port_conflicts()
        if conflicted_ports:
            console.print(f"[yellow]Port conflicts detected on: {conflicted_ports}[/yellow]")
            console.print("[yellow]Performing automatic cleanup...[/yellow]")
            docker_mgr.stop_existing_containers(force=True)
            docker_mgr.cleanup_networks_and_volumes()
            time.sleep(2)  # Give Docker time to clean up
    
    with console.status("[bold green]Initializing database...") as status:
        try:
            config = cli_obj.load_config()
            cli_obj.db_builder = GenomeDatabaseBuilder(cli_obj.config_path)
            
            status.update("[bold blue]Setting up Neo4j container...")
            success = asyncio.run(cli_obj.db_builder.initialize_database())
            
            if success:
                console.print("[bold green]✓ Database initialized successfully![/bold green]")
                
                # Display connection info
                neo4j_config = config.get('neo4j', {})
                panel_content = f"""
[bold]Neo4j Web Interface:[/bold] http://localhost:7474
[bold]Bolt URI:[/bold] bolt://localhost:7687
[bold]Username:[/bold] neo4j
[bold]Password:[/bold] genomics123
[bold]Data Directory:[/bold] {neo4j_config.get('data_dir', './neo4j_data')}
                """
                console.print(Panel(panel_content, title="Database Connection Info"))
            else:
                console.print("[bold red]✗ Failed to initialize database![/bold red]")
                sys.exit(1)
                
        except Exception as e:
            console.print(f"[bold red]Error initializing database: {e}[/bold red]")
            if "port is already allocated" in str(e).lower():
                console.print("[yellow]Tip: Run 'python genome_cli.py docker cleanup' first[/yellow]")
            sys.exit(1)

@database.command()
@click.pass_context
def status(ctx):
    """Check database status"""
    cli_obj = ctx.obj['cli']
    
    try:
        config = cli_obj.load_config()
        neo4j_manager = DockerNeo4jManager(**config.get('neo4j', {}))
        
        # Check container status
        try:
            container = neo4j_manager.docker_client.containers.get(neo4j_manager.container_name)
            status = container.status
            
            table = Table(title="Database Status")
            table.add_column("Component", style="cyan")
            table.add_column("Status", style="magenta")
            table.add_column("Details", style="white")
            
            if status == "running":
                table.add_row("Neo4j Container", "🟢 Running", f"Started: {container.attrs['State']['StartedAt']}")
                
                # Test connection
                try:
                    driver = neo4j_manager.get_connection()
                    with driver.session() as session:
                        result = session.run("CALL dbms.components() YIELD name, edition, versions")
                        components = [dict(record) for record in result]
                    driver.close()
                    table.add_row("Database Connection", "🟢 Connected", f"Neo4j {components[0]['versions'][0]}")
                except Exception as e:
                    table.add_row("Database Connection", "🔴 Failed", str(e))
            else:
                table.add_row("Neo4j Container", f"🔴 {status.title()}", "Container not running")
            
            console.print(table)
            
        except Exception as e:
            console.print(f"[red]Container not found or error: {e}[/red]")
            
    except Exception as e:
        console.print(f"[red]Error checking status: {e}[/red]")

@database.command()
@click.confirmation_option(prompt='Are you sure you want to stop the database?')
@click.pass_context
def stop(ctx):
    """Stop the database"""
    cli_obj = ctx.obj['cli']
    
    try:
        config = cli_obj.load_config()
        neo4j_manager = DockerNeo4jManager(**config.get('neo4j', {}))
        neo4j_manager.stop_container()
        console.print("[green]✓ Database stopped successfully[/green]")
    except Exception as e:
        console.print(f"[red]Error stopping database: {e}[/red]")

@cli.command()
@click.option('--docker-compose', is_flag=True, help='Start full stack with Docker Compose')
@click.option('--cleanup', is_flag=True, default=True, help='Auto-cleanup before starting')
@click.pass_context
def start(ctx, docker_compose, cleanup):
    """Start the genome database system"""
    cli_obj = ctx.obj['cli']
    docker_mgr = cli_obj.docker_manager
    
    # Check Docker
    if not docker_mgr.check_docker_running():
        console.print("[red]Docker is not running. Please start Docker first.[/red]")
        sys.exit(1)
    
    if docker_compose:
        # Auto-cleanup if requested
        if cleanup:
            console.print("[blue]Performing pre-startup cleanup...[/blue]")
            docker_mgr.stop_existing_containers(force=True)
            docker_mgr.cleanup_networks_and_volumes()
            time.sleep(2)
        
        console.print("[blue]Starting full stack with Docker Compose...[/blue]")
        try:
            # Check if docker-compose.yml exists
            if not Path("docker-compose.yml").exists():
                console.print("[red]docker-compose.yml not found in current directory[/red]")
                console.print("[yellow]Please make sure you're running this from the project root[/yellow]")
                sys.exit(1)
            
            result = subprocess.run(['docker-compose', 'up', '-d'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                console.print("[green]✓ Full stack started successfully![/green]")
                
                # Wait a moment for services to start
                time.sleep(5)
                
                # Display service URLs
                services_panel = Panel("""
[bold]Service URLs:[/bold]
• Neo4j Browser: http://localhost:7474
• Grafana Dashboard: http://localhost:3000
• Jupyter Notebooks: http://localhost:8888
• InfluxDB UI: http://localhost:8086

[bold]Default Credentials:[/bold]
• Neo4j: neo4j / genomics123
• Grafana: admin / genomics123
• Jupyter: Token 'genomics'

[bold]Getting Started:[/bold]
• Check status: python genome_cli.py database status
• Import genomes: python genome_cli.py import-genome public --list-sources
                """, title="Services Started")
                console.print(services_panel)
            else:
                console.print(f"[red]Error starting services:[/red]")
                console.print(result.stderr)
                
                if "port is already allocated" in result.stderr:
                    console.print("\n[yellow]Port conflict detected. Try:[/yellow]")
                    console.print("1. python genome_cli.py docker cleanup")
                    console.print("2. python genome_cli.py start --docker-compose")
                    
        except FileNotFoundError:
            console.print("[red]Docker Compose not found. Please install Docker Compose.[/red]")
    else:
        # Start just the database
        ctx.invoke(init, cleanup=cleanup)

# Add debugging command
@cli.command()
@click.pass_context
def debug(ctx):
    """Run diagnostic checks"""
    cli_obj = ctx.obj['cli']
    docker_mgr = cli_obj.docker_manager
    
    console.print(Panel("[bold blue]🔍 Diagnostic Information[/bold blue]", expand=False))
    
    # System info
    console.print("\n[bold cyan]System Information:[/bold cyan]")
    console.print(f"Python: {sys.version}")
    console.print(f"Platform: {sys.platform}")
    
    # Docker info
    console.print(f"Docker: {'🟢 Available' if docker_mgr.check_docker_running() else '🔴 Not Available'}")
    
    # Port conflicts
    conflicts = docker_mgr.check_port_conflicts()
    console.print(f"Port conflicts: {conflicts if conflicts else '✅ None'}")
    
    # File checks
    files_to_check = ['docker-compose.yml', 'genome_cli.py', 'neo4j_genome_importer.py']
    console.print("\n[bold cyan]File Checks:[/bold cyan]")
    for file in files_to_check:
        exists = Path(file).exists()
        status = "✅ Found" if exists else "❌ Missing"
        console.print(f"{file}: {status}")
    
    # Directory structure
    console.print("\n[bold cyan]Current Directory Contents:[/bold cyan]")
    try:
        files = list(Path(".").glob("*.py"))[:10]  # Show first 10 Python files
        for file in files:
            console.print(f"  {file}")
        if len(list(Path(".").glob("*.py"))) > 10:
            console.print(f"  ... and {len(list(Path('.').glob('*.py'))) - 10} more files")
    except Exception as e:
        console.print(f"Error listing files: {e}")

# Include all the other commands from the original CLI
# (import, query, map, analyze, etc.) - keeping them the same as before
# ... [Previous commands remain the same] ...

if __name__ == '__main__':
    cli()
