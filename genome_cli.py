"""
Command Line Interface for Fractal Pangenome Database
=====================================================
"""

import os
from dotenv import load_dotenv
import click
import asyncio
import logging
import sys
import json
import subprocess
from pathlib import Path
from typing import List, Optional, Dict
import yaml
from rich.console import Console
from rich.table import Table
from rich.progress import Progress, SpinnerColumn, TextColumn, BarColumn, TaskProgressColumn
from rich.panel import Panel
from rich.tree import Tree
import time
from datetime import datetime

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

# Load environment variables
load_dotenv()

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
    
    def check_existing_container(self, container_name: str = "neo4j-genomics") -> Optional[Dict[str, str]]:
        """Check if a container with the given name already exists"""
        try:
            result = subprocess.run(['docker', 'ps', '-a', '--filter', f'name={container_name}', 
                                   '--format', '{{.Names}}\t{{.Status}}\t{{.ID}}'], 
                                  capture_output=True, text=True)
            
            if result.returncode == 0 and result.stdout.strip():
                lines = result.stdout.strip().split('\n')
                for line in lines:
                    parts = line.split('\t')
                    if len(parts) >= 3 and parts[0] == container_name:
                        status = parts[1].lower()
                        container_id = parts[2]
                        
                        # Determine simple status
                        if 'up' in status:
                            simple_status = 'running'
                        elif 'exited' in status or 'stopped' in status:
                            simple_status = 'stopped'
                        else:
                            simple_status = 'unknown'
                        
                        return {
                            'name': container_name,
                            'status': simple_status,
                            'id': container_id,
                            'full_status': status
                        }
            
            return None
            
        except Exception as e:
            console.print(f"[yellow]Warning: Could not check for existing container: {e}[/yellow]")
            return None
    
    def start_existing_container(self, container_name: str = "neo4j-genomics") -> bool:
        """Start an existing stopped container"""
        try:
            result = subprocess.run(['docker', 'start', container_name], 
                                  capture_output=True, text=True)
            
            if result.returncode == 0:
                # Wait a moment for the container to start
                time.sleep(3)
                
                # Verify it's running
                container_info = self.check_existing_container(container_name)
                return container_info and container_info['status'] == 'running'
            else:
                console.print(f"[yellow]Failed to start container: {result.stderr}[/yellow]")
                return False
                
        except Exception as e:
            console.print(f"[yellow]Error starting container: {e}[/yellow]")
            return False

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
@click.option('--cleanup', is_flag=True, default=False, help='Auto-cleanup existing containers')
@click.option('--reuse-existing', is_flag=True, default=True, help='Reuse existing container if available')
@click.pass_context
def init(ctx, force, cleanup, reuse_existing):
    """Initialize the genome database"""
    cli_obj = ctx.obj['cli']
    docker_mgr = cli_obj.docker_manager
    
    # Check Docker
    if not docker_mgr.check_docker_running():
        console.print("[red]Docker is not running. Please start Docker first.[/red]")
        sys.exit(1)
    
    # Check if container already exists and handle appropriately
    if reuse_existing and not force:
        existing_container = docker_mgr.check_existing_container()
        if existing_container:
            if existing_container['status'] == 'running':
                console.print("[green]✓ Found running Neo4j container, using existing instance[/green]")
                config = cli_obj.load_config()
                neo4j_config = config.get('neo4j', {})
                panel_content = f"""
[bold]Neo4j Web Interface:[/bold] http://localhost:7474
[bold]Bolt URI:[/bold] bolt://localhost:7687
[bold]Username:[/bold] neo4j
[bold]Password:[/bold] genomics123
[bold]Data Directory:[/bold] {neo4j_config.get('data_dir', './neo4j_data')}
                """
                console.print(Panel(panel_content, title="Using Existing Database"))
                return
            elif existing_container['status'] in ['stopped', 'exited']:
                console.print("[yellow]Found stopped Neo4j container, starting it...[/yellow]")
                if docker_mgr.start_existing_container():
                    console.print("[green]✓ Successfully started existing container[/green]")
                    return
                else:
                    console.print("[yellow]Failed to start existing container, will create new one[/yellow]")
    
    # Auto-cleanup if requested or if force is specified
    if cleanup or force:
        conflicted_ports = docker_mgr.check_port_conflicts()
        if conflicted_ports or force:
            if conflicted_ports:
                console.print(f"[yellow]Port conflicts detected on: {conflicted_ports}[/yellow]")
            console.print("[yellow]Performing cleanup...[/yellow]")
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
            if "port is already allocated" in str(e).lower() or "already in use" in str(e).lower():
                console.print("[yellow]Container name conflict detected.[/yellow]")
                console.print("[yellow]Try: python genome_cli.py database init --force[/yellow]")
                console.print("[yellow]Or: python genome_cli.py docker cleanup[/yellow]")
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

# Import commands
@cli.group()
@click.pass_context
def import_cmd(ctx):
    """Import genome data"""
    pass

# Register the import command with the name 'import-genome'
cli.add_command(import_cmd, name='import-genome')

@import_cmd.command()
@click.option('--source', '-s', multiple=True, 
              help='Public genome sources to download (can specify multiple)')
@click.option('--list-sources', is_flag=True, help='List available sources')
@click.pass_context
def public(ctx, source, list_sources):
    """Import public reference genomes"""
    cli_obj = ctx.obj['cli']
    
    if list_sources:
        # Display available sources
        downloader = PublicGenomeDownloader()
        
        tree = Tree("[bold blue]Available Public Genomes[/bold blue]")
        
        for category, genomes in [
            ("Reference Genomes", ["GRCh38", "T2T-CHM13"]),
            ("Diverse Assemblies", ["HG002_maternal", "HG002_paternal", "HG00733_maternal", "HG00733_paternal"])
        ]:
            category_node = tree.add(f"[bold green]{category}[/bold green]")
            for genome_name in genomes:
                if hasattr(downloader, 'sources') and genome_name in downloader.sources:
                    info = downloader.sources[genome_name]
                    category_node.add(f"[cyan]{genome_name}[/cyan] - {info['description']}")
                else:
                    category_node.add(f"[cyan]{genome_name}[/cyan]")
        
        console.print(tree)
        return
    
    if not source:
        console.print("[red]Please specify at least one source with --source[/red]")
        console.print("Use --list-sources to see available genomes")
        return
    
    async def import_genomes():
        try:
            # Initialize database builder
            cli_obj.db_builder = GenomeDatabaseBuilder(cli_obj.config_path)
            await cli_obj.db_builder.initialize_database()
            
            # Import genomes with progress tracking
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                console=console
            ) as progress:
                
                task = progress.add_task("Importing genomes...", total=len(source))
                
                for genome_source in source:
                    progress.update(task, description=f"Processing {genome_source}...")
                    
                    # Download and import
                    success = await cli_obj.db_builder.import_public_genomes([genome_source])
                    
                    if success:
                        console.print(f"[green]✓ Successfully imported {genome_source}[/green]")
                    else:
                        console.print(f"[red]✗ Failed to import {genome_source}[/red]")
                    
                    progress.advance(task)
                
            console.print(f"[bold green]Import completed! Imported {len(source)} genomes.[/bold green]")
            
        except Exception as e:
            console.print(f"[red]Error during import: {e}[/red]")
    
    asyncio.run(import_genomes())

@import_cmd.command()
@click.option('--individual-id', '-i', required=True, help='Individual ID')
@click.option('--maternal', '-m', required=True, type=click.Path(exists=True),
              help='Maternal haplotype FASTA file')
@click.option('--paternal', '-p', required=True, type=click.Path(exists=True),
              help='Paternal haplotype FASTA file')
@click.option('--assembly-source', '-a', default='custom', 
              help='Assembly source name')
@click.pass_context
def diploid(ctx, individual_id, maternal, paternal, assembly_source):
    """Import a diploid genome from FASTA files"""
    
    async def import_diploid():
        try:
            cli_obj = ctx.obj['cli']
            cli_obj.db_builder = GenomeDatabaseBuilder(cli_obj.config_path)
            await cli_obj.db_builder.initialize_database()
            
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                console=console
            ) as progress:
                
                task = progress.add_task(f"Importing diploid genome for {individual_id}...")
                
                success = await cli_obj.db_builder.genome_importer.import_diploid_genome(
                    individual_id, maternal, paternal, assembly_source
                )
                
                if success:
                    console.print(f"[green]✓ Successfully imported diploid genome for {individual_id}[/green]")
                    
                    # Display import summary
                    summary_table = Table(title=f"Import Summary - {individual_id}")
                    summary_table.add_column("Property", style="cyan")
                    summary_table.add_column("Value", style="white")
                    
                    summary_table.add_row("Individual ID", individual_id)
                    summary_table.add_row("Assembly Source", assembly_source)
                    summary_table.add_row("Maternal FASTA", str(maternal))
                    summary_table.add_row("Paternal FASTA", str(paternal))
                    summary_table.add_row("Import Time", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
                    
                    console.print(summary_table)
                else:
                    console.print(f"[red]✗ Failed to import diploid genome for {individual_id}[/red]")
                
        except Exception as e:
            console.print(f"[red]Error during diploid import: {e}[/red]")
    
    asyncio.run(import_diploid())

# Query commands
@cli.group()
@click.pass_context
def query(ctx):
    """Query the genome database"""
    pass

@query.command()
@click.option('--individual', '-i', help='Individual ID to query')
@click.option('--chromosome', '-c', type=int, required=True, help='Chromosome number')
@click.option('--start', '-s', type=int, required=True, help='Start position')
@click.option('--end', '-e', type=int, required=True, help='End position')
@click.option('--scale-level', type=int, default=0, help='Scale level (0=nucleotide, 1=gene, etc.)')
@click.option('--output', '-o', type=click.Path(), help='Output file (JSON format)')
@click.pass_context
def region(ctx, individual, chromosome, start, end, scale_level, output):
    """Query a genomic region"""
    
    async def query_region():
        try:
            cli_obj = ctx.obj['cli']
            cli_obj.db_builder = GenomeDatabaseBuilder(cli_obj.config_path)
            await cli_obj.db_builder.initialize_database()
            
            with console.status(f"[bold green]Querying region chr{chromosome}:{start}-{end}..."):
                
                # Perform query using the fractal database
                result = await cli_obj.db_builder.comprehensive_genome_query(
                    individual or "any",
                    chromosome, start, end,
                    include_annotations=True,
                    optimization_level='balanced'
                )
                
                # Display results
                console.print(f"[green]✓ Query completed![/green]")
                
                # Create results table
                results_table = Table(title=f"Query Results - chr{chromosome}:{start}-{end}")
                results_table.add_column("Property", style="cyan")
                results_table.add_column("Value", style="white")
                
                results_table.add_row("Region", f"chr{chromosome}:{start:,}-{end:,}")
                results_table.add_row("Region Size", f"{end-start:,} bp")
                
                if 'diploid_routes' in result:
                    for route_type, path_nodes in result['diploid_routes'].items():
                        results_table.add_row(f"{route_type.title()} Path Length", str(len(path_nodes)))
                
                if 'metadata' in result:
                    metadata = result['metadata']
                    results_table.add_row("Query Time", metadata.get('query_time', 'Unknown'))
                    results_table.add_row("Optimization Level", metadata.get('optimization_level', 'Unknown'))
                
                console.print(results_table)
                
                # Save to file if requested
                if output:
                    with open(output, 'w') as f:
                        json.dump(result, f, indent=2, default=str)
                    console.print(f"[green]Results saved to {output}[/green]")
                
        except Exception as e:
            console.print(f"[red]Error during query: {e}[/red]")
    
    asyncio.run(query_region())

# Map command (flattened for easier use)
@cli.command()
@click.option('--fastq', '-f', multiple=True, required=True, type=click.Path(exists=True),
              help='FASTQ file(s) with genomic short reads (can specify multiple)')
@click.option('--individual', '-i', help='Target individual ID for mapping')
@click.option('--chromosome', '-c', type=int, required=True, help='Chromosome to analyze')
@click.option('--start', '-s', type=int, required=True, help='Start position')
@click.option('--end', '-e', type=int, required=True, help='End position')
@click.option('--output', '-o', type=click.Path(), help='Output file for results')
@click.option('--max-reads', type=int, default=10000, help='Maximum reads to process')
@click.pass_context
def map(ctx, fastq, individual, chromosome, start, end, output, max_reads):
    """Map genomic short reads to the pangenome and infer genomic paths"""
    
    async def map_reads_task():
        try:
            cli_obj = ctx.obj['cli']
            cli_obj.db_builder = GenomeDatabaseBuilder(cli_obj.config_path)
            await cli_obj.db_builder.initialize_database()
            
            console.print(f"[blue]🧬 Mapping {len(fastq)} FASTQ file(s) to pangenome[/blue]")
            
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                console=console
            ) as progress:
                
                # Parse reads from all FASTQ files
                task1 = progress.add_task("Parsing FASTQ files...", total=len(fastq))
                all_reads_data = []
                
                for fastq_file in fastq:
                    progress.update(task1, description=f"Parsing {Path(fastq_file).name}...")
                    reads_data = await cli_obj.db_builder._parse_fastq_file(fastq_file)
                    all_reads_data.extend(reads_data)
                    progress.advance(task1)
                
                # Limit reads if specified
                if len(all_reads_data) > max_reads:
                    all_reads_data = all_reads_data[:max_reads]
                    console.print(f"[yellow]Limited to {max_reads} reads for analysis[/yellow]")
                
                console.print(f"Parsed {len(all_reads_data):,} total reads from {len(fastq)} file(s)")
                
                # Map reads - use the first FASTQ file as reference for the method call
                task2 = progress.add_task("Mapping reads to pangenome...", total=1)
                result = await cli_obj.db_builder.analyze_short_reads(
                    list(fastq)[0], individual or "unknown", chromosome, start, end
                )
                progress.advance(task2)
                
                # Update result with actual read count
                result['total_reads'] = len(all_reads_data)
                result['input_files'] = list(fastq)
                
                # Display mapping results
                console.print(f"[green]✓ Genomic read mapping completed![/green]")
                
                mapping_table = Table(title="Genomic Read Mapping Results")
                mapping_table.add_column("Metric", style="cyan")
                mapping_table.add_column("Value", style="white")
                
                mapping_table.add_row("Input Files", str(len(fastq)))
                mapping_table.add_row("Total Reads", f"{result['total_reads']:,}")
                mapping_table.add_row("Mapped Reads", f"{result['mapped_reads']:,}")
                mapping_table.add_row("Mapping Rate", f"{result['mapping_rate']:.2%}")
                mapping_table.add_row("Target Region", f"chr{chromosome}:{start:,}-{end:,}")
                mapping_table.add_row("Individual", individual or "Any")
                
                # Display inferred paths
                inferred_paths = result.get('inferred_paths', {})
                if inferred_paths:
                    mapping_table.add_row("Inferred Paths", str(len(inferred_paths)))
                    
                    for path_type, path_nodes in inferred_paths.items():
                        console.print(f"[bold]{path_type}:[/bold] {len(path_nodes)} nodes")
                
                console.print(mapping_table)
                
                # Save results if requested
                if output:
                    with open(output, 'w') as f:
                        json.dump(result, f, indent=2, default=str)
                    console.print(f"[green]Results saved to {output}[/green]")
                
        except Exception as e:
            console.print(f"[red]Error during genomic read mapping: {e}[/red]")
            if ctx.obj.get('verbose'):
                import traceback
                console.print(f"[red]{traceback.format_exc()}[/red]")
    
    asyncio.run(map_reads_task())

# RNA-seq commands
@cli.group()
@click.pass_context
def rnaseq(ctx):
    """RNA-seq expression analysis (spliced reads, transcript quantification)"""
    pass

@rnaseq.command()
@click.option('--fastq', '-f', multiple=True, required=True,
              help='FASTQ file(s) with RNA-seq reads (can specify multiple)')
@click.option('--sample-id', '-s', required=True, help='Sample identifier')
@click.option('--gene-regions', '-g', type=click.Path(exists=True),
              help='JSON file with gene regions to analyze')
@click.option('--gene-id', help='Single gene ID (alternative to gene-regions file)')
@click.option('--chromosome', '-c', type=int, help='Chromosome (for single gene)')
@click.option('--start', type=int, help='Start position (for single gene)')
@click.option('--end', type=int, help='End position (for single gene)')
@click.option('--strand', type=click.Choice(['+', '-']), help='Gene strand (for single gene)')
@click.option('--individuals', '-i', multiple=True, 
              help='Individual IDs for eQTL analysis (specify multiple)')
@click.option('--output', '-o', type=click.Path(), default='rnaseq_results.json',
              help='Output file for results')
@click.option('--quantification-method', type=click.Choice(['em', 'unique_only', 'proportional']),
              default='em', help='Expression quantification method')
@click.option('--min-alignment-score', type=float, default=0.8,
              help='Minimum alignment score for read mapping')
@click.option('--eqtl-distance', type=int, default=1000000,
              help='Maximum distance for eQTL analysis (bp)')
@click.pass_context
def analyze(ctx, fastq, sample_id, gene_regions, gene_id, chromosome, start, end, strand,
           individuals, output, quantification_method, min_alignment_score, eqtl_distance):
    """Run complete RNA-seq expression analysis with allele-specific quantification and eQTL discovery"""
    
    async def rnaseq_analysis():
        try:
            # Import RNA-seq analyzer
            try:
                from rnaseq_expression_analyzer import RNASeqPipeline
            except ImportError:
                console.print("[red]Error: rnaseq_expression_analyzer.py not found or has import issues[/red]")
                console.print("[yellow]Make sure the RNA-seq analyzer module is available[/yellow]")
                return
            
            cli_obj = ctx.obj['cli']
            cli_obj.db_builder = GenomeDatabaseBuilder(cli_obj.config_path)
            await cli_obj.db_builder.initialize_database()
            
            # Initialize RNA-seq pipeline
            pipeline = RNASeqPipeline(cli_obj.db_builder)
            
            # Parse gene regions
            if gene_regions:
                with open(gene_regions, 'r') as f:
                    gene_regions_list = json.load(f)
            elif gene_id and chromosome and start and end:
                gene_regions_list = [{
                    'gene_id': gene_id,
                    'chromosome': chromosome,
                    'start': start,
                    'end': end,
                    'strand': strand or '+'
                }]
            else:
                console.print("[red]Error: Must specify either --gene-regions file or --gene-id with coordinates[/red]")
                return
            
            # Setup analysis parameters
            analysis_params = {
                'alignment': {
                    'min_alignment_score': min_alignment_score,
                    'max_mismatches': 3,
                    'allow_multimapping': True
                },
                'quantification': {
                    'method': quantification_method,
                    'min_expression': 0.1
                },
                'eqtl': {
                    'max_distance': eqtl_distance,
                    'min_samples': 20,
                    'significance_threshold': 0.05
                }
            }
            
            console.print(f"[blue]🧬 Starting RNA-seq analysis for sample: {sample_id}[/blue]")
            console.print(f"Processing {len(list(fastq))} FASTQ file(s)")
            console.print(f"Analyzing {len(gene_regions_list)} gene region(s)")
            if individuals:
                console.print(f"eQTL analysis for {len(individuals)} individuals")
            
            # Run analysis
            results = await pipeline.run_complete_analysis(
                fastq_files=list(fastq),
                sample_id=sample_id,
                gene_regions=gene_regions_list,
                individuals=list(individuals) if individuals else None,
                analysis_params=analysis_params
            )
            
            # Display results summary
            console.print(f"\n[green]✅ RNA-seq analysis completed for {sample_id}[/green]")
            
            summary_table = Table(title="RNA-seq Analysis Summary")
            summary_table.add_column("Metric", style="cyan")
            summary_table.add_column("Value", style="white")
            
            summary_table.add_row("Sample ID", sample_id)
            summary_table.add_row("Total Reads", f"{results.get('total_reads', 0):,}")
            summary_table.add_row("Total Alignments", f"{results.get('total_alignments', 0):,}")
            summary_table.add_row("Mapping Rate", f"{results.get('mapping_rate', 0):.1%}")
            
            if 'expression_results' in results:
                expr_count = len(results['expression_results'])
                expressed_count = sum(1 for expr in results['expression_results'].values()
                                    if expr.get('total_expression', 0) > 0)
                summary_table.add_row("Genes Analyzed", str(expr_count))
                summary_table.add_row("Expressed Genes", str(expressed_count))
            
            if 'eqtl_results' in results:
                eqtl_data = results['eqtl_results']
                summary_table.add_row("eQTL Tests", f"{eqtl_data.get('total_tests', 0):,}")
                summary_table.add_row("Significant eQTLs", str(len(eqtl_data.get('significant_eqtls', []))))
            
            console.print(summary_table)
            
            # Show top expressed genes
            if 'expression_results' in results:
                console.print(f"\n[bold]🧬 Top Expressed Genes:[/bold]")
                expr_results = results['expression_results']
                sorted_genes = sorted(expr_results.items(), 
                                    key=lambda x: x[1].get('total_expression', 0), reverse=True)
                
                for i, (gene_id, expr_data) in enumerate(sorted_genes[:5], 1):
                    total_expr = expr_data.get('total_expression', 0)
                    fpkm = expr_data.get('fpkm', 0)
                    console.print(f"  {i}. {gene_id}: Expression={total_expr:.3f}, FPKM={fpkm:.2f}")
                    
                    # Show allele-specific expression if available
                    allele_expr = expr_data.get('allele_specific', {})
                    if allele_expr:
                        for allele, expr_val in list(allele_expr.items())[:3]:  # Show top 3 alleles
                            console.print(f"     • {allele}: {expr_val:.3f}")
            
            # Show top eQTLs
            if 'eqtl_results' in results and results['eqtl_results'].get('top_associations'):
                console.print(f"\n[bold]🔗 Top eQTL Associations:[/bold]")
                for i, eqtl in enumerate(results['eqtl_results']['top_associations'][:5], 1):
                    console.print(f"  {i}. {eqtl['variant_id']} → {eqtl['gene_id']}")
                    console.print(f"     P-value: {eqtl['pvalue']:.2e}, β: {eqtl['beta']:.3f}")
            
            # Save results
            with open(output, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            
            console.print(f"\n[green]📁 Results saved to: {output}[/green]")
            
        except Exception as e:
            console.print(f"[red]❌ RNA-seq analysis error: {e}[/red]")
            if ctx.obj.get('verbose'):
                import traceback
                console.print(f"[red]{traceback.format_exc()}[/red]")
    
    asyncio.run(rnaseq_analysis())

@rnaseq.command()
@click.option('--fastq', '-f', multiple=True, required=True,
              help='FASTQ file(s) with RNA-seq reads')
@click.option('--gene-regions', '-g', type=click.Path(exists=True), required=True,
              help='JSON file with gene regions for transcript building')
@click.option('--output-dir', '-o', type=click.Path(), default='./transcripts',
              help='Output directory for transcript sequences')
@click.pass_context
def build_transcripts(ctx, fastq, gene_regions, output_dir):
    """Build transcript variants from pangenome for RNA-seq analysis"""
    
    async def build_transcripts_task():
        try:
            from rnaseq_expression_analyzer import PangenomeTranscriptBuilder
            
            cli_obj = ctx.obj['cli']
            cli_obj.db_builder = GenomeDatabaseBuilder(cli_obj.config_path)
            await cli_obj.db_builder.initialize_database()
            
            # Load gene regions
            with open(gene_regions, 'r') as f:
                gene_regions_list = json.load(f)
            
            # Initialize transcript builder
            builder = PangenomeTranscriptBuilder(cli_obj.db_builder)
            
            console.print(f"[blue]🧬 Building transcript variants for {len(gene_regions_list)} genes[/blue]")
            
            # Build transcripts
            with console.status("[bold green]Building transcript variants..."):
                transcript_variants = await builder.build_transcript_variants(gene_regions_list)
            
            # Create output directory
            output_path = Path(output_dir)
            output_path.mkdir(exist_ok=True)
            
            # Save transcript sequences
            for transcript_id, sequence in transcript_variants.items():
                transcript_file = output_path / f"{transcript_id}.fasta"
                with open(transcript_file, 'w') as f:
                    f.write(f">{transcript_id}\n")
                    f.write(f"{sequence}\n")
            
            console.print(f"[green]✅ Built {len(transcript_variants)} transcript variants[/green]")
            console.print(f"[green]📁 Transcripts saved to: {output_path}[/green]")
            
            # Display summary
            summary_table = Table(title="Transcript Building Summary")
            summary_table.add_column("Metric", style="cyan")
            summary_table.add_column("Value", style="white")
            
            summary_table.add_row("Input Genes", str(len(gene_regions_list)))
            summary_table.add_row("Transcript Variants", str(len(transcript_variants)))
            summary_table.add_row("Output Directory", str(output_path))
            
            # Calculate average transcript length
            if transcript_variants:
                avg_length = sum(len(seq) for seq in transcript_variants.values()) / len(transcript_variants)
                summary_table.add_row("Average Length", f"{avg_length:.0f} bp")
            
            console.print(summary_table)
            
        except Exception as e:
            console.print(f"[red]❌ Transcript building error: {e}[/red]")
    
    asyncio.run(build_transcripts_task())

# Analysis commands
@cli.group()
@click.pass_context
def analyze(ctx):
    """Analysis and statistics commands"""
    pass

@analyze.command()
@click.pass_context
def stats(ctx):
    """Display database statistics"""
    
    async def get_stats():
        try:
            cli_obj = ctx.obj['cli']
            config = cli_obj.load_config()
            neo4j_manager = DockerNeo4jManager(**config.get('neo4j', {}))
            driver = neo4j_manager.get_connection()
            
            with console.status("[bold green]Gathering database statistics..."):
                with driver.session() as session:
                    # Get node counts
                    node_result = session.run("""
                        MATCH (n:GenomeNode)
                        RETURN 
                            count(n) as total_nodes,
                            count(DISTINCT n.individual_id) as individuals,
                            count(DISTINCT n.chromosome) as chromosomes,
                            avg(n.segment_length) as avg_segment_length
                    """)
                    node_stats = node_result.single()
                    
                    # Get edge counts
                    edge_result = session.run("""
                        MATCH ()-[r:CONNECTS]->()
                        RETURN count(r) as total_edges,
                               count(DISTINCT r.edge_type) as edge_types
                    """)
                    edge_stats = edge_result.single()
                    
                    # Get scale level distribution
                    scale_result = session.run("""
                        MATCH (n:GenomeNode)
                        RETURN n.scale_level as scale_level, count(n) as count
                        ORDER BY n.scale_level
                    """)
                    scale_distribution = [(record['scale_level'], record['count']) 
                                        for record in scale_result]
            
            driver.close()
            
            # Display statistics
            stats_table = Table(title="Database Statistics")
            stats_table.add_column("Metric", style="cyan")
            stats_table.add_column("Value", style="white")
            
            stats_table.add_row("Total Nodes", f"{node_stats['total_nodes']:,}")
            stats_table.add_row("Total Edges", f"{edge_stats['total_edges']:,}")
            stats_table.add_row("Individuals", f"{node_stats['individuals']:,}")
            stats_table.add_row("Chromosomes", f"{node_stats['chromosomes']:,}")
            stats_table.add_row("Edge Types", f"{edge_stats['edge_types']:,}")
            
            if node_stats['avg_segment_length']:
                stats_table.add_row("Avg Segment Length", 
                                  f"{node_stats['avg_segment_length']:.0f} bp")
            
            console.print(stats_table)
            
            # Display scale level distribution
            if scale_distribution:
                scale_table = Table(title="Scale Level Distribution")
                scale_table.add_column("Scale Level", style="cyan")
                scale_table.add_column("Description", style="yellow")
                scale_table.add_column("Node Count", style="white")
                
                scale_descriptions = {
                    0: "Nucleotide Level",
                    1: "Gene Level", 
                    2: "Chromosome Segment",
                    3: "Chromosome Level"
                }
                
                for scale_level, count in scale_distribution:
                    description = scale_descriptions.get(scale_level, "Unknown")
                    scale_table.add_row(str(scale_level), description, f"{count:,}")
                
                console.print(scale_table)
                
        except Exception as e:
            console.print(f"[red]Error getting statistics: {e}[/red]")
    
    asyncio.run(get_stats())

@cli.command()
@click.option('--docker-compose', is_flag=True, help='Start full stack with Docker Compose')
@click.option('--cleanup', is_flag=True, default=False, help='Auto-cleanup before starting')
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
• RNA-seq analysis: python genome_cli.py rnaseq --help
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
        # Start just the database - use reuse-existing by default
        ctx.invoke(init, cleanup=cleanup, reuse_existing=True)

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
    files_to_check = ['docker-compose.yml', 'genome_cli.py', 'neo4j_genome_importer.py', 'rnaseq_expression_analyzer.py']
    console.print("\n[bold cyan]File Checks:[/bold cyan]")
    for file in files_to_check:
        exists = Path(file).exists()
        status = "✅ Found" if exists else "❌ Missing"
        console.print(f"{file}: {status}")
    
    # Available commands
    console.print("\n[bold cyan]Available Commands:[/bold cyan]")
    commands_info = {
        "Database": ["database init", "database status", "database stop"],
        "Import": ["import-genome public", "import-genome diploid"],
        "Query": ["query region"],
        "Genomic Mapping": ["map"],
        "RNA-seq Analysis": ["rnaseq analyze", "rnaseq build-transcripts"],
        "Analysis": ["analyze stats"],
        "Docker": ["docker status", "docker cleanup"],
        "Utilities": ["debug", "version", "start"]
    }
    
    for category, commands in commands_info.items():
        console.print(f"  [yellow]{category}:[/yellow] {', '.join(commands)}")
    
    console.print("\n[bold cyan]Quick Start:[/bold cyan]")
    console.print("1. python genome_cli.py database init")
    console.print("2. python genome_cli.py import-genome public --list-sources")
    console.print("3. python genome_cli.py map --fastq reads.fastq -c 1 -s 1000000 -e 2000000")
    console.print("4. python genome_cli.py rnaseq analyze --help  # for RNA-seq analysis")
    
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

@cli.command()
@click.pass_context  
def version(ctx):
    """Show version information"""
    version_info = {
        "Fractal Pangenome DB": "0.1.0",
        "Neo4j": "5.15-community",
        "Python": sys.version.split()[0]
    }
    
    version_table = Table(title="Version Information")
    version_table.add_column("Component", style="cyan")
    version_table.add_column("Version", style="white")
    
    for component, version in version_info.items():
        version_table.add_row(component, version)
    
    console.print(version_table)

if __name__ == '__main__':
    cli()
