"""
Command Line Interface for Fractal Pangenome Database
=====================================================
Provides easy-to-use commands for genome import, querying, and analysis
"""

import click
import asyncio
import logging
import sys
import json
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

# Import our modules
from .neo4j_genome_importer import (
    GenomeDatabaseBuilder, 
    DockerNeo4jManager,
    PublicGenomeDownloader,
    ShortReadMapper,
    GenomeImporter
)

console = Console()

class GenomeDBCLI:
    """Main CLI class for genome database operations"""
    
    def __init__(self, config_path: str = "genome_db_config.yaml"):
        self.config_path = config_path
        self.db_builder = None
        
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
    
    A command-line interface for managing genomic data in a fractal,
    scalable database system based on Hilbert curve indexing.
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
def database(ctx):
    """Database management commands"""
    pass

@database.command()
@click.option('--force', is_flag=True, help='Force recreate if exists')
@click.pass_context
def init(ctx, force):
    """Initialize the genome database"""
    cli_obj = ctx.obj['cli']
    
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
[bold]Password:[/bold] {neo4j_config.get('auth', {}).get('password', 'genomics123')}
[bold]Data Directory:[/bold] {neo4j_config.get('data_dir', './neo4j_data')}
                """
                console.print(Panel(panel_content, title="Database Connection Info"))
            else:
                console.print("[bold red]✗ Failed to initialize database![/bold red]")
                sys.exit(1)
                
        except Exception as e:
            console.print(f"[bold red]Error initializing database: {e}[/bold red]")
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

@cli.group()
@click.pass_context
def import_cmd(ctx):
    """Import genome data"""
    pass

# Rename the group to avoid Python keyword conflict
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
                if genome_name in downloader.sources:
                    info = downloader.sources[genome_name]
                    category_node.add(f"[cyan]{genome_name}[/cyan] - {info['description']}")
        
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

@cli.group()
@click.pass_context
def map_reads(ctx):
    """Map short reads to the pangenome"""
    pass

# Add the map command to cli
cli.add_command(map_reads, name='map')

@map_reads.command()
@click.option('--fastq', '-f', required=True, type=click.Path(exists=True),
              help='FASTQ file with short reads')
@click.option('--individual', '-i', help='Target individual ID for mapping')
@click.option('--chromosome', '-c', type=int, required=True, help='Chromosome to analyze')
@click.option('--start', '-s', type=int, required=True, help='Start position')
@click.option('--end', '-e', type=int, required=True, help='End position')
@click.option('--output', '-o', type=click.Path(), help='Output file for results')
@click.option('--max-reads', type=int, default=10000, help='Maximum reads to process')
@click.pass_context
def reads(ctx, fastq, individual, chromosome, start, end, output, max_reads):
    """Map short reads and infer genomic paths"""
    
    async def map_reads_task():
        try:
            cli_obj = ctx.obj['cli']
            cli_obj.db_builder = GenomeDatabaseBuilder(cli_obj.config_path)
            await cli_obj.db_builder.initialize_database()
            
            with Progress(
                SpinnerColumn(),
                TextColumn("[progress.description]{task.description}"),
                BarColumn(),
                TaskProgressColumn(),
                console=console
            ) as progress:
                
                # Parse reads
                task1 = progress.add_task("Parsing FASTQ file...", total=1)
                reads_data = await cli_obj.db_builder._parse_fastq_file(fastq)
                
                # Limit reads if specified
                if len(reads_data) > max_reads:
                    reads_data = reads_data[:max_reads]
                    console.print(f"[yellow]Limited to {max_reads} reads for analysis[/yellow]")
                
                progress.advance(task1)
                
                # Map reads
                task2 = progress.add_task("Mapping reads to pangenome...", total=1)
                result = await cli_obj.db_builder.analyze_short_reads(
                    fastq, individual or "unknown", chromosome, start, end
                )
                progress.advance(task2)
                
                # Display mapping results
                console.print(f"[green]✓ Read mapping completed![/green]")
                
                mapping_table = Table(title="Read Mapping Results")
                mapping_table.add_column("Metric", style="cyan")
                mapping_table.add_column("Value", style="white")
                
                mapping_table.add_row("Total Reads", f"{result['total_reads']:,}")
                mapping_table.add_row("Mapped Reads", f"{result['mapped_reads']:,}")
                mapping_table.add_row("Mapping Rate", f"{result['mapping_rate']:.2%}")
                mapping_table.add_row("Region", f"chr{chromosome}:{start:,}-{end:,}")
                
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
            console.print(f"[red]Error during read mapping: {e}[/red]")
    
    asyncio.run(map_reads_task())

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
@click.pass_context
def start(ctx, docker_compose):
    """Start the genome database system"""
    
    if docker_compose:
        console.print("[blue]Starting full stack with Docker Compose...[/blue]")
        try:
            import subprocess
            result = subprocess.run(['docker-compose', 'up', '-d'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                console.print("[green]✓ Full stack started successfully![/green]")
                
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
                """, title="Services Started")
                console.print(services_panel)
            else:
                console.print(f"[red]Error starting services: {result.stderr}[/red]")
        except FileNotFoundError:
            console.print("[red]Docker Compose not found. Please install Docker Compose.[/red]")
    else:
        # Start just the database
        asyncio.run(ctx.invoke(init))

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
