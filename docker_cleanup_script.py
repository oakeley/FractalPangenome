#!/usr/bin/env python3
"""
Docker Cleanup and Debugging Script for Fractal Pangenome
=========================================================
Provides utilities to clean up Docker containers and debug port conflicts
"""

import subprocess
import sys
import time
import json
from typing import List, Dict, Any
import click
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

console = Console()

class DockerCleaner:
    """Utility class for Docker cleanup operations"""
    
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
    
    def check_port_usage(self) -> Dict[int, str]:
        """Check which processes are using our required ports"""
        port_usage = {}
        
        for port in self.ports:
            try:
                # Use lsof to check port usage
                result = subprocess.run(['lsof', '-i', f':{port}'], 
                                      capture_output=True, text=True)
                if result.returncode == 0 and result.stdout.strip():
                    lines = result.stdout.strip().split('\n')[1:]  # Skip header
                    if lines:
                        process_info = lines[0].split()
                        port_usage[port] = f"{process_info[0]} (PID: {process_info[1]})"
                else:
                    port_usage[port] = "Available"
            except FileNotFoundError:
                # lsof not available, try netstat
                try:
                    result = subprocess.run(['netstat', '-tulpn'], 
                                          capture_output=True, text=True)
                    if f":{port}" in result.stdout:
                        port_usage[port] = "In use (unknown process)"
                    else:
                        port_usage[port] = "Available"
                except FileNotFoundError:
                    port_usage[port] = "Cannot check"
        
        return port_usage
    
    def list_genomics_containers(self) -> List[Dict[str, str]]:
        """List all containers related to our genomics stack"""
        try:
            result = subprocess.run(['docker', 'ps', '-a', '--format', 'json'], 
                                  capture_output=True, text=True)
            if result.returncode != 0:
                return []
            
            containers = []
            for line in result.stdout.strip().split('\n'):
                if line.strip():
                    try:
                        container = json.loads(line)
                        # Check if container name matches our patterns
                        if any(name in container.get('Names', '') for name in self.container_names):
                            containers.append(container)
                    except json.JSONDecodeError:
                        continue
            
            return containers
        except Exception:
            return []
    
    def stop_genomics_containers(self, force: bool = False) -> bool:
        """Stop all genomics-related containers"""
        try:
            containers = self.list_genomics_containers()
            if not containers:
                console.print("[green]No genomics containers found running[/green]")
                return True
            
            console.print(f"[yellow]Found {len(containers)} genomics containers[/yellow]")
            
            for container in containers:
                container_name = container.get('Names', 'unknown')
                status = container.get('State', 'unknown')
                
                console.print(f"Stopping container: {container_name} (Status: {status})")
                
                if status == 'running':
                    cmd = ['docker', 'stop', container_name]
                    if force:
                        cmd = ['docker', 'kill', container_name]
                    
                    result = subprocess.run(cmd, capture_output=True, text=True)
                    if result.returncode == 0:
                        console.print(f"[green]✓ Stopped {container_name}[/green]")
                    else:
                        console.print(f"[red]✗ Failed to stop {container_name}: {result.stderr}[/red]")
                
                # Remove the container
                result = subprocess.run(['docker', 'rm', container_name], 
                                      capture_output=True, text=True)
                if result.returncode == 0:
                    console.print(f"[green]✓ Removed {container_name}[/green]")
            
            return True
            
        except Exception as e:
            console.print(f"[red]Error stopping containers: {e}[/red]")
            return False
    
    def cleanup_docker_networks(self) -> bool:
        """Clean up Docker networks"""
        try:
            # Remove our specific network
            result = subprocess.run(['docker', 'network', 'rm', 'fractalpangenome_default'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                console.print("[green]✓ Removed fractalpangenome network[/green]")
            
            # Prune unused networks
            result = subprocess.run(['docker', 'network', 'prune', '-f'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                console.print("[green]✓ Pruned unused networks[/green]")
            
            return True
        except Exception as e:
            console.print(f"[red]Error cleaning networks: {e}[/red]")
            return False
    
    def cleanup_docker_volumes(self) -> bool:
        """Clean up Docker volumes"""
        try:
            # List volumes
            result = subprocess.run(['docker', 'volume', 'ls', '-q'], 
                                  capture_output=True, text=True)
            volumes = result.stdout.strip().split('\n') if result.stdout.strip() else []
            
            # Remove genomics-related volumes
            genomics_volumes = [v for v in volumes if 'genomics' in v or 'fractalpangenome' in v]
            
            for volume in genomics_volumes:
                result = subprocess.run(['docker', 'volume', 'rm', volume], 
                                      capture_output=True, text=True)
                if result.returncode == 0:
                    console.print(f"[green]✓ Removed volume {volume}[/green]")
            
            # Prune unused volumes
            result = subprocess.run(['docker', 'volume', 'prune', '-f'], 
                                  capture_output=True, text=True)
            if result.returncode == 0:
                console.print("[green]✓ Pruned unused volumes[/green]")
            
            return True
        except Exception as e:
            console.print(f"[red]Error cleaning volumes: {e}[/red]")
            return False
    
    def display_system_status(self):
        """Display comprehensive system status"""
        console.print(Panel("[bold blue]System Status Check[/bold blue]", expand=False))
        
        # Docker status
        docker_running = self.check_docker_running()
        docker_status = "🟢 Running" if docker_running else "🔴 Not Running"
        console.print(f"Docker Daemon: {docker_status}")
        
        if not docker_running:
            console.print("[red]Docker is not running. Please start Docker first.[/red]")
            return
        
        # Port usage
        console.print("\n[bold cyan]Port Usage:[/bold cyan]")
        port_usage = self.check_port_usage()
        
        port_table = Table()
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
        
        for port, usage in port_usage.items():
            service = port_services.get(port, "Unknown")
            status_color = "green" if usage == "Available" else "red"
            port_table.add_row(str(port), service, f"[{status_color}]{usage}[/{status_color}]")
        
        console.print(port_table)
        
        # Container status
        console.print("\n[bold cyan]Container Status:[/bold cyan]")
        containers = self.list_genomics_containers()
        
        if containers:
            container_table = Table()
            container_table.add_column("Container", style="cyan")
            container_table.add_column("Status", style="yellow")
            container_table.add_column("Image", style="white")
            
            for container in containers:
                name = container.get('Names', 'unknown')
                status = container.get('State', 'unknown')
                image = container.get('Image', 'unknown')
                
                status_color = "green" if status == "running" else "red"
                container_table.add_row(name, f"[{status_color}]{status}[/{status_color}]", image)
            
            console.print(container_table)
        else:
            console.print("[green]No genomics containers found[/green]")

@click.group()
def cli():
    """Docker cleanup and debugging utilities for Fractal Pangenome"""
    pass

@cli.command()
def status():
    """Check system status"""
    cleaner = DockerCleaner()
    cleaner.display_system_status()

@cli.command()
@click.option('--force', is_flag=True, help='Force kill containers')
def stop(force):
    """Stop all genomics containers"""
    cleaner = DockerCleaner()
    
    if not cleaner.check_docker_running():
        console.print("[red]Docker is not running[/red]")
        return
    
    cleaner.stop_genomics_containers(force)

@cli.command()
@click.option('--volumes', is_flag=True, help='Also clean volumes')
@click.option('--networks', is_flag=True, help='Also clean networks') 
@click.confirmation_option(prompt='This will remove all genomics containers and data. Continue?')
def cleanup(volumes, networks):
    """Complete cleanup of all genomics Docker resources"""
    cleaner = DockerCleaner()
    
    if not cleaner.check_docker_running():
        console.print("[red]Docker is not running[/red]")
        return
    
    console.print("[yellow]Starting complete cleanup...[/yellow]")
    
    # Stop containers
    cleaner.stop_genomics_containers(force=True)
    
    # Clean networks
    if networks:
        cleaner.cleanup_docker_networks()
    
    # Clean volumes
    if volumes:
        cleaner.cleanup_docker_volumes()
    
    console.print("[green]✅ Cleanup complete![/green]")

@cli.command()
def kill_ports():
    """Kill processes using required ports"""
    cleaner = DockerCleaner()
    port_usage = cleaner.check_port_usage()
    
    killed_any = False
    for port, usage in port_usage.items():
        if usage != "Available" and "PID:" in usage:
            try:
                # Extract PID
                pid = usage.split("PID: ")[1].split(")")[0]
                
                console.print(f"Killing process {pid} using port {port}")
                result = subprocess.run(['kill', '-9', pid], capture_output=True)
                
                if result.returncode == 0:
                    console.print(f"[green]✓ Killed process {pid}[/green]")
                    killed_any = True
                else:
                    console.print(f"[red]✗ Failed to kill process {pid}[/red]")
                    
            except Exception as e:
                console.print(f"[red]Error killing process: {e}[/red]")
    
    if not killed_any:
        console.print("[green]No processes needed to be killed[/green]")

if __name__ == '__main__':
    cli()
