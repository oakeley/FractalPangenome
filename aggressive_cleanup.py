#!/usr/bin/env python3
"""
Aggressive Cleanup Script for Port Conflicts
============================================
Finds and kills ANY process using our required ports
"""

import subprocess
import sys
import time
import signal
from typing import List, Dict, Tuple
from rich.console import Console
from rich.table import Table
from rich.panel import Panel

console = Console()

class AggressivePortCleaner:
    """Aggressively clean up port conflicts"""
    
    def __init__(self):
        self.required_ports = [7474, 7687, 6379, 8086, 3000, 8888]
        self.port_services = {
            7474: "Neo4j HTTP",
            7687: "Neo4j Bolt", 
            6379: "Redis",
            8086: "InfluxDB",
            3000: "Grafana",
            8888: "Jupyter"
        }
    
    def find_port_processes(self, port: int) -> List[Dict[str, str]]:
        """Find all processes using a specific port"""
        processes = []
        
        # Try multiple methods to find processes
        methods = [
            self._find_with_lsof,
            self._find_with_netstat,
            self._find_with_ss,
            self._find_with_fuser
        ]
        
        for method in methods:
            try:
                found = method(port)
                processes.extend(found)
            except Exception:
                continue
        
        # Remove duplicates
        seen_pids = set()
        unique_processes = []
        for proc in processes:
            if proc['pid'] not in seen_pids:
                unique_processes.append(proc)
                seen_pids.add(proc['pid'])
        
        return unique_processes
    
    def _find_with_lsof(self, port: int) -> List[Dict[str, str]]:
        """Find processes using lsof"""
        result = subprocess.run(['lsof', '-i', f':{port}'], 
                              capture_output=True, text=True)
        processes = []
        
        if result.returncode == 0 and result.stdout.strip():
            lines = result.stdout.strip().split('\n')[1:]  # Skip header
            for line in lines:
                parts = line.split()
                if len(parts) >= 2:
                    processes.append({
                        'command': parts[0],
                        'pid': parts[1],
                        'user': parts[2] if len(parts) > 2 else 'unknown',
                        'method': 'lsof'
                    })
        
        return processes
    
    def _find_with_netstat(self, port: int) -> List[Dict[str, str]]:
        """Find processes using netstat"""
        result = subprocess.run(['netstat', '-tulpn'], 
                              capture_output=True, text=True)
        processes = []
        
        if result.returncode == 0:
            for line in result.stdout.split('\n'):
                if f':{port}' in line and 'LISTEN' in line:
                    # Extract PID from the last column
                    parts = line.split()
                    if len(parts) > 6:
                        pid_part = parts[-1]
                        if '/' in pid_part:
                            pid = pid_part.split('/')[0]
                            command = pid_part.split('/')[1] if '/' in pid_part else 'unknown'
                            if pid.isdigit():
                                processes.append({
                                    'command': command,
                                    'pid': pid,
                                    'user': 'unknown',
                                    'method': 'netstat'
                                })
        
        return processes
    
    def _find_with_ss(self, port: int) -> List[Dict[str, str]]:
        """Find processes using ss (socket statistics)"""
        result = subprocess.run(['ss', '-tulpn'], 
                              capture_output=True, text=True)
        processes = []
        
        if result.returncode == 0:
            for line in result.stdout.split('\n'):
                if f':{port}' in line and 'LISTEN' in line:
                    # ss output format varies, try to extract PID
                    if 'pid=' in line:
                        pid_part = line.split('pid=')[1].split(',')[0].split(')')[0]
                        if pid_part.isdigit():
                            processes.append({
                                'command': 'unknown',
                                'pid': pid_part,
                                'user': 'unknown',
                                'method': 'ss'
                            })
        
        return processes
    
    def _find_with_fuser(self, port: int) -> List[Dict[str, str]]:
        """Find processes using fuser"""
        result = subprocess.run(['fuser', f'{port}/tcp'], 
                              capture_output=True, text=True)
        processes = []
        
        if result.returncode == 0 and result.stdout.strip():
            pids = result.stdout.strip().split()
            for pid in pids:
                if pid.isdigit():
                    processes.append({
                        'command': 'unknown',
                        'pid': pid,
                        'user': 'unknown',
                        'method': 'fuser'
                    })
        
        return processes
    
    def kill_process(self, pid: str, force: bool = False) -> bool:
        """Kill a process by PID"""
        try:
            signal_type = signal.SIGKILL if force else signal.SIGTERM
            
            # Try using kill command
            result = subprocess.run(['kill', f'-{signal_type}', pid], 
                                  capture_output=True, text=True)
            
            if result.returncode == 0:
                return True
            
            # If that fails, try with sudo
            result = subprocess.run(['sudo', 'kill', f'-{signal_type}', pid], 
                                  capture_output=True, text=True)
            
            return result.returncode == 0
            
        except Exception:
            return False
    
    def display_port_status(self):
        """Display current port usage"""
        console.print(Panel("[bold blue]🔍 Port Status Analysis[/bold blue]", expand=False))
        
        table = Table()
        table.add_column("Port", style="cyan")
        table.add_column("Service", style="yellow")
        table.add_column("Status", style="white")
        table.add_column("Process Details", style="green")
        
        for port in self.required_ports:
            service = self.port_services[port]
            processes = self.find_port_processes(port)
            
            if processes:
                status = "[red]OCCUPIED[/red]"
                details = []
                for proc in processes:
                    details.append(f"{proc['command']} (PID: {proc['pid']}, Method: {proc['method']})")
                detail_str = "; ".join(details)
            else:
                status = "[green]FREE[/green]"
                detail_str = "Available"
            
            table.add_row(str(port), service, status, detail_str)
        
        console.print(table)
    
    def aggressive_cleanup(self, confirm: bool = True) -> bool:
        """Aggressively clean up all port conflicts"""
        
        if confirm:
            console.print("\n[yellow]This will forcefully kill ALL processes using required ports.[/yellow]")
            console.print("[yellow]This includes any running Neo4j, Redis, InfluxDB, Grafana instances.[/yellow]")
            
            response = input("\nContinue? (y/N): ").lower().strip()
            if response != 'y':
                console.print("[yellow]Cleanup cancelled.[/yellow]")
                return False
        
        console.print("\n[bold red]🔥 Starting Aggressive Port Cleanup[/bold red]")
        
        killed_any = False
        
        for port in self.required_ports:
            console.print(f"\n[cyan]Checking port {port} ({self.port_services[port]})...[/cyan]")
            processes = self.find_port_processes(port)
            
            if not processes:
                console.print(f"  [green]✓ Port {port} is free[/green]")
                continue
            
            console.print(f"  [red]Found {len(processes)} process(es) using port {port}[/red]")
            
            for proc in processes:
                console.print(f"    Killing {proc['command']} (PID: {proc['pid']})")
                
                # Try gentle kill first
                if self.kill_process(proc['pid'], force=False):
                    console.print(f"      [green]✓ Terminated PID {proc['pid']}[/green]")
                    killed_any = True
                else:
                    # Force kill
                    console.print(f"      [yellow]Trying force kill...[/yellow]")
                    if self.kill_process(proc['pid'], force=True):
                        console.print(f"      [green]✓ Force killed PID {proc['pid']}[/green]")
                        killed_any = True
                    else:
                        console.print(f"      [red]✗ Failed to kill PID {proc['pid']}[/red]")
        
        if killed_any:
            console.print("\n[yellow]Waiting 5 seconds for processes to fully terminate...[/yellow]")
            time.sleep(5)
        
        return True
    
    def docker_cleanup(self):
        """Clean up Docker resources"""
        console.print("\n[cyan]Cleaning up Docker resources...[/cyan]")
        
        # Stop all containers with our names
        containers = ["neo4j-genomics", "redis-genomics", "influxdb-genomics", 
                     "grafana-genomics", "jupyter-genomics"]
        
        for container in containers:
            subprocess.run(['docker', 'kill', container], capture_output=True)
            subprocess.run(['docker', 'rm', container], capture_output=True)
        
        # Clean up networks
        networks = ["genomics-network", "fractalpangenome_default"]
        for network in networks:
            subprocess.run(['docker', 'network', 'rm', network], capture_output=True)
        
        # Prune everything
        subprocess.run(['docker', 'container', 'prune', '-f'], capture_output=True)
        subprocess.run(['docker', 'network', 'prune', '-f'], capture_output=True)
        subprocess.run(['docker', 'volume', 'prune', '-f'], capture_output=True)
        
        console.print("  [green]✓ Docker cleanup completed[/green]")

def main():
    """Main cleanup function"""
    cleaner = AggressivePortCleaner()
    
    console.print("[bold blue]🧹 Aggressive Port Cleanup Tool[/bold blue]")
    console.print("This tool will find and kill ANY process using the required ports.\n")
    
    # Show current status
    cleaner.display_port_status()
    
    # Docker cleanup first
    cleaner.docker_cleanup()
    
    # Port cleanup
    success = cleaner.aggressive_cleanup()
    
    if success:
        console.print("\n[green]✅ Aggressive cleanup completed![/green]")
        
        # Show status after cleanup
        console.print("\n[cyan]Post-cleanup status:[/cyan]")
        cleaner.display_port_status()
        
        console.print("\n[bold cyan]Next steps:[/bold cyan]")
        console.print("1. Wait 10 seconds for everything to settle")
        console.print("2. Run: [bold]python genome_cli.py start --docker-compose[/bold]")
        
        # Optional wait
        response = input("\nWait 10 seconds now? (Y/n): ").lower().strip()
        if response != 'n':
            console.print("Waiting 10 seconds...")
            time.sleep(10)
            console.print("[green]Ready to start the system![/green]")

if __name__ == "__main__":
    main()
