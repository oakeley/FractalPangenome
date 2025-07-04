"""
Performance Benchmarking and Optimization Tools
===============================================
Comprehensive performance analysis and optimization utilities for the fractal pangenome database
"""

import asyncio
import time
import psutil
import statistics
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from typing import Dict, List, Tuple, Optional, Callable, Any
from dataclasses import dataclass, field
from pathlib import Path
import json
import logging
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import multiprocessing as mp
from contextlib import contextmanager
import resource
import tracemalloc
from datetime import datetime, timedelta
import warnings
warnings.filterwarnings('ignore')

# Import our core modules
from neo4j import GraphDatabase
import redis
from hilbert_pangenome_architecture import HilbertCurve, GenomicCoordinate
from neo4j_genome_importer import GenomeDatabaseBuilder, DockerNeo4jManager

@dataclass
class PerformanceMetrics:
    """Container for performance measurement results"""
    operation_name: str
    execution_time: float
    memory_usage: float
    cpu_usage: float
    disk_io: Optional[Dict[str, float]] = None
    network_io: Optional[Dict[str, float]] = None
    query_complexity: Optional[str] = None
    data_size: Optional[int] = None
    success: bool = True
    error_message: Optional[str] = None
    timestamp: datetime = field(default_factory=datetime.now)

class SystemMonitor:
    """Monitor system resources during operations"""
    
    def __init__(self, sampling_interval: float = 0.1):
        self.sampling_interval = sampling_interval
        self.measurements = []
        self.is_monitoring = False
        
    def start_monitoring(self):
        """Start background monitoring"""
        self.is_monitoring = True
        self.measurements = []
        asyncio.create_task(self._monitor_loop())
    
    def stop_monitoring(self) -> Dict[str, List[float]]:
        """Stop monitoring and return aggregated results"""
        self.is_monitoring = False
        return self._aggregate_measurements()
    
    async def _monitor_loop(self):
        """Background monitoring loop"""
        while self.is_monitoring:
            measurement = {
                'timestamp': time.time(),
                'cpu_percent': psutil.cpu_percent(),
                'memory_percent': psutil.virtual_memory().percent,
                'memory_used': psutil.virtual_memory().used,
                'disk_io': psutil.disk_io_counters()._asdict() if psutil.disk_io_counters() else {},
                'network_io': psutil.net_io_counters()._asdict() if psutil.net_io_counters() else {}
            }
            self.measurements.append(measurement)
            await asyncio.sleep(self.sampling_interval)
    
    def _aggregate_measurements(self) -> Dict[str, Any]:
        """Aggregate monitoring measurements"""
        if not self.measurements:
            return {}
        
        cpu_values = [m['cpu_percent'] for m in self.measurements]
        memory_values = [m['memory_percent'] for m in self.measurements]
        
        return {
            'duration': self.measurements[-1]['timestamp'] - self.measurements[0]['timestamp'],
            'cpu_mean': statistics.mean(cpu_values),
            'cpu_max': max(cpu_values),
            'memory_mean': statistics.mean(memory_values),
            'memory_max': max(memory_values),
            'samples': len(self.measurements)
        }

@contextmanager
def performance_timer(operation_name: str, monitor_system: bool = True):
    """Context manager for timing operations and measuring resources"""
    monitor = SystemMonitor() if monitor_system else None
    
    # Start monitoring
    if monitor:
        monitor.start_monitoring()
    
    # Memory tracking
    tracemalloc.start()
    start_memory = tracemalloc.get_traced_memory()[0]
    start_time = time.perf_counter()
    
    try:
        yield
        success = True
        error_message = None
    except Exception as e:
        success = False
        error_message = str(e)
        raise
    finally:
        # Stop timing and monitoring
        end_time = time.perf_counter()
        end_memory = tracemalloc.get_traced_memory()[0]
        tracemalloc.stop()
        
        system_metrics = monitor.stop_monitoring() if monitor else {}
        
        # Create performance metrics
        metrics = PerformanceMetrics(
            operation_name=operation_name,
            execution_time=end_time - start_time,
            memory_usage=(end_memory - start_memory) / 1024 / 1024,  # MB
            cpu_usage=system_metrics.get('cpu_mean', 0),
            success=success,
            error_message=error_message
        )
        
        # Log the results
        logging.info(f"Performance: {operation_name} - {metrics.execution_time:.3f}s")

class DatabaseBenchmark:
    """Comprehensive database performance benchmarking"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        self.results = []
        
    async def run_comprehensive_benchmark(self) -> pd.DataFrame:
        """Run a comprehensive benchmark suite"""
        logging.info("Starting comprehensive database benchmark...")
        
        benchmark_suites = [
            self._benchmark_basic_queries(),
            self._benchmark_hilbert_queries(),
            self._benchmark_path_finding(),
            self._benchmark_import_performance(),
            self._benchmark_concurrent_access(),
            self._benchmark_memory_usage(),
            self._benchmark_scalability()
        ]
        
        for benchmark_suite in benchmark_suites:
            suite_results = await benchmark_suite
            self.results.extend(suite_results)
        
        return pd.DataFrame([vars(r) for r in self.results])
    
    async def _benchmark_basic_queries(self) -> List[PerformanceMetrics]:
        """Benchmark basic database operations"""
        results = []
        
        # Node count query
        with performance_timer("basic_node_count") as timer:
            with self.driver.session() as session:
                result = session.run("MATCH (n:GenomeNode) RETURN count(n) as count")
                count = result.single()['count']
        
        results.append(PerformanceMetrics(
            operation_name="basic_node_count",
            execution_time=timer.execution_time if hasattr(timer, 'execution_time') else 0,
            memory_usage=0,
            cpu_usage=0,
            data_size=count
        ))
        
        # Individual lookup
        test_queries = [
            ("individual_nodes", "MATCH (n:GenomeNode {individual_id: 'HG002'}) RETURN count(n)"),
            ("chromosome_nodes", "MATCH (n:GenomeNode {chromosome: 1}) RETURN count(n)"),
            ("frequency_filter", "MATCH (n:GenomeNode) WHERE n.frequency > 0.5 RETURN count(n)"),
            ("position_range", "MATCH (n:GenomeNode) WHERE n.start_pos BETWEEN 1000000 AND 2000000 RETURN count(n)")
        ]
        
        for query_name, query in test_queries:
            start_time = time.perf_counter()
            try:
                with self.driver.session() as session:
                    result = session.run(query)
                    count = result.single()[0] if result.single() else 0
                end_time = time.perf_counter()
                
                results.append(PerformanceMetrics(
                    operation_name=query_name,
                    execution_time=end_time - start_time,
                    memory_usage=0,
                    cpu_usage=0,
                    data_size=count,
                    query_complexity="simple"
                ))
            except Exception as e:
                results.append(PerformanceMetrics(
                    operation_name=query_name,
                    execution_time=0,
                    memory_usage=0,
                    cpu_usage=0,
                    success=False,
                    error_message=str(e)
                ))
        
        return results
    
    async def _benchmark_hilbert_queries(self) -> List[PerformanceMetrics]:
        """Benchmark Hilbert curve-based spatial queries"""
        results = []
        hilbert_curve = HilbertCurve(dimensions=4, order=16)
        
        # Test different region sizes
        test_regions = [
            (1, 1000000, 1010000),    # 10kb
            (1, 1000000, 1100000),    # 100kb
            (1, 1000000, 2000000),    # 1Mb
            (1, 1000000, 10000000),   # 10Mb
        ]
        
        for chromosome, start, end in test_regions:
            region_size = end - start
            
            # Calculate Hilbert range
            start_coord = GenomicCoordinate(chromosome, start, 0, 0)
            end_coord = GenomicCoordinate(chromosome, end, 0, 0)
            start_hilbert = hilbert_curve.encode(start_coord.to_hilbert_coords())
            end_hilbert = hilbert_curve.encode(end_coord.to_hilbert_coords())
            
            # Benchmark Hilbert-based query
            start_time = time.perf_counter()
            try:
                with self.driver.session() as session:
                    result = session.run("""
                        MATCH (n:GenomeNode)
                        WHERE n.hilbert_index >= $start_hilbert 
                          AND n.hilbert_index <= $end_hilbert
                          AND n.chromosome = $chromosome
                        RETURN count(n) as count
                    """, start_hilbert=start_hilbert, end_hilbert=end_hilbert, 
                        chromosome=chromosome)
                    count = result.single()['count']
                end_time = time.perf_counter()
                
                results.append(PerformanceMetrics(
                    operation_name=f"hilbert_query_{region_size//1000}kb",
                    execution_time=end_time - start_time,
                    memory_usage=0,
                    cpu_usage=0,
                    data_size=count,
                    query_complexity="spatial"
                ))
            except Exception as e:
                results.append(PerformanceMetrics(
                    operation_name=f"hilbert_query_{region_size//1000}kb",
                    execution_time=0,
                    memory_usage=0,
                    cpu_usage=0,
                    success=False,
                    error_message=str(e)
                ))
        
        return results
    
    async def _benchmark_path_finding(self) -> List[PerformanceMetrics]:
        """Benchmark graph pathfinding operations"""
        results = []
        
        path_queries = [
            ("short_path", "MATCH path = shortestPath((a:GenomeNode)-[:CONNECTS*1..5]->(b:GenomeNode)) WHERE a.chromosome = 1 RETURN length(path) LIMIT 10"),
            ("medium_path", "MATCH path = shortestPath((a:GenomeNode)-[:CONNECTS*1..10]->(b:GenomeNode)) WHERE a.chromosome = 1 RETURN length(path) LIMIT 10"),
            ("complex_path", "MATCH path = (a:GenomeNode)-[:CONNECTS*1..15]->(b:GenomeNode) WHERE a.chromosome = 1 AND b.chromosome = 1 RETURN length(path) LIMIT 5")
        ]
        
        for query_name, query in path_queries:
            start_time = time.perf_counter()
            try:
                with self.driver.session() as session:
                    result = session.run(query)
                    paths = list(result)
                end_time = time.perf_counter()
                
                results.append(PerformanceMetrics(
                    operation_name=query_name,
                    execution_time=end_time - start_time,
                    memory_usage=0,
                    cpu_usage=0,
                    data_size=len(paths),
                    query_complexity="graph_traversal"
                ))
            except Exception as e:
                results.append(PerformanceMetrics(
                    operation_name=query_name,
                    execution_time=0,
                    memory_usage=0,
                    cpu_usage=0,
                    success=False,
                    error_message=str(e),
                    query_complexity="graph_traversal"
                ))
        
        return results
    
    async def _benchmark_import_performance(self) -> List[PerformanceMetrics]:
        """Benchmark data import performance"""
        results = []
        
        # Simulate importing different sizes of genomic data
        import_sizes = [100, 500, 1000, 5000]  # Number of nodes
        
        for size in import_sizes:
            start_time = time.perf_counter()
            
            try:
                # Create simulated nodes
                nodes_data = []
                for i in range(size):
                    nodes_data.append({
                        'id': f'benchmark_node_{i}',
                        'sequence': 'ATCG' * 250,  # 1kb sequence
                        'chromosome': 1,
                        'start_pos': i * 1000,
                        'end_pos': (i + 1) * 1000,
                        'scale_level': 0,
                        'frequency': 0.5,
                        'hilbert_index': i
                    })
                
                # Batch import
                with self.driver.session() as session:
                    session.run("""
                        UNWIND $nodes as node
                        CREATE (n:BenchmarkNode {
                            id: node.id,
                            sequence: node.sequence,
                            chromosome: node.chromosome,
                            start_pos: node.start_pos,
                            end_pos: node.end_pos,
                            scale_level: node.scale_level,
                            frequency: node.frequency,
                            hilbert_index: node.hilbert_index
                        })
                    """, nodes=nodes_data)
                
                end_time = time.perf_counter()
                
                results.append(PerformanceMetrics(
                    operation_name=f"import_{size}_nodes",
                    execution_time=end_time - start_time,
                    memory_usage=0,
                    cpu_usage=0,
                    data_size=size,
                    query_complexity="batch_write"
                ))
                
                # Clean up benchmark nodes
                with self.driver.session() as session:
                    session.run("MATCH (n:BenchmarkNode) DELETE n")
                    
            except Exception as e:
                results.append(PerformanceMetrics(
                    operation_name=f"import_{size}_nodes",
                    execution_time=0,
                    memory_usage=0,
                    cpu_usage=0,
                    success=False,
                    error_message=str(e),
                    query_complexity="batch_write"
                ))
        
        return results
    
    async def _benchmark_concurrent_access(self) -> List[PerformanceMetrics]:
        """Benchmark concurrent database access"""
        results = []
        
        async def run_concurrent_queries(num_threads: int):
            """Run queries concurrently"""
            query = "MATCH (n:GenomeNode) WHERE n.chromosome = 1 RETURN count(n)"
            
            async def single_query():
                with self.driver.session() as session:
                    return session.run(query).single()[0]
            
            start_time = time.perf_counter()
            
            # Run concurrent queries
            tasks = [single_query() for _ in range(num_threads)]
            await asyncio.gather(*tasks)
            
            end_time = time.perf_counter()
            
            return end_time - start_time
        
        # Test different concurrency levels
        concurrency_levels = [1, 5, 10, 20]
        
        for num_threads in concurrency_levels:
            try:
                execution_time = await run_concurrent_queries(num_threads)
                
                results.append(PerformanceMetrics(
                    operation_name=f"concurrent_{num_threads}_threads",
                    execution_time=execution_time,
                    memory_usage=0,
                    cpu_usage=0,
                    data_size=num_threads,
                    query_complexity="concurrent"
                ))
            except Exception as e:
                results.append(PerformanceMetrics(
                    operation_name=f"concurrent_{num_threads}_threads",
                    execution_time=0,
                    memory_usage=0,
                    cpu_usage=0,
                    success=False,
                    error_message=str(e),
                    query_complexity="concurrent"
                ))
        
        return results
    
    async def _benchmark_memory_usage(self) -> List[PerformanceMetrics]:
        """Benchmark memory usage patterns"""
        results = []
        
        # Large result set query
        tracemalloc.start()
        start_memory = tracemalloc.get_traced_memory()[0]
        start_time = time.perf_counter()
        
        try:
            with self.driver.session() as session:
                result = session.run("""
                    MATCH (n:GenomeNode)
                    WHERE n.chromosome = 1
                    RETURN n.id, n.sequence, n.start_pos, n.end_pos
                    LIMIT 10000
                """)
                data = list(result)
            
            end_time = time.perf_counter()
            end_memory = tracemalloc.get_traced_memory()[0]
            
            results.append(PerformanceMetrics(
                operation_name="large_result_set",
                execution_time=end_time - start_time,
                memory_usage=(end_memory - start_memory) / 1024 / 1024,  # MB
                cpu_usage=0,
                data_size=len(data),
                query_complexity="memory_intensive"
            ))
        except Exception as e:
            results.append(PerformanceMetrics(
                operation_name="large_result_set",
                execution_time=0,
                memory_usage=0,
                cpu_usage=0,
                success=False,
                error_message=str(e),
                query_complexity="memory_intensive"
            ))
        finally:
            tracemalloc.stop()
        
        return results
    
    async def _benchmark_scalability(self) -> List[PerformanceMetrics]:
        """Benchmark scalability with different data sizes"""
        results = []
        
        # Test queries with different LIMIT values to simulate different data sizes
        data_sizes = [100, 1000, 10000, 50000]
        
        for limit in data_sizes:
            start_time = time.perf_counter()
            try:
                with self.driver.session() as session:
                    result = session.run(f"""
                        MATCH (n:GenomeNode)
                        RETURN n.id, n.chromosome, n.start_pos
                        LIMIT {limit}
                    """)
                    data = list(result)
                
                end_time = time.perf_counter()
                
                results.append(PerformanceMetrics(
                    operation_name=f"scalability_{limit}_records",
                    execution_time=end_time - start_time,
                    memory_usage=0,
                    cpu_usage=0,
                    data_size=len(data),
                    query_complexity="scalability"
                ))
            except Exception as e:
                results.append(PerformanceMetrics(
                    operation_name=f"scalability_{limit}_records",
                    execution_time=0,
                    memory_usage=0,
                    cpu_usage=0,
                    success=False,
                    error_message=str(e),
                    query_complexity="scalability"
                ))
        
        return results

class PerformanceOptimizer:
    """Automated performance optimization recommendations"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        
    async def analyze_and_optimize(self) -> Dict[str, Any]:
        """Run comprehensive performance analysis and provide optimization recommendations"""
        analysis = {
            'database_stats': await self._analyze_database_structure(),
            'index_analysis': await self._analyze_indices(),
            'query_patterns': await self._analyze_query_patterns(),
            'memory_usage': await self._analyze_memory_usage(),
            'recommendations': []
        }
        
        # Generate recommendations based on analysis
        analysis['recommendations'] = self._generate_recommendations(analysis)
        
        return analysis
    
    async def _analyze_database_structure(self) -> Dict[str, Any]:
        """Analyze database structure and size"""
        with self.driver.session() as session:
            # Get node counts by type
            node_stats = session.run("""
                MATCH (n)
                RETURN labels(n) as labels, count(n) as count
                ORDER BY count DESC
            """).data()
            
            # Get relationship counts
            rel_stats = session.run("""
                MATCH ()-[r]->()
                RETURN type(r) as type, count(r) as count
                ORDER BY count DESC
            """).data()
            
            # Get property statistics
            prop_stats = session.run("""
                MATCH (n:GenomeNode)
                RETURN 
                    count(n) as total_nodes,
                    avg(size(n.sequence)) as avg_sequence_length,
                    min(n.start_pos) as min_position,
                    max(n.end_pos) as max_position
            """).single()
            
            return {
                'node_statistics': node_stats,
                'relationship_statistics': rel_stats,
                'property_statistics': dict(prop_stats) if prop_stats else {}
            }
    
    async def _analyze_indices(self) -> Dict[str, Any]:
        """Analyze database indices and their usage"""
        with self.driver.session() as session:
            # Get index information
            indices = session.run("SHOW INDEXES").data()
            
            # Analyze index usage (simplified)
            index_analysis = []
            for index in indices:
                index_analysis.append({
                    'name': index.get('name', 'unknown'),
                    'type': index.get('type', 'unknown'),
                    'state': index.get('state', 'unknown'),
                    'properties': index.get('labelsOrTypes', [])
                })
            
            return {
                'existing_indices': index_analysis,
                'total_indices': len(indices)
            }
    
    async def _analyze_query_patterns(self) -> Dict[str, Any]:
        """Analyze common query patterns for optimization opportunities"""
        # This would analyze query logs in a production system
        # For now, we'll simulate common patterns
        
        common_patterns = [
            {
                'pattern': 'chromosome_range_queries',
                'frequency': 0.4,
                'optimization_potential': 'high',
                'description': 'Queries filtering by chromosome and position range'
            },
            {
                'pattern': 'individual_lookups',
                'frequency': 0.3,
                'optimization_potential': 'medium',
                'description': 'Queries filtering by individual_id'
            },
            {
                'pattern': 'frequency_filters',
                'frequency': 0.2,
                'optimization_potential': 'medium',
                'description': 'Queries filtering by variant frequency'
            },
            {
                'pattern': 'path_traversals',
                'frequency': 0.1,
                'optimization_potential': 'high',
                'description': 'Graph traversal queries for pathfinding'
            }
        ]
        
        return {
            'common_patterns': common_patterns,
            'analysis_date': datetime.now().isoformat()
        }
    
    async def _analyze_memory_usage(self) -> Dict[str, Any]:
        """Analyze memory usage patterns"""
        # Get system memory info
        memory = psutil.virtual_memory()
        
        # Check Neo4j memory configuration (would read from config in production)
        neo4j_memory = {
            'heap_size': '8G',  # Would read from actual config
            'page_cache': '4G',
            'available_memory': f"{memory.available // (1024**3)}G"
        }
        
        return {
            'system_memory': {
                'total': memory.total // (1024**3),
                'available': memory.available // (1024**3),
                'percent_used': memory.percent
            },
            'neo4j_memory': neo4j_memory
        }
    
    def _generate_recommendations(self, analysis: Dict[str, Any]) -> List[Dict[str, str]]:
        """Generate optimization recommendations based on analysis"""
        recommendations = []
        
        # Check node count for indexing recommendations
        total_nodes = 0
        for node_stat in analysis['database_stats']['node_statistics']:
            if 'GenomeNode' in node_stat['labels']:
                total_nodes = node_stat['count']
                break
        
        if total_nodes > 100000:
            recommendations.append({
                'type': 'indexing',
                'priority': 'high',
                'recommendation': 'Create composite index on (chromosome, start_pos, end_pos)',
                'impact': 'Significantly improve range queries',
                'implementation': 'CREATE INDEX composite_genomic_position FOR (n:GenomeNode) ON (n.chromosome, n.start_pos, n.end_pos)'
            })
        
        # Memory recommendations
        memory_analysis = analysis['memory_usage']
        if memory_analysis['system_memory']['percent_used'] > 80:
            recommendations.append({
                'type': 'memory',
                'priority': 'high',
                'recommendation': 'Increase system memory or optimize Neo4j memory settings',
                'impact': 'Reduce memory pressure and improve query performance',
                'implementation': 'Consider increasing heap size and page cache allocation'
            })
        
        # Query pattern recommendations
        for pattern in analysis['query_patterns']['common_patterns']:
            if pattern['optimization_potential'] == 'high':
                if pattern['pattern'] == 'chromosome_range_queries':
                    recommendations.append({
                        'type': 'indexing',
                        'priority': 'high',
                        'recommendation': 'Optimize chromosome range queries with Hilbert curve indexing',
                        'impact': 'Reduce query time from O(n) to O(log n)',
                        'implementation': 'Ensure Hilbert index is being used effectively'
                    })
                elif pattern['pattern'] == 'path_traversals':
                    recommendations.append({
                        'type': 'configuration',
                        'priority': 'medium',
                        'recommendation': 'Tune graph algorithms memory allocation',
                        'impact': 'Improve pathfinding query performance',
                        'implementation': 'Adjust dbms.memory.transaction.* settings'
                    })
        
        # Index recommendations
        existing_indices = analysis['index_analysis']['existing_indices']
        essential_indices = [
            'hilbert_index',
            'genomic_position',
            'individual_id_index'
        ]
        
        existing_index_names = [idx['name'] for idx in existing_indices]
        for essential_index in essential_indices:
            if essential_index not in existing_index_names:
                recommendations.append({
                    'type': 'indexing',
                    'priority': 'medium',
                    'recommendation': f'Create missing essential index: {essential_index}',
                    'impact': 'Improve query performance for common access patterns',
                    'implementation': f'Create appropriate index for {essential_index}'
                })
        
        return recommendations

class PerformanceReporter:
    """Generate comprehensive performance reports"""
    
    def __init__(self):
        self.report_data = {}
        
    def generate_report(self, benchmark_results: pd.DataFrame, 
                       optimization_analysis: Dict[str, Any]) -> Dict[str, Any]:
        """Generate a comprehensive performance report"""
        
        report = {
            'summary': self._generate_summary(benchmark_results),
            'detailed_results': self._process_benchmark_results(benchmark_results),
            'optimization_analysis': optimization_analysis,
            'performance_trends': self._analyze_trends(benchmark_results),
            'recommendations': optimization_analysis.get('recommendations', []),
            'report_metadata': {
                'generated_at': datetime.now().isoformat(),
                'total_tests': len(benchmark_results),
                'success_rate': (benchmark_results['success'].sum() / len(benchmark_results)) * 100
            }
        }
        
        return report
    
    def _generate_summary(self, results: pd.DataFrame) -> Dict[str, Any]:
        """Generate summary statistics"""
        if results.empty:
            return {'error': 'No benchmark results available'}
        
        successful_results = results[results['success'] == True]
        
        return {
            'total_operations': len(results),
            'successful_operations': len(successful_results),
            'success_rate': f"{(len(successful_results) / len(results)) * 100:.1f}%",
            'average_execution_time': f"{successful_results['execution_time'].mean():.3f}s",
            'fastest_operation': {
                'name': successful_results.loc[successful_results['execution_time'].idxmin(), 'operation_name'],
                'time': f"{successful_results['execution_time'].min():.3f}s"
            },
            'slowest_operation': {
                'name': successful_results.loc[successful_results['execution_time'].idxmax(), 'operation_name'],
                'time': f"{successful_results['execution_time'].max():.3f}s"
            }
        }
    
    def _process_benchmark_results(self, results: pd.DataFrame) -> Dict[str, Any]:
        """Process benchmark results by category"""
        if results.empty:
            return {}
        
        categories = {}
        
        # Group by query complexity
        for complexity in results['query_complexity'].dropna().unique():
            category_results = results[results['query_complexity'] == complexity]
            categories[complexity] = {
                'count': len(category_results),
                'avg_time': category_results['execution_time'].mean(),
                'operations': category_results['operation_name'].tolist()
            }
        
        return categories
    
    def _analyze_trends(self, results: pd.DataFrame) -> Dict[str, Any]:
        """Analyze performance trends"""
        if results.empty:
            return {}
        
        # Analyze scalability trends
        scalability_results = results[results['operation_name'].str.contains('scalability', na=False)]
        
        trends = {}
        
        if not scalability_results.empty:
            # Extract data sizes and execution times
            data_sizes = scalability_results['data_size'].values
            exec_times = scalability_results['execution_time'].values
            
            if len(data_sizes) > 1 and len(exec_times) > 1:
                # Calculate correlation between data size and execution time
                correlation = np.corrcoef(data_sizes, exec_times)[0, 1] if len(data_sizes) > 1 else 0
                
                trends['scalability'] = {
                    'correlation_coefficient': correlation,
                    'interpretation': 'Linear' if correlation > 0.8 else 'Sub-linear' if correlation > 0.5 else 'Optimal',
                    'data_points': len(scalability_results)
                }
        
        return trends
    
    def save_report(self, report: Dict[str, Any], filepath: str):
        """Save report to file"""
        with open(filepath, 'w') as f:
            json.dump(report, f, indent=2, default=str)
    
    def generate_visualizations(self, benchmark_results: pd.DataFrame, 
                              output_dir: str = "./performance_plots"):
        """Generate performance visualization plots"""
        Path(output_dir).mkdir(exist_ok=True)
        
        if benchmark_results.empty:
            print("No data available for visualization")
            return
        
        # Set up the plotting style
        plt.style.use('seaborn-v0_8')
        fig, axes = plt.subplots(2, 2, figsize=(16, 12))
        
        # 1. Execution time by operation
        successful_results = benchmark_results[benchmark_results['success'] == True]
        if not successful_results.empty:
            operations = successful_results['operation_name']
            times = successful_results['execution_time']
            
            axes[0, 0].barh(range(len(operations)), times)
            axes[0, 0].set_yticks(range(len(operations)))
            axes[0, 0].set_yticklabels(operations, fontsize=8)
            axes[0, 0].set_xlabel('Execution Time (seconds)')
            axes[0, 0].set_title('Execution Time by Operation')
        
        # 2. Success rate by query complexity
        if 'query_complexity' in benchmark_results.columns:
            complexity_success = benchmark_results.groupby('query_complexity')['success'].mean()
            axes[0, 1].bar(complexity_success.index, complexity_success.values)
            axes[0, 1].set_ylabel('Success Rate')
            axes[0, 1].set_title('Success Rate by Query Complexity')
            axes[0, 1].tick_params(axis='x', rotation=45)
        
        # 3. Scalability analysis
        scalability_data = benchmark_results[
            benchmark_results['operation_name'].str.contains('scalability', na=False)
        ]
        if not scalability_data.empty:
            axes[1, 0].scatter(scalability_data['data_size'], scalability_data['execution_time'])
            axes[1, 0].set_xlabel('Data Size')
            axes[1, 0].set_ylabel('Execution Time (seconds)')
            axes[1, 0].set_title('Scalability Analysis')
            
            # Add trend line
            if len(scalability_data) > 1:
                z = np.polyfit(scalability_data['data_size'], scalability_data['execution_time'], 1)
                p = np.poly1d(z)
                axes[1, 0].plot(scalability_data['data_size'], p(scalability_data['data_size']), "r--", alpha=0.8)
        
        # 4. Memory usage distribution
        if 'memory_usage' in benchmark_results.columns:
            memory_data = benchmark_results[benchmark_results['memory_usage'] > 0]['memory_usage']
            if not memory_data.empty:
                axes[1, 1].hist(memory_data, bins=20, alpha=0.7)
                axes[1, 1].set_xlabel('Memory Usage (MB)')
                axes[1, 1].set_ylabel('Frequency')
                axes[1, 1].set_title('Memory Usage Distribution')
        
        plt.tight_layout()
        plt.savefig(f"{output_dir}/performance_analysis.png", dpi=300, bbox_inches='tight')
        plt.close()
        
        print(f"Performance visualizations saved to {output_dir}")

# Main performance testing orchestrator
async def run_performance_suite():
    """Run the complete performance testing suite"""
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    
    print("🚀 Starting Fractal Pangenome Performance Suite")
    print("=" * 60)
    
    try:
        # Initialize database builder
        db_builder = GenomeDatabaseBuilder()
        await db_builder.initialize_database()
        
        # Run benchmarks
        print("📊 Running comprehensive benchmarks...")
        benchmark = DatabaseBenchmark(db_builder)
        results_df = await benchmark.run_comprehensive_benchmark()
        
        # Run optimization analysis
        print("🔧 Analyzing optimization opportunities...")
        optimizer = PerformanceOptimizer(db_builder)
        optimization_analysis = await optimizer.analyze_and_optimize()
        
        # Generate report
        print("📋 Generating performance report...")
        reporter = PerformanceReporter()
        report = reporter.generate_report(results_df, optimization_analysis)
        
        # Save results
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        results_df.to_csv(f"performance_results_{timestamp}.csv", index=False)
        reporter.save_report(report, f"performance_report_{timestamp}.json")
        
        # Generate visualizations
        print("📈 Creating performance visualizations...")
        reporter.generate_visualizations(results_df)
        
        # Print summary
        print("\n" + "=" * 60)
        print("🎯 Performance Suite Complete!")
        print("=" * 60)
        
        summary = report['summary']
        print(f"Total Operations: {summary['total_operations']}")
        print(f"Success Rate: {summary['success_rate']}")
        print(f"Average Execution Time: {summary['average_execution_time']}")
        print(f"Fastest Operation: {summary['fastest_operation']['name']} ({summary['fastest_operation']['time']})")
        print(f"Slowest Operation: {summary['slowest_operation']['name']} ({summary['slowest_operation']['time']})")
        
        print(f"\n📊 Recommendations ({len(report['recommendations'])} found):")
        for i, rec in enumerate(report['recommendations'][:5], 1):  # Show top 5
            print(f"{i}. [{rec['priority'].upper()}] {rec['recommendation']}")
        
        print(f"\n📁 Results saved:")
        print(f"   • CSV: performance_results_{timestamp}.csv")
        print(f"   • Report: performance_report_{timestamp}.json")
        print(f"   • Plots: ./performance_plots/")
        
    except Exception as e:
        print(f"❌ Error running performance suite: {e}")
        logging.error(f"Performance suite error: {e}")
    finally:
        print("\n🏁 Performance suite finished!")

if __name__ == "__main__":
    asyncio.run(run_performance_suite())
