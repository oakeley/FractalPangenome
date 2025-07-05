"""
Neo4j Genome Database Importer Module
=====================================
Docker-based Neo4j setup with FASTA genome import and short read mapping capabilities
"""

import os
import subprocess
import docker
import asyncio
import aiofiles
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set, AsyncGenerator
from dataclasses import dataclass, field
from neo4j import GraphDatabase
import hashlib
import gzip
from collections import defaultdict, deque
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import requests
import tempfile
import logging
import json
import aiohttp
import yaml
import time
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor
import multiprocessing as mp
from itertools import islice
import mmap

# Import from previous module
from hilbert_pangenome_architecture import HilbertCurve, GenomicCoordinate, PangenomeNode, PangenomeEdge

@dataclass
class GenomeAssembly:
    """Represents a complete genome assembly"""
    individual_id: str
    haplotype: int  # 0 or 1 for diploid
    assembly_source: str  # e.g., "T2T", "GRCh38", "personal_assembly"
    fasta_path: str
    metadata: Dict[str, any] = field(default_factory=dict)
    quality_metrics: Dict[str, float] = field(default_factory=dict)

@dataclass
class ShortRead:
    """Represents a short sequencing read"""
    read_id: str
    sequence: str
    quality: str
    mate_sequence: Optional[str] = None
    mate_quality: Optional[str] = None
    
class DockerNeo4jManager:
    """Manages Neo4j Docker container for genomic data"""
    
    def __init__(self, container_name: str = "neo4j-genomics", 
                 neo4j_version: str = "5.15-community",
                 data_dir: str = "./neo4j_data"):
        self.container_name = container_name
        self.neo4j_version = neo4j_version
        self.data_dir = Path(data_dir).absolute()
        self.docker_client = docker.from_env()
        self.neo4j_uri = "bolt://localhost:7687"
        self.neo4j_user = "neo4j"
        self.neo4j_password = os.getenv('NEO4J_PASSWORD', 'genomics123')
        self.container = None
        
    def setup_neo4j_container(self) -> bool:
        """Setup and start Neo4j Docker container"""
        try:
            # Create data directories
            self.data_dir.mkdir(parents=True, exist_ok=True)
            (self.data_dir / "data").mkdir(exist_ok=True)
            (self.data_dir / "logs").mkdir(exist_ok=True)
            (self.data_dir / "import").mkdir(exist_ok=True)
            (self.data_dir / "plugins").mkdir(exist_ok=True)
            
            # Stop existing container if running
            self.stop_container()
            
            # Configure Neo4j settings
            neo4j_conf = self._generate_neo4j_config()
            conf_path = self.data_dir / "conf"
            conf_path.mkdir(exist_ok=True)
            
            with open(conf_path / "neo4j.conf", "w") as f:
                f.write(neo4j_conf)
            
            # Run Neo4j container
            self.container = self.docker_client.containers.run(
                f"neo4j:{self.neo4j_version}",
                name=self.container_name,
                ports={'7474/tcp': 7474, '7687/tcp': 7687},
                volumes={
                    str(self.data_dir / "data"): {'bind': '/data', 'mode': 'rw'},
                    str(self.data_dir / "logs"): {'bind': '/logs', 'mode': 'rw'},
                    str(self.data_dir / "import"): {'bind': '/var/lib/neo4j/import', 'mode': 'rw'},
                    str(self.data_dir / "plugins"): {'bind': '/plugins', 'mode': 'rw'},
                    str(self.data_dir / "conf"): {'bind': '/conf', 'mode': 'rw'}
                },
                environment={
                    'NEO4J_AUTH': f'{self.neo4j_user}/{self.neo4j_password}',
                    'NEO4J_PLUGINS': '["graph-data-science", "apoc"]',
                    'NEO4J_dbms_memory_heap_initial__size': '2G',
                    'NEO4J_dbms_memory_heap_max__size': '4G',
                    'NEO4J_dbms_memory_pagecache_size': '2G'
                },
                detach=True,
                remove=True
            )
            
            # Wait for Neo4j to be ready
            self._wait_for_neo4j()
            logging.info(f"Neo4j container '{self.container_name}' started successfully")
            return True
            
        except Exception as e:
            logging.error(f"Failed to setup Neo4j container: {e}")
            return False
    
    def _generate_neo4j_config(self) -> str:
        """Generate optimized Neo4j configuration for genomics"""
        return """
# Genomics-optimized Neo4j configuration
dbms.default_listen_address=0.0.0.0
dbms.connector.bolt.listen_address=:7687
dbms.connector.http.listen_address=:7474

# Memory settings
dbms.memory.heap.initial_size=2G
dbms.memory.heap.max_size=4G
dbms.memory.pagecache.size=2G

# Transaction settings for large imports
dbms.transaction.timeout=300s
dbms.transaction.concurrent.maximum=1000

# Query settings
cypher.default_temporal_reuse=true
cypher.forbid_exhaustive_shortestpath=false

# Security
dbms.security.procedures.unrestricted=gds.*,apoc.*

# Performance
dbms.checkpoint.interval.time=300s
dbms.checkpoint.interval.tx=100000
"""
    
    def _wait_for_neo4j(self, max_wait: int = 120):
        """Wait for Neo4j to be ready"""
        start_time = time.time()
        while time.time() - start_time < max_wait:
            try:
                driver = GraphDatabase.driver(
                    self.neo4j_uri, 
                    auth=(self.neo4j_user, self.neo4j_password)
                )
                with driver.session() as session:
                    session.run("RETURN 1")
                driver.close()
                logging.info("Neo4j is ready")
                return True
            except Exception:
                time.sleep(2)
        
        raise TimeoutError("Neo4j failed to start within timeout period")
    
    def stop_container(self):
        """Stop the Neo4j container"""
        try:
            existing = self.docker_client.containers.get(self.container_name)
            existing.stop()
            existing.wait()
            logging.info(f"Stopped existing container '{self.container_name}'")
        except docker.errors.NotFound:
            pass
    
    def get_connection(self):
        """Get Neo4j driver connection"""
        return GraphDatabase.driver(
            self.neo4j_uri,
            auth=(self.neo4j_user, self.neo4j_password)
        )

class KmerIndex:
    """Efficient k-mer indexing for short read mapping"""
    
    def __init__(self, k: int = 31):
        self.k = k
        self.kmer_to_nodes = defaultdict(set)  # k-mer -> set of node_ids
        self.node_kmers = defaultdict(set)     # node_id -> set of k-mers
        self.kmer_positions = defaultdict(list)  # k-mer -> [(node_id, position)]
        
    def add_sequence(self, node_id: str, sequence: str, start_pos: int = 0):
        """Add sequence to k-mer index"""
        sequence = sequence.upper()
        for i in range(len(sequence) - self.k + 1):
            kmer = sequence[i:i + self.k]
            if 'N' not in kmer:  # Skip ambiguous k-mers
                self.kmer_to_nodes[kmer].add(node_id)
                self.node_kmers[node_id].add(kmer)
                self.kmer_positions[kmer].append((node_id, start_pos + i))
    
    def find_matching_nodes(self, read_sequence: str, min_kmers: int = 3) -> Dict[str, int]:
        """Find nodes that match k-mers from read sequence"""
        read_sequence = read_sequence.upper()
        node_matches = defaultdict(int)
        
        for i in range(len(read_sequence) - self.k + 1):
            kmer = read_sequence[i:i + self.k]
            if kmer in self.kmer_to_nodes:
                for node_id in self.kmer_to_nodes[kmer]:
                    node_matches[node_id] += 1
        
        # Filter nodes with sufficient k-mer matches
        return {node_id: count for node_id, count in node_matches.items() 
                if count >= min_kmers}
    
    def get_kmer_positions(self, kmer: str) -> List[Tuple[str, int]]:
        """Get all positions where k-mer occurs"""
        return self.kmer_positions.get(kmer, [])

class FastaGenomeProcessor:
    """Process FASTA files and build pangenome graph"""
    
    def __init__(self, hilbert_curve: HilbertCurve, segment_size: int = 10000):
        self.hilbert_curve = hilbert_curve
        self.segment_size = segment_size
        self.kmer_index = KmerIndex(k=31)
        self.reference_nodes = {}  # For linking assemblies
        
    async def process_fasta_file(self, fasta_path: str, individual_id: str, 
                               haplotype: int, assembly_source: str) -> List[PangenomeNode]:
        """Process a FASTA file and create pangenome nodes"""
        nodes = []
        
        async with aiofiles.open(fasta_path, 'r') as f:
            content = await f.read()
        
        # Parse FASTA sequences
        sequences = {}
        current_seq = ""
        current_header = ""
        
        for line in content.split('\n'):
            if line.startswith('>'):
                if current_header and current_seq:
                    sequences[current_header] = current_seq
                current_header = line[1:].strip()
                current_seq = ""
            else:
                current_seq += line.strip()
        
        # Add final sequence
        if current_header and current_seq:
            sequences[current_header] = current_seq
        
        # Process each chromosome/contig
        for header, sequence in sequences.items():
            chromosome = self._parse_chromosome(header)
            chromosome_nodes = await self._segment_sequence(
                sequence, chromosome, individual_id, haplotype, assembly_source
            )
            nodes.extend(chromosome_nodes)
            
            # Add to k-mer index
            for node in chromosome_nodes:
                self.kmer_index.add_sequence(node.node_id, node.sequence, node.start_pos)
        
        return nodes
    
    def _parse_chromosome(self, header: str) -> int:
        """Parse chromosome number from FASTA header"""
        # Handle different header formats
        header_lower = header.lower()
        
        # Standard formats: chr1, chromosome_1, etc.
        for i in range(1, 25):  # 1-22, X(23), Y(24)
            if f'chr{i}' in header_lower or f'chromosome_{i}' in header_lower:
                return i
        
        # Handle X and Y chromosomes
        if 'chrx' in header_lower or 'chromosome_x' in header_lower:
            return 23
        if 'chry' in header_lower or 'chromosome_y' in header_lower:
            return 24
        
        # Default to mitochondrial (25) for unrecognized
        return 25
    
    async def _segment_sequence(self, sequence: str, chromosome: int, 
                              individual_id: str, haplotype: int, 
                              assembly_source: str) -> List[PangenomeNode]:
        """Segment long sequence into manageable nodes"""
        nodes = []
        sequence_length = len(sequence)
        
        # Create multiple scale levels
        for scale_level in range(4):
            if scale_level == 0:  # Nucleotide level
                step_size = self.segment_size
            elif scale_level == 1:  # Gene level
                step_size = self.segment_size * 5
            elif scale_level == 2:  # Large segment level
                step_size = self.segment_size * 25
            else:  # Chromosome level
                step_size = sequence_length
            
            for start in range(0, sequence_length, step_size):
                end = min(start + step_size, sequence_length)
                
                # Create node
                node_id = f"{individual_id}_h{haplotype}_{assembly_source}_chr{chromosome}_{scale_level}_{start}_{end}"
                
                # Calculate Hilbert index
                coord = GenomicCoordinate(chromosome, start, scale_level, 0, haplotype)
                hilbert_index = self.hilbert_curve.encode(coord.to_hilbert_coords())
                
                node = PangenomeNode(
                    node_id=node_id,
                    sequence=sequence[start:end],
                    chromosome=chromosome,
                    start_pos=start,
                    end_pos=end,
                    scale_level=scale_level,
                    frequency=0.5,  # Will be updated based on population data
                    annotations={
                        'individual_id': individual_id,
                        'haplotype': haplotype,
                        'assembly_source': assembly_source,
                        'segment_length': end - start
                    },
                    hilbert_index=hilbert_index
                )
                
                nodes.append(node)
        
        return nodes

class GenomeImporter:
    """Main class for importing genomes into Neo4j"""
    
    def __init__(self, neo4j_manager: DockerNeo4jManager):
        self.neo4j_manager = neo4j_manager
        self.driver = neo4j_manager.get_connection()
        self.hilbert_curve = HilbertCurve(dimensions=4, order=20)
        self.fasta_processor = FastaGenomeProcessor(self.hilbert_curve)
        self.imported_genomes = set()
        
        # Initialize database schema
        self._initialize_database_schema()
    
    def _initialize_database_schema(self):
        """Create database constraints and indices"""
        with self.driver.session() as session:
            # Constraints
            session.run("""
                CREATE CONSTRAINT genome_node_id IF NOT EXISTS
                FOR (n:GenomeNode) REQUIRE n.id IS UNIQUE
            """)
            
            session.run("""
                CREATE CONSTRAINT individual_id IF NOT EXISTS
                FOR (i:Individual) REQUIRE i.id IS UNIQUE
            """)
            
            # Indices for efficient querying
            session.run("""
                CREATE INDEX hilbert_index IF NOT EXISTS
                FOR (n:GenomeNode) ON (n.hilbert_index)
            """)
            
            session.run("""
                CREATE INDEX genomic_position IF NOT EXISTS
                FOR (n:GenomeNode) ON (n.chromosome, n.start_pos, n.end_pos)
            """)
            
            session.run("""
                CREATE INDEX scale_level_index IF NOT EXISTS
                FOR (n:GenomeNode) ON (n.scale_level)
            """)
            
            session.run("""
                CREATE INDEX kmer_index IF NOT EXISTS
                FOR (k:Kmer) ON (k.sequence)
            """)
    
    async def import_diploid_genome(self, individual_id: str, 
                                  maternal_fasta: str, paternal_fasta: str,
                                  assembly_source: str = "custom") -> bool:
        """Import a complete diploid genome"""
        try:
            logging.info(f"Starting import of diploid genome for {individual_id}")
            
            # Process maternal haplotype
            maternal_nodes = await self.fasta_processor.process_fasta_file(
                maternal_fasta, individual_id, 0, assembly_source
            )
            
            # Process paternal haplotype  
            paternal_nodes = await self.fasta_processor.process_fasta_file(
                paternal_fasta, individual_id, 1, assembly_source
            )
            
            all_nodes = maternal_nodes + paternal_nodes
            
            # Import to Neo4j in batches
            await self._batch_import_nodes(all_nodes)
            await self._create_linear_edges(all_nodes)
            await self._create_individual_metadata(individual_id, assembly_source)
            
            # Update k-mer database
            await self._import_kmers_to_neo4j()
            
            self.imported_genomes.add(individual_id)
            logging.info(f"Successfully imported genome for {individual_id}")
            return True
            
        except Exception as e:
            logging.error(f"Failed to import genome for {individual_id}: {e}")
            return False
    
    async def _batch_import_nodes(self, nodes: List[PangenomeNode], batch_size: int = 1000):
        """Import nodes to Neo4j in batches"""
        for i in range(0, len(nodes), batch_size):
            batch = nodes[i:i + batch_size]
            
            with self.driver.session() as session:
                # Prepare batch data
                node_data = []
                for node in batch:
                    node_data.append({
                        'id': node.node_id,
                        'sequence': node.sequence,
                        'chromosome': node.chromosome,
                        'start_pos': node.start_pos,
                        'end_pos': node.end_pos,
                        'scale_level': node.scale_level,
                        'frequency': node.frequency,
                        'hilbert_index': node.hilbert_index,
                        'individual_id': node.annotations.get('individual_id'),
                        'haplotype': node.annotations.get('haplotype'),
                        'assembly_source': node.annotations.get('assembly_source'),
                        'segment_length': node.annotations.get('segment_length')
                    })
                
                # Batch insert
                session.run("""
                    UNWIND $nodes as node
                    CREATE (n:GenomeNode {
                        id: node.id,
                        sequence: node.sequence,
                        chromosome: node.chromosome,
                        start_pos: node.start_pos,
                        end_pos: node.end_pos,
                        scale_level: node.scale_level,
                        frequency: node.frequency,
                        hilbert_index: node.hilbert_index,
                        individual_id: node.individual_id,
                        haplotype: node.haplotype,
                        assembly_source: node.assembly_source,
                        segment_length: node.segment_length,
                        created_at: datetime()
                    })
                """, nodes=node_data)
            
            logging.info(f"Imported batch {i//batch_size + 1} ({len(batch)} nodes)")
    
    async def _create_linear_edges(self, nodes: List[PangenomeNode]):
        """Create linear progression edges between consecutive nodes"""
        # Group nodes by chromosome, haplotype, and scale level
        grouped_nodes = defaultdict(list)
        for node in nodes:
            key = (node.chromosome, node.annotations['haplotype'], node.scale_level)
            grouped_nodes[key].append(node)
        
        # Sort and create edges for each group
        for key, node_group in grouped_nodes.items():
            # Sort by start position
            node_group.sort(key=lambda x: x.start_pos)
            
            # Create edges between consecutive nodes
            edge_data = []
            for i in range(len(node_group) - 1):
                current_node = node_group[i]
                next_node = node_group[i + 1]
                
                # Only create edge if nodes are adjacent or overlapping
                if next_node.start_pos <= current_node.end_pos + 1000:  # Allow small gaps
                    edge_data.append({
                        'from_id': current_node.node_id,
                        'to_id': next_node.node_id,
                        'edge_type': 'linear',
                        'frequency': 1.0,
                        'distance': next_node.start_pos - current_node.end_pos
                    })
            
            # Batch create edges
            if edge_data:
                with self.driver.session() as session:
                    session.run("""
                        UNWIND $edges as edge
                        MATCH (a:GenomeNode {id: edge.from_id}),
                              (b:GenomeNode {id: edge.to_id})
                        CREATE (a)-[r:CONNECTS {
                            edge_type: edge.edge_type,
                            frequency: edge.frequency,
                            distance: edge.distance,
                            created_at: datetime()
                        }]->(b)
                    """, edges=edge_data)
    
    async def _create_individual_metadata(self, individual_id: str, assembly_source: str):
        """Create individual metadata node"""
        with self.driver.session() as session:
            session.run("""
                CREATE (i:Individual {
                    id: $individual_id,
                    assembly_source: $assembly_source,
                    import_date: datetime(),
                    genome_status: 'complete'
                })
            """, individual_id=individual_id, assembly_source=assembly_source)
            
            # Link individual to their genome nodes
            session.run("""
                MATCH (i:Individual {id: $individual_id}),
                      (n:GenomeNode {individual_id: $individual_id})
                CREATE (i)-[:HAS_GENOME]->(n)
            """, individual_id=individual_id)
    
    async def _import_kmers_to_neo4j(self):
        """Import k-mer index to Neo4j for efficient short read mapping"""
        kmer_data = []
        
        for kmer, positions in self.fasta_processor.kmer_index.kmer_positions.items():
            kmer_data.append({
                'sequence': kmer,
                'positions': [{'node_id': pos[0], 'position': pos[1]} for pos in positions],
                'node_count': len(set(pos[0] for pos in positions))
            })
        
        # Import in batches
        batch_size = 5000
        for i in range(0, len(kmer_data), batch_size):
            batch = kmer_data[i:i + batch_size]
            
            with self.driver.session() as session:
                session.run("""
                    UNWIND $kmers as kmer
                    CREATE (k:Kmer {
                        sequence: kmer.sequence,
                        node_count: kmer.node_count,
                        created_at: datetime()
                    })
                    WITH k, kmer
                    UNWIND kmer.positions as pos
                    MATCH (n:GenomeNode {id: pos.node_id})
                    CREATE (k)-[:FOUND_IN {position: pos.position}]->(n)
                """, kmers=batch)

class ShortReadMapper:
    """Map short reads to pangenome paths"""
    
    def __init__(self, driver, kmer_index: KmerIndex):
        self.driver = driver
        self.kmer_index = kmer_index
        self.min_kmer_matches = 5
        self.min_mapping_score = 0.8
    
    async def map_short_reads(self, reads: List[ShortRead], 
                            individual_id: Optional[str] = None) -> List[Dict]:
        """Map short reads to pangenome paths"""
        mapping_results = []
        
        for read in reads:
            # Find candidate nodes using k-mer matching
            candidate_nodes = self.kmer_index.find_matching_nodes(
                read.sequence, self.min_kmer_matches
            )
            
            if not candidate_nodes:
                continue
            
            # Score alignments and find best paths
            read_mappings = await self._score_read_alignments(
                read, candidate_nodes, individual_id
            )
            
            if read_mappings:
                mapping_results.append({
                    'read_id': read.read_id,
                    'mappings': read_mappings,
                    'best_score': max(m['score'] for m in read_mappings)
                })
        
        return mapping_results
    
    async def _score_read_alignments(self, read: ShortRead, 
                                   candidate_nodes: Dict[str, int],
                                   individual_id: Optional[str] = None) -> List[Dict]:
        """Score potential alignments for a read"""
        alignments = []
        
        # Query Neo4j for detailed node information
        with self.driver.session() as session:
            node_ids = list(candidate_nodes.keys())
            
            # Add individual-specific filtering if provided
            individual_filter = ""
            if individual_id:
                individual_filter = "AND n.individual_id = $individual_id"
            
            result = session.run(f"""
                MATCH (n:GenomeNode)
                WHERE n.id IN $node_ids {individual_filter}
                RETURN n.id as node_id, n.sequence as sequence, 
                       n.chromosome as chromosome, n.start_pos as start_pos,
                       n.end_pos as end_pos, n.frequency as frequency,
                       n.individual_id as individual_id, n.haplotype as haplotype
            """, node_ids=node_ids, individual_id=individual_id)
            
            for record in result:
                node_data = dict(record)
                
                # Calculate alignment score
                alignment_score = self._calculate_alignment_score(
                    read.sequence, node_data['sequence'], 
                    candidate_nodes[node_data['node_id']]
                )
                
                if alignment_score >= self.min_mapping_score:
                    alignments.append({
                        'node_id': node_data['node_id'],
                        'chromosome': node_data['chromosome'],
                        'start_pos': node_data['start_pos'],
                        'end_pos': node_data['end_pos'],
                        'score': alignment_score,
                        'kmer_matches': candidate_nodes[node_data['node_id']],
                        'frequency': node_data['frequency'],
                        'individual_id': node_data['individual_id'],
                        'haplotype': node_data['haplotype']
                    })
        
        # Sort by score
        alignments.sort(key=lambda x: x['score'], reverse=True)
        return alignments[:10]  # Return top 10 alignments
    
    def _calculate_alignment_score(self, read_seq: str, node_seq: str, kmer_matches: int) -> float:
        """Calculate alignment score between read and node"""
        read_len = len(read_seq)
        max_possible_kmers = read_len - self.kmer_index.k + 1
        
        if max_possible_kmers == 0:
            return 0.0
        
        # Base score from k-mer matches
        kmer_score = kmer_matches / max_possible_kmers
        
        # Bonus for exact sequence match (if read fits in node)
        if read_seq in node_seq:
            kmer_score += 0.2
        
        # Penalty for length mismatch
        if len(node_seq) > 0:
            length_ratio = min(read_len, len(node_seq)) / max(read_len, len(node_seq))
            kmer_score *= length_ratio
        
        return min(kmer_score, 1.0)
    
    async def infer_individual_paths(self, mapping_results: List[Dict], 
                                   chromosome: int, region_start: int, 
                                   region_end: int) -> Dict[str, List[str]]:
        """Infer likely genomic paths for an individual based on short read mappings"""
        # Collect all mapped nodes in the region
        region_nodes = defaultdict(list)  # haplotype -> [node_ids]
        node_scores = defaultdict(float)  # node_id -> total_score
        
        for result in mapping_results:
            for mapping in result['mappings']:
                if (mapping['chromosome'] == chromosome and 
                    mapping['start_pos'] >= region_start and 
                    mapping['end_pos'] <= region_end):
                    
                    node_id = mapping['node_id']
                    haplotype = mapping.get('haplotype', 0)
                    
                    region_nodes[haplotype].append(node_id)
                    node_scores[node_id] += mapping['score']
        
        # Find optimal paths for each haplotype
        inferred_paths = {}
        
        for haplotype, nodes in region_nodes.items():
            if nodes:
                # Remove duplicates and sort by genomic position
                unique_nodes = list(set(nodes))
                path = await self._construct_optimal_path(unique_nodes, node_scores)
                inferred_paths[f'haplotype_{haplotype}'] = path
        
        return inferred_paths
    
    async def _construct_optimal_path(self, node_ids: List[str], 
                                    node_scores: Dict[str, float]) -> List[str]:
        """Construct optimal path through nodes using Neo4j pathfinding"""
        if len(node_ids) <= 1:
            return node_ids
        
        # Query Neo4j to find connected path
        with self.driver.session() as session:
            # Sort nodes by genomic position first
            result = session.run("""
                MATCH (n:GenomeNode)
                WHERE n.id IN $node_ids
                RETURN n.id as node_id, n.start_pos as start_pos
                ORDER BY n.start_pos
            """, node_ids=node_ids)
            
            sorted_nodes = [record['node_id'] for record in result]
            
            if len(sorted_nodes) <= 1:
                return sorted_nodes
            
            # Find path connecting all nodes
            start_node = sorted_nodes[0]
            end_node = sorted_nodes[-1]
            
            path_result = session.run("""
                MATCH path = shortestPath((start:GenomeNode {id: $start_node})-
                    [r:CONNECTS*..20]->(end:GenomeNode {id: $end_node}))
                WHERE ALL(rel in relationships(path) WHERE rel.frequency > 0.01)
                RETURN [node in nodes(path) | node.id] as path_nodes
                LIMIT 1
            """, start_node=start_node, end_node=end_node)
            
            path_record = path_result.single()
            if path_record:
                return path_record['path_nodes']
            else:
                return sorted_nodes  # Fallback to sorted nodes

class PublicGenomeDownloader:
    """Download public reference genomes from various sources"""
    
    def __init__(self, download_dir: str = "./reference_genomes"):
        self.download_dir = Path(download_dir)
        self.download_dir.mkdir(parents=True, exist_ok=True)
        
        # Public genome sources
        self.sources = {
            'GRCh38': {
                'url': 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_genomic.fna.gz',
                'description': 'Human reference genome GRCh38'
            },
            'T2T-CHM13': {
                'url': 'https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz',
                'description': 'Telomere-to-telomere CHM13 assembly'
            },
            'HG002_maternal': {
                'url': 'https://s3-us-west-2.amazonaws.com/human-pangenomics/HPRC/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.maternal.f1_assembly_v2_genbank.fa.gz',
                'description': 'HG002 maternal haplotype assembly'
            },
            'HG002_paternal': {
                'url': 'https://s3-us-west-2.amazonaws.com/human-pangenomics/HPRC/HG002/assemblies/year1_f1_assembly_v2_genbank/HG002.paternal.f1_assembly_v2_genbank.fa.gz',
                'description': 'HG002 paternal haplotype assembly'
            }
        }
    
    async def download_genome(self, source_name: str) -> Optional[str]:
        """Download a public genome"""
        if source_name not in self.sources:
            logging.error(f"Unknown genome source: {source_name}")
            return None
        
        source_info = self.sources[source_name]
        url = source_info['url']
        
        # Generate local filename
        filename = f"{source_name}.fa.gz"
        local_path = self.download_dir / filename
        
        if local_path.exists():
            logging.info(f"Genome {source_name} already downloaded")
            return str(local_path)
        
        try:
            logging.info(f"Downloading {source_name} from {url}")
            
            async with aiohttp.ClientSession() as session:
                async with session.get(url) as response:
                    if response.status == 200:
                        with open(local_path, 'wb') as f:
                            async for chunk in response.content.iter_chunked(8192):
                                f.write(chunk)
                        
                        logging.info(f"Successfully downloaded {source_name}")
                        return str(local_path)
                    else:
                        logging.error(f"Failed to download {source_name}: HTTP {response.status}")
                        return None
        
        except Exception as e:
            logging.error(f"Error downloading {source_name}: {e}")
            return None
    
    async def decompress_genome(self, compressed_path: str) -> str:
        """Decompress a gzipped genome file"""
        input_path = Path(compressed_path)
        output_path = input_path.with_suffix('')  # Remove .gz extension
        
        if output_path.exists():
            return str(output_path)
        
        with gzip.open(input_path, 'rt') as f_in:
            with open(output_path, 'w') as f_out:
                f_out.write(f_in.read())
        
        logging.info(f"Decompressed {input_path} to {output_path}")
        return str(output_path)

# Main orchestration class
class GenomeDatabaseBuilder:
    """Main orchestrator for building the genome database"""
    
    def __init__(self, config_path: str = "genome_db_config.yaml"):
        self.config = self._load_config(config_path)
        self.neo4j_manager = DockerNeo4jManager(**self.config.get('neo4j', {}))
        self.genome_importer = None
        self.short_read_mapper = None
        self.downloader = PublicGenomeDownloader(**self.config.get('download', {}))
        
    def _load_config(self, config_path: str) -> Dict:
        """Load configuration from YAML file"""
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
        
        if Path(config_path).exists():
            with open(config_path, 'r') as f:
                user_config = yaml.safe_load(f)
            # Merge configurations
            for key, value in user_config.items():
                if key in default_config and isinstance(value, dict):
                    default_config[key].update(value)
                else:
                    default_config[key] = value
        
        return default_config
    
    async def initialize_database(self) -> bool:
        """Initialize the Neo4j database"""
        success = self.neo4j_manager.setup_neo4j_container()
        if success:
            self.genome_importer = GenomeImporter(self.neo4j_manager)
            self.short_read_mapper = ShortReadMapper(
                self.neo4j_manager.get_connection(),
                self.genome_importer.fasta_processor.kmer_index
            )
        return success
    
    async def import_public_genomes(self, genome_sources: List[str]) -> bool:
        """Download and import public reference genomes"""
        try:
            for source in genome_sources:
                # Download genome
                compressed_path = await self.downloader.download_genome(source)
                if not compressed_path:
                    continue
                
                # Decompress
                fasta_path = await self.downloader.decompress_genome(compressed_path)
                
                # Import to database
                if 'maternal' in source.lower():
                    haplotype = 0
                elif 'paternal' in source.lower():
                    haplotype = 1
                else:
                    haplotype = 0  # Default for reference genomes
                
                individual_id = source.replace('_maternal', '').replace('_paternal', '')
                
                # Process as single haplotype
                nodes = await self.genome_importer.fasta_processor.process_fasta_file(
                    fasta_path, individual_id, haplotype, source
                )
                
                await self.genome_importer._batch_import_nodes(nodes)
                await self.genome_importer._create_linear_edges(nodes)
                
                logging.info(f"Successfully imported {source}")
            
            return True
            
        except Exception as e:
            logging.error(f"Failed to import public genomes: {e}")
            return False
    
    async def analyze_short_reads(self, fastq_path: str, individual_id: str,
                                chromosome: int, region_start: int, 
                                region_end: int) -> Dict:
        """Analyze short reads and infer genomic paths"""
        # Parse FASTQ file (simplified - would use BioPython in production)
        reads = await self._parse_fastq_file(fastq_path)
        
        # Map reads to pangenome
        mapping_results = await self.short_read_mapper.map_short_reads(reads, individual_id)
        
        # Infer paths
        inferred_paths = await self.short_read_mapper.infer_individual_paths(
            mapping_results, chromosome, region_start, region_end
        )
        
        return {
            'individual_id': individual_id,
            'region': {'chromosome': chromosome, 'start': region_start, 'end': region_end},
            'total_reads': len(reads),
            'mapped_reads': len(mapping_results),
            'mapping_rate': len(mapping_results) / len(reads) if reads else 0,
            'inferred_paths': inferred_paths,
            'mapping_results': mapping_results[:100]  # Sample of results
        }
    
    async def _parse_fastq_file(self, fastq_path: str) -> List[ShortRead]:
        """Parse FASTQ file and return ShortRead objects"""
        reads = []
        
        async with aiofiles.open(fastq_path, 'r') as f:
            lines = []
            async for line in f:
                lines.append(line.strip())
                
                # Process every 4 lines (FASTQ format)
                if len(lines) == 4:
                    read_id = lines[0][1:]  # Remove @
                    sequence = lines[1]
                    quality = lines[3]
                    
                    reads.append(ShortRead(
                        read_id=read_id,
                        sequence=sequence,
                        quality=quality
                    ))
                    
                    lines = []
        
        return reads
    
    def cleanup(self):
        """Cleanup resources"""
        if self.neo4j_manager:
            self.neo4j_manager.stop_container()

# Example usage and demo
async def demo_genome_database():
    """Demonstrate the genome database system"""
    
    # Initialize system
    db_builder = GenomeDatabaseBuilder()
    
    try:
        # Setup database
        logging.info("Initializing Neo4j database...")
        success = await db_builder.initialize_database()
        if not success:
            logging.error("Failed to initialize database")
            return
        
        # Import public reference genomes
        logging.info("Importing public genomes...")
        public_genomes = ['T2T-CHM13', 'HG002_maternal', 'HG002_paternal']
        await db_builder.import_public_genomes(public_genomes)
        
        # Example: Import custom diploid genome
        # await db_builder.genome_importer.import_diploid_genome(
        #     "custom_individual_001",
        #     "/path/to/maternal.fa",
        #     "/path/to/paternal.fa",
        #     "custom_assembly"
        # )
        
        # Example: Analyze short reads
        # analysis_result = await db_builder.analyze_short_reads(
        #     "/path/to/reads.fastq",
        #     "individual_001",
        #     chromosome=1,
        #     region_start=1000000,
        #     region_end=2000000
        # )
        # print(f"Mapping rate: {analysis_result['mapping_rate']:.2%}")
        # print(f"Inferred paths: {analysis_result['inferred_paths']}")
        
        logging.info("Demo completed successfully!")
        
    finally:
        db_builder.cleanup()

if __name__ == "__main__":
    # Configure logging BEFORE importing other modules
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.StreamHandler(),
            logging.FileHandler('logs/genome_db.log', mode='a')
        ]
    )
    
    # Create logs directory if it doesn't exist
    Path('logs').mkdir(exist_ok=True)
    
    # Run demo
    asyncio.run(demo_genome_database())
