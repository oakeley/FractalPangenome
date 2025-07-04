"""
Hilbert Curve Fractal Pangenome Architecture
============================================
A scalable system for representing all possible human genomes as navigable paths
through a fractal data structure, analogous to Google Maps street navigation.
"""

import numpy as np
import networkx as nx
from dataclasses import dataclass, field
from typing import Dict, List, Tuple, Optional, Set, Union
from collections import defaultdict
from abc import ABC, abstractmethod
import hashlib
import json
from datetime import datetime, timezone
import logging

# Core Hilbert Curve Implementation
class HilbertCurve:
    """N-dimensional Hilbert curve for genomic coordinate mapping"""
    
    def __init__(self, dimensions: int = 4, order: int = 16):
        """
        Initialize Hilbert curve mapper
        Args:
            dimensions: Number of dimensions (chromosome, position, scale, time)
            order: Hilbert curve order (determines resolution)
        """
        self.dimensions = dimensions
        self.order = order
        self.max_coord = (1 << order) - 1
        
    def encode(self, coordinates: Tuple[int, ...]) -> int:
        """Convert multi-dimensional coordinates to Hilbert curve index"""
        # Normalize coordinates to fit in order bits
        normalized = []
        for coord in coordinates:
            normalized.append(min(coord, self.max_coord))
        
        return self._hilbert_encode(normalized, self.order)
    
    def decode(self, hilbert_index: int) -> Tuple[int, ...]:
        """Convert Hilbert curve index back to coordinates"""
        return self._hilbert_decode(hilbert_index, self.order, self.dimensions)
    
    def _hilbert_encode(self, coords: List[int], order: int) -> int:
        """Encode coordinates to Hilbert index (Gray code based)"""
        # Implementation of N-dimensional Hilbert encoding
        # Simplified version - production would use optimized algorithm
        result = 0
        for i in range(order):
            result <<= self.dimensions
            for dim in range(self.dimensions):
                bit = (coords[dim] >> (order - 1 - i)) & 1
                result |= bit << (self.dimensions - 1 - dim)
        return result
    
    def _hilbert_decode(self, index: int, order: int, dims: int) -> Tuple[int, ...]:
        """Decode Hilbert index back to coordinates"""
        coords = [0] * dims
        for i in range(order):
            for dim in range(dims):
                bit = (index >> (order * dims - 1 - i * dims - dim)) & 1
                coords[dim] |= bit << (order - 1 - i)
        return tuple(coords)

# Genomic Coordinate System
@dataclass
class GenomicCoordinate:
    """Represents a position in the genome with fractal level information"""
    chromosome: int
    position: int
    scale_level: int  # 0=nucleotide, 1=gene, 2=chromosome, etc.
    timestamp: int    # For temporal annotations
    haplotype: int = 0  # 0 or 1 for diploid genomes
    
    def to_hilbert_coords(self) -> Tuple[int, int, int, int]:
        """Convert to Hilbert curve coordinates"""
        return (self.chromosome, self.position, self.scale_level, self.timestamp)

# Pangenome Graph Node
@dataclass
class PangenomeNode:
    """Represents a segment of the pangenome graph"""
    node_id: str
    sequence: str
    chromosome: int
    start_pos: int
    end_pos: int
    scale_level: int
    frequency: float = 0.0  # Population frequency
    annotations: Dict[str, any] = field(default_factory=dict)
    hilbert_index: Optional[int] = None
    
    def __post_init__(self):
        if self.hilbert_index is None:
            coord = GenomicCoordinate(
                self.chromosome, self.start_pos, self.scale_level, 0
            )
            curve = HilbertCurve()
            self.hilbert_index = curve.encode(coord.to_hilbert_coords())

@dataclass
class PangenomeEdge:
    """Represents a transition between pangenome segments"""
    from_node: str
    to_node: str
    edge_type: str  # 'linear', 'variant', 'structural', 'inversion'
    frequency: float = 0.0
    metadata: Dict[str, any] = field(default_factory=dict)

# Individual Genome Path
@dataclass
class GenomePath:
    """Represents an individual's path through the pangenome"""
    individual_id: str
    haplotype: int  # 0 or 1
    path_nodes: List[str]
    variations: List[Dict[str, any]]
    metadata: Dict[str, any] = field(default_factory=dict)
    
    def get_route_signature(self) -> str:
        """Generate a unique signature for this genomic route"""
        path_str = "->".join(self.path_nodes)
        return hashlib.sha256(path_str.encode()).hexdigest()[:16]

# Fractal Database Architecture
class FractalPangenomeDB:
    """Main database class for fractal pangenome storage and querying"""
    
    def __init__(self, neo4j_uri: str = None, influx_config: Dict = None):
        self.hilbert_curve = HilbertCurve(dimensions=4, order=20)
        self.graph = nx.MultiDiGraph()  # Local graph for prototyping
        self.annotations_db = {}  # Placeholder for InfluxDB
        self.spatial_index = {}  # R-tree index for range queries
        self.individuals = {}  # Individual genome paths
        
        # Initialize external databases
        self._init_neo4j(neo4j_uri)
        self._init_influx(influx_config)
        
    def _init_neo4j(self, uri: str):
        """Initialize Neo4j connection for graph storage"""
        # Placeholder - would use neo4j driver in production
        logging.info(f"Neo4j connection initialized: {uri}")
        
    def _init_influx(self, config: Dict):
        """Initialize InfluxDB for temporal annotations"""
        # Placeholder - would use influxdb client in production
        logging.info("InfluxDB connection initialized")
    
    def add_pangenome_node(self, node: PangenomeNode):
        """Add a node to the pangenome graph"""
        # Add to NetworkX graph (local)
        self.graph.add_node(
            node.node_id,
            sequence=node.sequence,
            chromosome=node.chromosome,
            start_pos=node.start_pos,
            end_pos=node.end_pos,
            scale_level=node.scale_level,
            frequency=node.frequency,
            hilbert_index=node.hilbert_index,
            annotations=node.annotations
        )
        
        # Add to spatial index for efficient range queries
        hilbert_idx = node.hilbert_index
        if hilbert_idx not in self.spatial_index:
            self.spatial_index[hilbert_idx] = []
        self.spatial_index[hilbert_idx].append(node.node_id)
        
        # In production: store in Neo4j
        self._store_node_neo4j(node)
    
    def add_pangenome_edge(self, edge: PangenomeEdge):
        """Add an edge to the pangenome graph"""
        self.graph.add_edge(
            edge.from_node,
            edge.to_node,
            edge_type=edge.edge_type,
            frequency=edge.frequency,
            metadata=edge.metadata
        )
        
        # In production: store in Neo4j
        self._store_edge_neo4j(edge)
    
    def add_individual_genome(self, individual_id: str, 
                            maternal_path: GenomePath, 
                            paternal_path: GenomePath):
        """Add diploid genome for an individual"""
        self.individuals[individual_id] = {
            'maternal': maternal_path,
            'paternal': paternal_path,
            'metadata': {
                'added_timestamp': datetime.now(timezone.utc).isoformat(),
                'maternal_signature': maternal_path.get_route_signature(),
                'paternal_signature': paternal_path.get_route_signature()
            }
        }
    
    def query_genomic_region(self, chromosome: int, start: int, end: int, 
                           scale_level: int = 0) -> List[PangenomeNode]:
        """Query nodes in a genomic region using Hilbert curve indexing"""
        # Convert region to Hilbert coordinate range
        start_coord = GenomicCoordinate(chromosome, start, scale_level, 0)
        end_coord = GenomicCoordinate(chromosome, end, scale_level, 0)
        
        start_hilbert = self.hilbert_curve.encode(start_coord.to_hilbert_coords())
        end_hilbert = self.hilbert_curve.encode(end_coord.to_hilbert_coords())
        
        # Query spatial index
        result_nodes = []
        for hilbert_idx in range(start_hilbert, end_hilbert + 1):
            if hilbert_idx in self.spatial_index:
                for node_id in self.spatial_index[hilbert_idx]:
                    node_data = self.graph.nodes[node_id]
                    if (node_data['chromosome'] == chromosome and 
                        node_data['start_pos'] >= start and 
                        node_data['end_pos'] <= end and
                        node_data['scale_level'] == scale_level):
                        
                        result_nodes.append(PangenomeNode(
                            node_id=node_id,
                            sequence=node_data['sequence'],
                            chromosome=node_data['chromosome'],
                            start_pos=node_data['start_pos'],
                            end_pos=node_data['end_pos'],
                            scale_level=node_data['scale_level'],
                            frequency=node_data['frequency'],
                            annotations=node_data['annotations'],
                            hilbert_index=node_data['hilbert_index']
                        ))
        
        return result_nodes
    
    def find_individual_path(self, individual_id: str, chromosome: int, 
                           start: int, end: int) -> Dict[str, GenomePath]:
        """Get individual's diploid paths through a genomic region"""
        if individual_id not in self.individuals:
            return {}
        
        individual_data = self.individuals[individual_id]
        maternal = individual_data['maternal']
        paternal = individual_data['paternal']
        
        # Filter paths to region of interest
        region_nodes = self.query_genomic_region(chromosome, start, end)
        region_node_ids = {node.node_id for node in region_nodes}
        
        maternal_region = [node for node in maternal.path_nodes 
                          if node in region_node_ids]
        paternal_region = [node for node in paternal.path_nodes 
                          if node in region_node_ids]
        
        return {
            'maternal': GenomePath(
                individual_id=individual_id,
                haplotype=0,
                path_nodes=maternal_region,
                variations=maternal.variations,
                metadata=maternal.metadata
            ),
            'paternal': GenomePath(
                individual_id=individual_id,
                haplotype=1,
                path_nodes=paternal_region,
                variations=paternal.variations,
                metadata=paternal.metadata
            )
        }
    
    def adaptive_resolution_query(self, query_region: Tuple[int, int, int], 
                                computational_budget: int = 1000) -> Dict:
        """
        Adaptively select resolution based on query complexity and budget
        Returns data at optimal scale level
        """
        chromosome, start, end = query_region
        region_size = end - start
        
        # Determine optimal scale level based on region size and budget
        if region_size < 1000:  # Small region - nucleotide level
            scale_level = 0
        elif region_size < 100000:  # Medium region - gene level
            scale_level = 1
        elif region_size < 10000000:  # Large region - chromosome segment
            scale_level = 2
        else:  # Very large region - chromosome level
            scale_level = 3
        
        # Adjust based on computational budget
        if computational_budget < 100:
            scale_level = min(scale_level + 1, 3)
        
        nodes = self.query_genomic_region(chromosome, start, end, scale_level)
        
        return {
            'scale_level': scale_level,
            'nodes': nodes,
            'region_size': region_size,
            'node_count': len(nodes),
            'computational_cost': len(nodes) * (scale_level + 1)
        }
    
    def update_annotation(self, node_id: str, annotation_key: str, 
                         annotation_value: any, timestamp: datetime = None):
        """Update temporal annotation for a node"""
        if timestamp is None:
            timestamp = datetime.now(timezone.utc)
        
        # Update in graph
        if node_id in self.graph.nodes:
            if 'annotations' not in self.graph.nodes[node_id]:
                self.graph.nodes[node_id]['annotations'] = {}
            self.graph.nodes[node_id]['annotations'][annotation_key] = {
                'value': annotation_value,
                'timestamp': timestamp.isoformat()
            }
        
        # In production: store in InfluxDB for temporal queries
        self._store_annotation_influx(node_id, annotation_key, 
                                    annotation_value, timestamp)
    
    def _store_node_neo4j(self, node: PangenomeNode):
        """Store node in Neo4j (placeholder)"""
        # Production implementation would use Neo4j driver
        cypher_query = f"""
        CREATE (n:GenomeNode {{
            id: '{node.node_id}',
            sequence: '{node.sequence}',
            chromosome: {node.chromosome},
            start_pos: {node.start_pos},
            end_pos: {node.end_pos},
            scale_level: {node.scale_level},
            frequency: {node.frequency},
            hilbert_index: {node.hilbert_index}
        }})
        """
        logging.debug(f"Would execute: {cypher_query}")
    
    def _store_edge_neo4j(self, edge: PangenomeEdge):
        """Store edge in Neo4j (placeholder)"""
        cypher_query = f"""
        MATCH (a:GenomeNode {{id: '{edge.from_node}'}}),
              (b:GenomeNode {{id: '{edge.to_node}'}})
        CREATE (a)-[r:{edge.edge_type.upper()} {{
            frequency: {edge.frequency},
            metadata: '{json.dumps(edge.metadata)}'
        }}]->(b)
        """
        logging.debug(f"Would execute: {cypher_query}")
    
    def _store_annotation_influx(self, node_id: str, key: str, value: any, 
                               timestamp: datetime):
        """Store annotation in InfluxDB (placeholder)"""
        # Production implementation would use InfluxDB client
        influx_point = f"""
        annotations,node_id={node_id},annotation_key={key} 
        value="{value}" {int(timestamp.timestamp() * 1e9)}
        """
        logging.debug(f"Would write to InfluxDB: {influx_point}")

# Utility Functions for Pangenome Construction
class PangenomeBuilder:
    """Helper class for building pangenome from reference and variants"""
    
    def __init__(self, fractal_db: FractalPangenomeDB):
        self.db = fractal_db
        self.node_counter = 0
    
    def build_from_reference(self, reference_sequence: str, chromosome: int,
                           segment_size: int = 1000):
        """Build initial pangenome from reference sequence"""
        sequence_length = len(reference_sequence)
        
        for scale_level in range(4):  # Build multiple scale levels
            if scale_level == 0:  # Nucleotide level
                step_size = segment_size
            elif scale_level == 1:  # Gene level  
                step_size = segment_size * 10
            elif scale_level == 2:  # Chromosome segment
                step_size = segment_size * 100
            else:  # Chromosome level
                step_size = sequence_length
            
            for start in range(0, sequence_length, step_size):
                end = min(start + step_size, sequence_length)
                node_id = f"ref_chr{chromosome}_{scale_level}_{self.node_counter}"
                self.node_counter += 1
                
                node = PangenomeNode(
                    node_id=node_id,
                    sequence=reference_sequence[start:end],
                    chromosome=chromosome,
                    start_pos=start,
                    end_pos=end,
                    scale_level=scale_level,
                    frequency=1.0  # Reference frequency
                )
                
                self.db.add_pangenome_node(node)
                
                # Add linear progression edges within same scale
                if start > 0:
                    prev_node_id = f"ref_chr{chromosome}_{scale_level}_{self.node_counter-2}"
                    if prev_node_id in self.db.graph.nodes:
                        edge = PangenomeEdge(
                            from_node=prev_node_id,
                            to_node=node_id,
                            edge_type='linear',
                            frequency=1.0
                        )
                        self.db.add_pangenome_edge(edge)
    
    def add_structural_variant(self, chromosome: int, start: int, end: int,
                             variant_sequence: str, variant_type: str,
                             frequency: float = 0.1):
        """Add structural variant to pangenome"""
        variant_node_id = f"sv_{variant_type}_chr{chromosome}_{start}_{self.node_counter}"
        self.node_counter += 1
        
        # Create variant node
        variant_node = PangenomeNode(
            node_id=variant_node_id,
            sequence=variant_sequence,
            chromosome=chromosome,
            start_pos=start,
            end_pos=end,
            scale_level=0,  # Start at nucleotide level
            frequency=frequency,
            annotations={'variant_type': variant_type}
        )
        
        self.db.add_pangenome_node(variant_node)
        
        # Find reference nodes that this variant replaces
        ref_nodes = self.db.query_genomic_region(chromosome, start, end, 0)
        
        # Create alternative path through variant
        for ref_node in ref_nodes:
            # Edge from previous reference to variant
            edge_to_variant = PangenomeEdge(
                from_node=ref_node.node_id,
                to_node=variant_node_id,
                edge_type='variant',
                frequency=frequency,
                metadata={'variant_type': variant_type}
            )
            self.db.add_pangenome_edge(edge_to_variant)

# Example Usage and Testing
def demo_fractal_pangenome():
    """Demonstrate the fractal pangenome architecture"""
    
    # Initialize database
    db = FractalPangenomeDB()
    builder = PangenomeBuilder(db)
    
    # Build reference pangenome for chromosome 1 (simplified)
    reference_seq = "ATCGATCGATCG" * 1000  # 12kb sequence
    builder.build_from_reference(reference_seq, chromosome=1, segment_size=100)
    
    # Add some structural variants
    builder.add_structural_variant(
        chromosome=1, start=500, end=600,
        variant_sequence="GGGGGGGGGG",  # Insertion
        variant_type="insertion",
        frequency=0.15
    )
    
    # Add individual genomes
    maternal_path = GenomePath(
        individual_id="individual_001",
        haplotype=0,
        path_nodes=["ref_chr1_0_1", "ref_chr1_0_2", "sv_insertion_chr1_500_100"],
        variations=[{"type": "SNP", "position": 750, "allele": "T"}]
    )
    
    paternal_path = GenomePath(
        individual_id="individual_001", 
        haplotype=1,
        path_nodes=["ref_chr1_0_1", "ref_chr1_0_2", "ref_chr1_0_3"],
        variations=[{"type": "SNP", "position": 850, "allele": "G"}]
    )
    
    db.add_individual_genome("individual_001", maternal_path, paternal_path)
    
    # Demonstrate adaptive resolution querying
    print("=== Adaptive Resolution Query Demo ===")
    
    # Small region query
    small_query = db.adaptive_resolution_query((1, 500, 600), computational_budget=1000)
    print(f"Small region query - Scale level: {small_query['scale_level']}, Nodes: {small_query['node_count']}")
    
    # Large region query with limited budget
    large_query = db.adaptive_resolution_query((1, 0, 10000), computational_budget=50)
    print(f"Large region query - Scale level: {large_query['scale_level']}, Nodes: {large_query['node_count']}")
    
    # Query individual paths
    individual_paths = db.find_individual_path("individual_001", 1, 400, 700)
    print(f"Individual paths found: {len(individual_paths)}")
    
    # Update annotations
    db.update_annotation("ref_chr1_0_1", "gene_name", "GENE_ABC123")
    db.update_annotation("ref_chr1_0_1", "function", "DNA_repair")
    
    print("Demo completed successfully!")

if __name__ == "__main__":
    demo_fractal_pangenome()
