# %% [markdown]
"""
# Fractal Pangenome Analysis - Getting Started

This notebook demonstrates the core functionality of the fractal pangenome database system.
We'll explore how genomic data is represented as navigable "street maps" and perform
various analysis tasks.

## Table of Contents
1. Database Connection and Setup
2. Exploring the Pangenome Street Map
3. Individual Genome Routes Analysis
4. Short Read Mapping and Path Inference
5. Population-Level Analysis
6. Visualization of Genomic Paths
7. Advanced Queries and Custom Analysis
"""

# %% [markdown]
"""
## 1. Database Connection and Setup

First, let's connect to our fractal pangenome database and explore its structure.
"""

# %%
import sys
import os
import asyncio
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import networkx as nx
from neo4j import GraphDatabase
import redis
from collections import defaultdict, Counter
import json
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

# Set up plotting
plt.style.use('seaborn-v0_8')
sns.set_palette("husl")

# %%
# Database connection configuration
NEO4J_URI = "bolt://localhost:7687"
NEO4J_USER = "neo4j"
NEO4J_PASSWORD = "genomics123"
REDIS_HOST = "localhost"
REDIS_PORT = 6379

# Connect to databases
def get_neo4j_driver():
    return GraphDatabase.driver(NEO4J_URI, auth=(NEO4J_USER, NEO4J_PASSWORD))

def get_redis_client():
    return redis.Redis(host=REDIS_HOST, port=REDIS_PORT, decode_responses=True)

# Test connections
try:
    driver = get_neo4j_driver()
    with driver.session() as session:
        result = session.run("RETURN 'Connected to Neo4j!' as message")
        print(result.single()['message'])
    
    redis_client = get_redis_client()
    redis_client.ping()
    print("Connected to Redis!")
    
except Exception as e:
    print(f"Connection error: {e}")
    print("Make sure the database is running: python -m genome_cli database init")

# %% [markdown]
"""
## 2. Exploring the Pangenome Street Map

Let's explore how our genomic data is structured as a navigable street map.
"""

# %%
def get_database_overview():
    """Get an overview of the pangenome database structure"""
    with driver.session() as session:
        # Basic statistics
        stats_query = """
        MATCH (n:GenomeNode)
        RETURN 
            count(n) as total_nodes,
            count(DISTINCT n.individual_id) as individuals,
            count(DISTINCT n.chromosome) as chromosomes,
            min(n.start_pos) as min_position,
            max(n.end_pos) as max_position,
            avg(n.frequency) as avg_frequency
        """
        stats = session.run(stats_query).single()
        
        # Scale level distribution  
        scale_query = """
        MATCH (n:GenomeNode)
        RETURN n.scale_level as level, count(n) as count
        ORDER BY n.scale_level
        """
        scale_dist = [(r['level'], r['count']) for r in session.run(scale_query)]
        
        # Edge types
        edge_query = """
        MATCH ()-[r:CONNECTS]->()
        RETURN r.edge_type as type, count(r) as count
        ORDER BY count DESC
        """
        edge_types = [(r['type'], r['count']) for r in session.run(edge_query)]
        
        return dict(stats), scale_dist, edge_types

# Get database overview
stats, scale_distribution, edge_types = get_database_overview()

print("=== Pangenome Database Overview ===")
print(f"Total genome nodes (street segments): {stats['total_nodes']:,}")
print(f"Individuals (unique routes): {stats['individuals']:,}")
print(f"Chromosomes: {stats['chromosomes']:,}")
print(f"Genomic span: {stats['min_position']:,} - {stats['max_position']:,} bp")
print(f"Average variant frequency: {stats['avg_frequency']:.3f}")

# %%
# Visualize scale level distribution
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Scale level distribution
scale_levels, scale_counts = zip(*scale_distribution) if scale_distribution else ([], [])
scale_labels = {0: 'Nucleotide', 1: 'Gene', 2: 'Segment', 3: 'Chromosome'}

ax1.bar([scale_labels.get(l, f'Level {l}') for l in scale_levels], scale_counts)
ax1.set_title('Genomic Resolution Levels\n(Street Map Detail Levels)')
ax1.set_ylabel('Number of Nodes')
ax1.tick_params(axis='x', rotation=45)

# Edge type distribution
if edge_types:
    edge_labels, edge_counts = zip(*edge_types)
    ax2.pie(edge_counts, labels=edge_labels, autopct='%1.1f%%')
    ax2.set_title('Connection Types\n(Street Connections)')

plt.tight_layout()
plt.show()

# %% [markdown]
"""
## 3. Individual Genome Routes Analysis

Now let's analyze specific individual genome routes - how different people navigate
through the same genomic region using different "streets" (variants).
"""

# %%
def analyze_individual_routes(chromosome=1, start_pos=1000000, end_pos=2000000):
    """Analyze individual genomic routes through a specific region"""
    with driver.session() as session:
        # Find all individuals with data in this region
        individuals_query = """
        MATCH (n:GenomeNode)
        WHERE n.chromosome = $chr AND n.start_pos >= $start AND n.end_pos <= $end
        RETURN DISTINCT n.individual_id as individual, n.haplotype as haplotype,
               count(n) as nodes_in_region
        ORDER BY individual, haplotype
        """
        
        individuals = session.run(individuals_query, 
                                chr=chromosome, start=start_pos, end=end_pos).data()
        
        # Get detailed path information for each individual
        routes_data = []
        for ind_data in individuals:
            individual = ind_data['individual']
            haplotype = ind_data['haplotype']
            
            # Get the actual path through this region
            path_query = """
            MATCH (n:GenomeNode)
            WHERE n.chromosome = $chr AND n.individual_id = $individual 
              AND n.haplotype = $haplotype
              AND n.start_pos >= $start AND n.end_pos <= $end
              AND n.scale_level = 0
            RETURN n.node_id as node_id, n.start_pos as pos, n.frequency as freq,
                   n.sequence as sequence
            ORDER BY n.start_pos
            """
            
            path_nodes = session.run(path_query, 
                                   chr=chromosome, individual=individual,
                                   haplotype=haplotype, start=start_pos, end=end_pos).data()
            
            if path_nodes:
                # Calculate route characteristics
                frequencies = [node['freq'] for node in path_nodes]
                route_rarity = 1 - np.mean(frequencies)  # How rare is this route?
                
                routes_data.append({
                    'individual': individual,
                    'haplotype': f'h{haplotype}',
                    'route_type': 'maternal' if haplotype == 0 else 'paternal',
                    'nodes_count': len(path_nodes),
                    'route_rarity': route_rarity,
                    'path_nodes': path_nodes
                })
        
        return routes_data

# Analyze routes for a specific region
region_chr, region_start, region_end = 1, 1000000, 1500000
routes = analyze_individual_routes(region_chr, region_start, region_end)

print(f"=== Individual Routes Analysis ===")
print(f"Region: chr{region_chr}:{region_start:,}-{region_end:,}")
print(f"Found {len(routes)} individual routes")

if routes:
    routes_df = pd.DataFrame([{
        'Individual': r['individual'],
        'Haplotype': r['haplotype'], 
        'Route Type': r['route_type'],
        'Nodes': r['nodes_count'],
        'Route Rarity': r['route_rarity']
    } for r in routes])
    
    print("\nRoute Summary:")
    print(routes_df.to_string(index=False))

# %%
# Visualize route characteristics
if routes:
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    # Route rarity distribution
    rarities = [r['route_rarity'] for r in routes]
    axes[0,0].hist(rarities, bins=20, alpha=0.7, edgecolor='black')
    axes[0,0].set_title('Route Rarity Distribution\n(0=common path, 1=unique path)')
    axes[0,0].set_xlabel('Route Rarity Score')
    axes[0,0].set_ylabel('Number of Routes')
    
    # Nodes per route
    node_counts = [r['nodes_count'] for r in routes]
    axes[0,1].hist(node_counts, bins=15, alpha=0.7, edgecolor='black', color='orange')
    axes[0,1].set_title('Route Complexity\n(Number of Genomic Segments)')
    axes[0,1].set_xlabel('Nodes per Route')
    axes[0,1].set_ylabel('Number of Routes')
    
    # Maternal vs Paternal comparison
    route_types = [r['route_type'] for r in routes]
    type_counts = Counter(route_types)
    axes[1,0].bar(type_counts.keys(), type_counts.values(), color=['lightblue', 'lightcoral'])
    axes[1,0].set_title('Maternal vs Paternal Routes')
    axes[1,0].set_ylabel('Number of Routes')
    
    # Route rarity by haplotype
    maternal_rarity = [r['route_rarity'] for r in routes if r['route_type'] == 'maternal']
    paternal_rarity = [r['route_rarity'] for r in routes if r['route_type'] == 'paternal']
    
    axes[1,1].boxplot([maternal_rarity, paternal_rarity], labels=['Maternal', 'Paternal'])
    axes[1,1].set_title('Route Rarity by Haplotype')
    axes[1,1].set_ylabel('Route Rarity Score')
    
    plt.tight_layout()
    plt.show()

# %% [markdown]
"""
## 4. Short Read Mapping and Path Inference

This section demonstrates how short sequencing reads are mapped to the pangenome
to infer the most likely genomic paths an individual took.
"""

# %%
def simulate_short_reads(routes, num_reads=1000, read_length=150):
    """Simulate short reads from known routes for demonstration"""
    simulated_reads = []
    
    for route in routes[:3]:  # Use first 3 routes for simulation
        nodes = route['path_nodes']
        individual = route['individual']
        haplotype = route['haplotype']
        
        reads_per_route = num_reads // len(routes[:3])
        
        for i in range(reads_per_route):
            # Select random node
            if nodes:
                node = np.random.choice(nodes)
                sequence = node['sequence']
                
                if len(sequence) >= read_length:
                    # Extract random subsequence
                    start = np.random.randint(0, len(sequence) - read_length + 1)
                    read_seq = sequence[start:start + read_length]
                    
                    simulated_reads.append({
                        'read_id': f"{individual}_{haplotype}_read_{i}",
                        'sequence': read_seq,
                        'true_source': {
                            'individual': individual,
                            'haplotype': haplotype,
                            'node_id': node['node_id'],
                            'position': node['pos'] + start
                        }
                    })
    
    return simulated_reads

def map_reads_to_pangenome(reads, kmer_size=31):
    """Map reads to pangenome using k-mer matching"""
    # This is a simplified version - the full implementation is in the main module
    mapping_results = []
    
    with driver.session() as session:
        for read in reads[:20]:  # Limit for demo
            read_seq = read['sequence']
            
            # Extract k-mers from read
            kmers = [read_seq[i:i+kmer_size] for i in range(len(read_seq) - kmer_size + 1)]
            
            # Find matching nodes for each k-mer
            kmer_matches = defaultdict(int)
            
            for kmer in kmers[:5]:  # Check first 5 k-mers for speed
                kmer_query = """
                MATCH (n:GenomeNode)
                WHERE n.sequence CONTAINS $kmer
                RETURN n.node_id as node_id, n.individual_id as individual,
                       n.chromosome as chr, n.start_pos as pos
                LIMIT 10
                """
                
                matches = session.run(kmer_query, kmer=kmer).data()
                for match in matches:
                    kmer_matches[match['node_id']] += 1
            
            # Find best matches
            best_matches = sorted(kmer_matches.items(), key=lambda x: x[1], reverse=True)[:5]
            
            mapping_results.append({
                'read_id': read['read_id'],
                'true_source': read['true_source'],
                'predicted_matches': best_matches,
                'mapping_quality': len(best_matches)
            })
    
    return mapping_results

# Simulate and map reads
if routes:
    print("=== Short Read Mapping Simulation ===")
    
    # Generate simulated reads
    sim_reads = simulate_short_reads(routes, num_reads=100)
    print(f"Generated {len(sim_reads)} simulated reads")
    
    # Map reads to pangenome
    print("Mapping reads to pangenome...")
    mapping_results = map_reads_to_pangenome(sim_reads)
    
    # Analyze mapping accuracy
    correct_mappings = 0
    for result in mapping_results:
        true_individual = result['true_source']['individual']
        
        # Check if any of the top matches are from the correct individual
        for node_id, score in result['predicted_matches']:
            if true_individual in node_id:
                correct_mappings += 1
                break
    
    accuracy = correct_mappings / len(mapping_results) if mapping_results else 0
    print(f"Mapping accuracy: {accuracy:.2%}")
    
    # Show sample results
    print("\nSample Mapping Results:")
    for i, result in enumerate(mapping_results[:3]):
        print(f"\nRead {i+1}: {result['read_id']}")
        print(f"True source: {result['true_source']['individual']}")
        print(f"Top matches: {result['predicted_matches'][:2]}")

# %% [markdown]
"""
## 5. Population-Level Analysis

Let's analyze patterns across the entire population to understand genomic diversity
and identify common vs. rare genomic "routes".
"""

# %%
def analyze_population_diversity():
    """Analyze genomic diversity at the population level"""
    with driver.session() as session:
        # Variant frequency distribution
        freq_query = """
        MATCH (n:GenomeNode)
        WHERE n.scale_level = 0
        RETURN n.frequency as freq, count(n) as count
        ORDER BY n.frequency
        """
        freq_data = session.run(freq_query).data()
        
        # Chromosome-level diversity
        chr_query = """
        MATCH (n:GenomeNode)
        WHERE n.scale_level = 0
        RETURN n.chromosome as chr,
               count(DISTINCT n.individual_id) as individuals,
               avg(n.frequency) as avg_freq,
               stddev(n.frequency) as freq_std,
               count(n) as total_nodes
        ORDER BY n.chromosome
        """
        chr_diversity = session.run(chr_query).data()
        
        # Find rare variants (population-specific routes)
        rare_query = """
        MATCH (n:GenomeNode)
        WHERE n.frequency < 0.05 AND n.scale_level = 0
        RETURN n.chromosome as chr, count(n) as rare_variants,
               collect(DISTINCT n.individual_id)[0..5] as sample_individuals
        ORDER BY chr
        """
        rare_variants = session.run(rare_query).data()
        
        return freq_data, chr_diversity, rare_variants

# Analyze population diversity
freq_data, chr_diversity, rare_variants = analyze_population_diversity()

print("=== Population-Level Genomic Diversity ===")

# Display chromosome diversity
if chr_diversity:
    chr_df = pd.DataFrame(chr_diversity)
    print("\nChromosome Diversity Summary:")
    print(chr_df.to_string(index=False))

# %%
# Visualize population diversity
if freq_data and chr_diversity:
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=('Variant Frequency Distribution', 'Diversity by Chromosome',
                       'Rare Variants by Chromosome', 'Population Structure'),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # Variant frequency distribution
    frequencies = [d['freq'] for d in freq_data]
    fig.add_trace(
        go.Histogram(x=frequencies, nbinsx=50, name='Variant Frequency'),
        row=1, col=1
    )
    
    # Diversity by chromosome
    chr_nums = [d['chr'] for d in chr_diversity]
    avg_freqs = [d['avg_freq'] for d in chr_diversity]
    fig.add_trace(
        go.Scatter(x=chr_nums, y=avg_freqs, mode='lines+markers', 
                  name='Avg Frequency'),
        row=1, col=2
    )
    
    # Rare variants by chromosome
    rare_chr = [d['chr'] for d in rare_variants]
    rare_counts = [d['rare_variants'] for d in rare_variants]
    fig.add_trace(
        go.Bar(x=rare_chr, y=rare_counts, name='Rare Variants'),
        row=2, col=1
    )
    
    # Population structure (simulated)
    # In a real analysis, this would show population genetic structure
    pop_data = np.random.multivariate_normal([0, 0], [[1, 0.5], [0.5, 1]], 100)
    fig.add_trace(
        go.Scatter(x=pop_data[:, 0], y=pop_data[:, 1], mode='markers',
                  name='Population Structure'),
        row=2, col=2
    )
    
    fig.update_layout(height=800, title_text="Population Genomic Analysis")
    fig.show()

# %% [markdown]
"""
## 6. Visualization of Genomic Paths

Let's create interactive visualizations of how individuals navigate through
the genomic landscape - their unique "routes" through the pangenome.
"""

# %%
def visualize_genomic_network(chromosome=1, start_pos=1000000, end_pos=1200000):
    """Create a network visualization of genomic paths"""
    with driver.session() as session:
        # Get nodes and edges in region
        nodes_query = """
        MATCH (n:GenomeNode)
        WHERE n.chromosome = $chr AND n.start_pos >= $start AND n.end_pos <= $end
          AND n.scale_level = 0
        RETURN n.node_id as id, n.start_pos as pos, n.frequency as freq,
               n.individual_id as individual, n.haplotype as haplotype
        ORDER BY n.start_pos
        """
        
        edges_query = """
        MATCH (a:GenomeNode)-[r:CONNECTS]->(b:GenomeNode)
        WHERE a.chromosome = $chr AND a.start_pos >= $start AND a.end_pos <= $end
          AND b.chromosome = $chr AND b.start_pos >= $start AND b.end_pos <= $end
          AND a.scale_level = 0 AND b.scale_level = 0
        RETURN a.node_id as source, b.node_id as target, r.frequency as weight
        """
        
        nodes = session.run(nodes_query, chr=chromosome, start=start_pos, end=end_pos).data()
        edges = session.run(edges_query, chr=chromosome, start=start_pos, end=end_pos).data()
        
        if not nodes:
            print("No nodes found in specified region")
            return None
        
        # Create NetworkX graph
        G = nx.DiGraph()
        
        # Add nodes
        for node in nodes:
            G.add_node(node['id'], 
                      pos=node['pos'],
                      freq=node['freq'],
                      individual=node['individual'],
                      haplotype=node['haplotype'])
        
        # Add edges
        for edge in edges:
            if edge['source'] in G.nodes and edge['target'] in G.nodes:
                G.add_edge(edge['source'], edge['target'], weight=edge['weight'])
        
        return G, nodes, edges

# Create genomic network visualization
print("=== Genomic Network Visualization ===")
G, nodes, edges = visualize_genomic_network(1, 1000000, 1100000)

if G and len(G.nodes) > 0:
    print(f"Network has {len(G.nodes)} nodes and {len(G.edges)} edges")
    
    # Create layout based on genomic position
    pos = {}
    for node_id in G.nodes:
        node_data = G.nodes[node_id]
        pos[node_id] = (node_data['pos'], np.random.random() * 0.1)  # X=position, Y=random jitter
    
    # Plot network
    plt.figure(figsize=(16, 8))
    
    # Color nodes by individual
    individuals = list(set([G.nodes[n]['individual'] for n in G.nodes]))
    colors = plt.cm.Set3(np.linspace(0, 1, len(individuals)))
    individual_colors = {ind: colors[i] for i, ind in enumerate(individuals)}
    
    node_colors = [individual_colors[G.nodes[n]['individual']] for n in G.nodes]
    node_sizes = [G.nodes[n]['freq'] * 1000 + 100 for n in G.nodes]  # Size by frequency
    
    nx.draw(G, pos, node_color=node_colors, node_size=node_sizes,
            with_labels=False, edge_color='gray', alpha=0.7,
            arrows=True, arrowsize=10)
    
    plt.title('Genomic Street Map Network\n(Nodes=genomic segments, Edges=connections, Colors=individuals)')
    plt.xlabel('Genomic Position')
    plt.ylabel('Individual Routes')
    
    # Add legend
    legend_elements = [plt.Line2D([0], [0], marker='o', color='w', 
                                 markerfacecolor=individual_colors[ind], 
                                 markersize=8, label=ind) 
                      for ind in individuals[:5]]  # Show first 5 individuals
    plt.legend(handles=legend_elements, title='Individuals', bbox_to_anchor=(1.05, 1), loc='upper left')
    
    plt.tight_layout()
    plt.show()

# %%
# Interactive Plotly visualization
if G and len(G.nodes) > 0:
    # Prepare data for Plotly
    edge_x = []
    edge_y = []
    
    for edge in G.edges():
        x0, y0 = pos[edge[0]]
        x1, y1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
    
    # Create edge trace
    edge_trace = go.Scatter(x=edge_x, y=edge_y,
                           line=dict(width=0.5, color='#888'),
                           hoverinfo='none',
                           mode='lines')
    
    # Create node trace
    node_x = [pos[node][0] for node in G.nodes()]
    node_y = [pos[node][1] for node in G.nodes()]
    
    node_trace = go.Scatter(x=node_x, y=node_y,
                           mode='markers',
                           hoverinfo='text',
                           marker=dict(size=[G.nodes[n]['freq'] * 50 + 10 for n in G.nodes()],
                                     color=[hash(G.nodes[n]['individual']) % 10 for n in G.nodes()],
                                     colorscale='Viridis',
                                     line=dict(width=2)))
    
    # Add hover text
    node_text = []
    for node in G.nodes():
        node_info = G.nodes[node]
        text = f"Position: {node_info['pos']:,}<br>"
        text += f"Individual: {node_info['individual']}<br>"
        text += f"Frequency: {node_info['freq']:.3f}<br>"
        text += f"Connections: {len(list(G.neighbors(node)))}"
        node_text.append(text)
    
    node_trace.text = node_text
    
    # Create figure
    fig = go.Figure(data=[edge_trace, node_trace],
                   layout=go.Layout(
                       title='Interactive Genomic Street Map',
                       titlefont_size=16,
                       showlegend=False,
                       hovermode='closest',
                       margin=dict(b=20,l=5,r=5,t=40),
                       annotations=[ dict(
                           text="Genomic Position",
                           showarrow=False,
                           xref="paper", yref="paper",
                           x=0.005, y=-0.002 ) ],
                       xaxis=dict(showgrid=False, zeroline=False, showticklabels=True),
                       yaxis=dict(showgrid=False, zeroline=False, showticklabels=False)))
    
    fig.show()

# %% [markdown]
"""
## 7. Advanced Queries and Custom Analysis

Finally, let's demonstrate some advanced analysis capabilities including
custom Cypher queries and specialized genomic analyses.
"""

# %%
def find_complex_variants():
    """Find complex structural variants - like highway interchanges in our street map"""
    with driver.session() as session:
        # Find nodes with many connections (complex regions)
        complex_query = """
        MATCH (n:GenomeNode)-[r:CONNECTS]->(m:GenomeNode)
        WHERE n.scale_level = 0
        WITH n, count(r) as connections
        WHERE connections > 3
        RETURN n.node_id as node_id, n.chromosome as chr, n.start_pos as pos,
               n.individual_id as individual, connections
        ORDER BY connections DESC
        LIMIT 20
        """
        complex_nodes = session.run(complex_query).data()
        
        # Find rare paths (private variants)
        rare_paths_query = """
        MATCH path = (a:GenomeNode)-[:CONNECTS*2..5]->(b:GenomeNode)
        WHERE ALL(n in nodes(path) WHERE n.frequency < 0.01)
        AND length(path) >= 3
        RETURN nodes(path)[0].individual_id as individual,
               nodes(path)[0].chromosome as chr,
               [n in nodes(path) | n.start_pos] as positions,
               length(path) as path_length
        LIMIT 10
        """
        rare_paths = session.run(rare_paths_query).data()
        
        return complex_nodes, rare_paths

def calculate_genomic_distances():
    """Calculate genomic distances between individuals"""
    with driver.session() as session:
        # Find shared vs unique nodes between individuals
        comparison_query = """
        MATCH (n1:GenomeNode), (n2:GenomeNode)
        WHERE n1.individual_id < n2.individual_id
          AND n1.chromosome = n2.chromosome
          AND abs(n1.start_pos - n2.start_pos) < 1000
          AND n1.scale_level = 0 AND n2.scale_level = 0
        RETURN n1.individual_id as ind1, n2.individual_id as ind2,
               count(*) as shared_regions,
               n1.chromosome as chr
        ORDER BY shared_regions DESC
        LIMIT 20
        """
        similarities = session.run(comparison_query).data()
        
        return similarities

# Run advanced analyses
print("=== Advanced Genomic Analysis ===")

complex_nodes, rare_paths = find_complex_variants()
similarities = calculate_genomic_distances()

print(f"\nFound {len(complex_nodes)} complex genomic regions (highway interchanges)")
if complex_nodes:
    print("Top complex regions:")
    for node in complex_nodes[:5]:
        print(f"  chr{node['chr']}:{node['pos']:,} - {node['connections']} connections ({node['individual']})")

print(f"\nFound {len(rare_paths)} rare genomic paths (private routes)")
if rare_paths:
    print("Sample rare paths:")
    for path in rare_paths[:3]:
        positions = path['positions']
        print(f"  {path['individual']} chr{path['chr']}: {positions[0]:,} → {positions[-1]:,} ({path['path_length']} segments)")

print(f"\nGenome similarity analysis:")
if similarities:
    sim_df = pd.DataFrame(similarities[:10])
    print(sim_df.to_string(index=False))

# %%
# Create a genomic distance matrix visualization
if similarities:
    # Convert to distance matrix format
    individuals = list(set([s['ind1'] for s in similarities] + [s['ind2'] for s in similarities]))
    
    # Create similarity matrix
    n_ind = len(individuals)
    similarity_matrix = np.zeros((n_ind, n_ind))
    
    ind_to_idx = {ind: i for i, ind in enumerate(individuals)}
    
    for sim in similarities:
        i = ind_to_idx[sim['ind1']]
        j = ind_to_idx[sim['ind2']]
        similarity_matrix[i, j] = sim['shared_regions']
        similarity_matrix[j, i] = sim['shared_regions']  # Make symmetric
    
    # Plot heatmap
    plt.figure(figsize=(12, 10))
    sns.heatmap(similarity_matrix, 
                xticklabels=individuals, 
                yticklabels=individuals,
                annot=True, 
                fmt='.0f',
                cmap='viridis',
                cbar_kws={'label': 'Shared Genomic Regions'})
    plt.title('Genomic Similarity Between Individuals\n(Street Map Route Overlap)')
    plt.xlabel('Individual')
    plt.ylabel('Individual')
    plt.xticks(rotation=45)
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.show()

# %%
# Performance analysis
def analyze_query_performance():
    """Analyze database query performance"""
    with driver.session() as session:
        # Simple performance test
        import time
        
        test_queries = [
            ("Node count", "MATCH (n:GenomeNode) RETURN count(n)"),
            ("Region query", "MATCH (n:GenomeNode) WHERE n.chromosome = 1 AND n.start_pos BETWEEN 1000000 AND 2000000 RETURN count(n)"),
            ("Path finding", "MATCH path = (a:GenomeNode)-[:CONNECTS*1..3]->(b:GenomeNode) WHERE a.chromosome = 1 RETURN count(path) LIMIT 1000"),
            ("Frequency filter", "MATCH (n:GenomeNode) WHERE n.frequency > 0.5 RETURN count(n)")
        ]
        
        performance_results = []
        
        for query_name, query in test_queries:
            start_time = time.time()
            result = session.run(query)
            list(result)  # Force execution
            end_time = time.time()
            
            performance_results.append({
                'Query': query_name,
                'Time (seconds)': round(end_time - start_time, 3)
            })
        
        return performance_results

# Test query performance
print("\n=== Database Performance Analysis ===")
perf_results = analyze_query_performance()

perf_df = pd.DataFrame(perf_results)
print(perf_df.to_string(index=False))

# Visualize performance
plt.figure(figsize=(10, 6))
plt.bar(perf_df['Query'], perf_df['Time (seconds)'])
plt.title('Query Performance Analysis')
plt.xlabel('Query Type')
plt.ylabel('Execution Time (seconds)')
plt.xticks(rotation=45)
plt.tight_layout()
plt.show()

# %% [markdown]
"""
## Summary and Next Steps

This notebook demonstrated the core capabilities of the fractal pangenome database:

### Key Insights
1. **Genomic Street Map**: Individual genomes are represented as navigable routes through a shared street network
2. **Multi-Scale Analysis**: The system efficiently handles queries from nucleotide to population level
3. **Route Diversity**: Different individuals take different paths through the same genomic regions
4. **Short Read Mapping**: Reads can be mapped to infer likely genomic paths
5. **Population Structure**: The system reveals population-level genomic diversity patterns

### Performance Characteristics
- Query times scale logarithmically with data size thanks to Hilbert indexing
- Complex path-finding queries complete in seconds
- Multi-resolution queries automatically optimize for available computational budget

### Next Steps for Analysis
1. **Import more genomes** to build a comprehensive pangenome
2. **Analyze specific genes** or disease-associated regions
3. **Compare populations** to identify ancestry-specific variants
4. **Integrate with phenotype data** for association studies
5. **Develop custom visualization** tools for your specific use cases

### Custom Analysis Templates
Use this notebook as a template for your own genomic analyses. The modular design allows you to:
- Modify queries for your specific research questions
- Integrate with existing genomics pipelines
- Scale analysis to larger datasets
- Develop custom visualization tools

Happy genomic exploration! 🧬🗺️
"""

print("=== Analysis Complete ===")
print("Notebook execution finished successfully!")
print("\nTo continue exploring:")
print("1. Modify the region coordinates to analyze different genomic areas")
print("2. Import your own genome assemblies")
print("3. Try the command-line interface: python -m genome_cli --help")
print("4. Explore the Neo4j browser at http://localhost:7474")
