# Fractal Pangenome Database - Quick Start Guide

## Overview

This system transforms genomic data into a navigable "street map" where individual genomes are represented as routes through a pangenome graph. It uses Hilbert curve indexing for efficient multi-scale queries and supports both reference genomes and custom assemblies.

## Key Features

- 🗺️ **Street Map Analogy**: Individual genomes as routes, variants as alternative paths
- 📏 **Multi-Scale Querying**: Zoom from nucleotides to populations seamlessly  
- 🧬 **Diploid Support**: Separate maternal/paternal routes with structural variants
- 🔍 **Short Read Mapping**: Map Illumina reads to infer genomic paths
- 📊 **Fractal Architecture**: Efficient storage and querying at any resolution
- 🐳 **Docker Integration**: Complete containerized deployment

## Quick Installation

### Prerequisites

- Docker and Docker Compose
- Python 3.9+
- 8GB+ RAM recommended

### 1. Clone and Setup

```bash
# Clone the repository
git clone https://github.com/genomics/fractal-pangenome-db
cd fractal-pangenome-db

# Install Python dependencies
pip install -r requirements.txt

# Copy environment template
cp .env.example .env
```

### 2. Start the Database

```bash
# Option A: Start full stack (recommended)
python -m genome_cli start --docker-compose

# Option B: Start just Neo4j
python -m genome_cli database init
```

This will start:
- **Neo4j** (graph database) on port 7474/7687
- **Redis** (caching) on port 6379  
- **InfluxDB** (temporal annotations) on port 8086
- **Grafana** (monitoring) on port 3000
- **Jupyter** (analysis) on port 8888

### 3. Import Reference Genomes

```bash
# List available public genomes
python -m genome_cli import public --list-sources

# Import reference genomes
python -m genome_cli import public \
  --source T2T-CHM13 \
  --source HG002_maternal \
  --source HG002_paternal
```

### 4. Import Your Custom Genomes

```bash
# Import diploid genome from FASTA files
python -m genome_cli import diploid \
  --individual-id "patient_001" \
  --maternal /path/to/maternal.fa \
  --paternal /path/to/paternal.fa \
  --assembly-source "custom_assembly"
```

## Basic Usage Examples

### Query a Genomic Region

```bash
# Query a specific region
python -m genome_cli query region \
  --chromosome 1 \
  --start 1000000 \
  --end 2000000 \
  --individual patient_001 \
  --output results.json
```

### Map Short Reads

```bash
# Map Illumina reads to find genomic paths
python -m genome_cli map reads \
  --fastq /path/to/reads.fastq \
  --chromosome 1 \
  --start 1000000 \
  --end 2000000 \
  --individual patient_001 \
  --output mapping_results.json
```

### Database Statistics

```bash
# View database statistics
python -m genome_cli analyze stats

# Check database status
python -m genome_cli database status
```

## Understanding the Street Map Analogy

### Genomic "Streets" (Pangenome Nodes)
- **Reference paths**: Main highways (high frequency)
- **Common variants**: Side streets (medium frequency) 
- **Rare variants**: Back alleys (low frequency)
- **Structural variants**: Bridges, tunnels, detours

### Individual "Routes" (Genome Paths)
- **Maternal route**: "Route to work" - one path through the street network
- **Paternal route**: "Route home" - different path, may use one-way streets (inversions)
- **Both routes**: Must navigate the same start/end points but can take different paths

### Short Read "GPS" (Read Mapping)
- **Illumina reads**: GPS coordinates of where you are
- **K-mer matching**: Finding which streets contain your coordinates
- **Path inference**: Reconstructing the most likely route taken

## Advanced Features

### Multi-Resolution Queries

The system automatically selects the appropriate resolution level:

```python
# Small region → nucleotide level detail
result = query_region(chr=1, start=1000, end=2000)

# Large region → gene level aggregation  
result = query_region(chr=1, start=1000000, end=10000000)

# Whole chromosome → chromosome-level summary
result = query_region(chr=1, start=1, end=248000000)
```

### Adaptive Path Finding

```python
# Find optimal diploid routes with constraints
routes = find_diploid_routes(
    individual="patient_001",
    chromosome=1, 
    start=1000000, 
    end=2000000
)

print(routes['maternal'])  # High-frequency path
print(routes['paternal'])  # Alternative path with different variants
```

### Temporal Annotations

```python
# Annotations evolve over time (like shops changing purpose)
update_annotation(
    node_id="chr1_gene_BRCA1",
    annotation="function", 
    value="tumor_suppressor_updated_2024"
)
```

## API Usage (Python)

```python
import asyncio
from genome_processor import GenomeDatabaseBuilder

async def example_usage():
    # Initialize system
    db = GenomeDatabaseBuilder()
    await db.initialize_database()
    
    # Import a genome
    success = await db.genome_importer.import_diploid_genome(
        individual_id="test_patient",
        maternal_fasta="maternal.fa",
        paternal_fasta="paternal.fa"
    )
    
    # Query genomic region
    result = await db.comprehensive_genome_query(
        individual_id="test_patient",
        chromosome=1,
        start=1000000,
        end=2000000
    )
    
    print(f"Found {len(result['diploid_routes'])} paths")
    
    # Map short reads
    analysis = await db.analyze_short_reads(
        fastq_path="reads.fastq",
        individual_id="test_patient", 
        chromosome=1,
        start=1000000,
        end=2000000
    )
    
    print(f"Mapping rate: {analysis['mapping_rate']:.2%}")

# Run example
asyncio.run(example_usage())
```

## Configuration

Edit `genome_db_config.yaml` to customize:

```yaml
# Memory settings
neo4j:
  memory:
    heap_max: "16G"      # Increase for large datasets
    pagecache: "8G"

# Processing settings  
processing:
  segment_size: 5000     # Smaller segments = higher resolution
  kmer_size: 31          # K-mer size for read mapping
  batch_size: 2000       # Larger batches = faster imports

# Performance tuning
mapping:
  min_mapping_score: 0.9  # Stricter mapping threshold
  max_alignments_per_read: 5  # Fewer alignments = faster
```

## Performance Guidelines

### Small Dataset (< 10 individuals)
- Default settings work well
- 8GB RAM sufficient
- Import time: ~30 minutes per genome

### Medium Dataset (10-100 individuals) 
- Increase Neo4j heap to 16GB
- Use SSD storage
- Consider clustering common variants

### Large Dataset (100+ individuals)
- 32GB+ RAM recommended
- Distribute across multiple servers
- Use pre-computed k-mer indices

## Troubleshooting

### Common Issues

**Database won't start:**
```bash
# Check Docker status
docker ps

# View logs
docker logs neo4j-genomics

# Restart with fresh data
python -m genome_cli database stop
rm -rf neo4j_data
python -m genome_cli database init
```

**Import fails:**
```bash
# Check file formats
file genome.fa  # Should be: ASCII text

# Verify FASTA format
head -20 genome.fa

# Check available disk space
df -h
```

**Slow queries:**
```bash
# Check database statistics
python -m genome_cli analyze stats

# Monitor resource usage
docker stats neo4j-genomics
```

## Web Interfaces

After starting the system, access:

- **Neo4j Browser**: http://localhost:7474
  - Username: `neo4j` 
  - Password: `genomics123`
  - Run Cypher queries directly

- **Grafana Dashboard**: http://localhost:3000
  - Username: `admin`
  - Password: `genomics123` 
  - Monitor system performance

- **Jupyter Notebooks**: http://localhost:8888
  - Token: `genomics`
  - Interactive analysis environment

## Example Cypher Queries

```cypher
// Find all nodes for an individual
MATCH (n:GenomeNode {individual_id: 'patient_001'})
RETURN n.chromosome, count(n) as node_count
ORDER BY n.chromosome

// Find high-frequency variants
MATCH (n:GenomeNode)
WHERE n.frequency > 0.1
RETURN n.chromosome, n.start_pos, n.frequency
ORDER BY n.frequency DESC

// Find paths between two nodes
MATCH path = shortestPath(
  (start:GenomeNode {id: 'node1'})-[:CONNECTS*..10]->(end:GenomeNode {id: 'node2'})
)
RETURN path
```

## Next Steps

1. **Explore Data**: Use the Jupyter notebooks to analyze your imported genomes
2. **Custom Analysis**: Write Python scripts using the API
3. **Scale Up**: Import more diverse genomes to build your pangenome
4. **Optimize**: Tune configuration for your specific use case
5. **Extend**: Add custom annotations and analysis pipelines

## Support

- **Documentation**: See `/docs` directory for detailed API reference
- **Examples**: Check `/examples` for common use cases  
- **Issues**: Report bugs on GitHub
- **Community**: Join our Discord for discussions

## Citation

If you use this system in your research, please cite:

```bibtex
@software{fractal_pangenome_db,
  title={Fractal Pangenome Database: A Scalable System for Genomic Variation Analysis},
  author={Edward J. Oakeley, PhD},
  year={2025},
  url={[https://github.com/oakeley/FractalPangenome](https://github.com/oakeley/FractalPangenome)}
}
```
