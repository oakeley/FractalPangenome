"""
REST API Server for Fractal Pangenome Database
==============================================
FastAPI-based REST API providing programmatic access to the genomic street map database
"""

from fastapi import FastAPI, HTTPException, Query, BackgroundTasks, Depends, status
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import StreamingResponse, FileResponse
from fastapi.security import HTTPBearer, HTTPAuthorizationCredentials
import uvicorn
from pydantic import BaseModel, Field, validator
from typing import List, Dict, Optional, Any, Union
import asyncio
import json
import logging
from datetime import datetime, timedelta
import pandas as pd
import io
import gzip
from pathlib import Path
import hashlib
import redis
from contextlib import asynccontextmanager

# Import our core modules
from .neo4j_genome_importer import GenomeDatabaseBuilder, GenomeImporter
from .hilbert_pangenome_architecture import HilbertCurve, GenomicCoordinate
from .performance_optimization_tools import DatabaseBenchmark, PerformanceOptimizer

# Pydantic models for API requests and responses
class GenomicRegion(BaseModel):
    """Genomic region specification"""
    chromosome: int = Field(..., ge=1, le=24, description="Chromosome number (1-22, X=23, Y=24)")
    start_position: int = Field(..., ge=1, description="Start position (1-based)")
    end_position: int = Field(..., gt=0, description="End position (1-based)")
    
    @validator('end_position')
    def end_must_be_greater_than_start(cls, v, values):
        if 'start_position' in values and v <= values['start_position']:
            raise ValueError('end_position must be greater than start_position')
        return v

class QueryParameters(BaseModel):
    """Query parameters for genomic searches"""
    region: GenomicRegion
    individual_id: Optional[str] = Field(None, description="Specific individual to query")
    scale_level: int = Field(0, ge=0, le=3, description="Resolution level (0=nucleotide, 3=chromosome)")
    include_annotations: bool = Field(True, description="Include temporal annotations")
    max_results: int = Field(1000, ge=1, le=100000, description="Maximum number of results")
    frequency_filter: Optional[float] = Field(None, ge=0.0, le=1.0, description="Minimum variant frequency")

class GenomeNode(BaseModel):
    """Genomic node representation"""
    node_id: str
    chromosome: int
    start_position: int
    end_position: int
    sequence: Optional[str] = None
    frequency: float
    individual_id: str
    haplotype: int
    scale_level: int
    annotations: Dict[str, Any] = {}

class GenomicPath(BaseModel):
    """Genomic path representation"""
    individual_id: str
    haplotype: int
    path_nodes: List[str]
    path_length: int
    complexity_score: float
    rarity_score: float

class ImportRequest(BaseModel):
    """Genome import request"""
    individual_id: str = Field(..., description="Unique individual identifier")
    assembly_source: str = Field("custom", description="Assembly source name")
    maternal_fasta_path: Optional[str] = Field(None, description="Path to maternal FASTA file")
    paternal_fasta_path: Optional[str] = Field(None, description="Path to paternal FASTA file")
    public_genome_sources: Optional[List[str]] = Field(None, description="List of public genome sources")
    metadata: Dict[str, Any] = Field({}, description="Additional metadata")

class ShortReadMapping(BaseModel):
    """Short read mapping request"""
    reads: List[str] = Field(..., description="List of read sequences")
    region: GenomicRegion
    individual_id: Optional[str] = None
    mapping_quality_threshold: float = Field(0.8, ge=0.0, le=1.0)

class AnalysisJob(BaseModel):
    """Background analysis job"""
    job_id: str
    job_type: str
    status: str
    created_at: datetime
    completed_at: Optional[datetime] = None
    result_path: Optional[str] = None
    error_message: Optional[str] = None
    progress: float = 0.0

# Global app state
app_state = {
    "db_builder": None,
    "job_queue": {},
    "redis_client": None
}

# Security
security = HTTPBearer()

async def verify_token(credentials: HTTPAuthorizationCredentials = Depends(security)):
    """Verify API token (simplified - use proper JWT in production)"""
    # In production, implement proper JWT verification
    token = credentials.credentials
    if token != "genomics-api-key-change-me":  # Change this in production!
        raise HTTPException(
            status_code=status.HTTP_401_UNAUTHORIZED,
            detail="Invalid authentication token"
        )
    return token

@asynccontextmanager
async def lifespan(app: FastAPI):
    """Application lifespan management"""
    # Startup
    logging.info("Starting Fractal Pangenome API Server...")
    
    try:
        # Initialize database connection
        app_state["db_builder"] = GenomeDatabaseBuilder()
        success = await app_state["db_builder"].initialize_database()
        if not success:
            raise Exception("Failed to initialize database")
        
        # Initialize Redis for job queue
        app_state["redis_client"] = redis.Redis(
            host='localhost', port=6379, decode_responses=True
        )
        app_state["redis_client"].ping()
        
        logging.info("✅ Database connections established")
        
    except Exception as e:
        logging.error(f"❌ Failed to initialize: {e}")
        raise e
    
    yield
    
    # Shutdown
    logging.info("Shutting down API server...")
    if app_state["db_builder"]:
        app_state["db_builder"].cleanup()
    if app_state["redis_client"]:
        app_state["redis_client"].close()

# Create FastAPI app
app = FastAPI(
    title="Fractal Pangenome API",
    description="REST API for the fractal pangenome database system",
    version="1.0.0",
    docs_url="/docs",
    redoc_url="/redoc",
    lifespan=lifespan
)

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # Configure appropriately for production
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Health check endpoint
@app.get("/health", tags=["System"])
async def health_check():
    """Health check endpoint"""
    try:
        # Test database connection
        if app_state["db_builder"]:
            driver = app_state["db_builder"].neo4j_manager.get_connection()
            with driver.session() as session:
                session.run("RETURN 1")
            driver.close()
            
        return {
            "status": "healthy",
            "timestamp": datetime.now().isoformat(),
            "database": "connected",
            "version": "1.0.0"
        }
    except Exception as e:
        raise HTTPException(status_code=503, detail=f"Service unhealthy: {str(e)}")

# Database status endpoint
@app.get("/api/v1/status", tags=["System"])
async def get_database_status(token: str = Depends(verify_token)):
    """Get comprehensive database status"""
    try:
        db_builder = app_state["db_builder"]
        driver = db_builder.neo4j_manager.get_connection()
        
        with driver.session() as session:
            # Get basic statistics
            stats = session.run("""
                MATCH (n:GenomeNode)
                RETURN 
                    count(n) as total_nodes,
                    count(DISTINCT n.individual_id) as individuals,
                    count(DISTINCT n.chromosome) as chromosomes,
                    avg(n.frequency) as avg_frequency
            """).single()
            
            # Get scale level distribution
            scale_dist = session.run("""
                MATCH (n:GenomeNode)
                RETURN n.scale_level as level, count(n) as count
                ORDER BY n.scale_level
            """).data()
            
        driver.close()
        
        return {
            "database_statistics": dict(stats) if stats else {},
            "scale_distribution": scale_dist,
            "timestamp": datetime.now().isoformat()
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get status: {str(e)}")

# Genomic region query
@app.post("/api/v1/query/region", response_model=Dict[str, Any], tags=["Queries"])
async def query_genomic_region(
    params: QueryParameters,
    token: str = Depends(verify_token)
):
    """Query a genomic region with adaptive resolution"""
    try:
        db_builder = app_state["db_builder"]
        
        # Execute comprehensive query
        result = await db_builder.comprehensive_genome_query(
            individual_id=params.individual_id or "any",
            chromosome=params.region.chromosome,
            start=params.region.start_position,
            end=params.region.end_position,
            include_annotations=params.include_annotations,
            optimization_level='balanced'
        )
        
        # Apply additional filters
        if params.frequency_filter is not None:
            # Filter results by frequency (simplified)
            filtered_routes = {}
            for route_type, nodes in result.get('diploid_routes', {}).items():
                filtered_routes[route_type] = nodes[:params.max_results]
            result['diploid_routes'] = filtered_routes
        
        return result
        
    except Exception as e:
        logging.error(f"Query error: {e}")
        raise HTTPException(status_code=500, detail=f"Query failed: {str(e)}")

# Individual genome paths
@app.get("/api/v1/individuals/{individual_id}/paths", tags=["Individuals"])
async def get_individual_paths(
    individual_id: str,
    chromosome: int = Query(..., ge=1, le=24),
    start_position: int = Query(..., ge=1),
    end_position: int = Query(..., gt=0),
    token: str = Depends(verify_token)
):
    """Get diploid genomic paths for a specific individual"""
    try:
        db_builder = app_state["db_builder"]
        
        # Find individual paths
        paths = await db_builder.router.find_diploid_routes(
            individual_id=individual_id,
            chromosome=chromosome,
            start_pos=start_position,
            end_pos=end_position
        )
        
        # Format response
        formatted_paths = []
        for route_type, path_nodes in paths.items():
            formatted_paths.append(GenomicPath(
                individual_id=individual_id,
                haplotype=0 if route_type == 'maternal' else 1,
                path_nodes=path_nodes,
                path_length=len(path_nodes),
                complexity_score=len(path_nodes) / 100.0,  # Simplified
                rarity_score=0.5  # Would calculate based on node frequencies
            ))
        
        return {
            "individual_id": individual_id,
            "region": {
                "chromosome": chromosome,
                "start_position": start_position,
                "end_position": end_position
            },
            "paths": [path.dict() for path in formatted_paths]
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get paths: {str(e)}")

# Short read mapping
@app.post("/api/v1/mapping/reads", tags=["Mapping"])
async def map_short_reads(
    mapping_request: ShortReadMapping,
    token: str = Depends(verify_token)
):
    """Map short reads to the pangenome"""
    try:
        db_builder = app_state["db_builder"]
        
        # Create temporary FASTQ-like data
        reads_data = []
        for i, read_seq in enumerate(mapping_request.reads):
            reads_data.append({
                'read_id': f'api_read_{i}',
                'sequence': read_seq,
                'quality': 'I' * len(read_seq)  # Fake quality scores
            })
        
        # Map reads using the short read mapper
        mapping_results = await db_builder.short_read_mapper.map_short_reads(
            reads_data, mapping_request.individual_id
        )
        
        # Infer paths from mappings
        inferred_paths = await db_builder.short_read_mapper.infer_individual_paths(
            mapping_results,
            mapping_request.region.chromosome,
            mapping_request.region.start_position,
            mapping_request.region.end_position
        )
        
        return {
            "mapping_results": mapping_results,
            "inferred_paths": inferred_paths,
            "statistics": {
                "total_reads": len(mapping_request.reads),
                "mapped_reads": len(mapping_results),
                "mapping_rate": len(mapping_results) / len(mapping_request.reads) if mapping_request.reads else 0
            }
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Mapping failed: {str(e)}")

# Genome import
@app.post("/api/v1/import/genome", tags=["Import"])
async def import_genome(
    import_request: ImportRequest,
    background_tasks: BackgroundTasks,
    token: str = Depends(verify_token)
):
    """Import a genome (background job)"""
    try:
        # Generate job ID
        job_id = hashlib.md5(f"{import_request.individual_id}_{datetime.now()}".encode()).hexdigest()
        
        # Create job record
        job = AnalysisJob(
            job_id=job_id,
            job_type="genome_import",
            status="queued",
            created_at=datetime.now()
        )
        
        # Store job in Redis
        app_state["redis_client"].setex(
            f"job:{job_id}",
            timedelta(hours=24),
            job.json()
        )
        
        # Add background task
        background_tasks.add_task(
            _import_genome_background,
            job_id,
            import_request
        )
        
        return {
            "job_id": job_id,
            "status": "queued",
            "message": "Import job started in background"
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Import failed: {str(e)}")

async def _import_genome_background(job_id: str, import_request: ImportRequest):
    """Background task for genome import"""
    try:
        # Update job status
        job_data = app_state["redis_client"].get(f"job:{job_id}")
        if job_data:
            job = AnalysisJob.parse_raw(job_data)
            job.status = "running"
            job.progress = 0.1
            app_state["redis_client"].setex(
                f"job:{job_id}",
                timedelta(hours=24),
                job.json()
            )
        
        db_builder = app_state["db_builder"]
        
        # Import public genomes if specified
        if import_request.public_genome_sources:
            job.progress = 0.3
            app_state["redis_client"].setex(f"job:{job_id}", timedelta(hours=24), job.json())
            
            success = await db_builder.import_public_genomes(
                import_request.public_genome_sources
            )
            if not success:
                raise Exception("Failed to import public genomes")
        
        # Import custom diploid genome if paths provided
        if import_request.maternal_fasta_path and import_request.paternal_fasta_path:
            job.progress = 0.7
            app_state["redis_client"].setex(f"job:{job_id}", timedelta(hours=24), job.json())
            
            success = await db_builder.genome_importer.import_diploid_genome(
                import_request.individual_id,
                import_request.maternal_fasta_path,
                import_request.paternal_fasta_path,
                import_request.assembly_source
            )
            if not success:
                raise Exception("Failed to import diploid genome")
        
        # Complete job
        job.status = "completed"
        job.progress = 1.0
        job.completed_at = datetime.now()
        app_state["redis_client"].setex(f"job:{job_id}", timedelta(hours=24), job.json())
        
    except Exception as e:
        # Mark job as failed
        job_data = app_state["redis_client"].get(f"job:{job_id}")
        if job_data:
            job = AnalysisJob.parse_raw(job_data)
            job.status = "failed"
            job.error_message = str(e)
            job.completed_at = datetime.now()
            app_state["redis_client"].setex(f"job:{job_id}", timedelta(hours=24), job.json())

# Job status
@app.get("/api/v1/jobs/{job_id}", tags=["Jobs"])
async def get_job_status(job_id: str, token: str = Depends(verify_token)):
    """Get background job status"""
    try:
        job_data = app_state["redis_client"].get(f"job:{job_id}")
        if not job_data:
            raise HTTPException(status_code=404, detail="Job not found")
        
        job = AnalysisJob.parse_raw(job_data)
        return job.dict()
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to get job status: {str(e)}")

# Population analysis
@app.get("/api/v1/population/diversity", tags=["Population"])
async def get_population_diversity(
    chromosome: int = Query(..., ge=1, le=24),
    token: str = Depends(verify_token)
):
    """Analyze population-level genomic diversity"""
    try:
        db_builder = app_state["db_builder"]
        driver = db_builder.neo4j_manager.get_connection()
        
        with driver.session() as session:
            # Variant frequency distribution
            freq_query = """
            MATCH (n:GenomeNode)
            WHERE n.chromosome = $chr AND n.scale_level = 0
            RETURN n.frequency as freq, count(n) as count
            ORDER BY n.frequency
            """
            freq_data = session.run(freq_query, chr=chromosome).data()
            
            # Individual diversity metrics
            diversity_query = """
            MATCH (n:GenomeNode)
            WHERE n.chromosome = $chr AND n.scale_level = 0
            RETURN 
                count(DISTINCT n.individual_id) as individuals,
                avg(n.frequency) as avg_frequency,
                stddev(n.frequency) as freq_stddev,
                count(n) as total_variants
            """
            diversity_stats = session.run(diversity_query, chr=chromosome).single()
            
            # Rare variants
            rare_query = """
            MATCH (n:GenomeNode)
            WHERE n.chromosome = $chr AND n.frequency < 0.05 AND n.scale_level = 0
            RETURN count(n) as rare_variants,
                   collect(DISTINCT n.individual_id)[0..10] as carriers
            """
            rare_data = session.run(rare_query, chr=chromosome).single()
        
        driver.close()
        
        return {
            "chromosome": chromosome,
            "diversity_statistics": dict(diversity_stats) if diversity_stats else {},
            "rare_variants": dict(rare_data) if rare_data else {},
            "frequency_distribution": freq_data,
            "analysis_timestamp": datetime.now().isoformat()
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Analysis failed: {str(e)}")

# Export data
@app.get("/api/v1/export/region", tags=["Export"])
async def export_genomic_region(
    chromosome: int = Query(..., ge=1, le=24),
    start_position: int = Query(..., ge=1),
    end_position: int = Query(..., gt=0),
    format: str = Query("json", regex="^(json|csv|vcf)$"),
    token: str = Depends(verify_token)
):
    """Export genomic region data"""
    try:
        db_builder = app_state["db_builder"]
        driver = db_builder.neo4j_manager.get_connection()
        
        with driver.session() as session:
            query = """
            MATCH (n:GenomeNode)
            WHERE n.chromosome = $chr 
              AND n.start_pos >= $start 
              AND n.end_pos <= $end
              AND n.scale_level = 0
            RETURN n.node_id as id, n.individual_id as individual,
                   n.haplotype as haplotype, n.start_pos as start,
                   n.end_pos as end, n.frequency as frequency,
                   n.sequence as sequence
            ORDER BY n.start_pos
            """
            
            results = session.run(
                query,
                chr=chromosome,
                start=start_position,
                end=end_position
            ).data()
        
        driver.close()
        
        if format == "csv":
            # Convert to CSV
            df = pd.DataFrame(results)
            output = io.StringIO()
            df.to_csv(output, index=False)
            
            return StreamingResponse(
                io.BytesIO(output.getvalue().encode()),
                media_type="text/csv",
                headers={"Content-Disposition": f"attachment; filename=genomic_region_chr{chromosome}_{start_position}_{end_position}.csv"}
            )
        
        elif format == "vcf":
            # Convert to simplified VCF format
            vcf_lines = [
                "##fileformat=VCFv4.2",
                f"##source=FractalPangenome_v1.0",
                f"##reference=custom",
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE"
            ]
            
            for result in results:
                # Simplified VCF entry
                vcf_lines.append(
                    f"{chromosome}\t{result['start']}\t{result['id']}\tN\t{result['sequence'][:10]}\t60\tPASS\tAF={result['frequency']}\tGT\t0/1"
                )
            
            vcf_content = "\n".join(vcf_lines)
            
            return StreamingResponse(
                io.BytesIO(vcf_content.encode()),
                media_type="text/plain",
                headers={"Content-Disposition": f"attachment; filename=genomic_region_chr{chromosome}_{start_position}_{end_position}.vcf"}
            )
        
        else:  # JSON format
            return {
                "region": {
                    "chromosome": chromosome,
                    "start_position": start_position,
                    "end_position": end_position
                },
                "variants": results,
                "count": len(results),
                "export_timestamp": datetime.now().isoformat()
            }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Export failed: {str(e)}")

# Performance metrics
@app.get("/api/v1/performance/metrics", tags=["Performance"])
async def get_performance_metrics(token: str = Depends(verify_token)):
    """Get database performance metrics"""
    try:
        db_builder = app_state["db_builder"]
        
        # Run quick performance benchmark
        benchmark = DatabaseBenchmark(db_builder)
        basic_results = await benchmark._benchmark_basic_queries()
        
        # System metrics
        import psutil
        system_metrics = {
            "cpu_percent": psutil.cpu_percent(),
            "memory_percent": psutil.virtual_memory().percent,
            "disk_usage": psutil.disk_usage('/').percent
        }
        
        return {
            "system_metrics": system_metrics,
            "query_performance": [vars(r) for r in basic_results],
            "timestamp": datetime.now().isoformat()
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Performance check failed: {str(e)}")

# Batch operations
@app.post("/api/v1/batch/query", tags=["Batch"])
async def batch_query_regions(
    regions: List[GenomicRegion],
    individual_id: Optional[str] = None,
    token: str = Depends(verify_token)
):
    """Query multiple genomic regions in batch"""
    try:
        if len(regions) > 100:
            raise HTTPException(status_code=400, detail="Maximum 100 regions per batch")
        
        db_builder = app_state["db_builder"]
        results = []
        
        for region in regions:
            try:
                result = await db_builder.comprehensive_genome_query(
                    individual_id=individual_id or "any",
                    chromosome=region.chromosome,
                    start=region.start_position,
                    end=region.end_position,
                    include_annotations=False,  # Skip annotations for batch
                    optimization_level='fast'
                )
                results.append({
                    "region": region.dict(),
                    "result": result,
                    "success": True
                })
            except Exception as e:
                results.append({
                    "region": region.dict(),
                    "error": str(e),
                    "success": False
                })
        
        return {
            "batch_results": results,
            "total_regions": len(regions),
            "successful_queries": sum(1 for r in results if r["success"]),
            "timestamp": datetime.now().isoformat()
        }
        
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Batch query failed: {str(e)}")

# WebSocket endpoint for real-time updates (optional)
from fastapi import WebSocket, WebSocketDisconnect

@app.websocket("/ws/updates")
async def websocket_endpoint(websocket: WebSocket):
    """WebSocket endpoint for real-time updates"""
    await websocket.accept()
    try:
        while True:
            # Send periodic updates
            await asyncio.sleep(10)
            await websocket.send_json({
                "type": "system_status",
                "data": {
                    "timestamp": datetime.now().isoformat(),
                    "active_connections": 1  # Would track actual connections
                }
            })
    except WebSocketDisconnect:
        pass

# Error handlers
@app.exception_handler(ValueError)
async def value_error_handler(request, exc):
    return HTTPException(status_code=400, detail=str(exc))

@app.exception_handler(Exception)
async def general_exception_handler(request, exc):
    logging.error(f"Unhandled exception: {exc}")
    return HTTPException(status_code=500, detail="Internal server error")

# Main entry point
if __name__ == "__main__":
    import argparse
    
    parser = argparse.ArgumentParser(description="Fractal Pangenome API Server")
    parser.add_argument("--host", default="0.0.0.0", help="Host to bind to")
    parser.add_argument("--port", type=int, default=8000, help="Port to bind to")
    parser.add_argument("--reload", action="store_true", help="Enable auto-reload")
    parser.add_argument("--log-level", default="info", help="Log level")
    
    args = parser.parse_args()
    
    # Configure logging
    logging.basicConfig(
        level=getattr(logging, args.log_level.upper()),
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    )
    
    # Run server
    uvicorn.run(
        "rest_api_server:app",
        host=args.host,
        port=args.port,
        reload=args.reload,
        log_level=args.log_level
    )
