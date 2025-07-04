"""
RNA-seq Expression Analysis for Fractal Pangenome Database
==========================================================
Comprehensive RNA-seq alignment and expression quantification system that:
- Aligns short reads to all genome variants in the fractal pangenome
- Calculates allele-specific expression values
- Links expression to specific pangenome variants
- Enables eQTL discovery and personalized expression profiling
"""

import asyncio
import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Set, Any
from dataclasses import dataclass, field
from pathlib import Path
import tempfile
import subprocess
import logging
from collections import defaultdict, Counter
import json
import gzip
from datetime import datetime
import hashlib
import pickle
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots

# Import our core modules
from neo4j_genome_importer import GenomeDatabaseBuilder
from hilbert_pangenome_architecture import HilbertCurve, GenomicCoordinate

@dataclass
class RNASeqRead:
    """RNA-seq short read representation"""
    read_id: str
    sequence: str
    quality: str
    mate_sequence: Optional[str] = None
    mate_quality: Optional[str] = None
    is_paired: bool = False

@dataclass
class TranscriptAlignment:
    """Alignment of read to transcript variant"""
    read_id: str
    transcript_variant_id: str
    gene_id: str
    allele_id: str
    alignment_score: float
    mapping_quality: int
    start_pos: int
    end_pos: int
    is_unique: bool
    mismatch_count: int
    
@dataclass
class ExpressionResult:
    """Gene expression quantification result"""
    gene_id: str
    transcript_variants: Dict[str, float]  # variant_id -> expression
    total_expression: float
    allele_specific_expression: Dict[str, float]  # allele_id -> expression
    confidence_interval: Tuple[float, float]
    fpkm: float
    tpm: float
    raw_counts: int
    effective_length: float

@dataclass
class eQTLResult:
    """Expression quantitative trait loci result"""
    variant_id: str
    gene_id: str
    chromosome: int
    variant_position: int
    gene_position: int
    distance_to_tss: int
    beta: float
    pvalue: float
    qvalue: float
    r_squared: float
    allele_frequencies: Dict[str, float]
    effect_direction: str  # "positive" or "negative"

class PangenomeTranscriptBuilder:
    """Build transcript references from pangenome variants"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        self.transcript_cache = {}
        
    async def build_transcript_variants(self, gene_regions: List[Dict[str, Any]]) -> Dict[str, str]:
        """Build all possible transcript variants from pangenome for given genes"""
        
        transcript_variants = {}
        
        with self.driver.session() as session:
            for gene_region in gene_regions:
                gene_id = gene_region['gene_id']
                chromosome = gene_region['chromosome']
                start = gene_region['start']
                end = gene_region['end']
                strand = gene_region.get('strand', '+')
                
                # Get all genomic variants in this gene region
                variant_query = """
                MATCH (n:GenomeNode)
                WHERE n.chromosome = $chromosome
                  AND n.start_pos >= $start
                  AND n.end_pos <= $end
                  AND n.scale_level = 0
                RETURN n.node_id as variant_id, n.individual_id as individual,
                       n.haplotype as haplotype, n.start_pos as pos,
                       n.end_pos as end_pos, n.sequence as sequence,
                       n.frequency as frequency, n.variant_type as var_type
                ORDER BY n.start_pos, n.individual_id, n.haplotype
                """
                
                variants = session.run(variant_query, 
                                     chromosome=chromosome, start=start, end=end).data()
                
                # Group variants by individual and haplotype
                individual_variants = defaultdict(lambda: defaultdict(list))
                for variant in variants:
                    individual_variants[variant['individual']][variant['haplotype']].append(variant)
                
                # Build transcript sequences for each individual/haplotype combination
                for individual, haplotypes in individual_variants.items():
                    for haplotype, hap_variants in haplotypes.items():
                        transcript_id = f"{gene_id}_{individual}_h{haplotype}"
                        
                        # Construct transcript sequence from variants
                        transcript_seq = await self._construct_transcript_sequence(
                            hap_variants, gene_region, strand
                        )
                        
                        if transcript_seq:
                            transcript_variants[transcript_id] = transcript_seq
                            
                            # Store in database for future reference
                            session.run("""
                                CREATE (t:TranscriptVariant {
                                    transcript_id: $transcript_id,
                                    gene_id: $gene_id,
                                    individual_id: $individual,
                                    haplotype: $haplotype,
                                    sequence: $sequence,
                                    length: $length,
                                    chromosome: $chromosome,
                                    created_at: datetime()
                                })
                            """, transcript_id=transcript_id, gene_id=gene_id,
                                individual=individual, haplotype=haplotype,
                                sequence=transcript_seq, length=len(transcript_seq),
                                chromosome=chromosome)
        
        logging.info(f"Built {len(transcript_variants)} transcript variants for {len(gene_regions)} genes")
        return transcript_variants
    
    async def _construct_transcript_sequence(self, variants: List[Dict], 
                                           gene_region: Dict, strand: str) -> str:
        """Construct transcript sequence from genomic variants"""
        
        # Sort variants by position
        variants.sort(key=lambda x: x['pos'])
        
        # Start with reference sequence or reconstruct from variants
        transcript_sequence = ""
        last_pos = gene_region['start']
        
        for variant in variants:
            # Add gap sequence (would normally come from reference)
            gap_size = variant['pos'] - last_pos
            if gap_size > 0:
                # In a real implementation, this would come from reference genome
                gap_seq = 'N' * min(gap_size, 1000)  # Limit gap size for demo
                transcript_sequence += gap_seq
            
            # Add variant sequence
            transcript_sequence += variant['sequence']
            last_pos = variant['end_pos']
        
        # Reverse complement if on negative strand
        if strand == '-':
            transcript_sequence = self._reverse_complement(transcript_sequence)
        
        return transcript_sequence
    
    def _reverse_complement(self, sequence: str) -> str:
        """Generate reverse complement of DNA sequence"""
        complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
        return ''.join(complement.get(base, base) for base in reversed(sequence))

class RNASeqAligner:
    """RNA-seq read alignment to pangenome transcript variants"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        self.transcript_builder = PangenomeTranscriptBuilder(db_builder)
        self.alignment_cache = {}
        
    async def align_reads_to_transcripts(self, 
                                       reads: List[RNASeqRead],
                                       gene_regions: List[Dict[str, Any]],
                                       alignment_params: Dict[str, Any] = None) -> List[TranscriptAlignment]:
        """Align RNA-seq reads to all transcript variants in the pangenome"""
        
        if alignment_params is None:
            alignment_params = {
                'min_alignment_score': 0.8,
                'max_mismatches': 3,
                'allow_multimapping': True,
                'max_alignments_per_read': 10
            }
        
        # Build transcript variants for target genes
        transcript_variants = await self.transcript_builder.build_transcript_variants(gene_regions)
        
        if not transcript_variants:
            logging.warning("No transcript variants found for alignment")
            return []
        
        # Create k-mer index for fast alignment
        kmer_index = self._build_transcript_kmer_index(transcript_variants)
        
        alignments = []
        
        for read in reads:
            # Find candidate transcripts using k-mer matching
            candidates = self._find_candidate_transcripts(read.sequence, kmer_index, transcript_variants)
            
            # Perform detailed alignment for each candidate
            read_alignments = []
            for transcript_id, transcript_seq in candidates.items():
                alignment = self._align_read_to_transcript(
                    read, transcript_id, transcript_seq, alignment_params
                )
                if alignment and alignment.alignment_score >= alignment_params['min_alignment_score']:
                    read_alignments.append(alignment)
            
            # Filter and rank alignments
            read_alignments.sort(key=lambda x: x.alignment_score, reverse=True)
            
            # Keep top alignments
            max_alignments = alignment_params['max_alignments_per_read']
            read_alignments = read_alignments[:max_alignments]
            
            # Mark unique vs multi-mapping reads
            if len(read_alignments) == 1:
                read_alignments[0].is_unique = True
            
            alignments.extend(read_alignments)
        
        logging.info(f"Generated {len(alignments)} alignments for {len(reads)} reads")
        return alignments
    
    def _build_transcript_kmer_index(self, transcript_variants: Dict[str, str], k: int = 31) -> Dict[str, Set[str]]:
        """Build k-mer index for fast transcript matching"""
        
        kmer_index = defaultdict(set)
        
        for transcript_id, sequence in transcript_variants.items():
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                if 'N' not in kmer:  # Skip k-mers with ambiguous bases
                    kmer_index[kmer].add(transcript_id)
        
        return kmer_index
    
    def _find_candidate_transcripts(self, read_sequence: str, 
                                  kmer_index: Dict[str, Set[str]], 
                                  transcript_variants: Dict[str, str],
                                  k: int = 31) -> Dict[str, str]:
        """Find candidate transcripts for read using k-mer matching"""
        
        transcript_scores = defaultdict(int)
        
        # Extract k-mers from read
        for i in range(len(read_sequence) - k + 1):
            kmer = read_sequence[i:i+k]
            if kmer in kmer_index:
                for transcript_id in kmer_index[kmer]:
                    transcript_scores[transcript_id] += 1
        
        # Filter candidates with sufficient k-mer matches
        min_kmers = max(3, (len(read_sequence) - k + 1) // 10)  # At least 10% k-mer match
        candidates = {
            transcript_id: transcript_variants[transcript_id]
            for transcript_id, score in transcript_scores.items()
            if score >= min_kmers
        }
        
        return candidates
    
    def _align_read_to_transcript(self, read: RNASeqRead, transcript_id: str,
                                transcript_seq: str, params: Dict[str, Any]) -> Optional[TranscriptAlignment]:
        """Perform detailed alignment of read to transcript"""
        
        # Simplified alignment - in production, use proper alignment algorithms
        best_score = 0
        best_pos = 0
        best_mismatches = float('inf')
        
        read_seq = read.sequence
        read_len = len(read_seq)
        
        # Sliding window alignment
        for start_pos in range(len(transcript_seq) - read_len + 1):
            transcript_window = transcript_seq[start_pos:start_pos + read_len]
            
            # Calculate alignment score
            matches = sum(1 for r, t in zip(read_seq, transcript_window) if r == t)
            mismatches = read_len - matches
            score = matches / read_len
            
            if score > best_score and mismatches <= params['max_mismatches']:
                best_score = score
                best_pos = start_pos
                best_mismatches = mismatches
        
        if best_score >= params['min_alignment_score']:
            # Parse gene and allele info from transcript ID
            parts = transcript_id.split('_')
            gene_id = parts[0]
            individual_id = '_'.join(parts[1:-1])
            haplotype = parts[-1]
            allele_id = f"{individual_id}_{haplotype}"
            
            return TranscriptAlignment(
                read_id=read.read_id,
                transcript_variant_id=transcript_id,
                gene_id=gene_id,
                allele_id=allele_id,
                alignment_score=best_score,
                mapping_quality=int(best_score * 60),  # Convert to MAPQ-like score
                start_pos=best_pos,
                end_pos=best_pos + read_len,
                is_unique=False,  # Will be determined later
                mismatch_count=int(best_mismatches)
            )
        
        return None

class ExpressionQuantifier:
    """Quantify gene expression from RNA-seq alignments"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        
    async def quantify_expression(self, alignments: List[TranscriptAlignment],
                                transcript_lengths: Dict[str, int],
                                total_mapped_reads: int,
                                quantification_method: str = "em") -> Dict[str, ExpressionResult]:
        """Quantify gene expression with allele-specific resolution"""
        
        # Group alignments by gene
        gene_alignments = defaultdict(list)
        for alignment in alignments:
            gene_alignments[alignment.gene_id].append(alignment)
        
        expression_results = {}
        
        for gene_id, gene_alns in gene_alignments.items():
            if quantification_method == "em":
                result = await self._quantify_with_em(gene_id, gene_alns, transcript_lengths, total_mapped_reads)
            elif quantification_method == "unique_only":
                result = await self._quantify_unique_only(gene_id, gene_alns, transcript_lengths, total_mapped_reads)
            else:
                result = await self._quantify_proportional(gene_id, gene_alns, transcript_lengths, total_mapped_reads)
            
            expression_results[gene_id] = result
        
        # Store results in database
        await self._store_expression_results(expression_results)
        
        return expression_results
    
    async def _quantify_with_em(self, gene_id: str, alignments: List[TranscriptAlignment],
                              transcript_lengths: Dict[str, int], total_reads: int) -> ExpressionResult:
        """Quantify expression using Expectation-Maximization algorithm"""
        
        # Get unique transcript variants for this gene
        transcript_variants = list(set(aln.transcript_variant_id for aln in alignments))
        
        if not transcript_variants:
            return self._create_empty_result(gene_id)
        
        # Initialize transcript expression estimates
        transcript_expr = {tv: 1.0 for tv in transcript_variants}
        
        # EM algorithm iterations
        for iteration in range(50):  # Max 50 iterations
            old_expr = transcript_expr.copy()
            
            # E-step: Calculate expected counts
            expected_counts = defaultdict(float)
            
            for alignment in alignments:
                tv = alignment.transcript_variant_id
                
                # Find all alignments for this read
                read_alignments = [a for a in alignments if a.read_id == alignment.read_id]
                
                # Calculate probability of alignment to this transcript
                total_prob = sum(transcript_expr[a.transcript_variant_id] * a.alignment_score 
                               for a in read_alignments)
                
                if total_prob > 0:
                    prob = (transcript_expr[tv] * alignment.alignment_score) / total_prob
                    expected_counts[tv] += prob
            
            # M-step: Update expression estimates
            total_expected = sum(expected_counts.values())
            if total_expected > 0:
                for tv in transcript_variants:
                    transcript_expr[tv] = expected_counts[tv] / total_expected
            
            # Check convergence
            if self._check_convergence(old_expr, transcript_expr, threshold=1e-6):
                break
        
        # Calculate final expression metrics
        return self._calculate_expression_metrics(gene_id, transcript_expr, alignments, 
                                                transcript_lengths, total_reads)
    
    async def _quantify_unique_only(self, gene_id: str, alignments: List[TranscriptAlignment],
                                  transcript_lengths: Dict[str, int], total_reads: int) -> ExpressionResult:
        """Quantify expression using only uniquely mapped reads"""
        
        unique_alignments = [aln for aln in alignments if aln.is_unique]
        
        # Count reads per transcript variant
        transcript_counts = defaultdict(int)
        for alignment in unique_alignments:
            transcript_counts[alignment.transcript_variant_id] += 1
        
        # Convert counts to expression levels
        transcript_expr = {}
        total_counts = sum(transcript_counts.values())
        
        if total_counts > 0:
            for tv, count in transcript_counts.items():
                transcript_expr[tv] = count / total_counts
        
        return self._calculate_expression_metrics(gene_id, transcript_expr, unique_alignments,
                                                transcript_lengths, total_reads)
    
    async def _quantify_proportional(self, gene_id: str, alignments: List[TranscriptAlignment],
                                   transcript_lengths: Dict[str, int], total_reads: int) -> ExpressionResult:
        """Quantify expression by proportionally distributing multi-mapped reads"""
        
        # First, count unique alignments
        unique_counts = defaultdict(int)
        multi_alignments = []
        
        read_alignment_counts = defaultdict(int)
        for alignment in alignments:
            read_alignment_counts[alignment.read_id] += 1
        
        for alignment in alignments:
            if read_alignment_counts[alignment.read_id] == 1:
                unique_counts[alignment.transcript_variant_id] += 1
            else:
                multi_alignments.append(alignment)
        
        # Distribute multi-mapped reads proportionally
        total_unique = sum(unique_counts.values())
        
        transcript_expr = {}
        if total_unique > 0:
            # Calculate proportions from unique reads
            proportions = {tv: count / total_unique for tv, count in unique_counts.items()}
            
            # Distribute multi-mapped reads
            for alignment in multi_alignments:
                tv = alignment.transcript_variant_id
                read_id = alignment.read_id
                read_alignments = [a for a in multi_alignments if a.read_id == read_id]
                
                # Distribute read proportionally
                if tv in proportions:
                    transcript_expr[tv] = proportions[tv] + (proportions[tv] / len(read_alignments))
        
        return self._calculate_expression_metrics(gene_id, transcript_expr, alignments,
                                                transcript_lengths, total_reads)
    
    def _calculate_expression_metrics(self, gene_id: str, transcript_expr: Dict[str, float],
                                    alignments: List[TranscriptAlignment],
                                    transcript_lengths: Dict[str, int], total_reads: int) -> ExpressionResult:
        """Calculate comprehensive expression metrics"""
        
        # Total expression across all variants
        total_expression = sum(transcript_expr.values())
        
        # Raw read counts
        raw_counts = len(alignments)
        
        # Calculate effective length (average of transcript lengths weighted by expression)
        effective_length = 0
        if transcript_expr:
            total_expr = sum(transcript_expr.values())
            if total_expr > 0:
                for tv, expr in transcript_expr.items():
                    if tv in transcript_lengths:
                        effective_length += (expr / total_expr) * transcript_lengths[tv]
        
        if effective_length == 0:
            effective_length = 1000  # Default length
        
        # Calculate FPKM (Fragments Per Kilobase Million)
        fpkm = 0
        if total_reads > 0 and effective_length > 0:
            fpkm = (raw_counts * 1e9) / (effective_length * total_reads)
        
        # Calculate TPM (Transcripts Per Million) - simplified
        tpm = fpkm  # In practice, TPM requires normalization across all genes
        
        # Allele-specific expression
        allele_expression = defaultdict(float)
        for tv, expr in transcript_expr.items():
            # Extract allele info from transcript variant ID
            parts = tv.split('_')
            if len(parts) >= 3:
                allele_id = '_'.join(parts[-2:])  # individual_haplotype
                allele_expression[allele_id] += expr
        
        # Confidence interval (simplified)
        confidence_interval = (total_expression * 0.9, total_expression * 1.1)
        
        return ExpressionResult(
            gene_id=gene_id,
            transcript_variants=transcript_expr,
            total_expression=total_expression,
            allele_specific_expression=dict(allele_expression),
            confidence_interval=confidence_interval,
            fpkm=fpkm,
            tpm=tpm,
            raw_counts=raw_counts,
            effective_length=effective_length
        )
    
    def _create_empty_result(self, gene_id: str) -> ExpressionResult:
        """Create empty expression result for genes with no alignments"""
        return ExpressionResult(
            gene_id=gene_id,
            transcript_variants={},
            total_expression=0.0,
            allele_specific_expression={},
            confidence_interval=(0.0, 0.0),
            fpkm=0.0,
            tpm=0.0,
            raw_counts=0,
            effective_length=0.0
        )
    
    def _check_convergence(self, old_expr: Dict[str, float], new_expr: Dict[str, float], 
                          threshold: float = 1e-6) -> bool:
        """Check if EM algorithm has converged"""
        for tv in old_expr:
            if abs(old_expr[tv] - new_expr.get(tv, 0)) > threshold:
                return False
        return True
    
    async def _store_expression_results(self, results: Dict[str, ExpressionResult]):
        """Store expression results in the database"""
        
        with self.driver.session() as session:
            for gene_id, result in results.items():
                # Store gene-level expression
                session.run("""
                    MERGE (g:Gene {gene_id: $gene_id})
                    SET g.total_expression = $total_expr,
                        g.fpkm = $fpkm,
                        g.tpm = $tpm,
                        g.raw_counts = $raw_counts,
                        g.effective_length = $effective_length,
                        g.analysis_timestamp = datetime()
                """, gene_id=gene_id, total_expr=result.total_expression,
                    fpkm=result.fpkm, tpm=result.tpm, raw_counts=result.raw_counts,
                    effective_length=result.effective_length)
                
                # Store transcript variant expressions
                for tv, expr in result.transcript_variants.items():
                    session.run("""
                        MATCH (t:TranscriptVariant {transcript_id: $tv})
                        SET t.expression_level = $expr,
                            t.analysis_timestamp = datetime()
                    """, tv=tv, expr=expr)
                
                # Store allele-specific expression
                for allele_id, expr in result.allele_specific_expression.items():
                    session.run("""
                        MERGE (a:AlleleExpression {
                            gene_id: $gene_id,
                            allele_id: $allele_id
                        })
                        SET a.expression_level = $expr,
                            a.analysis_timestamp = datetime()
                    """, gene_id=gene_id, allele_id=allele_id, expr=expr)

class eQTLAnalyzer:
    """Expression Quantitative Trait Loci (eQTL) analysis"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        
    async def discover_eqtls(self, expression_results: Dict[str, ExpressionResult],
                           individuals: List[str], chromosome: int = None,
                           max_distance: int = 1000000) -> List[eQTLResult]:
        """Discover expression quantitative trait loci"""
        
        # Get variant data for individuals
        variant_data = await self._get_variant_matrix(individuals, chromosome)
        if not variant_data:
            return []
        
        # Get expression data
        expression_matrix = self._build_expression_matrix(expression_results, individuals)
        if expression_matrix.empty:
            return []
        
        eqtl_results = []
        
        # Test each gene against nearby variants
        for gene_id in expression_matrix.columns:
            gene_expression = expression_matrix[gene_id].values
            
            # Get gene position
            gene_pos = await self._get_gene_position(gene_id)
            if not gene_pos:
                continue
            
            # Test variants within max_distance
            for variant_id, variant_info in variant_data.items():
                var_pos = variant_info['position']
                distance = abs(var_pos - gene_pos['tss'])
                
                if distance <= max_distance:
                    # Get genotype data for this variant
                    genotypes = self._get_variant_genotypes(variant_id, individuals, variant_data)
                    
                    if len(set(genotypes)) > 1:  # Only test polymorphic variants
                        # Perform association test
                        eqtl_result = self._test_variant_gene_association(
                            variant_id, variant_info, gene_id, gene_pos,
                            genotypes, gene_expression, distance
                        )
                        
                        if eqtl_result and eqtl_result.pvalue < 0.05:
                            eqtl_results.append(eqtl_result)
        
        # Multiple testing correction
        if eqtl_results:
            eqtl_results = self._apply_multiple_testing_correction(eqtl_results)
        
        # Store results
        await self._store_eqtl_results(eqtl_results)
        
        return eqtl_results
    
    async def _get_variant_matrix(self, individuals: List[str], 
                                chromosome: int = None) -> Dict[str, Dict]:
        """Get variant data for individuals"""
        
        with self.driver.session() as session:
            query = """
            MATCH (n:GenomeNode)
            WHERE n.individual_id IN $individuals
              AND n.scale_level = 0
              AND n.variant_type IS NOT NULL
            """
            
            if chromosome:
                query += " AND n.chromosome = $chromosome"
            
            query += """
            RETURN n.node_id as variant_id, n.chromosome as chr,
                   n.start_pos as position, n.individual_id as individual,
                   n.haplotype as haplotype, n.frequency as frequency,
                   n.variant_type as var_type
            """
            
            params = {"individuals": individuals}
            if chromosome:
                params["chromosome"] = chromosome
            
            results = session.run(query, **params).data()
        
        # Organize variant data
        variant_data = {}
        for result in results:
            variant_id = result['variant_id']
            if variant_id not in variant_data:
                variant_data[variant_id] = {
                    'chromosome': result['chr'],
                    'position': result['position'],
                    'frequency': result['frequency'],
                    'variant_type': result['var_type'],
                    'genotypes': {}
                }
            
            # Store genotype for this individual
            individual = result['individual']
            haplotype = result['haplotype']
            if individual not in variant_data[variant_id]['genotypes']:
                variant_data[variant_id]['genotypes'][individual] = [0, 0]  # Diploid
            
            variant_data[variant_id]['genotypes'][individual][haplotype] = 1
        
        return variant_data
    
    def _build_expression_matrix(self, expression_results: Dict[str, ExpressionResult],
                               individuals: List[str]) -> pd.DataFrame:
        """Build expression matrix for individuals"""
        
        expression_data = []
        
        for individual in individuals:
            individual_expression = {}
            
            for gene_id, result in expression_results.items():
                # Sum expression across alleles for this individual
                total_expr = 0
                for allele_id, expr in result.allele_specific_expression.items():
                    if individual in allele_id:
                        total_expr += expr
                
                individual_expression[gene_id] = total_expr
            
            expression_data.append(individual_expression)
        
        return pd.DataFrame(expression_data, index=individuals)
    
    async def _get_gene_position(self, gene_id: str) -> Optional[Dict[str, int]]:
        """Get gene position information"""
        
        with self.driver.session() as session:
            result = session.run("""
                MATCH (g:Gene {gene_id: $gene_id})
                RETURN g.chromosome as chr, g.start_pos as start,
                       g.end_pos as end, g.tss as tss
            """, gene_id=gene_id).single()
            
            if result:
                return dict(result)
            
            # If not found, try to infer from transcript variants
            result = session.run("""
                MATCH (t:TranscriptVariant {gene_id: $gene_id})
                RETURN t.chromosome as chr, 
                       min(t.start_pos) as start,
                       max(t.end_pos) as end,
                       min(t.start_pos) as tss
            """, gene_id=gene_id).single()
            
            return dict(result) if result else None
    
    def _get_variant_genotypes(self, variant_id: str, individuals: List[str],
                             variant_data: Dict[str, Dict]) -> List[int]:
        """Get genotype vector for variant across individuals"""
        
        genotypes = []
        variant_info = variant_data[variant_id]
        
        for individual in individuals:
            if individual in variant_info['genotypes']:
                # Sum across haplotypes (0, 1, or 2 copies of variant)
                genotype = sum(variant_info['genotypes'][individual])
            else:
                genotype = 0  # No variant
            
            genotypes.append(genotype)
        
        return genotypes
    
    def _test_variant_gene_association(self, variant_id: str, variant_info: Dict,
                                     gene_id: str, gene_pos: Dict,
                                     genotypes: List[int], expression: List[float],
                                     distance: int) -> Optional[eQTLResult]:
        """Test association between variant and gene expression"""
        
        try:
            # Remove missing data
            valid_indices = [i for i, (g, e) in enumerate(zip(genotypes, expression))
                           if g is not None and e is not None and not np.isnan(e)]
            
            if len(valid_indices) < 10:  # Need sufficient samples
                return None
            
            clean_genotypes = [genotypes[i] for i in valid_indices]
            clean_expression = [expression[i] for i in valid_indices]
            
            # Perform linear regression
            slope, intercept, r_value, p_value, std_err = stats.linregress(clean_genotypes, clean_expression)
            
            # Calculate allele frequencies
            total_alleles = len(clean_genotypes) * 2
            alt_alleles = sum(clean_genotypes)
            allele_frequencies = {
                'ref': (total_alleles - alt_alleles) / total_alleles,
                'alt': alt_alleles / total_alleles
            }
            
            # Determine effect direction
            effect_direction = "positive" if slope > 0 else "negative"
            
            return eQTLResult(
                variant_id=variant_id,
                gene_id=gene_id,
                chromosome=variant_info['chromosome'],
                variant_position=variant_info['position'],
                gene_position=gene_pos['tss'],
                distance_to_tss=distance,
                beta=slope,
                pvalue=p_value,
                qvalue=p_value,  # Will be corrected later
                r_squared=r_value**2,
                allele_frequencies=allele_frequencies,
                effect_direction=effect_direction
            )
            
        except Exception as e:
            logging.warning(f"Failed to test association for {variant_id}-{gene_id}: {e}")
            return None
    
    def _apply_multiple_testing_correction(self, eqtl_results: List[eQTLResult]) -> List[eQTLResult]:
        """Apply Benjamini-Hochberg correction for multiple testing"""
        
        if not eqtl_results:
            return eqtl_results
        
        # Extract p-values
        pvalues = [result.pvalue for result in eqtl_results]
        
        # Benjamini-Hochberg correction
        from scipy.stats import false_discovery_control
        qvalues = false_discovery_control(pvalues, method='bh')
        
        # Update q-values
        for i, result in enumerate(eqtl_results):
            result.qvalue = qvalues[i]
        
        return eqtl_results
    
    async def _store_eqtl_results(self, eqtl_results: List[eQTLResult]):
        """Store eQTL results in database"""
        
        with self.driver.session() as session:
            for result in eqtl_results:
                session.run("""
                    CREATE (e:eQTL {
                        variant_id: $variant_id,
                        gene_id: $gene_id,
                        chromosome: $chromosome,
                        variant_position: $var_pos,
                        gene_position: $gene_pos,
                        distance_to_tss: $distance,
                        beta: $beta,
                        pvalue: $pvalue,
                        qvalue: $qvalue,
                        r_squared: $r_squared,
                        effect_direction: $effect_direction,
                        allele_frequencies: $allele_freq,
                        analysis_timestamp: datetime()
                    })
                """, 
                variant_id=result.variant_id,
                gene_id=result.gene_id,
                chromosome=result.chromosome,
                var_pos=result.variant_position,
                gene_pos=result.gene_position,
                distance=result.distance_to_tss,
                beta=result.beta,
                pvalue=result.pvalue,
                qvalue=result.qvalue,
                r_squared=result.r_squared,
                effect_direction=result.effect_direction,
                allele_freq=json.dumps(result.allele_frequencies))

class RNASeqPipeline:
    """Complete RNA-seq analysis pipeline"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.aligner = RNASeqAligner(db_builder)
        self.quantifier = ExpressionQuantifier(db_builder)
        self.eqtl_analyzer = eQTLAnalyzer(db_builder)
        
    async def run_complete_analysis(self, 
                                  fastq_files: List[str],
                                  sample_id: str,
                                  gene_regions: List[Dict[str, Any]],
                                  individuals: List[str] = None,
                                  analysis_params: Dict[str, Any] = None) -> Dict[str, Any]:
        """Run complete RNA-seq analysis pipeline"""
        
        if analysis_params is None:
            analysis_params = {
                'alignment': {
                    'min_alignment_score': 0.8,
                    'max_mismatches': 3,
                    'allow_multimapping': True
                },
                'quantification': {
                    'method': 'em',  # 'em', 'unique_only', 'proportional'
                    'min_expression': 0.1
                },
                'eqtl': {
                    'max_distance': 1000000,
                    'min_samples': 20,
                    'significance_threshold': 0.05
                }
            }
        
        results = {
            'sample_id': sample_id,
            'analysis_timestamp': datetime.now().isoformat(),
            'parameters': analysis_params
        }
        
        try:
            # Step 1: Parse RNA-seq reads
            logging.info(f"Parsing RNA-seq reads from {len(fastq_files)} files")
            reads = await self._parse_fastq_files(fastq_files)
            results['total_reads'] = len(reads)
            
            # Step 2: Align reads to transcript variants
            logging.info("Aligning reads to pangenome transcript variants")
            alignments = await self.aligner.align_reads_to_transcripts(
                reads, gene_regions, analysis_params['alignment']
            )
            results['total_alignments'] = len(alignments)
            results['mapping_rate'] = len(alignments) / len(reads) if reads else 0
            
            # Step 3: Quantify expression
            logging.info("Quantifying gene expression")
            transcript_lengths = await self._get_transcript_lengths(gene_regions)
            expression_results = await self.quantifier.quantify_expression(
                alignments, transcript_lengths, len(reads), 
                analysis_params['quantification']['method']
            )
            results['expression_results'] = {
                gene_id: {
                    'total_expression': expr.total_expression,
                    'fpkm': expr.fpkm,
                    'tpm': expr.tpm,
                    'allele_specific': expr.allele_specific_expression
                }
                for gene_id, expr in expression_results.items()
            }
            
            # Step 4: eQTL analysis (if individuals provided)
            if individuals and len(individuals) >= analysis_params['eqtl']['min_samples']:
                logging.info("Performing eQTL analysis")
                eqtl_results = await self.eqtl_analyzer.discover_eqtls(
                    expression_results, individuals,
                    max_distance=analysis_params['eqtl']['max_distance']
                )
                
                significant_eqtls = [
                    {
                        'variant_id': eqtl.variant_id,
                        'gene_id': eqtl.gene_id,
                        'pvalue': eqtl.pvalue,
                        'qvalue': eqtl.qvalue,
                        'beta': eqtl.beta,
                        'r_squared': eqtl.r_squared
                    }
                    for eqtl in eqtl_results 
                    if eqtl.qvalue < analysis_params['eqtl']['significance_threshold']
                ]
                
                results['eqtl_results'] = {
                    'total_tests': len(eqtl_results),
                    'significant_eqtls': significant_eqtls,
                    'top_associations': sorted(significant_eqtls, key=lambda x: x['pvalue'])[:10]
                }
            
            # Step 5: Generate summary statistics
            results['summary'] = await self._generate_summary_stats(expression_results, alignments)
            
            logging.info(f"RNA-seq analysis completed for sample {sample_id}")
            return results
            
        except Exception as e:
            logging.error(f"RNA-seq analysis failed for sample {sample_id}: {e}")
            results['error'] = str(e)
            return results
    
    async def _parse_fastq_files(self, fastq_files: List[str]) -> List[RNASeqRead]:
        """Parse FASTQ files into RNASeqRead objects"""
        
        reads = []
        
        for fastq_file in fastq_files:
            try:
                with open(fastq_file, 'r') as f:
                    lines = []
                    for line in f:
                        lines.append(line.strip())
                        
                        if len(lines) == 4:  # Complete FASTQ record
                            read_id = lines[0][1:]  # Remove @
                            sequence = lines[1]
                            quality = lines[3]
                            
                            reads.append(RNASeqRead(
                                read_id=read_id,
                                sequence=sequence,
                                quality=quality
                            ))
                            
                            lines = []
                            
            except Exception as e:
                logging.error(f"Error parsing FASTQ file {fastq_file}: {e}")
        
        return reads
    
    async def _get_transcript_lengths(self, gene_regions: List[Dict[str, Any]]) -> Dict[str, int]:
        """Get transcript lengths for expression normalization"""
        
        with self.db_builder.driver.session() as session:
            result = session.run("""
                MATCH (t:TranscriptVariant)
                WHERE t.gene_id IN $gene_ids
                RETURN t.transcript_id as transcript_id, t.length as length
            """, gene_ids=[gene['gene_id'] for gene in gene_regions])
            
            return {record['transcript_id']: record['length'] for record in result}
    
    async def _generate_summary_stats(self, expression_results: Dict[str, ExpressionResult],
                                    alignments: List[TranscriptAlignment]) -> Dict[str, Any]:
        """Generate summary statistics for the analysis"""
        
        # Expression statistics
        total_genes = len(expression_results)
        expressed_genes = sum(1 for expr in expression_results.values() if expr.total_expression > 0)
        
        expression_values = [expr.total_expression for expr in expression_results.values()]
        
        # Alignment statistics
        unique_alignments = sum(1 for aln in alignments if aln.is_unique)
        multi_alignments = len(alignments) - unique_alignments
        
        # Gene-level statistics
        avg_expression = np.mean(expression_values) if expression_values else 0
        median_expression = np.median(expression_values) if expression_values else 0
        
        return {
            'total_genes_tested': total_genes,
            'expressed_genes': expressed_genes,
            'expression_rate': expressed_genes / total_genes if total_genes > 0 else 0,
            'average_expression': avg_expression,
            'median_expression': median_expression,
            'unique_alignments': unique_alignments,
            'multi_alignments': multi_alignments,
            'unique_alignment_rate': unique_alignments / len(alignments) if alignments else 0
        }

# Visualization and reporting
class RNASeqVisualizer:
    """Create visualizations for RNA-seq analysis results"""
    
    @staticmethod
    def create_expression_heatmap(expression_data: Dict[str, Dict[str, float]], 
                                title: str = "Gene Expression Heatmap") -> go.Figure:
        """Create heatmap of gene expression across samples/alleles"""
        
        # Convert to matrix format
        genes = list(expression_data.keys())
        samples = list(set(sample for gene_data in expression_data.values() 
                          for sample in gene_data.keys()))
        
        matrix = []
        for gene in genes:
            row = []
            for sample in samples:
                value = expression_data[gene].get(sample, 0)
                row.append(value)
            matrix.append(row)
        
        fig = go.Figure(data=go.Heatmap(
            z=matrix,
            x=samples,
            y=genes,
            colorscale='Viridis',
            hoverongaps=False
        ))
        
        fig.update_layout(
            title=title,
            xaxis_title="Samples/Alleles",
            yaxis_title="Genes",
            width=800,
            height=600
        )
        
        return fig
    
    @staticmethod
    def create_eqtl_manhattan_plot(eqtl_results: List[eQTLResult],
                                 title: str = "eQTL Manhattan Plot") -> go.Figure:
        """Create Manhattan plot for eQTL results"""
        
        chromosomes = [eqtl.chromosome for eqtl in eqtl_results]
        positions = [eqtl.variant_position for eqtl in eqtl_results]
        pvalues = [-np.log10(eqtl.pvalue) for eqtl in eqtl_results]
        
        fig = go.Figure()
        
        # Plot points colored by chromosome
        colors = px.colors.qualitative.Set1
        for chr_num in sorted(set(chromosomes)):
            chr_mask = [i for i, c in enumerate(chromosomes) if c == chr_num]
            chr_positions = [positions[i] for i in chr_mask]
            chr_pvalues = [pvalues[i] for i in chr_mask]
            
            fig.add_trace(go.Scatter(
                x=chr_positions,
                y=chr_pvalues,
                mode='markers',
                marker=dict(size=6, color=colors[chr_num % len(colors)]),
                name=f'Chr {chr_num}',
                hovertemplate='Chr %{text}<br>Position: %{x}<br>-log10(p): %{y}<extra></extra>',
                text=[chr_num] * len(chr_positions)
            ))
        
        # Add significance threshold
        fig.add_hline(y=-np.log10(0.05), line_dash="dash", line_color="red",
                     annotation_text="Significance threshold (p=0.05)")
        
        fig.update_layout(
            title=title,
            xaxis_title="Genomic Position",
            yaxis_title="-log10(p-value)",
            showlegend=False,
            width=1200,
            height=600
        )
        
        return fig
    
    @staticmethod
    def create_allele_expression_plot(allele_expression: Dict[str, float],
                                    title: str = "Allele-Specific Expression") -> go.Figure:
        """Create plot showing allele-specific expression"""
        
        alleles = list(allele_expression.keys())
        expressions = list(allele_expression.values())
        
        fig = go.Figure(data=[
            go.Bar(x=alleles, y=expressions, marker_color='lightblue')
        ])
        
        fig.update_layout(
            title=title,
            xaxis_title="Alleles",
            yaxis_title="Expression Level",
            width=600,
            height=400
        )
        
        return fig

# Example usage and demonstration
async def demonstrate_rnaseq_analysis():
    """Demonstrate RNA-seq analysis with fractal pangenome"""
    
    print("🧬 RNA-seq Expression Analysis Demo")
    print("=" * 50)
    
    try:
        # Initialize system
        db_builder = GenomeDatabaseBuilder()
        await db_builder.initialize_database()
        
        pipeline = RNASeqPipeline(db_builder)
        
        # Define gene regions of interest
        gene_regions = [
            {
                'gene_id': 'BRCA1',
                'chromosome': 17,
                'start': 43044295,
                'end': 43125364,
                'strand': '-'
            },
            {
                'gene_id': 'TP53',
                'chromosome': 17,
                'start': 7661779,
                'end': 7687550,
                'strand': '-'
            },
            {
                'gene_id': 'APOE',
                'chromosome': 19,
                'start': 45409011,
                'end': 45412650,
                'strand': '+'
            }
        ]
        
        # Simulate FASTQ files
        temp_dir = Path(tempfile.mkdtemp())
        fastq_files = []
        
        # Create demo FASTQ file
        demo_fastq = temp_dir / "demo_rnaseq.fastq"
        with open(demo_fastq, 'w') as f:
            for i in range(1000):  # 1000 demo reads
                f.write(f"@read_{i}\n")
                f.write("ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG\n")
                f.write("+\n")
                f.write("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n")
        
        fastq_files.append(str(demo_fastq))
        
        # Define individuals for eQTL analysis
        individuals = ["HG002", "HG00733", "HG00096", "HG00097", "HG00419"]
        
        # Run complete analysis
        print("🔬 Running RNA-seq analysis pipeline...")
        results = await pipeline.run_complete_analysis(
            fastq_files=fastq_files,
            sample_id="DEMO_SAMPLE_001",
            gene_regions=gene_regions,
            individuals=individuals
        )
        
        # Display results
        print(f"\n📊 Analysis Results for {results['sample_id']}")
        print("-" * 40)
        print(f"Total reads processed: {results['total_reads']:,}")
        print(f"Total alignments: {results['total_alignments']:,}")
        print(f"Mapping rate: {results['mapping_rate']:.1%}")
        
        # Expression results
        if 'expression_results' in results:
            expr_results = results['expression_results']
            print(f"\n🧬 Gene Expression Results")
            print("-" * 30)
            
            for gene_id, expr_data in expr_results.items():
                print(f"{gene_id}:")
                print(f"  Total expression: {expr_data['total_expression']:.3f}")
                print(f"  FPKM: {expr_data['fpkm']:.2f}")
                print(f"  TPM: {expr_data['tpm']:.2f}")
                
                # Allele-specific expression
                if expr_data['allele_specific']:
                    print(f"  Allele-specific expression:")
                    for allele, expr in expr_data['allele_specific'].items():
                        print(f"    {allele}: {expr:.3f}")
                print()
        
        # eQTL results
        if 'eqtl_results' in results:
            eqtl_data = results['eqtl_results']
            print(f"🔗 eQTL Analysis Results")
            print("-" * 25)
            print(f"Total tests performed: {eqtl_data['total_tests']:,}")
            print(f"Significant eQTLs found: {len(eqtl_data['significant_eqtls'])}")
            
            if eqtl_data['top_associations']:
                print("\nTop eQTL associations:")
                for i, eqtl in enumerate(eqtl_data['top_associations'][:5], 1):
                    print(f"  {i}. {eqtl['variant_id']} → {eqtl['gene_id']}")
                    print(f"     P-value: {eqtl['pvalue']:.2e}")
                    print(f"     Effect size (β): {eqtl['beta']:.3f}")
                    print(f"     R²: {eqtl['r_squared']:.3f}")
        
        # Summary statistics
        if 'summary' in results:
            summary = results['summary']
            print(f"\n📈 Summary Statistics")
            print("-" * 20)
            print(f"Genes tested: {summary['total_genes_tested']}")
            print(f"Expressed genes: {summary['expressed_genes']}")
            print(f"Expression rate: {summary['expression_rate']:.1%}")
            print(f"Average expression: {summary['average_expression']:.3f}")
            print(f"Unique alignment rate: {summary['unique_alignment_rate']:.1%}")
        
        # Generate visualizations
        print(f"\n📊 Generating visualizations...")
        visualizer = RNASeqVisualizer()
        
        if 'expression_results' in results:
            # Expression heatmap
            expr_data_for_viz = {}
            for gene_id, expr_data in results['expression_results'].items():
                expr_data_for_viz[gene_id] = expr_data['allele_specific']
            
            if expr_data_for_viz:
                heatmap = visualizer.create_expression_heatmap(
                    expr_data_for_viz, "Demo Gene Expression"
                )
                print("  ✓ Expression heatmap created")
        
        # Save results
        results_file = temp_dir / "rnaseq_analysis_results.json"
        with open(results_file, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        
        print(f"\n📁 Results saved to: {results_file}")
        print("✅ RNA-seq analysis demonstration completed!")
        
    except Exception as e:
        print(f"❌ RNA-seq analysis demo error: {e}")
        logging.error(f"RNA-seq demo error: {e}")

# CLI integration
async def run_rnaseq_cli(args):
    """Command-line interface for RNA-seq analysis"""
    
    db_builder = GenomeDatabaseBuilder()
    await db_builder.initialize_database()
    
    pipeline = RNASeqPipeline(db_builder)
    
    # Parse gene regions from file or command line
    if args.gene_regions_file:
        with open(args.gene_regions_file, 'r') as f:
            gene_regions = json.load(f)
    else:
        gene_regions = [{
            'gene_id': args.gene_id,
            'chromosome': args.chromosome,
            'start': args.start,
            'end': args.end,
            'strand': args.strand or '+'
        }]
    
    # Run analysis
    results = await pipeline.run_complete_analysis(
        fastq_files=args.fastq_files,
        sample_id=args.sample_id,
        gene_regions=gene_regions,
        individuals=args.individuals
    )
    
    # Save results
    with open(args.output, 'w') as f:
        json.dump(results, f, indent=2, default=str)
    
    print(f"RNA-seq analysis completed. Results saved to {args.output}")

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Run demonstration
    asyncio.run(demonstrate_rnaseq_analysis())
