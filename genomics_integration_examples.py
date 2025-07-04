"""
Integration Examples with Existing Genomics Tools
=================================================
Demonstrations of how to integrate the fractal pangenome database 
with popular genomics tools and workflows
"""

import pandas as pd
import numpy as np
import asyncio
import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Union, Any
import subprocess
import tempfile
import json
import gzip
from datetime import datetime
import requests
from io import StringIO

# Bioinformatics libraries
try:
    import pysam
    import cyvcf2
    from pybedtools import BedTool
    import pyBigWig
    HAS_BIO_LIBS = True
except ImportError:
    HAS_BIO_LIBS = False
    print("Warning: Some bioinformatics libraries not available. Install with:")
    print("pip install pysam cyvcf2 pybedtools pyBigWig")

# Import our core modules
from rest_api_server import app_state
from neo4j_genome_importer import GenomeDatabaseBuilder
from hilbert_pangenome_architecture import HilbertCurve

class VCFIntegrator:
    """Integration with VCF (Variant Call Format) files"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        
    async def import_vcf_variants(self, vcf_path: str, sample_mapping: Dict[str, str] = None) -> Dict[str, Any]:
        """Import variants from VCF file into the pangenome database"""
        if not HAS_BIO_LIBS:
            raise ImportError("pysam and cyvcf2 required for VCF integration")
        
        results = {
            'imported_variants': 0,
            'skipped_variants': 0,
            'samples_processed': 0,
            'errors': []
        }
        
        try:
            # Read VCF file
            vcf = cyvcf2.VCF(vcf_path)
            samples = vcf.samples
            results['samples_processed'] = len(samples)
            
            # Process variants in batches
            batch_size = 1000
            variant_batch = []
            
            for variant in vcf:
                try:
                    # Extract variant information
                    chrom = self._parse_chromosome(variant.CHROM)
                    pos = variant.POS
                    ref = variant.REF
                    alts = variant.ALT
                    
                    if not alts:  # Skip reference-only variants
                        results['skipped_variants'] += 1
                        continue
                    
                    # Process each sample's genotype
                    for i, sample in enumerate(samples):
                        gt = variant.genotypes[i]
                        
                        # Map sample name to individual ID
                        individual_id = sample_mapping.get(sample, sample) if sample_mapping else sample
                        
                        # Create variant nodes for each allele
                        for haplotype in [0, 1]:  # Diploid
                            allele_idx = gt[haplotype] if len(gt) > haplotype else 0
                            
                            if allele_idx > 0:  # Non-reference allele
                                sequence = alts[allele_idx - 1] if allele_idx <= len(alts) else ref
                                frequency = self._calculate_allele_frequency(variant, allele_idx)
                                
                                variant_data = {
                                    'node_id': f"vcf_{chrom}_{pos}_{individual_id}_h{haplotype}_{allele_idx}",
                                    'individual_id': individual_id,
                                    'haplotype': haplotype,
                                    'chromosome': chrom,
                                    'start_pos': pos,
                                    'end_pos': pos + len(sequence) - 1,
                                    'sequence': sequence,
                                    'frequency': frequency,
                                    'scale_level': 0,
                                    'variant_type': self._classify_variant(ref, sequence),
                                    'vcf_info': {
                                        'qual': variant.QUAL,
                                        'filter': variant.FILTER,
                                        'ref': ref,
                                        'alt': sequence,
                                        'genotype': f"{gt[0]}/{gt[1]}" if len(gt) >= 2 else str(gt[0])
                                    }
                                }
                                
                                variant_batch.append(variant_data)
                    
                    # Import batch when full
                    if len(variant_batch) >= batch_size:
                        await self._import_variant_batch(variant_batch)
                        results['imported_variants'] += len(variant_batch)
                        variant_batch = []
                        
                except Exception as e:
                    results['errors'].append(f"Error processing variant {variant.CHROM}:{variant.POS}: {str(e)}")
                    results['skipped_variants'] += 1
            
            # Import remaining variants
            if variant_batch:
                await self._import_variant_batch(variant_batch)
                results['imported_variants'] += len(variant_batch)
            
        except Exception as e:
            results['errors'].append(f"VCF import error: {str(e)}")
        
        return results
    
    async def export_to_vcf(self, output_path: str, chromosome: int, 
                          start_pos: int, end_pos: int) -> str:
        """Export genomic region to VCF format"""
        
        with self.driver.session() as session:
            # Query variants in region
            query = """
            MATCH (n:GenomeNode)
            WHERE n.chromosome = $chr 
              AND n.start_pos >= $start 
              AND n.end_pos <= $end
              AND n.scale_level = 0
            RETURN n.individual_id as sample, n.chromosome as chrom,
                   n.start_pos as pos, n.sequence as alt,
                   n.frequency as af, n.haplotype as haplotype,
                   n.vcf_info as info
            ORDER BY n.start_pos, n.individual_id
            """
            
            results = session.run(query, chr=chromosome, start=start_pos, end=end_pos).data()
        
        # Group by position and create VCF entries
        variants_by_pos = {}
        for result in results:
            pos = result['pos']
            if pos not in variants_by_pos:
                variants_by_pos[pos] = {
                    'chrom': result['chrom'],
                    'pos': pos,
                    'samples': {},
                    'alts': set(),
                    'af': result['af']
                }
            
            sample = result['sample']
            if sample not in variants_by_pos[pos]['samples']:
                variants_by_pos[pos]['samples'][sample] = ['0', '0']  # Default genotype
            
            # Set haplotype allele
            haplotype = result['haplotype']
            if haplotype < 2:
                variants_by_pos[pos]['samples'][sample][haplotype] = '1'
                variants_by_pos[pos]['alts'].add(result['alt'])
        
        # Write VCF file
        with open(output_path, 'w') as f:
            # Write header
            f.write("##fileformat=VCFv4.2\n")
            f.write(f"##source=FractalPangenome\n")
            f.write(f"##reference=custom\n")
            f.write("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency\">\n")
            f.write("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n")
            
            # Sample names
            samples = set()
            for var_data in variants_by_pos.values():
                samples.update(var_data['samples'].keys())
            samples = sorted(samples)
            
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples) + "\n")
            
            # Write variants
            for pos in sorted(variants_by_pos.keys()):
                var_data = variants_by_pos[pos]
                alts = list(var_data['alts'])
                
                f.write(f"{var_data['chrom']}\t{pos}\t.\tN\t{','.join(alts)}\t60\tPASS\t")
                f.write(f"AF={var_data['af']:.3f}\tGT")
                
                for sample in samples:
                    gt = var_data['samples'].get(sample, ['0', '0'])
                    f.write(f"\t{gt[0]}/{gt[1]}")
                f.write("\n")
        
        return output_path
    
    def _parse_chromosome(self, chrom_str: str) -> int:
        """Parse chromosome string to integer"""
        chrom_str = chrom_str.replace('chr', '').replace('CHR', '')
        if chrom_str == 'X':
            return 23
        elif chrom_str == 'Y':
            return 24
        elif chrom_str == 'M' or chrom_str == 'MT':
            return 25
        else:
            try:
                return int(chrom_str)
            except ValueError:
                return 0
    
    def _calculate_allele_frequency(self, variant, allele_idx: int) -> float:
        """Calculate allele frequency from VCF variant"""
        # Simplified calculation - would use INFO field in real implementation
        total_alleles = len(variant.genotypes) * 2
        allele_count = sum(1 for gt in variant.genotypes 
                          for allele in gt[:2] if allele == allele_idx)
        return allele_count / total_alleles if total_alleles > 0 else 0.0
    
    def _classify_variant(self, ref: str, alt: str) -> str:
        """Classify variant type"""
        if len(ref) == 1 and len(alt) == 1:
            return 'SNP'
        elif len(ref) < len(alt):
            return 'insertion'
        elif len(ref) > len(alt):
            return 'deletion'
        else:
            return 'complex'
    
    async def _import_variant_batch(self, variant_batch: List[Dict]) -> None:
        """Import a batch of variants to the database"""
        with self.driver.session() as session:
            session.run("""
                UNWIND $variants as variant
                CREATE (n:GenomeNode {
                    id: variant.node_id,
                    individual_id: variant.individual_id,
                    haplotype: variant.haplotype,
                    chromosome: variant.chromosome,
                    start_pos: variant.start_pos,
                    end_pos: variant.end_pos,
                    sequence: variant.sequence,
                    frequency: variant.frequency,
                    scale_level: variant.scale_level,
                    variant_type: variant.variant_type,
                    vcf_info: variant.vcf_info,
                    import_timestamp: datetime()
                })
            """, variants=variant_batch)

class BEDIntegrator:
    """Integration with BED (Browser Extensible Data) format"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
    
    async def annotate_with_bed(self, bed_path: str, annotation_name: str) -> Dict[str, Any]:
        """Annotate genomic nodes with BED file annotations"""
        if not HAS_BIO_LIBS:
            raise ImportError("pybedtools required for BED integration")
        
        results = {'annotated_nodes': 0, 'total_annotations': 0}
        
        try:
            # Load BED file
            bed = BedTool(bed_path)
            
            # Process each BED interval
            for interval in bed:
                chrom = self._parse_chromosome(interval.chrom)
                start = interval.start
                end = interval.end
                name = interval.name if len(interval.fields) > 3 else f"{annotation_name}_{results['total_annotations']}"
                score = interval.score if len(interval.fields) > 4 else 1000
                
                # Find overlapping nodes in database
                with self.driver.session() as session:
                    overlapping_nodes = session.run("""
                        MATCH (n:GenomeNode)
                        WHERE n.chromosome = $chrom
                          AND n.start_pos <= $end
                          AND n.end_pos >= $start
                        RETURN n.id as node_id, n.start_pos as start_pos, n.end_pos as end_pos
                    """, chrom=chrom, start=start, end=end).data()
                    
                    # Add annotations
                    for node in overlapping_nodes:
                        overlap_start = max(start, node['start_pos'])
                        overlap_end = min(end, node['end_pos'])
                        overlap_length = overlap_end - overlap_start + 1
                        
                        session.run("""
                            MATCH (n:GenomeNode {id: $node_id})
                            SET n.annotations = coalesce(n.annotations, {}) + {
                                $annotation_key: {
                                    name: $name,
                                    score: $score,
                                    overlap_start: $overlap_start,
                                    overlap_end: $overlap_end,
                                    overlap_length: $overlap_length,
                                    annotation_type: $annotation_name
                                }
                            }
                        """, node_id=node['node_id'], annotation_key=f"{annotation_name}_{results['total_annotations']}",
                            name=name, score=score, overlap_start=overlap_start,
                            overlap_end=overlap_end, overlap_length=overlap_length,
                            annotation_name=annotation_name)
                        
                        results['annotated_nodes'] += 1
                
                results['total_annotations'] += 1
        
        except Exception as e:
            logging.error(f"BED annotation error: {e}")
            raise e
        
        return results
    
    async def export_to_bed(self, output_path: str, chromosome: int, 
                          annotation_filter: str = None) -> str:
        """Export genomic regions to BED format"""
        
        with self.driver.session() as session:
            query = """
            MATCH (n:GenomeNode)
            WHERE n.chromosome = $chr AND n.scale_level = 0
            """
            
            if annotation_filter:
                query += " AND any(key in keys(n.annotations) WHERE key CONTAINS $filter)"
            
            query += """
            RETURN n.chromosome as chrom, n.start_pos as start, n.end_pos as end,
                   n.individual_id as name, n.frequency as score, n.annotations as annotations
            ORDER BY n.start_pos
            """
            
            params = {'chr': chromosome}
            if annotation_filter:
                params['filter'] = annotation_filter
            
            results = session.run(query, **params).data()
        
        # Write BED file
        with open(output_path, 'w') as f:
            for result in results:
                chrom = f"chr{result['chrom']}"
                start = result['start'] - 1  # BED is 0-based
                end = result['end']
                name = result['name']
                score = int(result['score'] * 1000)  # Scale to 0-1000
                
                f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\n")
        
        return output_path
    
    def _parse_chromosome(self, chrom_str: str) -> int:
        """Parse chromosome string to integer"""
        chrom_str = chrom_str.replace('chr', '').replace('CHR', '')
        if chrom_str == 'X':
            return 23
        elif chrom_str == 'Y':
            return 24
        elif chrom_str == 'M' or chrom_str == 'MT':
            return 25
        else:
            try:
                return int(chrom_str)
            except ValueError:
                return 0

class GATKIntegrator:
    """Integration with GATK (Genome Analysis Toolkit) workflows"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        
    async def create_reference_from_pangenome(self, individual_id: str, 
                                            output_fasta: str, 
                                            chromosome: int = None) -> str:
        """Create a FASTA reference from an individual's pangenome path"""
        
        driver = self.db_builder.neo4j_manager.get_connection()
        
        with driver.session() as session:
            # Query for individual's genomic nodes
            query = """
            MATCH (n:GenomeNode {individual_id: $individual_id})
            WHERE n.scale_level = 0
            """
            
            if chromosome:
                query += " AND n.chromosome = $chromosome"
            
            query += """
            RETURN n.chromosome as chr, n.start_pos as start, n.end_pos as end,
                   n.sequence as sequence, n.haplotype as haplotype
            ORDER BY n.chromosome, n.haplotype, n.start_pos
            """
            
            params = {'individual_id': individual_id}
            if chromosome:
                params['chromosome'] = chromosome
            
            results = session.run(query, **params).data()
        
        driver.close()
        
        # Group by chromosome and haplotype
        sequences = {}
        for result in results:
            chr_hap = f"chr{result['chr']}_h{result['haplotype']}"
            if chr_hap not in sequences:
                sequences[chr_hap] = []
            sequences[chr_hap].append({
                'start': result['start'],
                'end': result['end'],
                'sequence': result['sequence']
            })
        
        # Write FASTA file
        with open(output_fasta, 'w') as f:
            for chr_hap, segments in sequences.items():
                # Sort segments by position
                segments.sort(key=lambda x: x['start'])
                
                # Concatenate sequences
                full_sequence = ""
                last_end = 0
                
                for segment in segments:
                    # Fill gaps with Ns
                    gap_size = segment['start'] - last_end - 1
                    if gap_size > 0:
                        full_sequence += 'N' * gap_size
                    
                    full_sequence += segment['sequence']
                    last_end = segment['end']
                
                # Write FASTA entry
                f.write(f">{chr_hap}_{individual_id}\n")
                
                # Write sequence in 80-character lines
                for i in range(0, len(full_sequence), 80):
                    f.write(full_sequence[i:i+80] + "\n")
        
        return output_fasta
    
    async def run_gatk_haplotypecaller(self, bam_file: str, reference_fasta: str,
                                     output_vcf: str, region: str = None) -> Dict[str, Any]:
        """Run GATK HaplotypeCaller and import results"""
        
        # Construct GATK command
        cmd = [
            'gatk', 'HaplotypeCaller',
            '-R', reference_fasta,
            '-I', bam_file,
            '-O', output_vcf
        ]
        
        if region:
            cmd.extend(['-L', region])
        
        try:
            # Run GATK
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # Import resulting VCF
            vcf_integrator = VCFIntegrator(self.db_builder)
            import_results = await vcf_integrator.import_vcf_variants(output_vcf)
            
            return {
                'gatk_success': True,
                'gatk_output': result.stdout,
                'import_results': import_results
            }
            
        except subprocess.CalledProcessError as e:
            return {
                'gatk_success': False,
                'error': e.stderr,
                'return_code': e.returncode
            }
        except Exception as e:
            return {
                'gatk_success': False,
                'error': str(e)
            }

class EnsemblIntegrator:
    """Integration with Ensembl REST API for gene annotations"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        self.ensembl_base_url = "https://rest.ensembl.org"
    
    async def annotate_genes(self, chromosome: int, start: int, end: int) -> Dict[str, Any]:
        """Fetch gene annotations from Ensembl and add to database"""
        
        results = {'genes_annotated': 0, 'annotations_added': 0}
        
        try:
            # Query Ensembl for genes in region
            chr_name = f"chr{chromosome}" if chromosome <= 22 else ("chrX" if chromosome == 23 else "chrY")
            url = f"{self.ensembl_base_url}/overlap/region/human/{chr_name}:{start}-{end}"
            
            params = {
                'feature': 'gene',
                'content-type': 'application/json'
            }
            
            response = requests.get(url, params=params)
            response.raise_for_status()
            
            genes = response.json()
            
            # Process each gene
            for gene in genes:
                gene_start = gene.get('start', 0)
                gene_end = gene.get('end', 0)
                gene_id = gene.get('id', '')
                gene_name = gene.get('external_name', gene_id)
                gene_type = gene.get('biotype', 'unknown')
                strand = gene.get('strand', 0)
                
                # Find overlapping nodes in database
                with self.driver.session() as session:
                    overlapping_nodes = session.run("""
                        MATCH (n:GenomeNode)
                        WHERE n.chromosome = $chrom
                          AND n.start_pos <= $gene_end
                          AND n.end_pos >= $gene_start
                          AND n.scale_level = 0
                        RETURN n.id as node_id
                    """, chrom=chromosome, gene_start=gene_start, gene_end=gene_end).data()
                    
                    # Add gene annotations
                    for node in overlapping_nodes:
                        session.run("""
                            MATCH (n:GenomeNode {id: $node_id})
                            SET n.gene_annotations = coalesce(n.gene_annotations, []) + [{
                                gene_id: $gene_id,
                                gene_name: $gene_name,
                                gene_type: $gene_type,
                                gene_start: $gene_start,
                                gene_end: $gene_end,
                                strand: $strand,
                                source: 'Ensembl'
                            }]
                        """, node_id=node['node_id'], gene_id=gene_id, gene_name=gene_name,
                            gene_type=gene_type, gene_start=gene_start, gene_end=gene_end,
                            strand=strand)
                        
                        results['annotations_added'] += 1
                
                results['genes_annotated'] += 1
        
        except Exception as e:
            logging.error(f"Ensembl annotation error: {e}")
            raise e
        
        return results

class ClinVarIntegrator:
    """Integration with ClinVar for clinical variant annotations"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
    
    async def annotate_clinical_variants(self, chromosome: int, start: int, end: int) -> Dict[str, Any]:
        """Fetch clinical variant annotations from ClinVar"""
        
        results = {'clinical_variants': 0, 'pathogenic_variants': 0}
        
        try:
            # This would typically involve downloading ClinVar VCF and processing
            # For demo purposes, we'll simulate clinical annotations
            
            with self.driver.session() as session:
                # Find variants in region
                variants = session.run("""
                    MATCH (n:GenomeNode)
                    WHERE n.chromosome = $chrom
                      AND n.start_pos >= $start
                      AND n.end_pos <= $end
                      AND n.variant_type IS NOT NULL
                    RETURN n.id as node_id, n.start_pos as pos, n.sequence as alt
                """, chrom=chromosome, start=start, end=end).data()
                
                # Simulate clinical annotations for some variants
                for i, variant in enumerate(variants):
                    if i % 10 == 0:  # Annotate every 10th variant for demo
                        clinical_significance = np.random.choice([
                            'Pathogenic', 'Likely_pathogenic', 'Uncertain_significance',
                            'Likely_benign', 'Benign'
                        ])
                        
                        session.run("""
                            MATCH (n:GenomeNode {id: $node_id})
                            SET n.clinical_annotations = {
                                clinical_significance: $significance,
                                review_status: 'criteria_provided_single_submitter',
                                condition: 'Simulated_condition',
                                source: 'ClinVar_simulation',
                                last_evaluated: date()
                            }
                        """, node_id=variant['node_id'], significance=clinical_significance)
                        
                        results['clinical_variants'] += 1
                        
                        if clinical_significance in ['Pathogenic', 'Likely_pathogenic']:
                            results['pathogenic_variants'] += 1
        
        except Exception as e:
            logging.error(f"ClinVar annotation error: {e}")
            raise e
        
        return results

class PipelineOrchestrator:
    """Orchestrate complex genomics workflows"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.vcf_integrator = VCFIntegrator(db_builder)
        self.bed_integrator = BEDIntegrator(db_builder)
        self.gatk_integrator = GATKIntegrator(db_builder)
        self.ensembl_integrator = EnsemblIntegrator(db_builder)
        self.clinvar_integrator = ClinVarIntegrator(db_builder)
    
    async def comprehensive_annotation_pipeline(self, 
                                              vcf_file: str,
                                              bed_files: Dict[str, str],
                                              chromosome: int,
                                              start: int,
                                              end: int) -> Dict[str, Any]:
        """Run comprehensive annotation pipeline"""
        
        pipeline_results = {
            'pipeline_start': datetime.now().isoformat(),
            'steps_completed': [],
            'total_annotations': 0
        }
        
        try:
            # Step 1: Import VCF variants
            print("Step 1: Importing VCF variants...")
            vcf_results = await self.vcf_integrator.import_vcf_variants(vcf_file)
            pipeline_results['vcf_import'] = vcf_results
            pipeline_results['steps_completed'].append('vcf_import')
            
            # Step 2: Annotate with BED files
            print("Step 2: Adding BED annotations...")
            bed_results = {}
            for annotation_name, bed_file in bed_files.items():
                result = await self.bed_integrator.annotate_with_bed(bed_file, annotation_name)
                bed_results[annotation_name] = result
                pipeline_results['total_annotations'] += result['annotated_nodes']
            
            pipeline_results['bed_annotations'] = bed_results
            pipeline_results['steps_completed'].append('bed_annotations')
            
            # Step 3: Add gene annotations from Ensembl
            print("Step 3: Fetching gene annotations...")
            gene_results = await self.ensembl_integrator.annotate_genes(chromosome, start, end)
            pipeline_results['gene_annotations'] = gene_results
            pipeline_results['total_annotations'] += gene_results['annotations_added']
            pipeline_results['steps_completed'].append('gene_annotations')
            
            # Step 4: Add clinical annotations
            print("Step 4: Adding clinical annotations...")
            clinical_results = await self.clinvar_integrator.annotate_clinical_variants(chromosome, start, end)
            pipeline_results['clinical_annotations'] = clinical_results
            pipeline_results['steps_completed'].append('clinical_annotations')
            
            # Step 5: Generate summary report
            pipeline_results['summary'] = await self._generate_annotation_summary(chromosome, start, end)
            pipeline_results['steps_completed'].append('summary_generation')
            
            pipeline_results['pipeline_end'] = datetime.now().isoformat()
            pipeline_results['success'] = True
            
        except Exception as e:
            pipeline_results['error'] = str(e)
            pipeline_results['success'] = False
            pipeline_results['pipeline_end'] = datetime.now().isoformat()
        
        return pipeline_results
    
    async def _generate_annotation_summary(self, chromosome: int, start: int, end: int) -> Dict[str, Any]:
        """Generate summary of annotations in region"""
        
        driver = self.db_builder.neo4j_manager.get_connection()
        
        with driver.session() as session:
            # Count different types of annotations
            summary_query = """
            MATCH (n:GenomeNode)
            WHERE n.chromosome = $chrom AND n.start_pos >= $start AND n.end_pos <= $end
            RETURN 
                count(n) as total_nodes,
                count(CASE WHEN n.gene_annotations IS NOT NULL THEN 1 END) as nodes_with_genes,
                count(CASE WHEN n.clinical_annotations IS NOT NULL THEN 1 END) as nodes_with_clinical,
                count(CASE WHEN n.annotations IS NOT NULL THEN 1 END) as nodes_with_bed_annotations,
                count(CASE WHEN n.variant_type IS NOT NULL THEN 1 END) as variant_nodes
            """
            
            summary = session.run(summary_query, chrom=chromosome, start=start, end=end).single()
            
            # Get variant type distribution
            variant_types = session.run("""
                MATCH (n:GenomeNode)
                WHERE n.chromosome = $chrom AND n.start_pos >= $start AND n.end_pos <= $end
                  AND n.variant_type IS NOT NULL
                RETURN n.variant_type as type, count(n) as count
            """, chrom=chromosome, start=start, end=end).data()
        
        driver.close()
        
        return {
            'region': f"chr{chromosome}:{start}-{end}",
            'statistics': dict(summary) if summary else {},
            'variant_types': {vt['type']: vt['count'] for vt in variant_types}
        }

# Example usage and demonstration
async def demonstrate_integrations():
    """Demonstrate various tool integrations"""
    
    print("🧬 Fractal Pangenome Database - Tool Integration Demo")
    print("=" * 60)
    
    try:
        # Initialize database
        db_builder = GenomeDatabaseBuilder()
        await db_builder.initialize_database()
        
        # Create demo data files
        demo_dir = Path("./integration_demo")
        demo_dir.mkdir(exist_ok=True)
        
        # 1. VCF Integration Demo
        print("\n📄 VCF Integration Demo")
        print("-" * 30)
        
        # Create sample VCF
        sample_vcf = demo_dir / "sample.vcf"
        with open(sample_vcf, 'w') as f:
            f.write("""##fileformat=VCFv4.2
##source=DemoData
##reference=demo
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Sample1	Sample2
1	1000000	.	A	T	60	PASS	AF=0.3	GT	0/1	1/1
1	1000010	.	G	C	70	PASS	AF=0.2	GT	0/0	0/1
1	1000020	.	C	CTTT	80	PASS	AF=0.1	GT	0/1	0/0
""")
        
        vcf_integrator = VCFIntegrator(db_builder)
        vcf_results = await vcf_integrator.import_vcf_variants(str(sample_vcf))
        print(f"✅ VCF Import: {vcf_results['imported_variants']} variants imported")
        
        # 2. BED Integration Demo
        print("\n🛏️ BED Integration Demo")
        print("-" * 30)
        
        # Create sample BED file
        sample_bed = demo_dir / "genes.bed"
        with open(sample_bed, 'w') as f:
            f.write("""chr1	999000	1001000	GENE1	800
chr1	1000500	1002000	GENE2	900
""")
        
        bed_integrator = BEDIntegrator(db_builder)
        bed_results = await bed_integrator.annotate_with_bed(str(sample_bed), "test_genes")
        print(f"✅ BED Annotation: {bed_results['annotated_nodes']} nodes annotated")
        
        # 3. Ensembl Integration Demo (if internet available)
        print("\n🧬 Ensembl Integration Demo")
        print("-" * 30)
        
        try:
            ensembl_integrator = EnsemblIntegrator(db_builder)
            ensembl_results = await ensembl_integrator.annotate_genes(1, 1000000, 1100000)
            print(f"✅ Ensembl Annotation: {ensembl_results['genes_annotated']} genes annotated")
        except Exception as e:
            print(f"⚠️ Ensembl annotation skipped: {e}")
        
        # 4. Comprehensive Pipeline Demo
        print("\n🔄 Comprehensive Pipeline Demo")
        print("-" * 30)
        
        orchestrator = PipelineOrchestrator(db_builder)
        pipeline_results = await orchestrator.comprehensive_annotation_pipeline(
            vcf_file=str(sample_vcf),
            bed_files={"genes": str(sample_bed)},
            chromosome=1,
            start=999000,
            end=1003000
        )
        
        print(f"✅ Pipeline completed: {len(pipeline_results['steps_completed'])} steps")
        print(f"   Total annotations: {pipeline_results['total_annotations']}")
        
        # Save results
        results_file = demo_dir / "integration_results.json"
        with open(results_file, 'w') as f:
            json.dump({
                'vcf_results': vcf_results,
                'bed_results': bed_results,
                'pipeline_results': pipeline_results
            }, f, indent=2, default=str)
        
        print(f"\n📁 Results saved to: {results_file}")
        
    except Exception as e:
        print(f"❌ Integration demo error: {e}")
        logging.error(f"Integration demo error: {e}")

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Run demonstration
    asyncio.run(demonstrate_integrations())
