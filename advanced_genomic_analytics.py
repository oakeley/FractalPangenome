"""
Advanced Genomic Analytics and Population Studies
=================================================
Sophisticated analytics leveraging the fractal pangenome database for:
- Population genomics analysis
- Genome-wide association studies (GWAS)
- Pharmacogenomics
- Evolutionary genomics
- Clinical decision support
"""

import numpy as np
import pandas as pd
import scipy.stats as stats
from scipy.spatial.distance import pdist, squareform
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
import matplotlib.pyplot as plt
import seaborn as sns
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import asyncio
import logging
from typing import Dict, List, Tuple, Optional, Union, Any
from dataclasses import dataclass, field
from datetime import datetime, timedelta
import json
from pathlib import Path
import warnings
warnings.filterwarnings('ignore')

# Import our core modules
from neo4j_genome_importer import GenomeDatabaseBuilder
from hilbert_pangenome_architecture import HilbertCurve

@dataclass
class PopulationMetrics:
    """Population genetics metrics"""
    allele_frequencies: Dict[str, float]
    heterozygosity: float
    inbreeding_coefficient: float
    tajimas_d: float
    nucleotide_diversity: float
    population_size: int
    
@dataclass
class GWASResult:
    """GWAS analysis result"""
    variant_id: str
    chromosome: int
    position: int
    p_value: float
    odds_ratio: float
    beta: float
    standard_error: float
    allele_frequency: float
    phenotype_association: str

@dataclass
class PharmacogenomicProfile:
    """Pharmacogenomic profile for an individual"""
    individual_id: str
    drug_responses: Dict[str, str]  # drug -> response prediction
    metabolizer_status: Dict[str, str]  # enzyme -> status
    risk_alleles: List[str]
    recommendations: List[str]

class PopulationGenomicsAnalyzer:
    """Advanced population genomics analysis"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        
    async def calculate_population_structure(self, 
                                           individuals: List[str],
                                           chromosome: Optional[int] = None) -> Dict[str, Any]:
        """Calculate population genetic structure using PCA and clustering"""
        
        # Extract variant matrix
        variant_matrix, variant_info = await self._extract_variant_matrix(individuals, chromosome)
        
        if variant_matrix.shape[0] == 0:
            return {"error": "No variants found for analysis"}
        
        # Calculate allele frequencies
        allele_frequencies = np.mean(variant_matrix, axis=0)
        
        # Filter variants by MAF (Minor Allele Frequency)
        maf_threshold = 0.05
        maf_filter = (allele_frequencies >= maf_threshold) & (allele_frequencies <= 1 - maf_threshold)
        filtered_matrix = variant_matrix[:, maf_filter]
        filtered_variants = [v for i, v in enumerate(variant_info) if maf_filter[i]]
        
        if filtered_matrix.shape[1] == 0:
            return {"error": "No variants passed MAF filter"}
        
        # Principal Component Analysis
        pca = PCA(n_components=min(10, filtered_matrix.shape[0] - 1))
        pca_coords = pca.fit_transform(filtered_matrix)
        
        # Population clustering
        n_clusters = min(5, len(individuals))
        kmeans = KMeans(n_clusters=n_clusters, random_state=42)
        cluster_labels = kmeans.fit_predict(filtered_matrix)
        
        # Calculate population metrics
        pop_metrics = self._calculate_population_metrics(filtered_matrix, allele_frequencies[maf_filter])
        
        # Pairwise genetic distances
        genetic_distances = pdist(filtered_matrix, metric='hamming')
        distance_matrix = squareform(genetic_distances)
        
        return {
            "individuals": individuals,
            "pca_coordinates": pca_coords.tolist(),
            "pca_explained_variance": pca.explained_variance_ratio_.tolist(),
            "cluster_assignments": cluster_labels.tolist(),
            "population_metrics": pop_metrics,
            "genetic_distance_matrix": distance_matrix.tolist(),
            "variant_count": filtered_matrix.shape[1],
            "individual_count": len(individuals)
        }
    
    async def analyze_population_diversity(self, 
                                         populations: Dict[str, List[str]],
                                         chromosome: int) -> Dict[str, Any]:
        """Analyze genetic diversity within and between populations"""
        
        diversity_results = {}
        
        for pop_name, individuals in populations.items():
            # Calculate within-population diversity
            variant_matrix, variant_info = await self._extract_variant_matrix(individuals, chromosome)
            
            if variant_matrix.shape[0] > 0:
                allele_frequencies = np.mean(variant_matrix, axis=0)
                diversity_results[pop_name] = self._calculate_population_metrics(
                    variant_matrix, allele_frequencies
                )
        
        # Calculate between-population differentiation (Fst)
        fst_matrix = await self._calculate_fst_matrix(populations, chromosome)
        
        # Find population-specific variants
        pop_specific_variants = await self._find_population_specific_variants(populations, chromosome)
        
        return {
            "within_population_diversity": diversity_results,
            "fst_matrix": fst_matrix,
            "population_specific_variants": pop_specific_variants,
            "analysis_timestamp": datetime.now().isoformat()
        }
    
    async def _extract_variant_matrix(self, individuals: List[str], 
                                    chromosome: Optional[int] = None) -> Tuple[np.ndarray, List[Dict]]:
        """Extract variant matrix for individuals"""
        
        with self.driver.session() as session:
            # Build query
            query = """
            MATCH (n:GenomeNode)
            WHERE n.individual_id IN $individuals
              AND n.scale_level = 0
              AND n.variant_type IS NOT NULL
            """
            
            if chromosome:
                query += " AND n.chromosome = $chromosome"
            
            query += """
            RETURN n.individual_id as individual, n.chromosome as chr,
                   n.start_pos as pos, n.haplotype as haplotype,
                   n.frequency as freq, n.variant_type as type,
                   n.node_id as variant_id
            ORDER BY n.chromosome, n.start_pos, n.individual_id, n.haplotype
            """
            
            params = {"individuals": individuals}
            if chromosome:
                params["chromosome"] = chromosome
            
            results = session.run(query, **params).data()
        
        if not results:
            return np.array([]), []
        
        # Convert to variant matrix (individuals x variants)
        df = pd.DataFrame(results)
        
        # Create variant identifier
        df['variant_pos'] = df['chr'].astype(str) + '_' + df['pos'].astype(str)
        
        # Pivot to create matrix (individual-haplotype x variant)
        df['individual_haplotype'] = df['individual'] + '_h' + df['haplotype'].astype(str)
        
        # Create presence/absence matrix
        variant_matrix = df.pivot_table(
            index='individual_haplotype',
            columns='variant_pos',
            values='freq',
            fill_value=0,
            aggfunc='mean'
        ).values
        
        # Convert to binary (presence/absence)
        variant_matrix = (variant_matrix > 0).astype(int)
        
        # Get variant information
        variant_info = df.groupby('variant_pos').first()[['chr', 'pos', 'type']].to_dict('records')
        
        return variant_matrix, variant_info
    
    def _calculate_population_metrics(self, variant_matrix: np.ndarray, 
                                    allele_frequencies: np.ndarray) -> PopulationMetrics:
        """Calculate population genetics metrics"""
        
        n_individuals = variant_matrix.shape[0]
        n_variants = variant_matrix.shape[1]
        
        # Heterozygosity (expected vs observed)
        expected_het = 2 * allele_frequencies * (1 - allele_frequencies)
        observed_het = np.mean([np.mean(variant_matrix[i::2] != variant_matrix[1+i::2]) 
                               for i in range(0, n_individuals-1, 2)])
        
        # Inbreeding coefficient (Fis)
        mean_expected_het = np.mean(expected_het)
        fis = (mean_expected_het - observed_het) / mean_expected_het if mean_expected_het > 0 else 0
        
        # Nucleotide diversity (π)
        nucleotide_diversity = np.mean(expected_het)
        
        # Tajima's D (simplified calculation)
        tajimas_d = self._calculate_tajimas_d(variant_matrix)
        
        return PopulationMetrics(
            allele_frequencies={f"var_{i}": freq for i, freq in enumerate(allele_frequencies)},
            heterozygosity=observed_het,
            inbreeding_coefficient=fis,
            tajimas_d=tajimas_d,
            nucleotide_diversity=nucleotide_diversity,
            population_size=n_individuals
        )
    
    def _calculate_tajimas_d(self, variant_matrix: np.ndarray) -> float:
        """Calculate Tajima's D statistic (simplified)"""
        n = variant_matrix.shape[0]
        if n < 3:
            return 0.0
        
        # Number of segregating sites
        S = np.sum(np.any(variant_matrix == 1, axis=0) & np.any(variant_matrix == 0, axis=0))
        
        if S == 0:
            return 0.0
        
        # Average number of pairwise differences
        pi = np.mean([np.sum(variant_matrix[i] != variant_matrix[j]) 
                     for i in range(n) for j in range(i+1, n)])
        
        # Expected value under neutrality
        a1 = sum(1/i for i in range(1, n))
        expected_pi = S / a1
        
        # Tajima's D (simplified)
        if expected_pi > 0:
            return (pi - expected_pi) / np.sqrt(expected_pi)
        else:
            return 0.0
    
    async def _calculate_fst_matrix(self, populations: Dict[str, List[str]], 
                                  chromosome: int) -> Dict[str, Any]:
        """Calculate pairwise Fst between populations"""
        
        pop_names = list(populations.keys())
        n_pops = len(pop_names)
        fst_matrix = np.zeros((n_pops, n_pops))
        
        # Calculate Fst for each pair of populations
        for i in range(n_pops):
            for j in range(i+1, n_pops):
                pop1_individuals = populations[pop_names[i]]
                pop2_individuals = populations[pop_names[j]]
                
                # Get variant data for both populations
                matrix1, _ = await self._extract_variant_matrix(pop1_individuals, chromosome)
                matrix2, _ = await self._extract_variant_matrix(pop2_individuals, chromosome)
                
                if matrix1.shape[1] > 0 and matrix2.shape[1] > 0:
                    # Calculate Fst (simplified Hudson's Fst)
                    fst = self._calculate_hudson_fst(matrix1, matrix2)
                    fst_matrix[i, j] = fst
                    fst_matrix[j, i] = fst
        
        return {
            "population_names": pop_names,
            "fst_matrix": fst_matrix.tolist(),
            "mean_fst": np.mean(fst_matrix[fst_matrix > 0])
        }
    
    def _calculate_hudson_fst(self, matrix1: np.ndarray, matrix2: np.ndarray) -> float:
        """Calculate Hudson's Fst between two populations"""
        
        # Align matrices to common variants (simplified)
        min_variants = min(matrix1.shape[1], matrix2.shape[1])
        if min_variants == 0:
            return 0.0
        
        matrix1 = matrix1[:, :min_variants]
        matrix2 = matrix2[:, :min_variants]
        
        # Calculate allele frequencies
        p1 = np.mean(matrix1, axis=0)
        p2 = np.mean(matrix2, axis=0)
        
        # Hudson's Fst calculation
        numerator = np.var([p1, p2], axis=0)
        denominator = p1 * (1 - p1) + p2 * (1 - p2)
        
        # Avoid division by zero
        valid_sites = denominator > 0
        if np.sum(valid_sites) == 0:
            return 0.0
        
        fst_values = numerator[valid_sites] / denominator[valid_sites]
        return np.mean(fst_values)
    
    async def _find_population_specific_variants(self, populations: Dict[str, List[str]], 
                                               chromosome: int) -> Dict[str, List[Dict]]:
        """Find variants specific to each population"""
        
        pop_specific = {}
        
        for target_pop, target_individuals in populations.items():
            # Get variants for target population
            with self.driver.session() as session:
                target_variants = session.run("""
                    MATCH (n:GenomeNode)
                    WHERE n.individual_id IN $individuals
                      AND n.chromosome = $chromosome
                      AND n.scale_level = 0
                      AND n.variant_type IS NOT NULL
                    RETURN DISTINCT n.chromosome as chr, n.start_pos as pos, 
                           n.variant_type as type, count(n) as count
                """, individuals=target_individuals, chromosome=chromosome).data()
                
                # Get variants in other populations
                other_individuals = [ind for pop, inds in populations.items() 
                                   if pop != target_pop for ind in inds]
                
                if other_individuals:
                    other_variants = session.run("""
                        MATCH (n:GenomeNode)
                        WHERE n.individual_id IN $individuals
                          AND n.chromosome = $chromosome
                          AND n.scale_level = 0
                          AND n.variant_type IS NOT NULL
                        RETURN DISTINCT n.chromosome as chr, n.start_pos as pos
                    """, individuals=other_individuals, chromosome=chromosome).data()
                    
                    other_positions = {(v['chr'], v['pos']) for v in other_variants}
                else:
                    other_positions = set()
                
                # Find population-specific variants
                specific_variants = []
                for variant in target_variants:
                    pos_key = (variant['chr'], variant['pos'])
                    if pos_key not in other_positions:
                        specific_variants.append(variant)
                
                pop_specific[target_pop] = specific_variants
        
        return pop_specific

class GWASAnalyzer:
    """Genome-Wide Association Study analysis"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
    
    async def run_gwas_analysis(self, 
                              cases: List[str],
                              controls: List[str],
                              phenotype_name: str,
                              chromosome: Optional[int] = None,
                              p_value_threshold: float = 5e-8) -> Dict[str, Any]:
        """Run GWAS analysis comparing cases vs controls"""
        
        logging.info(f"Starting GWAS analysis for {phenotype_name}")
        
        # Extract variant data
        all_individuals = cases + controls
        variant_matrix, variant_info = await self._extract_variant_data(all_individuals, chromosome)
        
        if variant_matrix.shape[0] == 0:
            return {"error": "No variant data found"}
        
        # Create phenotype vector (1 = case, 0 = control)
        phenotype = np.array([1] * len(cases) + [0] * len(controls))
        
        # Run association tests
        gwas_results = []
        
        for i, variant in enumerate(variant_info):
            genotypes = variant_matrix[:, i]
            
            # Skip monomorphic variants
            if np.var(genotypes) == 0:
                continue
            
            # Calculate association statistics
            try:
                # Chi-square test for association
                contingency_table = self._create_contingency_table(genotypes, phenotype)
                chi2, p_value = stats.chi2_contingency(contingency_table)[:2]
                
                # Odds ratio calculation
                odds_ratio = self._calculate_odds_ratio(contingency_table)
                
                # Allele frequency
                allele_freq = np.mean(genotypes)
                
                # Effect size (beta)
                beta = np.log(odds_ratio) if odds_ratio > 0 else 0
                
                gwas_result = GWASResult(
                    variant_id=variant.get('variant_id', f"var_{variant['chr']}_{variant['pos']}"),
                    chromosome=variant['chr'],
                    position=variant['pos'],
                    p_value=p_value,
                    odds_ratio=odds_ratio,
                    beta=beta,
                    standard_error=0.1,  # Simplified
                    allele_frequency=allele_freq,
                    phenotype_association=phenotype_name
                )
                
                gwas_results.append(gwas_result)
                
            except Exception as e:
                logging.warning(f"Error processing variant {variant}: {e}")
                continue
        
        # Sort by p-value
        gwas_results.sort(key=lambda x: x.p_value)
        
        # Find significant variants
        significant_variants = [r for r in gwas_results if r.p_value < p_value_threshold]
        
        # Calculate genomic inflation factor
        observed_chi2 = [-2 * np.log(r.p_value) for r in gwas_results if r.p_value > 0]
        lambda_gc = np.median(observed_chi2) / 0.456 if observed_chi2 else 1.0
        
        return {
            "phenotype": phenotype_name,
            "total_variants": len(gwas_results),
            "significant_variants": len(significant_variants),
            "lambda_gc": lambda_gc,
            "top_associations": [vars(r) for r in gwas_results[:20]],
            "significant_associations": [vars(r) for r in significant_variants],
            "case_count": len(cases),
            "control_count": len(controls)
        }
    
    async def _extract_variant_data(self, individuals: List[str], 
                                  chromosome: Optional[int] = None) -> Tuple[np.ndarray, List[Dict]]:
        """Extract variant data for GWAS analysis"""
        
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
            RETURN n.individual_id as individual, n.chromosome as chr,
                   n.start_pos as pos, n.haplotype as haplotype,
                   n.variant_type as type, n.node_id as variant_id
            ORDER BY n.chromosome, n.start_pos
            """
            
            params = {"individuals": individuals}
            if chromosome:
                params["chromosome"] = chromosome
            
            results = session.run(query, **params).data()
        
        if not results:
            return np.array([]), []
        
        # Process results into genotype matrix
        df = pd.DataFrame(results)
        
        # Create variant position identifier
        df['variant_pos'] = df['chr'].astype(str) + '_' + df['pos'].astype(str)
        
        # Group by individual and variant position to get diploid genotypes
        genotype_data = []
        variant_positions = df['variant_pos'].unique()
        
        for individual in individuals:
            individual_genotypes = []
            for var_pos in variant_positions:
                # Get haplotypes for this individual and variant
                ind_var_data = df[(df['individual'] == individual) & (df['variant_pos'] == var_pos)]
                
                # Count alleles (simplified: presence = 1, absence = 0)
                allele_count = len(ind_var_data)
                individual_genotypes.append(min(allele_count, 2))  # Cap at 2 for diploid
            
            genotype_data.append(individual_genotypes)
        
        variant_matrix = np.array(genotype_data)
        
        # Get variant information
        variant_info = df.groupby('variant_pos').first()[['chr', 'pos', 'type', 'variant_id']].to_dict('records')
        
        return variant_matrix, variant_info
    
    def _create_contingency_table(self, genotypes: np.ndarray, phenotype: np.ndarray) -> np.ndarray:
        """Create contingency table for chi-square test"""
        
        # Convert genotypes to binary (0 = no variant, 1 = has variant)
        has_variant = (genotypes > 0).astype(int)
        
        # Create 2x2 contingency table
        case_with_variant = np.sum((phenotype == 1) & (has_variant == 1))
        case_without_variant = np.sum((phenotype == 1) & (has_variant == 0))
        control_with_variant = np.sum((phenotype == 0) & (has_variant == 1))
        control_without_variant = np.sum((phenotype == 0) & (has_variant == 0))
        
        return np.array([[case_with_variant, case_without_variant],
                        [control_with_variant, control_without_variant]])
    
    def _calculate_odds_ratio(self, contingency_table: np.ndarray) -> float:
        """Calculate odds ratio from contingency table"""
        
        a, b = contingency_table[0]  # cases with/without variant
        c, d = contingency_table[1]  # controls with/without variant
        
        # Add pseudocount to avoid division by zero
        if b == 0 or c == 0:
            a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5
        
        return (a * d) / (b * c) if (b * c) > 0 else 1.0

class PharmacogenomicsAnalyzer:
    """Pharmacogenomics analysis for personalized medicine"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.driver = db_builder.neo4j_manager.get_connection()
        
        # Pharmacogenomic knowledge base (simplified)
        self.pgx_knowledge = {
            "CYP2D6": {
                "drugs": ["codeine", "tramadol", "metoprolol", "risperidone"],
                "variants": {
                    "CYP2D6*1": "normal_metabolizer",
                    "CYP2D6*2": "normal_metabolizer", 
                    "CYP2D6*4": "poor_metabolizer",
                    "CYP2D6*10": "intermediate_metabolizer",
                    "CYP2D6*17": "ultrarapid_metabolizer"
                }
            },
            "CYP2C19": {
                "drugs": ["clopidogrel", "omeprazole", "escitalopram"],
                "variants": {
                    "CYP2C19*1": "normal_metabolizer",
                    "CYP2C19*2": "poor_metabolizer",
                    "CYP2C19*3": "poor_metabolizer",
                    "CYP2C19*17": "ultrarapid_metabolizer"
                }
            },
            "DPYD": {
                "drugs": ["5-fluorouracil", "capecitabine"],
                "variants": {
                    "DPYD*2A": "poor_metabolizer",
                    "DPYD*13": "poor_metabolizer"
                }
            }
        }
    
    async def generate_pharmacogenomic_profile(self, individual_id: str) -> PharmacogenomicProfile:
        """Generate comprehensive pharmacogenomic profile for an individual"""
        
        # Extract relevant variants for this individual
        pgx_variants = await self._extract_pgx_variants(individual_id)
        
        # Determine metabolizer status for each enzyme
        metabolizer_status = {}
        for enzyme, info in self.pgx_knowledge.items():
            status = self._determine_metabolizer_status(pgx_variants, enzyme, info["variants"])
            metabolizer_status[enzyme] = status
        
        # Predict drug responses
        drug_responses = {}
        recommendations = []
        
        for enzyme, info in self.pgx_knowledge.items():
            status = metabolizer_status[enzyme]
            
            for drug in info["drugs"]:
                response = self._predict_drug_response(drug, enzyme, status)
                drug_responses[drug] = response
                
                if response in ["poor_response", "toxicity_risk"]:
                    recommendations.append(f"Consider alternative to {drug} due to {enzyme} {status}")
                elif response == "reduced_efficacy":
                    recommendations.append(f"May need dose adjustment for {drug} due to {enzyme} {status}")
        
        # Find risk alleles
        risk_alleles = []
        for variant, annotation in pgx_variants.items():
            if annotation.get("clinical_significance") in ["pathogenic", "likely_pathogenic"]:
                risk_alleles.append(variant)
        
        return PharmacogenomicProfile(
            individual_id=individual_id,
            drug_responses=drug_responses,
            metabolizer_status=metabolizer_status,
            risk_alleles=risk_alleles,
            recommendations=recommendations
        )
    
    async def _extract_pgx_variants(self, individual_id: str) -> Dict[str, Dict]:
        """Extract pharmacogenomically relevant variants for an individual"""
        
        # Known pharmacogenes (simplified list)
        pharmacogenes = ["CYP2D6", "CYP2C19", "CYP2C9", "DPYD", "TPMT", "UGT1A1"]
        
        pgx_variants = {}
        
        with self.driver.session() as session:
            # This would typically involve querying known pharmacogene positions
            # For demo, we'll simulate finding relevant variants
            
            query = """
            MATCH (n:GenomeNode)
            WHERE n.individual_id = $individual_id
              AND n.scale_level = 0
              AND n.variant_type IS NOT NULL
            RETURN n.chromosome as chr, n.start_pos as pos, n.sequence as alt,
                   n.variant_type as type, n.gene_annotations as genes
            """
            
            results = session.run(query, individual_id=individual_id).data()
            
            # Simulate pharmacogene variant identification
            for i, result in enumerate(results[:20]):  # Limit for demo
                if i % 5 == 0:  # Simulate finding PGx variants
                    gene = np.random.choice(pharmacogenes)
                    variant_id = f"{gene}*{np.random.randint(1, 20)}"
                    
                    pgx_variants[variant_id] = {
                        "chromosome": result["chr"],
                        "position": result["pos"],
                        "gene": gene,
                        "variant_type": result["type"],
                        "clinical_significance": np.random.choice([
                            "normal", "reduced_function", "no_function", "increased_function"
                        ])
                    }
        
        return pgx_variants
    
    def _determine_metabolizer_status(self, pgx_variants: Dict, enzyme: str, 
                                    variant_status: Dict[str, str]) -> str:
        """Determine metabolizer status based on detected variants"""
        
        detected_variants = [v for v in pgx_variants.keys() if v.startswith(enzyme)]
        
        if not detected_variants:
            return "normal_metabolizer"  # Default assumption
        
        # Simplified logic: use the most severe status found
        statuses = []
        for variant in detected_variants:
            if variant in variant_status:
                statuses.append(variant_status[variant])
        
        if "poor_metabolizer" in statuses:
            return "poor_metabolizer"
        elif "intermediate_metabolizer" in statuses:
            return "intermediate_metabolizer"
        elif "ultrarapid_metabolizer" in statuses:
            return "ultrarapid_metabolizer"
        else:
            return "normal_metabolizer"
    
    def _predict_drug_response(self, drug: str, enzyme: str, status: str) -> str:
        """Predict drug response based on metabolizer status"""
        
        # Simplified prediction rules
        response_rules = {
            "poor_metabolizer": {
                "codeine": "poor_response",  # Can't convert to active metabolite
                "tramadol": "poor_response",
                "clopidogrel": "reduced_efficacy",  # Can't activate prodrug
                "5-fluorouracil": "toxicity_risk"  # Can't clear drug
            },
            "ultrarapid_metabolizer": {
                "codeine": "toxicity_risk",  # Too much active metabolite
                "metoprolol": "reduced_efficacy",  # Clears drug too quickly
                "omeprazole": "reduced_efficacy"
            },
            "intermediate_metabolizer": {
                "clopidogrel": "reduced_efficacy",
                "omeprazole": "normal_response"
            }
        }
        
        if status in response_rules and drug in response_rules[status]:
            return response_rules[status][drug]
        else:
            return "normal_response"

class ClinicalDecisionSupport:
    """Clinical decision support system"""
    
    def __init__(self, db_builder: GenomeDatabaseBuilder):
        self.db_builder = db_builder
        self.pgx_analyzer = PharmacogenomicsAnalyzer(db_builder)
        self.gwas_analyzer = GWASAnalyzer(db_builder)
    
    async def generate_clinical_report(self, individual_id: str, 
                                     indication: str = None) -> Dict[str, Any]:
        """Generate comprehensive clinical genomics report"""
        
        report = {
            "individual_id": individual_id,
            "report_date": datetime.now().isoformat(),
            "indication": indication,
            "sections": {}
        }
        
        try:
            # 1. Pharmacogenomic profile
            pgx_profile = await self.pgx_analyzer.generate_pharmacogenomic_profile(individual_id)
            report["sections"]["pharmacogenomics"] = {
                "metabolizer_status": pgx_profile.metabolizer_status,
                "drug_responses": pgx_profile.drug_responses,
                "recommendations": pgx_profile.recommendations,
                "risk_alleles": pgx_profile.risk_alleles
            }
            
            # 2. Disease risk assessment
            disease_risks = await self._assess_disease_risks(individual_id)
            report["sections"]["disease_risk"] = disease_risks
            
            # 3. Variant classification
            variant_summary = await self._classify_variants(individual_id)
            report["sections"]["variant_summary"] = variant_summary
            
            # 4. Clinical recommendations
            recommendations = await self._generate_clinical_recommendations(
                pgx_profile, disease_risks, indication
            )
            report["sections"]["clinical_recommendations"] = recommendations
            
            report["status"] = "complete"
            
        except Exception as e:
            report["status"] = "error"
            report["error"] = str(e)
            logging.error(f"Clinical report generation error: {e}")
        
        return report
    
    async def _assess_disease_risks(self, individual_id: str) -> Dict[str, Any]:
        """Assess disease risks based on known variants"""
        
        # This would typically involve:
        # 1. Querying known disease-associated variants
        # 2. Calculating polygenic risk scores
        # 3. Integrating family history
        
        # Simplified simulation
        diseases = ["coronary_artery_disease", "type_2_diabetes", "alzheimers", "breast_cancer"]
        
        disease_risks = {}
        for disease in diseases:
            # Simulate risk calculation
            risk_score = np.random.beta(2, 5)  # Skewed toward lower risk
            risk_category = "low" if risk_score < 0.3 else "moderate" if risk_score < 0.7 else "high"
            
            disease_risks[disease] = {
                "risk_score": round(risk_score, 3),
                "risk_category": risk_category,
                "population_percentile": round(risk_score * 100, 1),
                "contributing_variants": [],  # Would list specific variants
                "recommendations": self._get_disease_recommendations(disease, risk_category)
            }
        
        return disease_risks
    
    async def _classify_variants(self, individual_id: str) -> Dict[str, Any]:
        """Classify variants by clinical significance"""
        
        with self.driver.session() as session:
            query = """
            MATCH (n:GenomeNode)
            WHERE n.individual_id = $individual_id
              AND n.scale_level = 0
              AND n.variant_type IS NOT NULL
            RETURN n.variant_type as type, n.clinical_annotations as clinical,
                   count(n) as count
            """
            
            results = session.run(query, individual_id=individual_id).data()
        
        classification = {
            "total_variants": sum(r["count"] for r in results),
            "by_type": {r["type"]: r["count"] for r in results},
            "clinical_significance": {
                "pathogenic": 0,
                "likely_pathogenic": 0,
                "uncertain_significance": 0,
                "likely_benign": 0,
                "benign": 0
            }
        }
        
        # Simulate clinical significance distribution
        for result in results:
            if result.get("clinical"):
                # Would parse actual clinical annotations
                pass
        
        # Simulate some pathogenic variants
        classification["clinical_significance"]["pathogenic"] = np.random.poisson(2)
        classification["clinical_significance"]["likely_pathogenic"] = np.random.poisson(5)
        classification["clinical_significance"]["uncertain_significance"] = np.random.poisson(20)
        
        return classification
    
    async def _generate_clinical_recommendations(self, pgx_profile: PharmacogenomicProfile,
                                               disease_risks: Dict, indication: str) -> List[str]:
        """Generate clinical recommendations"""
        
        recommendations = []
        
        # Pharmacogenomic recommendations
        recommendations.extend(pgx_profile.recommendations)
        
        # Disease risk recommendations
        for disease, risk_data in disease_risks.items():
            if risk_data["risk_category"] == "high":
                recommendations.extend(risk_data["recommendations"])
        
        # Indication-specific recommendations
        if indication:
            recommendations.extend(self._get_indication_recommendations(indication, pgx_profile))
        
        # General recommendations
        recommendations.extend([
            "Consider genetic counseling for family planning",
            "Regular monitoring recommended for high-risk conditions",
            "Maintain updated medication list for drug interactions"
        ])
        
        return list(set(recommendations))  # Remove duplicates
    
    def _get_disease_recommendations(self, disease: str, risk_category: str) -> List[str]:
        """Get disease-specific recommendations"""
        
        recommendations_map = {
            "coronary_artery_disease": {
                "high": ["Cardiology consultation", "Lipid monitoring", "Lifestyle modifications"],
                "moderate": ["Regular BP monitoring", "Exercise program"],
                "low": ["Standard preventive care"]
            },
            "type_2_diabetes": {
                "high": ["Endocrinology referral", "HbA1c monitoring", "Diet counseling"],
                "moderate": ["Annual glucose screening", "Weight management"],
                "low": ["Standard screening"]
            }
        }
        
        return recommendations_map.get(disease, {}).get(risk_category, [])
    
    def _get_indication_recommendations(self, indication: str, 
                                      pgx_profile: PharmacogenomicProfile) -> List[str]:
        """Get indication-specific recommendations"""
        
        indication_map = {
            "depression": [
                "Consider pharmacogenomic testing before prescribing antidepressants",
                "Monitor for drug interactions with CYP2D6 substrates"
            ],
            "cancer": [
                "Tumor genomic profiling recommended",
                "Consider hereditary cancer syndrome testing"
            ],
            "cardiovascular": [
                "Clopidogrel effectiveness testing if CYP2C19 variants present",
                "Warfarin dosing guided by CYP2C9/VKORC1 status"
            ]
        }
        
        return indication_map.get(indication, [])

# Visualization utilities
class GenomicsVisualizer:
    """Create publication-quality genomics visualizations"""
    
    @staticmethod
    def create_manhattan_plot(gwas_results: List[GWASResult], 
                            title: str = "GWAS Manhattan Plot") -> go.Figure:
        """Create Manhattan plot for GWAS results"""
        
        # Prepare data
        chromosomes = [r.chromosome for r in gwas_results]
        positions = [r.position for r in gwas_results]
        p_values = [-np.log10(r.p_value) if r.p_value > 0 else 50 for r in gwas_results]
        
        # Create cumulative positions for x-axis
        cum_pos = []
        cum_offset = 0
        chr_positions = {}
        
        for chr_num in sorted(set(chromosomes)):
            chr_positions[chr_num] = cum_offset
            chr_data = [(pos, -np.log10(p)) for pos, p, chr in 
                       zip(positions, [r.p_value for r in gwas_results], chromosomes) 
                       if chr == chr_num]
            
            if chr_data:
                for pos, log_p in chr_data:
                    cum_pos.append(pos + cum_offset)
                cum_offset += max(pos for pos, _ in chr_data) + 1000000
        
        # Create plot
        fig = go.Figure()
        
        # Add points colored by chromosome
        colors = px.colors.qualitative.Set3
        for i, chr_num in enumerate(sorted(set(chromosomes))):
            chr_indices = [j for j, c in enumerate(chromosomes) if c == chr_num]
            chr_p_values = [p_values[j] for j in chr_indices]
            chr_cum_pos = [cum_pos[j] for j in chr_indices]
            
            fig.add_trace(go.Scatter(
                x=chr_cum_pos,
                y=chr_p_values,
                mode='markers',
                marker=dict(size=4, color=colors[i % len(colors)]),
                name=f'Chr {chr_num}',
                hovertemplate='Chr %{text}<br>Position: %{x}<br>-log10(p): %{y}<extra></extra>',
                text=[chr_num] * len(chr_indices)
            ))
        
        # Add significance threshold
        fig.add_hline(y=-np.log10(5e-8), line_dash="dash", line_color="red",
                     annotation_text="Genome-wide significance (p=5e-8)")
        
        fig.update_layout(
            title=title,
            xaxis_title="Chromosome",
            yaxis_title="-log10(p-value)",
            showlegend=False,
            width=1200,
            height=600
        )
        
        return fig
    
    @staticmethod
    def create_pca_plot(pca_coords: np.ndarray, 
                       cluster_labels: np.ndarray,
                       individuals: List[str],
                       title: str = "Population Structure PCA") -> go.Figure:
        """Create PCA plot for population structure"""
        
        fig = go.Figure()
        
        # Color by cluster
        colors = px.colors.qualitative.Set1
        for cluster in set(cluster_labels):
            cluster_mask = cluster_labels == cluster
            fig.add_trace(go.Scatter(
                x=pca_coords[cluster_mask, 0],
                y=pca_coords[cluster_mask, 1],
                mode='markers',
                marker=dict(size=8, color=colors[cluster % len(colors)]),
                name=f'Cluster {cluster}',
                text=[individuals[i] for i in range(len(individuals)) if cluster_mask[i]],
                hovertemplate='%{text}<br>PC1: %{x:.3f}<br>PC2: %{y:.3f}<extra></extra>'
            ))
        
        fig.update_layout(
            title=title,
            xaxis_title="Principal Component 1",
            yaxis_title="Principal Component 2",
            width=800,
            height=600
        )
        
        return fig

# Example usage and demonstration
async def demonstrate_advanced_analytics():
    """Demonstrate advanced genomics analytics"""
    
    print("🧬 Advanced Genomics Analytics Demo")
    print("=" * 50)
    
    try:
        # Initialize components
        db_builder = GenomeDatabaseBuilder()
        await db_builder.initialize_database()
        
        # 1. Population genomics analysis
        print("\n📊 Population Genomics Analysis")
        print("-" * 30)
        
        pop_analyzer = PopulationGenomicsAnalyzer(db_builder)
        
        # Simulate population data
        populations = {
            "European": ["HG002", "HG00096", "HG00097"],
            "African": ["HG00733", "NA19240", "NA18502"],
            "Asian": ["HG00514", "HG00419", "HG00403"]
        }
        
        diversity_results = await pop_analyzer.analyze_population_diversity(populations, chromosome=1)
        print(f"✅ Population diversity analysis completed")
        print(f"   Populations analyzed: {len(populations)}")
        
        # 2. GWAS analysis
        print("\n🔬 GWAS Analysis")
        print("-" * 30)
        
        gwas_analyzer = GWASAnalyzer(db_builder)
        
        # Simulate case-control study
        cases = ["HG002", "HG00096"]  # Simulated cases
        controls = ["HG00733", "HG00097"]  # Simulated controls
        
        gwas_results = await gwas_analyzer.run_gwas_analysis(
            cases, controls, "hypertension", chromosome=1
        )
        print(f"✅ GWAS analysis completed")
        print(f"   Variants tested: {gwas_results['total_variants']}")
        print(f"   Significant hits: {gwas_results['significant_variants']}")
        
        # 3. Pharmacogenomics analysis
        print("\n💊 Pharmacogenomics Analysis")
        print("-" * 30)
        
        pgx_analyzer = PharmacogenomicsAnalyzer(db_builder)
        pgx_profile = await pgx_analyzer.generate_pharmacogenomic_profile("HG002")
        
        print(f"✅ Pharmacogenomic profile generated for HG002")
        print(f"   Drug responses: {len(pgx_profile.drug_responses)}")
        print(f"   Recommendations: {len(pgx_profile.recommendations)}")
        
        # 4. Clinical decision support
        print("\n🏥 Clinical Decision Support")
        print("-" * 30)
        
        clinical_support = ClinicalDecisionSupport(db_builder)
        clinical_report = await clinical_support.generate_clinical_report(
            "HG002", indication="cardiovascular"
        )
        
        print(f"✅ Clinical report generated")
        print(f"   Report sections: {len(clinical_report['sections'])}")
        print(f"   Status: {clinical_report['status']}")
        
        # 5. Create visualizations
        print("\n📈 Creating Visualizations")
        print("-" * 30)
        
        # Save results
        results_dir = Path("./analytics_results")
        results_dir.mkdir(exist_ok=True)
        
        # Save comprehensive results
        with open(results_dir / "analytics_results.json", 'w') as f:
            json.dump({
                "population_diversity": diversity_results,
                "gwas_results": gwas_results,
                "pharmacogenomics": vars(pgx_profile),
                "clinical_report": clinical_report
            }, f, indent=2, default=str)
        
        print(f"📁 Results saved to: {results_dir}")
        
    except Exception as e:
        print(f"❌ Analytics demo error: {e}")
        logging.error(f"Analytics demo error: {e}")

if __name__ == "__main__":
    # Configure logging
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s'
    )
    
    # Run demonstration
    asyncio.run(demonstrate_advanced_analytics())
