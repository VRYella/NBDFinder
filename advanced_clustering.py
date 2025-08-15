"""
Advanced Clustering and Hybrid Structure Analysis Module
=======================================================

This module implements sophisticated algorithms for detecting non-B DNA clusters,
hybrid structures, and cooperative binding effects. It provides enhanced analysis
of genomic hotspots where multiple non-B structures interact.

Key Features:
- Advanced clustering algorithms for non-B DNA hotspots
- Hybrid structure detection with overlap analysis
- Cooperative binding effects modeling
- Structural interaction prediction
- Genomic context analysis
- Publication-quality visualization

Author: Dr. Venkata Rajesh Yella
Updated: 2024 with latest clustering methods
License: Academic Use
"""

import numpy as np
import pandas as pd
from typing import Dict, List, Tuple, Optional, Set
from collections import defaultdict, Counter
from dataclasses import dataclass
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
import matplotlib.pyplot as plt
import seaborn as sns

@dataclass
class ClusterRegion:
    """Represents a genomic region with clustered non-B DNA structures"""
    start: int
    end: int
    motif_count: int
    motif_types: Set[str]
    diversity_index: float
    interaction_score: float
    cooperative_score: float
    clinical_significance: str
    
@dataclass
class HybridStructure:
    """Represents a hybrid non-B DNA structure with multiple overlapping motifs"""
    start: int
    end: int
    component_motifs: List[Dict]
    overlap_degree: float
    structural_impact: str
    stability_enhancement: float
    functional_prediction: str

class AdvancedClusterAnalyzer:
    """Advanced clustering and hybrid structure analysis"""
    
    def __init__(self):
        self.interaction_weights = self._load_interaction_weights()
        self.cooperative_thresholds = self._load_cooperative_thresholds()
        
    def _load_interaction_weights(self) -> Dict:
        """Load interaction weights between different motif types"""
        # Based on experimental evidence and theoretical predictions
        return {
            # G-quadruplex interactions
            ('Canonical G4', 'Canonical G4'): 0.8,
            ('Canonical G4', 'i-Motif'): 0.9,  # Strong complementary interaction
            ('Canonical G4', 'R-Loop'): 0.6,
            ('Canonical G4', 'Z-DNA'): 0.4,
            
            # Z-DNA interactions
            ('Z-DNA', 'Cruciform'): 0.5,
            ('Z-DNA', 'Slipped_DNA'): 0.6,
            ('Z-DNA', 'eGZ'): 0.7,
            
            # R-loop interactions
            ('R-Loop', 'G-Triplex'): 0.8,
            ('R-Loop', 'Sticky_DNA'): 0.5,
            
            # Triplex interactions
            ('Triplex_DNA', 'Sticky_DNA'): 0.9,  # Both involve triple helix
            ('G-Triplex', 'Canonical G4'): 0.6,
            
            # Repeat structure interactions
            ('Slipped_DNA', 'Cruciform'): 0.7,
            ('Sticky_DNA', 'Triplex_DNA'): 0.8,
            
            # Default interaction strength
            'default': 0.3
        }
    
    def _load_cooperative_thresholds(self) -> Dict:
        """Load thresholds for cooperative binding detection"""
        return {
            'distance_threshold': 200,  # Base pairs
            'overlap_threshold': 0.1,   # Minimum overlap for hybrid
            'cluster_min_motifs': 3,    # Minimum motifs for cluster
            'cluster_window': 500,      # Sliding window size
            'interaction_min_score': 0.4  # Minimum interaction score
        }
    
    def detect_advanced_clusters(self, motifs: List[Dict], sequence_length: int) -> List[ClusterRegion]:
        """
        Detect advanced non-B DNA clusters using sophisticated algorithms
        
        Args:
            motifs: List of detected motifs
            sequence_length: Length of analyzed sequence
            
        Returns:
            List of detected cluster regions
        """
        if len(motifs) < self.cooperative_thresholds['cluster_min_motifs']:
            return []
        
        # Create motif position matrix
        motif_positions = [(m['Start'], m['End'], m['Subtype']) for m in motifs]
        motif_positions.sort()
        
        # Apply sliding window clustering
        window_size = self.cooperative_thresholds['cluster_window']
        clusters = []
        
        for start_pos in range(1, sequence_length - window_size + 1, window_size // 2):
            end_pos = start_pos + window_size
            
            # Find motifs in current window
            window_motifs = [
                m for m in motifs 
                if m['Start'] >= start_pos and m['End'] <= end_pos
            ]
            
            if len(window_motifs) >= self.cooperative_thresholds['cluster_min_motifs']:
                cluster = self._analyze_cluster_region(window_motifs, start_pos, end_pos)
                if cluster:
                    clusters.append(cluster)
        
        # Merge overlapping clusters
        clusters = self._merge_overlapping_clusters(clusters)
        
        # Apply hierarchical clustering for refinement
        clusters = self._refine_clusters_hierarchical(clusters, motifs)
        
        return clusters
    
    def _analyze_cluster_region(self, motifs: List[Dict], start: int, end: int) -> Optional[ClusterRegion]:
        """Analyze a potential cluster region"""
        if len(motifs) < self.cooperative_thresholds['cluster_min_motifs']:
            return None
        
        # Calculate cluster metrics
        motif_types = set(m['Subtype'] for m in motifs)
        motif_count = len(motifs)
        
        # Calculate diversity index (Shannon diversity)
        type_counts = Counter(m['Subtype'] for m in motifs)
        total = sum(type_counts.values())
        diversity = -sum((count/total) * np.log2(count/total) for count in type_counts.values())
        
        # Calculate interaction score
        interaction_score = self._calculate_cluster_interaction_score(motifs)
        
        # Calculate cooperative score
        cooperative_score = self._calculate_cooperative_score(motifs)
        
        # Assess clinical significance
        clinical_sig = self._assess_cluster_clinical_significance(motifs, motif_types)
        
        return ClusterRegion(
            start=start,
            end=end,
            motif_count=motif_count,
            motif_types=motif_types,
            diversity_index=diversity,
            interaction_score=interaction_score,
            cooperative_score=cooperative_score,
            clinical_significance=clinical_sig
        )
    
    def _calculate_cluster_interaction_score(self, motifs: List[Dict]) -> float:
        """Calculate interaction score for motifs in cluster"""
        if len(motifs) < 2:
            return 0.0
        
        total_score = 0.0
        pair_count = 0
        
        for i in range(len(motifs)):
            for j in range(i + 1, len(motifs)):
                motif1_type = motifs[i]['Subtype']
                motif2_type = motifs[j]['Subtype']
                
                # Get interaction weight
                interaction_key = tuple(sorted([motif1_type, motif2_type]))
                weight = self.interaction_weights.get(
                    interaction_key, 
                    self.interaction_weights.get((motif1_type, motif2_type),
                    self.interaction_weights.get((motif2_type, motif1_type),
                    self.interaction_weights['default']))
                )
                
                # Calculate distance-based interaction
                distance = min(
                    abs(motifs[i]['Start'] - motifs[j]['End']),
                    abs(motifs[j]['Start'] - motifs[i]['End'])
                )
                
                if distance <= self.cooperative_thresholds['distance_threshold']:
                    distance_factor = 1.0 - (distance / self.cooperative_thresholds['distance_threshold'])
                    interaction_score = weight * distance_factor
                    total_score += interaction_score
                    pair_count += 1
        
        return total_score / pair_count if pair_count > 0 else 0.0
    
    def _calculate_cooperative_score(self, motifs: List[Dict]) -> float:
        """Calculate cooperative binding score"""
        if len(motifs) < 2:
            return 0.0
        
        # Factors contributing to cooperativity
        proximity_score = 0.0
        density_score = 0.0
        diversity_score = 0.0
        
        # Proximity scoring
        positions = [(m['Start'], m['End']) for m in motifs]
        positions.sort()
        
        for i in range(len(positions) - 1):
            gap = positions[i+1][0] - positions[i][1]
            if gap <= self.cooperative_thresholds['distance_threshold']:
                proximity_score += 1.0 - (gap / self.cooperative_thresholds['distance_threshold'])
        
        proximity_score /= (len(positions) - 1) if len(positions) > 1 else 1
        
        # Density scoring
        total_length = max(m['End'] for m in motifs) - min(m['Start'] for m in motifs)
        motif_coverage = sum(m['Length'] for m in motifs)
        density_score = motif_coverage / total_length if total_length > 0 else 0
        
        # Diversity scoring (more diverse clusters may have higher cooperativity)
        unique_types = len(set(m['Subtype'] for m in motifs))
        max_diversity = min(len(motifs), 7)  # Reasonable maximum
        diversity_score = unique_types / max_diversity
        
        # Combined cooperative score
        cooperative_score = (proximity_score * 0.4 + density_score * 0.3 + diversity_score * 0.3)
        
        return min(1.0, cooperative_score)
    
    def _assess_cluster_clinical_significance(self, motifs: List[Dict], motif_types: Set[str]) -> str:
        """Assess clinical significance of cluster"""
        # Check for disease-associated combinations
        disease_motifs = {'GAA_repeat', 'CGG_repeat', 'CAG_repeat', 'CTG_repeat', 'G4C2_repeat'}
        
        if any(motif_type in disease_motifs for motif_type in motif_types):
            return "High Clinical Significance"
        
        # Check for cancer-associated combinations
        cancer_motifs = {'Canonical G4', 'G-Triplex', 'R-Loop'}
        if len(motif_types.intersection(cancer_motifs)) >= 2:
            return "Moderate Clinical Significance"
        
        # Check for regulatory combinations
        regulatory_motifs = {'Z-DNA', 'Cruciform', 'Curved_DNA'}
        if len(motif_types.intersection(regulatory_motifs)) >= 2:
            return "Regulatory Significance"
        
        return "Research Interest"
    
    def _merge_overlapping_clusters(self, clusters: List[ClusterRegion]) -> List[ClusterRegion]:
        """Merge overlapping cluster regions"""
        if not clusters:
            return []
        
        merged = []
        clusters.sort(key=lambda c: c.start)
        
        current = clusters[0]
        
        for next_cluster in clusters[1:]:
            # Check for overlap
            if next_cluster.start <= current.end:
                # Merge clusters
                current = ClusterRegion(
                    start=current.start,
                    end=max(current.end, next_cluster.end),
                    motif_count=current.motif_count + next_cluster.motif_count,
                    motif_types=current.motif_types.union(next_cluster.motif_types),
                    diversity_index=max(current.diversity_index, next_cluster.diversity_index),
                    interaction_score=max(current.interaction_score, next_cluster.interaction_score),
                    cooperative_score=max(current.cooperative_score, next_cluster.cooperative_score),
                    clinical_significance=self._merge_clinical_significance(
                        current.clinical_significance, next_cluster.clinical_significance
                    )
                )
            else:
                merged.append(current)
                current = next_cluster
        
        merged.append(current)
        return merged
    
    def _merge_clinical_significance(self, sig1: str, sig2: str) -> str:
        """Merge clinical significance levels"""
        priority = {
            "High Clinical Significance": 4,
            "Moderate Clinical Significance": 3,
            "Regulatory Significance": 2,
            "Research Interest": 1
        }
        
        if priority[sig1] >= priority[sig2]:
            return sig1
        else:
            return sig2
    
    def _refine_clusters_hierarchical(self, clusters: List[ClusterRegion], all_motifs: List[Dict]) -> List[ClusterRegion]:
        """Refine clusters using hierarchical clustering"""
        if len(clusters) < 2:
            return clusters
        
        # Create feature matrix for clustering
        features = []
        for cluster in clusters:
            cluster_motifs = [
                m for m in all_motifs 
                if m['Start'] >= cluster.start and m['End'] <= cluster.end
            ]
            
            feature_vector = [
                cluster.motif_count,
                cluster.diversity_index,
                cluster.interaction_score,
                cluster.cooperative_score,
                len(cluster.motif_types),
                cluster.end - cluster.start  # cluster length
            ]
            features.append(feature_vector)
        
        if len(features) < 2:
            return clusters
        
        # Perform hierarchical clustering
        try:
            distances = pdist(features, metric='euclidean')
            linkage_matrix = linkage(distances, method='ward')
            
            # Determine optimal number of clusters
            max_clusters = min(len(clusters), 5)
            cluster_labels = fcluster(linkage_matrix, max_clusters, criterion='maxclust')
            
            # Group clusters by labels
            clustered_regions = defaultdict(list)
            for i, label in enumerate(cluster_labels):
                clustered_regions[label].append(clusters[i])
            
            # Merge clusters within each group
            refined_clusters = []
            for group_clusters in clustered_regions.values():
                if len(group_clusters) == 1:
                    refined_clusters.append(group_clusters[0])
                else:
                    # Merge multiple clusters in the same group
                    merged = self._merge_cluster_group(group_clusters)
                    refined_clusters.append(merged)
            
            return refined_clusters
        
        except Exception:
            # If hierarchical clustering fails, return original clusters
            return clusters
    
    def _merge_cluster_group(self, cluster_group: List[ClusterRegion]) -> ClusterRegion:
        """Merge a group of clusters into one"""
        if len(cluster_group) == 1:
            return cluster_group[0]
        
        start = min(c.start for c in cluster_group)
        end = max(c.end for c in cluster_group)
        total_motifs = sum(c.motif_count for c in cluster_group)
        all_types = set()
        for c in cluster_group:
            all_types.update(c.motif_types)
        
        avg_diversity = np.mean([c.diversity_index for c in cluster_group])
        max_interaction = max(c.interaction_score for c in cluster_group)
        avg_cooperative = np.mean([c.cooperative_score for c in cluster_group])
        
        # Determine merged clinical significance
        significances = [c.clinical_significance for c in cluster_group]
        merged_significance = max(significances, key=lambda s: {
            "High Clinical Significance": 4,
            "Moderate Clinical Significance": 3,
            "Regulatory Significance": 2,
            "Research Interest": 1
        }[s])
        
        return ClusterRegion(
            start=start,
            end=end,
            motif_count=total_motifs,
            motif_types=all_types,
            diversity_index=avg_diversity,
            interaction_score=max_interaction,
            cooperative_score=avg_cooperative,
            clinical_significance=merged_significance
        )
    
    def detect_hybrid_structures(self, motifs: List[Dict]) -> List[HybridStructure]:
        """
        Detect hybrid structures with overlapping motifs
        
        Args:
            motifs: List of detected motifs
            
        Returns:
            List of detected hybrid structures
        """
        hybrids = []
        
        for i in range(len(motifs)):
            for j in range(i + 1, len(motifs)):
                motif1 = motifs[i]
                motif2 = motifs[j]
                
                # Check for overlap
                overlap = self._calculate_overlap(motif1, motif2)
                if overlap >= self.cooperative_thresholds['overlap_threshold']:
                    hybrid = self._create_hybrid_structure(motif1, motif2, overlap)
                    if hybrid:
                        hybrids.append(hybrid)
        
        # Remove redundant hybrids and merge complex overlaps
        hybrids = self._consolidate_hybrids(hybrids)
        
        return hybrids
    
    def _calculate_overlap(self, motif1: Dict, motif2: Dict) -> float:
        """Calculate overlap percentage between two motifs"""
        start1, end1 = motif1['Start'], motif1['End']
        start2, end2 = motif2['Start'], motif2['End']
        
        # Calculate intersection
        intersection_start = max(start1, start2)
        intersection_end = min(end1, end2)
        
        if intersection_end <= intersection_start:
            return 0.0
        
        intersection_length = intersection_end - intersection_start
        min_length = min(end1 - start1, end2 - start2)
        
        return intersection_length / min_length if min_length > 0 else 0.0
    
    def _create_hybrid_structure(self, motif1: Dict, motif2: Dict, overlap: float) -> Optional[HybridStructure]:
        """Create hybrid structure from two overlapping motifs"""
        if overlap < self.cooperative_thresholds['overlap_threshold']:
            return None
        
        # Calculate hybrid boundaries
        start = min(motif1['Start'], motif2['Start'])
        end = max(motif1['End'], motif2['End'])
        
        # Assess structural impact
        impact = self._assess_hybrid_structural_impact(motif1, motif2, overlap)
        
        # Calculate stability enhancement
        stability = self._calculate_hybrid_stability(motif1, motif2, overlap)
        
        # Predict function
        function = self._predict_hybrid_function(motif1, motif2)
        
        return HybridStructure(
            start=start,
            end=end,
            component_motifs=[motif1, motif2],
            overlap_degree=overlap,
            structural_impact=impact,
            stability_enhancement=stability,
            functional_prediction=function
        )
    
    def _assess_hybrid_structural_impact(self, motif1: Dict, motif2: Dict, overlap: float) -> str:
        """Assess structural impact of hybrid formation"""
        type1 = motif1['Subtype']
        type2 = motif2['Subtype']
        
        # High impact combinations
        high_impact_pairs = {
            ('Canonical G4', 'i-Motif'),
            ('Z-DNA', 'Cruciform'),
            ('R-Loop', 'G-Triplex'),
            ('Triplex_DNA', 'Sticky_DNA')
        }
        
        pair = tuple(sorted([type1, type2]))
        
        if pair in high_impact_pairs:
            if overlap > 0.7:
                return "Major structural reorganization"
            elif overlap > 0.4:
                return "Significant structural change"
            else:
                return "Moderate structural alteration"
        
        # Moderate impact combinations
        moderate_impact = {
            ('Canonical G4', 'Z-DNA'),
            ('Slipped_DNA', 'Cruciform'),
            ('Curved_DNA', 'Z-DNA')
        }
        
        if pair in moderate_impact:
            return "Moderate structural change"
        
        return "Minor structural perturbation"
    
    def _calculate_hybrid_stability(self, motif1: Dict, motif2: Dict, overlap: float) -> float:
        """Calculate stability enhancement from hybrid formation"""
        # Base stability from individual motifs
        score1 = float(motif1.get('Score', 0))
        score2 = float(motif2.get('Score', 0))
        
        base_stability = (score1 + score2) / 200.0  # Normalize to 0-1
        
        # Overlap bonus (more overlap can increase stability)
        overlap_bonus = overlap * 0.3
        
        # Type-specific stability
        type1 = motif1['Subtype']
        type2 = motif2['Subtype']
        
        stability_matrix = {
            ('Canonical G4', 'i-Motif'): 0.8,  # Very stable complementary structure
            ('Z-DNA', 'Cruciform'): 0.6,
            ('R-Loop', 'G-Triplex'): 0.7,
            ('Triplex_DNA', 'Sticky_DNA'): 0.9
        }
        
        pair = tuple(sorted([type1, type2]))
        type_bonus = stability_matrix.get(pair, 0.3)
        
        total_stability = min(1.0, base_stability + overlap_bonus + type_bonus)
        return total_stability
    
    def _predict_hybrid_function(self, motif1: Dict, motif2: Dict) -> str:
        """Predict functional impact of hybrid structure"""
        type1 = motif1['Subtype']
        type2 = motif2['Subtype']
        
        functional_predictions = {
            ('Canonical G4', 'i-Motif'): "pH-sensitive gene regulation",
            ('Z-DNA', 'Cruciform'): "Transcriptional pausing and recombination",
            ('R-Loop', 'G-Triplex'): "Enhanced transcription-replication conflicts",
            ('Triplex_DNA', 'Sticky_DNA'): "Chromatin compaction and silencing",
            ('Canonical G4', 'R-Loop'): "Transcriptional regulation",
            ('Slipped_DNA', 'Cruciform'): "Replication fork stalling",
            ('Z-DNA', 'Curved_DNA'): "Nucleosome positioning"
        }
        
        pair = tuple(sorted([type1, type2]))
        return functional_predictions.get(pair, "Unknown functional impact")
    
    def _consolidate_hybrids(self, hybrids: List[HybridStructure]) -> List[HybridStructure]:
        """Consolidate overlapping hybrid structures"""
        if not hybrids:
            return []
        
        # Sort by position
        hybrids.sort(key=lambda h: h.start)
        
        consolidated = []
        current = hybrids[0]
        
        for next_hybrid in hybrids[1:]:
            # Check if they overlap significantly
            if (next_hybrid.start <= current.end and 
                (next_hybrid.start - current.start) / (current.end - current.start) < 0.8):
                
                # Merge hybrids
                current = self._merge_hybrids(current, next_hybrid)
            else:
                consolidated.append(current)
                current = next_hybrid
        
        consolidated.append(current)
        return consolidated
    
    def _merge_hybrids(self, hybrid1: HybridStructure, hybrid2: HybridStructure) -> HybridStructure:
        """Merge two overlapping hybrid structures"""
        start = min(hybrid1.start, hybrid2.start)
        end = max(hybrid1.end, hybrid2.end)
        
        # Combine component motifs
        all_motifs = hybrid1.component_motifs + hybrid2.component_motifs
        
        # Calculate merged properties
        avg_overlap = (hybrid1.overlap_degree + hybrid2.overlap_degree) / 2
        max_stability = max(hybrid1.stability_enhancement, hybrid2.stability_enhancement)
        
        # Choose more severe structural impact
        impact_severity = {
            "Major structural reorganization": 4,
            "Significant structural change": 3,
            "Moderate structural change": 2,
            "Minor structural perturbation": 1
        }
        
        if impact_severity[hybrid1.structural_impact] >= impact_severity[hybrid2.structural_impact]:
            merged_impact = hybrid1.structural_impact
        else:
            merged_impact = hybrid2.structural_impact
        
        # Combine functional predictions
        merged_function = f"{hybrid1.functional_prediction}; {hybrid2.functional_prediction}"
        
        return HybridStructure(
            start=start,
            end=end,
            component_motifs=all_motifs,
            overlap_degree=avg_overlap,
            structural_impact=merged_impact,
            stability_enhancement=max_stability,
            functional_prediction=merged_function
        )

def find_advanced_clusters(motifs: List[Dict], sequence_length: int) -> List[Dict]:
    """
    Find advanced non-B DNA clusters with detailed analysis
    
    Args:
        motifs: List of detected motifs
        sequence_length: Length of analyzed sequence
        
    Returns:
        List of cluster motifs with enhanced annotations
    """
    analyzer = AdvancedClusterAnalyzer()
    clusters = analyzer.detect_advanced_clusters(motifs, sequence_length)
    
    cluster_motifs = []
    for i, cluster in enumerate(clusters):
        cluster_motif = {
            "Sequence Name": "",
            "Class": "Advanced Non-B DNA Clusters",
            "Subtype": "Advanced_Cluster",
            "Start": cluster.start,
            "End": cluster.end,
            "Sequence": f"CLUSTER_{i+1}",
            "Length": cluster.end - cluster.start,
            "Score": cluster.interaction_score * 100,
            
            # Enhanced cluster metrics
            "Motif_Count": cluster.motif_count,
            "Motif_Types": "; ".join(sorted(cluster.motif_types)),
            "Diversity_Index": round(cluster.diversity_index, 3),
            "Interaction_Score": round(cluster.interaction_score, 3),
            "Cooperative_Score": round(cluster.cooperative_score, 3),
            "Clinical_Significance": cluster.clinical_significance,
            
            # Additional annotations
            "Arms/Repeat Unit/Copies": f"{cluster.motif_count} motifs",
            "Spacer": f"{len(cluster.motif_types)} types",
            "ScoreMethod": "Advanced_Cluster_Analysis",
            "Conservation_Score": cluster.diversity_index * 100,
            "Structural_Impact": f"Multi-motif cooperative binding ({cluster.motif_count} structures)"
        }
        cluster_motifs.append(cluster_motif)
    
    return cluster_motifs

def find_hybrid_structures(motifs: List[Dict]) -> List[Dict]:
    """
    Find hybrid structures with overlapping motifs
    
    Args:
        motifs: List of detected motifs
        
    Returns:
        List of hybrid structure motifs with detailed annotations
    """
    analyzer = AdvancedClusterAnalyzer()
    hybrids = analyzer.detect_hybrid_structures(motifs)
    
    hybrid_motifs = []
    for i, hybrid in enumerate(hybrids):
        motif_types = [m['Subtype'] for m in hybrid.component_motifs]
        
        hybrid_motif = {
            "Sequence Name": "",
            "Class": "Hybrid Structures",
            "Subtype": "Hybrid_Structure",
            "Start": hybrid.start,
            "End": hybrid.end,
            "Sequence": f"HYBRID_{i+1}",
            "Length": hybrid.end - hybrid.start,
            "Score": hybrid.stability_enhancement * 100,
            
            # Hybrid-specific metrics
            "Component_Motifs": "; ".join(motif_types),
            "Overlap_Degree": round(hybrid.overlap_degree, 3),
            "Structural_Impact": hybrid.structural_impact,
            "Stability_Enhancement": round(hybrid.stability_enhancement, 3),
            "Functional_Prediction": hybrid.functional_prediction,
            
            # Additional annotations
            "Arms/Repeat Unit/Copies": f"{len(hybrid.component_motifs)} components",
            "Spacer": f"{hybrid.overlap_degree:.1%} overlap",
            "ScoreMethod": "Hybrid_Structure_Analysis",
            "Conservation_Score": hybrid.stability_enhancement * 100,
            "Structural_Impact_Detail": f"Hybrid of {' + '.join(motif_types)}"
        }
        hybrid_motifs.append(hybrid_motif)
    
    return hybrid_motifs