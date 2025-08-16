"""
Modular Scoring System for Non-B DNA Motifs
==========================================

This module provides a flexible, combinable scoring system that integrates:
1. Sequence pattern scores
2. Structural feature scores  
3. Evolutionary conservation metrics
4. Thermodynamic stability estimates

Scientific Basis:
Scoring functions are based on experimental validation studies and 
structural biology data for each motif type.

References:
- Bedrat et al. (2016) NAR (G4Hunter scoring)
- Ho et al. (1986) Nucleic Acids Res (Z-DNA energetics)
- Aguilera & García-Muse (2012) Mol Cell (R-loop thermodynamics)
- Huppert & Balasubramanian (2005) NAR (conservation scoring)

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import math
import random
from typing import Dict, List, Tuple, Optional, Any, Callable
from dataclasses import dataclass
from abc import ABC, abstractmethod

@dataclass
class ScoreComponents:
    """Container for individual scoring components."""
    pattern_score: float = 0.0
    structural_score: float = 0.0
    conservation_score: float = 0.0
    thermodynamic_score: float = 0.0
    context_score: float = 0.0
    
    def total(self, weights: Optional[Dict[str, float]] = None) -> float:
        """Calculate weighted total score."""
        if weights is None:
            weights = {
                'pattern': 1.0,
                'structural': 1.0, 
                'conservation': 0.5,
                'thermodynamic': 0.8,
                'context': 0.3
            }
        
        return (
            self.pattern_score * weights.get('pattern', 1.0) +
            self.structural_score * weights.get('structural', 1.0) +
            self.conservation_score * weights.get('conservation', 0.5) +
            self.thermodynamic_score * weights.get('thermodynamic', 0.8) +
            self.context_score * weights.get('context', 0.3)
        )

class BaseScorer(ABC):
    """Abstract base class for all scoring functions."""
    
    @abstractmethod
    def score(self, sequence: str, **kwargs) -> float:
        """Calculate score for given sequence."""
        pass
    
    @abstractmethod
    def get_score_range(self) -> Tuple[float, float]:
        """Get theoretical min/max score range."""
        pass

class G4HunterScorer(BaseScorer):
    """
    G4Hunter scoring algorithm with structural factors.
    
    Scientific Basis: G4Hunter assigns scores based on G/C run lengths
    and averages over sequence. Structural factors account for loop
    constraints and G-tetrad stacking energies.
    
    Reference: Bedrat et al. (2016) NAR
    """
    
    def __init__(self, max_score: int = 4, include_structural: bool = True):
        self.max_score = max_score
        self.include_structural = include_structural
    
    def score(self, sequence: str, **kwargs) -> float:
        """Calculate G4Hunter score with optional structural factors."""
        if not sequence:
            return 0.0
        
        sequence = sequence.upper()
        scores = []
        i = 0
        
        while i < len(sequence):
            if sequence[i] == 'G':
                # Count consecutive Gs
                run_length = 0
                while i < len(sequence) and sequence[i] == 'G':
                    run_length += 1
                    i += 1
                # Score each G in the run
                run_score = min(run_length, self.max_score)
                scores.extend([run_score] * run_length)
                
            elif sequence[i] == 'C':
                # Count consecutive Cs (negative scores)
                run_length = 0
                while i < len(sequence) and sequence[i] == 'C':
                    run_length += 1
                    i += 1
                # Score each C in the run (negative)
                run_score = -min(run_length, self.max_score)
                scores.extend([run_score] * run_length)
                
            else:
                # A/T get score 0
                scores.append(0.0)
                i += 1
        
        base_score = sum(scores) / len(scores) if scores else 0.0
        
        if self.include_structural:
            structural_factor = self._calculate_structural_factor(sequence)
            return base_score * structural_factor
        
        return base_score
    
    def _calculate_structural_factor(self, sequence: str) -> float:
        """
        Calculate structural factor based on G4 folding constraints.
        
        Factors considered:
        - Loop length distribution
        - G-tetrad stacking potential
        - Sequence context effects
        """
        # Simple structural factor based on loop constraints
        import re
        
        # Find G-runs and loops
        g_runs = [(m.start(), m.end()) for m in re.finditer(r'G{3,}', sequence)]
        
        if len(g_runs) < 4:
            return 0.8  # Reduced factor for insufficient G-runs
        
        # Calculate loop lengths between first 4 G-runs
        loops = []
        for i in range(min(3, len(g_runs) - 1)):
            loop_start = g_runs[i][1]
            loop_end = g_runs[i + 1][0]
            loop_length = loop_end - loop_start
            loops.append(loop_length)
        
        # Optimal loop lengths: 1-7 nt (penalty for very long loops)
        loop_factor = 1.0
        for loop_len in loops:
            if loop_len > 7:
                loop_factor *= 0.9  # Penalty for long loops
            elif loop_len == 0:
                loop_factor *= 0.7  # Penalty for no loop
        
        return max(0.5, min(1.2, loop_factor))
    
    def get_score_range(self) -> Tuple[float, float]:
        """G4Hunter scores typically range from -4 to +4."""
        return (-self.max_score, self.max_score)

class iMotifScorer(BaseScorer):
    """
    i-motif scoring using reverse G4Hunter logic.
    
    Scientific Basis: i-motifs are the C-rich complement to G4s,
    forming at acidic pH through C-C+ base pairs.
    """
    
    def __init__(self, max_score: int = 4):
        self.max_score = max_score
    
    def score(self, sequence: str, **kwargs) -> float:
        """Calculate i-motif score (reverse G4Hunter logic)."""
        if not sequence:
            return 0.0
        
        sequence = sequence.upper()
        scores = []
        i = 0
        
        while i < len(sequence):
            if sequence[i] == 'C':
                # Count consecutive Cs (positive for i-motif)
                run_length = 0
                while i < len(sequence) and sequence[i] == 'C':
                    run_length += 1
                    i += 1
                run_score = min(run_length, self.max_score)
                scores.extend([run_score] * run_length)
                
            elif sequence[i] == 'G':
                # Count consecutive Gs (negative for i-motif)
                run_length = 0
                while i < len(sequence) and sequence[i] == 'G':
                    run_length += 1
                    i += 1
                run_score = -min(run_length, self.max_score)
                scores.extend([run_score] * run_length)
                
            else:
                scores.append(0.0)
                i += 1
        
        return sum(scores) / len(scores) if scores else 0.0
    
    def get_score_range(self) -> Tuple[float, float]:
        return (-self.max_score, self.max_score)

class ZDNAScorer(BaseScorer):
    """
    Z-DNA scoring using Kadane's maximum subarray algorithm.
    
    Scientific Basis: Z-DNA forms from alternating purine-pyrimidine
    sequences with specific dinucleotide preferences. Uses weighted
    dinucleotide scoring based on experimental data.
    
    Reference: Ho et al. (1986) Nucleic Acids Res
    """
    
    def __init__(self):
        # Dinucleotide weights based on Z-DNA formation propensity
        self.weights = {
            'CG': 3.0, 'GC': 3.0,  # Strong Z-DNA formers
            'CA': 2.0, 'TG': 2.0,  # Moderate Z-DNA formers
            'CT': 1.5, 'AG': 1.5,  # Weak Z-DNA formers
            'AT': -1.0, 'TA': -1.0,  # Z-DNA inhibitors
            'AA': -2.0, 'TT': -2.0, 'GG': -2.0, 'CC': -2.0  # Strong inhibitors
        }
    
    def score(self, sequence: str, **kwargs) -> float:
        """Calculate Z-DNA score using Kadane's algorithm."""
        if len(sequence) < 4:
            return 0.0
        
        sequence = sequence.upper()
        
        # Convert to dinucleotide scores
        dinuc_scores = []
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            score = self.weights.get(dinuc, 0.0)
            dinuc_scores.append(score)
        
        if not dinuc_scores:
            return 0.0
        
        # Apply Kadane's maximum subarray algorithm
        max_score = current_score = dinuc_scores[0]
        
        for score in dinuc_scores[1:]:
            current_score = max(score, current_score + score)
            max_score = max(max_score, current_score)
        
        return max_score
    
    def get_score_range(self) -> Tuple[float, float]:
        """Z-DNA scores can range widely based on sequence length."""
        return (-100.0, 1000.0)

class ConservationScorer(BaseScorer):
    """
    Evolutionary conservation scoring using composition-preserving shuffles.
    
    Scientific Basis: Functionally important sequences show higher
    conservation than expected by chance. Uses k-mer enrichment
    analysis with null model from shuffled sequences.
    """
    
    def __init__(self, num_shuffles: int = 100):
        self.num_shuffles = num_shuffles
    
    def score(self, sequence: str, motif_type: str = "general", **kwargs) -> float:
        """Calculate conservation score using k-mer enrichment."""
        if not sequence or len(sequence) < 4:
            return 0.0
        
        sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        # Choose k-mer length based on sequence and motif type
        k, core_kmers = self._get_motif_kmers(sequence, motif_type)
        
        if not core_kmers or k > len(sequence):
            return 0.0
        
        # Count observed k-mer occurrences
        observed_count = self._count_kmers(sequence, core_kmers, k)
        
        # Generate null distribution from shuffled sequences
        shuffle_counts = []
        seq_list = list(sequence)
        
        for _ in range(self.num_shuffles):
            random.shuffle(seq_list)
            shuffled_seq = ''.join(seq_list)
            shuffle_count = self._count_kmers(shuffled_seq, core_kmers, k)
            shuffle_counts.append(shuffle_count)
        
        # Calculate enrichment score and p-value
        if not shuffle_counts:
            return 0.0
        
        mean_shuffle = sum(shuffle_counts) / len(shuffle_counts)
        
        if mean_shuffle == 0:
            enrichment = 0.0
        else:
            enrichment = math.log2(max(1, observed_count) / max(1, mean_shuffle))
        
        # Calculate empirical p-value
        extreme_count = sum(1 for c in shuffle_counts if c >= observed_count)
        p_value = extreme_count / len(shuffle_counts)
        
        # Convert to conservation score (higher = more conserved)
        if p_value == 0:
            conservation_score = 10.0  # Highly significant
        else:
            conservation_score = -math.log10(p_value)
        
        return enrichment * conservation_score
    
    def _get_motif_kmers(self, sequence: str, motif_type: str) -> Tuple[int, List[str]]:
        """Get motif-specific k-mers for conservation analysis."""
        seq_len = len(sequence)
        
        if motif_type == "g4":
            k = 3
            core_kmers = ["GGG", "CCC"]
        elif motif_type == "zdna":
            k = 2
            core_kmers = ["CG", "GC"]
        elif motif_type == "triplex":
            k = 3
            core_kmers = ["GAA", "TTC", "AGG", "CCT"]
        elif motif_type == "curved":
            k = 3
            core_kmers = ["AAA", "TTT"]
        else:  # general
            k = min(4, seq_len // 4)
            # Use most frequent k-mers as core elements
            kmer_counts = {}
            for i in range(len(sequence) - k + 1):
                kmer = sequence[i:i+k]
                kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            
            # Select top k-mers as core elements
            sorted_kmers = sorted(kmer_counts.items(), key=lambda x: x[1], reverse=True)
            core_kmers = [kmer for kmer, count in sorted_kmers[:5]]
        
        return k, core_kmers
    
    def _count_kmers(self, sequence: str, core_kmers: List[str], k: int) -> int:
        """Count occurrences of core k-mers in sequence."""
        count = 0
        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer in core_kmers:
                count += 1
        return count
    
    def get_score_range(self) -> Tuple[float, float]:
        """Conservation scores can range widely."""
        return (-10.0, 50.0)

class ThermodynamicScorer(BaseScorer):
    """
    Thermodynamic stability scoring for DNA structures.
    
    Scientific Basis: Structural stability correlates with biological
    function. Uses nearest-neighbor parameters and loop penalty models.
    """
    
    def __init__(self, temperature: float = 310.15):  # 37°C
        self.temperature = temperature
        # Simplified nearest-neighbor parameters (kcal/mol)
        self.stacking_energies = {
            'AA': -1.0, 'AT': -0.88, 'AG': -1.3, 'AC': -2.3,
            'TA': -0.58, 'TT': -1.0, 'TG': -1.45, 'TC': -1.3,
            'GA': -1.29, 'GT': -1.45, 'GG': -1.84, 'GC': -2.17,
            'CA': -1.45, 'CT': -1.29, 'CG': -2.36, 'CC': -1.84
        }
    
    def score(self, sequence: str, structure_type: str = "hairpin", **kwargs) -> float:
        """Calculate thermodynamic stability score."""
        if len(sequence) < 4:
            return 0.0
        
        sequence = sequence.upper()
        
        # Calculate stacking energy
        stacking_energy = 0.0
        for i in range(len(sequence) - 1):
            dinuc = sequence[i:i+2]
            stacking_energy += self.stacking_energies.get(dinuc, 0.0)
        
        # Add structure-specific penalties
        loop_penalty = self._calculate_loop_penalty(sequence, structure_type)
        
        total_energy = stacking_energy + loop_penalty
        
        # Convert to stability score (more negative = more stable)
        return -total_energy
    
    def _calculate_loop_penalty(self, sequence: str, structure_type: str) -> float:
        """Calculate loop penalty based on structure type."""
        seq_len = len(sequence)
        
        if structure_type == "hairpin":
            # Hairpin loop penalty
            return 4.0 + 0.5 * seq_len
        elif structure_type == "g4":
            # G4 loop penalty (simplified)
            return 2.0 + 0.3 * seq_len
        elif structure_type == "triplex":
            # Triplex stability depends on length
            return -0.1 * seq_len  # Longer = more stable
        else:
            return 3.0  # Default penalty
    
    def get_score_range(self) -> Tuple[float, float]:
        """Thermodynamic scores in kcal/mol units."""
        return (-50.0, 20.0)

class CombinatorialScorer:
    """
    Combines multiple scoring functions with flexible weighting.
    
    Enables creation of custom scoring schemes by combining
    pattern, structural, conservation, and thermodynamic scores.
    """
    
    def __init__(self):
        self.scorers = {}
        self.weights = {}
    
    def add_scorer(self, name: str, scorer: BaseScorer, weight: float = 1.0):
        """Add a scoring function with specified weight."""
        self.scorers[name] = scorer
        self.weights[name] = weight
    
    def score(self, sequence: str, **kwargs) -> ScoreComponents:
        """Calculate combined score using all registered scorers."""
        components = ScoreComponents()
        
        for name, scorer in self.scorers.items():
            try:
                score_value = scorer.score(sequence, **kwargs)
                
                # Map to appropriate component
                if isinstance(scorer, (G4HunterScorer, iMotifScorer, ZDNAScorer)):
                    components.pattern_score += score_value * self.weights[name]
                elif isinstance(scorer, ConservationScorer):
                    components.conservation_score += score_value * self.weights[name]
                elif isinstance(scorer, ThermodynamicScorer):
                    components.thermodynamic_score += score_value * self.weights[name]
                else:
                    components.structural_score += score_value * self.weights[name]
                    
            except Exception as e:
                # Log error but continue with other scorers
                print(f"Warning: Scorer {name} failed: {e}")
                continue
        
        return components
    
    def get_total_score(self, sequence: str, score_weights: Optional[Dict[str, float]] = None, **kwargs) -> float:
        """Get weighted total score."""
        components = self.score(sequence, **kwargs)
        return components.total(score_weights)

# Pre-configured scoring schemes for different motif types
def create_g4_scorer() -> CombinatorialScorer:
    """Create optimized G4 scoring scheme."""
    scorer = CombinatorialScorer()
    scorer.add_scorer('g4hunter', G4HunterScorer(include_structural=True), weight=1.0)
    scorer.add_scorer('conservation', ConservationScorer(num_shuffles=50), weight=0.3)
    scorer.add_scorer('thermodynamic', ThermodynamicScorer(), weight=0.2)
    return scorer

def create_zdna_scorer() -> CombinatorialScorer:
    """Create optimized Z-DNA scoring scheme."""
    scorer = CombinatorialScorer()
    scorer.add_scorer('zdna', ZDNAScorer(), weight=1.0)
    scorer.add_scorer('conservation', ConservationScorer(num_shuffles=50), weight=0.4)
    return scorer

def create_imotif_scorer() -> CombinatorialScorer:
    """Create optimized i-motif scoring scheme."""
    scorer = CombinatorialScorer()
    scorer.add_scorer('imotif', iMotifScorer(), weight=1.0)
    scorer.add_scorer('conservation', ConservationScorer(num_shuffles=50), weight=0.3)
    scorer.add_scorer('thermodynamic', ThermodynamicScorer(), weight=0.2)
    return scorer