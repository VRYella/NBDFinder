"""
Configurable Regex Engine for Non-B DNA Motif Detection
=====================================================

This module provides a centralized, configurable regex pattern system for
detecting non-B DNA motifs with adjustable parameters.

Features:
- Pre-compiled regex patterns for performance
- Configurable motif parameters (length, spacing, loops)
- Pattern caching and optimization
- Batch processing support
- Scientific parameter validation

Scientific Basis:
Pattern parameters are based on experimental structure determination and
molecular dynamics studies of non-B DNA forms.

References:
- Wells et al. (1988) J Biol Chem (curved DNA spacing)
- Frank-Kamenetskii & Mirkin (1995) Annu Rev Biochem (DNA structures)
- Champ et al. (2004) Nucleic Acids Res (G4 parameters)

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

import re
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass, field
from functools import lru_cache
import threading

@dataclass
class MotifConfig:
    """Configuration for a specific motif type."""
    name: str
    description: str
    min_length: int = 10
    max_length: int = 1000
    min_spacing: int = 1
    max_spacing: int = 50
    min_loop_size: int = 1
    max_loop_size: int = 30
    min_repeats: int = 2
    max_repeats: int = 10
    gc_content_range: Tuple[float, float] = (0.0, 100.0)
    context_requirements: Dict[str, Any] = field(default_factory=dict)
    score_threshold: float = 0.0
    overlap_tolerance: float = 0.3
    
    def validate(self) -> bool:
        """Validate configuration parameters."""
        if self.min_length > self.max_length:
            raise ValueError(f"min_length ({self.min_length}) > max_length ({self.max_length})")
        if self.min_spacing > self.max_spacing:
            raise ValueError(f"min_spacing ({self.min_spacing}) > max_spacing ({self.max_spacing})")
        if not (0 <= self.gc_content_range[0] <= self.gc_content_range[1] <= 100):
            raise ValueError(f"Invalid GC content range: {self.gc_content_range}")
        return True

class RegexPatternEngine:
    """
    High-performance regex engine for non-B DNA motif detection.
    
    This engine provides:
    1. Pre-compiled pattern caching
    2. Configurable motif parameters
    3. Batch processing capabilities
    4. Thread-safe operations
    """
    
    def __init__(self):
        self._pattern_cache: Dict[str, re.Pattern] = {}
        self._config_cache: Dict[str, MotifConfig] = {}
        self._lock = threading.RLock()
        self._initialize_default_configs()
    
    def _initialize_default_configs(self):
        """Initialize scientifically validated default configurations."""
        
        # G-Quadruplex configurations based on experimental data
        self._config_cache['g4_canonical'] = MotifConfig(
            name="G4 Canonical",
            description="Canonical G-quadruplex with standard loop constraints",
            min_length=15,
            max_length=200,
            min_spacing=1,
            max_spacing=7,
            min_loop_size=1,
            max_loop_size=7,
            min_repeats=4,
            score_threshold=1.0
        )
        
        self._config_cache['g4_relaxed'] = MotifConfig(
            name="G4 Relaxed",
            description="Relaxed G-quadruplex with extended loop lengths",
            min_length=15,
            max_length=300,
            min_spacing=1,
            max_spacing=12,
            min_loop_size=1,
            max_loop_size=12,
            min_repeats=4,
            score_threshold=0.8
        )
        
        # Curved DNA configurations
        self._config_cache['curved_periodic'] = MotifConfig(
            name="Curved DNA Periodic",
            description="Periodically curved DNA with A/T tracts",
            min_length=21,
            max_length=500,
            min_spacing=8,
            max_spacing=12,
            min_repeats=3,
            context_requirements={'helical_period': 10.5}
        )
        
        # Triplex DNA configurations
        self._config_cache['triplex_purine'] = MotifConfig(
            name="Triplex Purine",
            description="Purine-rich triplex-forming sequences",
            min_length=15,
            max_length=1000,
            min_repeats=15,
            score_threshold=25.0,
            context_requirements={'purine_fraction': 0.7}
        )
        
        # Z-DNA configurations
        self._config_cache['zdna'] = MotifConfig(
            name="Z-DNA",
            description="Left-handed Z-DNA forming sequences",
            min_length=6,
            max_length=1000,
            context_requirements={'alternating_pyrimidine_purine': True},
            score_threshold=500.0
        )
        
        # R-loop configurations
        self._config_cache['rloop'] = MotifConfig(
            name="R-loop",
            description="R-loop forming sequences with RIZ and REZ",
            min_length=50,
            max_length=5000,
            gc_content_range=(60.0, 100.0),
            score_threshold=50.0
        )
        
        # i-motif configurations
        self._config_cache['imotif'] = MotifConfig(
            name="i-motif",
            description="Cytosine-rich i-motif structures",
            min_length=12,
            max_length=200,
            min_spacing=1,
            max_spacing=15,
            min_repeats=4,
            score_threshold=1.0
        )
    
    @lru_cache(maxsize=1000)
    def _compile_pattern(self, pattern: str, flags: int = re.IGNORECASE) -> re.Pattern:
        """Compile and cache regex patterns for performance."""
        return re.compile(pattern, flags)
    
    def get_config(self, motif_type: str) -> MotifConfig:
        """Get configuration for a motif type."""
        with self._lock:
            if motif_type not in self._config_cache:
                raise ValueError(f"Unknown motif type: {motif_type}")
            return self._config_cache[motif_type]
    
    def set_config(self, motif_type: str, config: MotifConfig):
        """Set configuration for a motif type."""
        config.validate()
        with self._lock:
            self._config_cache[motif_type] = config
            # Clear related pattern cache
            keys_to_remove = [k for k in self._pattern_cache.keys() if k.startswith(motif_type)]
            for key in keys_to_remove:
                del self._pattern_cache[key]
    
    def build_g4_pattern(self, config: Optional[MotifConfig] = None) -> str:
        """
        Build G-quadruplex detection pattern with configurable parameters.
        
        Scientific Basis: G4 motifs require 4 G-runs of ≥3 Gs separated by loops.
        Pattern captures all possible G4 topologies including parallel, antiparallel.
        
        Returns optimized regex pattern for G4 detection.
        """
        if config is None:
            config = self.get_config('g4_canonical')
        
        # Base G-run pattern (minimum 3 Gs, configurable)
        g_run = f"G{{3,}}"
        
        # Loop pattern with configurable size constraints
        loop_pattern = f"[ATGC]{{{config.min_loop_size},{config.max_loop_size}}}"
        
        # Complete G4 pattern: 4 G-runs separated by loops
        pattern = f"(?=({g_run}{loop_pattern}{g_run}{loop_pattern}{g_run}{loop_pattern}{g_run}))"
        
        return pattern
    
    def build_imotif_pattern(self, config: Optional[MotifConfig] = None) -> str:
        """Build i-motif detection pattern (C-rich counterpart to G4)."""
        if config is None:
            config = self.get_config('imotif')
        
        c_run = f"C{{3,}}"
        loop_pattern = f"[ATGC]{{{config.min_loop_size},{config.max_loop_size}}}"
        pattern = f"(?=({c_run}{loop_pattern}{c_run}{loop_pattern}{c_run}{loop_pattern}{c_run}))"
        
        return pattern
    
    def build_curved_dna_pattern(self, config: Optional[MotifConfig] = None) -> str:
        """
        Build curved DNA detection pattern.
        
        Scientific Basis: Curved DNA results from phased A-tracts or T-tracts 
        with ~10.5 bp spacing (helical periodicity). Intrinsic curvature depends
        on sequence context and bendability.
        """
        if config is None:
            config = self.get_config('curved_periodic')
        
        # A-tract pattern
        a_tract = f"A{{{config.min_repeats},}}"
        t_tract = f"T{{{config.min_repeats},}}"
        
        # Spacing based on helical periodicity
        spacing = f"[ATGC]{{{config.min_spacing},{config.max_spacing}}}"
        
        # Periodic pattern for multiple tracts
        a_periodic = f"({a_tract}{spacing}){{3,}}{a_tract}"
        t_periodic = f"({t_tract}{spacing}){{3,}}{t_tract}"
        
        return f"({a_periodic}|{t_periodic})"
    
    def build_triplex_pattern(self, config: Optional[MotifConfig] = None) -> str:
        """
        Build triplex DNA (H-DNA) detection pattern.
        
        Scientific Basis: Triplex DNA forms from purine-rich or pyrimidine-rich
        homopurine/homopyrimidine mirror repeats. Requires minimum length for stability.
        """
        if config is None:
            config = self.get_config('triplex_purine')
        
        # Purine-rich patterns (A/G repeats)
        purine_repeat = f"[AG]{{{config.min_repeats},}}"
        
        # Pyrimidine-rich patterns (C/T repeats) 
        pyrimidine_repeat = f"[CT]{{{config.min_repeats},}}"
        
        return f"({purine_repeat}|{pyrimidine_repeat})"
    
    def build_zdna_pattern(self, config: Optional[MotifConfig] = None) -> str:
        """
        Build Z-DNA detection pattern.
        
        Scientific Basis: Z-DNA forms from alternating purine-pyrimidine sequences,
        particularly CG dinucleotides. Left-handed helix with distinctive structure.
        """
        if config is None:
            config = self.get_config('zdna')
        
        # Alternating CG pattern (classic Z-DNA)
        cg_alternating = f"(?:CG){{3,}}"
        
        # More general alternating purine-pyrimidine
        alt_pattern = f"(?:[CG][AT]|[AT][CG]){{6,}}"
        
        return f"({cg_alternating}|{alt_pattern})"
    
    def find_motifs(self, sequence: str, motif_type: str, 
                   config: Optional[MotifConfig] = None) -> List[Dict[str, Any]]:
        """
        Find motifs in sequence using configurable patterns.
        
        Parameters:
        sequence (str): DNA sequence to search
        motif_type (str): Type of motif to detect
        config (MotifConfig): Optional configuration override
        
        Returns:
        List of motif matches with positions and scores
        """
        sequence = sequence.upper().replace('\n', '').replace(' ', '')
        
        # Get pattern based on motif type
        pattern_builders = {
            'g4': self.build_g4_pattern,
            'imotif': self.build_imotif_pattern,
            'curved': self.build_curved_dna_pattern,
            'triplex': self.build_triplex_pattern,
            'zdna': self.build_zdna_pattern
        }
        
        if motif_type not in pattern_builders:
            raise ValueError(f"Unsupported motif type: {motif_type}")
        
        pattern_str = pattern_builders[motif_type](config)
        pattern = self._compile_pattern(pattern_str)
        
        motifs = []
        for match in pattern.finditer(sequence):
            motif_seq = match.group(1) if match.groups() else match.group(0)
            
            motif_data = {
                'sequence': motif_seq,
                'start': match.start(),
                'end': match.end(),
                'length': len(motif_seq),
                'motif_type': motif_type,
                'pattern': pattern_str
            }
            
            motifs.append(motif_data)
        
        return motifs
    
    def batch_find_motifs(self, sequences: List[str], motif_types: List[str],
                         configs: Optional[Dict[str, MotifConfig]] = None) -> Dict[str, List[Dict[str, Any]]]:
        """
        Batch process multiple sequences and motif types.
        
        Optimized for high-throughput analysis with minimal overhead.
        """
        results = {}
        
        for seq_idx, sequence in enumerate(sequences):
            seq_key = f"sequence_{seq_idx}"
            results[seq_key] = []
            
            for motif_type in motif_types:
                config = configs.get(motif_type) if configs else None
                motifs = self.find_motifs(sequence, motif_type, config)
                results[seq_key].extend(motifs)
        
        return results
    
    def get_pattern_stats(self) -> Dict[str, Any]:
        """Get statistics about cached patterns and configurations."""
        with self._lock:
            return {
                'cached_patterns': len(self._pattern_cache),
                'cached_configs': len(self._config_cache),
                'available_motif_types': list(self._config_cache.keys()),
                'cache_info': self._compile_pattern.cache_info()._asdict()
            }

# Global instance for shared use
regex_engine = RegexPatternEngine()

# Convenience functions for backward compatibility
def get_g4_pattern(min_loop=1, max_loop=7):
    """Get G4 pattern with specified loop constraints."""
    config = MotifConfig(
        name="custom_g4",
        min_loop_size=min_loop,
        max_loop_size=max_loop
    )
    return regex_engine.build_g4_pattern(config)

def get_curved_pattern(min_spacing=8, max_spacing=12, min_repeats=3):
    """Get curved DNA pattern with specified spacing."""
    config = MotifConfig(
        name="custom_curved",
        min_spacing=min_spacing,
        max_spacing=max_spacing,
        min_repeats=min_repeats
    )
    return regex_engine.build_curved_dna_pattern(config)