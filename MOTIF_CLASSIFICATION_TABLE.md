# NBDFinder - Complete Non-B DNA Motif Classification Table

This comprehensive table documents all non-B DNA motif classes and subclasses implemented in NBDFinder, following the classification framework requested in the issue.

## Complete Motif Classification System

| Class | Subclass/Motif | Detection Logic / Pattern Description | Regex/Pattern Example | Scoring System | Score Range |
|-------|----------------|--------------------------------------|----------------------|----------------|-------------|
| **G-Quadruplex (G4)** | Canonical PQS (G3+L1–7) | Four runs of ≥3 Gs, loops 1–7 nt; "gold standard" PQS | `([G]{3,}N{1,7}){3,}[G]{3,}` | G4Hunter: avg +1 per G, -1 per C, capped ±4/run, mean per base | <1.0 (Low), 1.0–1.5 (Mod), ≥1.5 (High) |
| | Extended Canonical (G3+L1–12) | Four runs of ≥3 Gs, loops 1–12 nt | `([G]{3,}N{1,12}){3,}[G]{3,}` | G4Hunter (as above), lower threshold | ≥0.8 (Mod–High) |
| | Long-Loop (G3+L8–12) | At least one loop 8–12 nt; subset of above | Subset: ≥1 loop 8–12 nt | G4Hunter | Moderate–Low |
| | Two-Tetrad PQS (G2L1–12) | Four G-runs, each 2 Gs, loops 1–12 nt | `([G]{2}N{1,12}){3,}[G]{2}` | G4Hunter | Moderate–Low |
| | Generalized PQS (G2+L1–12) | Four G-runs, ≥2 Gs per run, loops 1–12 nt | `([G]{2,}N{1,12}){3,}[G]{2,}` | G4Hunter | Mixed |
| | Bulged G4s | G-runs interrupted by 1–3 nt bulges | `([G]{2,}N{0,3}[G]{1,}N{1,12}){3,}[G]{2,}` | G4Hunter, bulge penalty | Variable |
| | Snapback G4s | One G-run provided by strand foldback | Not explicit (secondary structure) | G4Hunter logic, but topological annotation | Moderate |
| | Multimeric/Tandem G4s | Two+ canonical G4s in proximity | `(([G]{3,}N{1,7}){3,}[G]{3,}){2,}` | Stacking: G4Hunter + stacking bonus | Very High |
| | Split/Bipartite G4s | Distant G-runs with long loops (20–50+ nt) | `[G]{3,}N{20,50}[G]{3,}` | G4Hunter, penalty for long loop | Moderate–Low |
| | G-Triplex | Three G-runs, loops ≤15 nt | `G{3,}N{1,15}G{3,}N{1,15}G{3,}` | Length × G-content × #runs | Moderate–Low |
| **i-Motif** | Canonical i-motif | Four runs of ≥3 C, loops 1–7 nt | `C{3,}N{1,7}C{3,}N{1,7}C{3,}N{1,7}C{3,}` | Hunter-style: +1/C, -1/G, capped, mean per base | High (acidic, crowded) |
| | Relaxed i-motif | Four runs of ≥3 C, loops 1–12 nt | `C{3,}N{1,12}C{3,}N{1,12}C{3,}N{1,12}C{3,}` | Hunter-style (lower threshold) | Moderate |
| | AC-motif | Four C-runs, A/C-rich spacers | See Hur et al. 2021 (use N for loop, C/A enrichment) | C/A content, loop filter | Moderate |
| **R-loop** | RLFS (m1) | G{3,}, 1–10 nt, G{3,}, at least two runs | `G{3,}N{1,10}?G{3,}(?:N{1,10}?G{3,}){1,}` | (GC%×50+G-runs×10)×len^0.25 | >10 (high) |
| | RLFS (m2) | ≥4 Gs/run, at least two runs | `G{4,}(?:N{1,10}?G{4,}){1,}` | As above | High–Moderate |
| | RLFS+REZ | Downstream GC-rich region after RLFS (window, min. GC%) | Sliding window, GC%>40 | Combined RLFS+REZ score | High |
| **Z-DNA** | Z-DNA (Kadane) | Alternating purine-pyrimidine, dinucleotide scoring: GC/CG(+7), GT/TG/AC/CA(+1.25), AT/TA(0), others(–3) | Not regex; dinucleotide scan | Max subarray sum (Kadane); region score | ≥50 (threshold) |
| | eGZ (Extruded-G) | (CGG) repeat expansions, ≥3 repeats | `(CGG){3,}` | n_repeats×3×(1+2G_frac) | >10 (pathological) |
| **Triplex (H-DNA)** | Mirror repeat triplex | Mirror A/G or T/C runs, min 10–15 nt; polypurine/polypyrimidine tracts | e.g., `(A{n,}N{0,5}A{n,})` | Repeat count, symmetry | Variable |
| | Sticky DNA | (GAA)n or similar long direct repeats | `(GAA){n,}` | Repeat count | High (expanded) |
| **Cruciform** | Palindromic cruciform | Inverted repeats, stem-loop, stem/loop min lengths | `([ACTG]{n,})N{0,10}\1` | Stem/loop threshold | Variable |
| **Slipped DNA (STRs)** | Mono/di/tri/tetra STRs | Tandem repeats; motif-specific min repeat | `(AT){n,}`, `(GAA){n,}`, `(CGG){n,}` | Repeat count | Variable |
| **Curved DNA** | PolyA/PolyT tracts | AT-rich tracts, polyA/polyT ≥5 nt | `A{5,}`, `T{5,}` | Length, AT content | Low |
| **Hybrids** | G4/R-loop, G4/Z-DNA | Overlap of G4 and R-loop/Z-DNA motifs | Both motifs overlap | Combined, presence of both motifs | High |
| **Cluster** | Non-B DNA cluster | ≥2 motif types in region/window (aggregation) | Sliding window, motif types≥2 | Motif density | High |

## Implementation Notes

### Performance Optimizations Applied:
1. **Consolidated G4 Detection**: Merged duplicate G4 modules for better performance
2. **Optimized Regex Patterns**: Pre-compiled patterns with priority ordering
3. **Caching**: Sequence preprocessing and score caching
4. **Memory Efficiency**: Reduced redundant data structures

### Scientific Accuracy:
- All patterns follow literature-standard definitions
- G4Hunter scoring uses exact algorithm from Bedrat et al. NAR 2016
- Priority-based selection prevents overlapping predictions
- Maintains biological relevance while optimizing computational efficiency

### Quality Assurance:
- Comprehensive test suite for all motif types
- Validation against known non-B DNA databases
- Performance benchmarking for large sequences
- Error handling and edge case management

## Usage Examples

```python
from motifs.g4_related import find_all_g4_motifs
from motifs.r_loop import find_rlfs
from motifs.imotif_ac import find_imotifs

# Find all G4 motifs with non-overlapping selection
g4_results = find_all_g4_motifs(sequence, use_non_overlapping=True)

# Find R-loop forming sequences
rloop_results = find_rlfs(sequence, models=["m1", "m2"])

# Find i-motifs
imotif_results = find_imotifs(sequence)
```

## References

1. Bedrat et al. (2016) "Re-evaluation of G-quadruplex propensity with G4Hunter" Nucleic Acids Research
2. Hänsel-Hertsch et al. (2017) "G-quadruplex structures mark human regulatory chromatin" Nature Reviews Molecular Cell Biology  
3. Chambers et al. (2015) "High-throughput sequencing of DNA G-quadruplex structures in the human genome" Nature Biotechnology
4. QuadBase2 database for comprehensive G4 classification
5. Hur et al. (2021) "AC-motifs in mammalian genomes" Various references for AC-motif detection

---
*NBDFinder v2.0 - Advanced Non-B DNA Analysis Platform*  
*Author: Dr. Venkata Rajesh Yella*