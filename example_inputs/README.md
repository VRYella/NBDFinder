# NBDFinder Example Input Files

This directory contains example FASTA files for testing and demonstrating the capabilities of NBDFinder.

## File Descriptions

### 1. g4_rich_sequence.fasta
Contains sequences rich in G-quadruplex forming motifs:
- **Human telomere sequence**: Classic TTAGGG repeats
- **MYC oncogene promoter**: Cancer-associated G4 sequences
- **Immunoglobulin switch region**: Antibody diversity sequences

**Expected motifs**: G4, Relaxed G4, G-Triplex, potentially i-Motif

### 2. disease_repeats.fasta  
Contains pathogenic repeat sequences associated with genetic diseases:
- **Friedreich's ataxia**: GAA repeats (>59 repeats = pathogenic)
- **Fragile X syndrome**: CGG repeats (>200 repeats = full mutation)
- **Huntington's disease**: CAG repeats (>36 repeats = pathogenic)

**Expected motifs**: Sticky DNA, eGZ (Extruded-G), Slipped DNA, Triplex DNA

### 3. structural_motifs.fasta
Contains sequences designed to form specific structural motifs:
- **Z-DNA sequences**: Alternating purine-pyrimidine patterns
- **Curved DNA**: Phased A-tracts causing DNA bending
- **Cruciform structures**: Palindromic sequences forming four-way junctions

**Expected motifs**: Z-DNA, Curved DNA, Cruciform

### 4. comprehensive_example.fasta
Contains mixed sequences with multiple motif types for comprehensive analysis:
- **Comprehensive mixed motifs**: Multiple structural elements in one sequence
- **R-loop prone sequences**: G-rich regions for RNA-DNA hybrid formation
- **Cancer oncogene region**: Real oncogene sequences with regulatory motifs

**Expected motifs**: Multiple classes including G4, Z-DNA, R-Loop, Hybrid, Non-B DNA Clusters

## Usage Instructions

1. **Upload via web interface**: Use the "Upload & Analyze" tab in NBDFinder
2. **Select appropriate motif classes**: Choose specific motifs or analyze all
3. **Review results**: Check the Results tab for detailed analysis
4. **Download data**: Export findings via the Download tab

## Scientific Validation

These sequences are based on:
- Experimentally validated non-B DNA forming sequences
- Clinically relevant pathogenic repeats  
- Published genomic regions with known structural properties
- Synthetic sequences designed for educational purposes

## Expected Analysis Results

| File | Primary Motifs | Secondary Motifs | Biological Relevance |
|------|---------------|------------------|---------------------|
| g4_rich_sequence.fasta | G4, G-Triplex | i-Motif, Relaxed G4 | Telomere biology, cancer |
| disease_repeats.fasta | Sticky DNA, Slipped DNA | eGZ, Triplex DNA | Genetic diseases |
| structural_motifs.fasta | Z-DNA, Curved DNA, Cruciform | AC-Motif | Structural biology |
| comprehensive_example.fasta | Multiple classes | Hybrid, Clusters | Genomic complexity |

For questions about these examples, please refer to the NBDFinder documentation or contact the development team.