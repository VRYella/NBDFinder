# NBDFinder: Non-B DNA Analysis Platform

[![Python 3.8+](https://img.shields.io/badge/python-3.8+-blue.svg)](https://www.python.org/downloads/)
[![Streamlit](https://img.shields.io/badge/Built%20with-Streamlit-red.svg)](https://streamlit.io/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-10.1038%2Fs41467--023--42156--1-blue.svg)](https://doi.org/10.1038/s41467-023-42156-1)

## Non-B DNA Structure Detection Platform

**NBDFinder** is a computational framework for genome-wide detection and analysis of Non-B DNA structural motifs. The platform provides various algorithms for identifying different types of non-canonical DNA structures.

### Key Features

- **Multiple Detection Algorithms**: Implements various algorithms for different motif types
- **Sequence Processing**: Handles FASTA format sequences and NCBI database queries
- **Interactive Visualizations**: Generates plots and statistical analysis of detected motifs
- **Export Options**: Results can be exported in CSV and Excel formats
- **Web Interface**: User-friendly Streamlit-based interface for analysis

## Comprehensive Motif Detection Suite

NBDFinder detects and analyzes **19 distinct Non-B DNA motif classes** across all major structural categories:

### G-Quadruplex Family
- **Canonical G4**: Standard four-stranded structures with Hoogsteen base pairing
- **Relaxed G4**: Extended loop regions maintaining quadruplex stability
- **Bulged G4**: Nucleotide bulges within G-tracts preserving topology
- **Bipartite G4**: Two G4-forming sequences enabling long-range DNA looping
- **Multimeric G4**: Tandem arrays creating higher-order chromatin structures
- **Imperfect G4**: Alternative structures with imperfect G-tracts

### Triplex & i-Motif Structures
- **G-Triplex**: Three-stranded parallel structures with G·G·G base triads
- **Triplex DNA**: Watson-Crick and Hoogsteen base pairing combinations
- **i-Motif**: pH-dependent cytosine tetraplex structures

### Helix Deviations
- **Z-DNA**: Left-handed double helix with alternating purine-pyrimidine sequences
- **eGZ (Extruded-G)**: CGG repeat expansions causing chromatin silencing
- **Curved DNA**: Intrinsic DNA curvature from phased A/T tracts
- **AC-Motif**: Alternating purine-pyrimidine patterns with enhanced bendability

### Repeat & Junction Structures
- **Slipped DNA**: Mispairing during replication creating hairpin structures
- **Cruciform**: Four-way Holliday junctions from palindromic sequences
- **Sticky DNA**: GAA/TTC expansions associated with neurodegeneration
- **R-Loop**: RNA-DNA hybrids in transcriptionally active regions

### Advanced Analysis
- **Hybrid Motifs**: Superposition of multiple non-B structures
- **Non-B DNA Clusters**: Genomic hotspots with multiple motif types

## Quick Start

### Installation

```bash
git clone https://github.com/VRYella/NBDFinder.git
cd NBDFinder
pip install -r requirements.txt
```

### Launch the Application

```bash
streamlit run app.py
```

The application will open in your browser at `http://localhost:8501`

## Key Features

### Detection Algorithms
- **G4Hunter**: G-quadruplex detection with structural factors
- **Kadane's Algorithm**: Z-DNA detection with dinucleotide weighting
- **RLFS+REZ Method**: R-loop prediction based on thermodynamic principles
- **Various Methods**: Additional algorithms for different motif types

### Analysis Features
- **Sequence Input**: Supports FASTA format and NCBI database queries
- **Interactive Visualizations**: Plotly-based charts and graphs
- **Statistical Analysis**: Basic motif distribution analysis
- **Export Options**: Results available in CSV and Excel formats
- **Batch Processing**: Multi-FASTA file support

### Data Integration
- **NCBI Access**: Direct GenBank sequence retrieval
- **Example Datasets**: Curated sequences for testing
- **Format Support**: FASTA, Multi-FASTA, plain text
- **Progress Tracking**: Basic progress monitoring

## Scientific Documentation

### Core Algorithm References

1. **Bedrat, A., Lacroix, L., & Mergny, J.L. (2016)** "Re-evaluation of G-quadruplex propensity with G4Hunter." *Nucleic Acids Research* 44(4): 1746-1759. [DOI: 10.1093/nar/gkw006](https://doi.org/10.1093/nar/gkw006)

2. **Hänsel-Hertsch, R. et al. (2017)** "G-quadruplex structures mark human regulatory chromatin." *Nature Genetics* 49: 1212-1221. [DOI: 10.1038/ng.3917](https://doi.org/10.1038/ng.3917)

3. **Santos-Pereira, J.M. & Aguilera, A. (2015)** "R loops: new modulators of genome dynamics and function." *Nature Reviews Genetics* 16(10): 583-597. [DOI: 10.1038/nrg3961](https://doi.org/10.1038/nrg3961)

### Implementation Notes

The algorithms implemented in NBDFinder are based on published methods with various optimizations for computational efficiency. Performance may vary depending on sequence length and complexity.

## Clinical Applications

### Disease Associations
- **Neurological Disorders**: Huntington's, Friedreich's ataxia, ALS
- **Cancer Biology**: Oncogene regulation, telomere maintenance
- **Genetic Instability**: Microsatellite disorders, repeat expansions
- **Therapeutic Targeting**: Drug design, antigene therapy

### Research Applications
- **Genomics**: Whole genome structure analysis
- **Epigenetics**: Chromatin organization studies
- **Evolution**: Phylogenetic conservation analysis
- **Biotechnology**: Synthetic biology applications

## 🤝 Contributing

We welcome contributions from the scientific community:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 Citation

If you use NBDFinder in your research, please cite:

```bibtex
@article{yella2024nbdfinder,
  title={NBDFinder: A Comprehensive Computational Framework for Non-B DNA Structure Detection},
  author={Yella, Venkata Rajesh},
  journal={Nature Communications},
  volume={15},
  pages={8234},
  year={2024},
  publisher={Nature Publishing Group},
  doi={10.1038/s41467-023-42156-1}
}
```

## 📞 Contact & Support

**Developer**: Dr. Venkata Rajesh Yella  
**Email**: [yvrajesh_bt@kluniversity.in](mailto:yvrajesh_bt@kluniversity.in)  
**GitHub**: [@VRYella](https://github.com/VRYella)  
**Institution**: KL University, India

### 🆘 Support Resources
- **Documentation**: Comprehensive guides in the application
- **Issues**: Report bugs via GitHub Issues
- **Discussions**: Scientific discussions via GitHub Discussions
- **Updates**: Follow releases for latest features

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Acknowledgments

- **Scientific Community**: For published algorithms and methods
- **Streamlit Team**: For the web framework
- **Plotly**: For visualization capabilities
- **Bioinformatics Community**: For algorithm development

---

**NBDFinder** - A computational tool for Non-B DNA structure detection