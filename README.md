# NBDFinder - Advanced Non-B DNA Analysis Platform

NBDFinder is a comprehensive computational platform for detecting and analyzing non-canonical DNA structures (non-B DNA motifs) in genomic sequences. The platform has been completely reorganized and optimized for better performance and maintainability.

## 🚀 Recent Improvements (2024)

### Code Organization & Performance
- **Consolidated duplicate modules**: Merged redundant G4 detection modules for better performance
- **Streamlined architecture**: Replaced large monolithic files with modular design
- **Optimized algorithms**: Improved speed and memory efficiency of core detection functions
- **Archive management**: Moved obsolete files to organized archive structure

### Comprehensive Documentation
- **Complete motif classification table**: Detailed documentation of all 22+ subclasses across 10 major motif families
- **Scientific accuracy**: All patterns follow literature-standard definitions with proper citations
- **Performance benchmarks**: Documented optimization improvements and computational complexity

### Advanced Features
- **Enhanced G4 detection**: Exact G4Hunter algorithm with 9 distinct subclasses
- **Optimized R-loop detection**: Improved REZ detection with reduced computational complexity  
- **Machine learning integration**: Optional ML-enhanced motif predictions
- **Hybrid structure detection**: Dynamic detection of overlapping motif combinations

## 📊 Supported Non-B DNA Motifs

| Motif Class | Subclasses | Key Features |
|-------------|------------|--------------|
| **G-Quadruplex** | 10 subclasses | Canonical, Extended, Bulged, Multimeric, Split, G-Triplex |
| **i-Motif** | 3 subclasses | Canonical, Relaxed, AC-motif variants |
| **R-loop** | 3 subclasses | RLFS (m1, m2), RLFS+REZ combined |
| **Z-DNA** | 2 subclasses | Kadane algorithm, eGZ expansions |
| **Triplex/H-DNA** | 2 subclasses | Mirror repeats, Sticky DNA |
| **Cruciform** | 1 class | Palindromic inverted repeats |
| **Slipped DNA** | 4 subclasses | Mono/di/tri/tetra STRs |
| **Curved DNA** | 1 class | PolyA/PolyT tracts |
| **Hybrids** | Dynamic | G4/R-loop, G4/Z-DNA overlaps |
| **Clusters** | Dynamic | Multi-motif dense regions |

See [MOTIF_CLASSIFICATION_TABLE.md](MOTIF_CLASSIFICATION_TABLE.md) for complete technical specifications.

## 🔬 Scientific Accuracy

All detection algorithms are based on peer-reviewed literature and validated methodologies:
- **G4Hunter algorithm** (Bedrat et al. NAR 2016) for G-quadruplex scoring
- **Kadane's maximum subarray** for Z-DNA region detection
- **RLFS+REZ methodology** for R-loop formation prediction
- **Priority-based selection** to prevent overlapping predictions

## 🛠 Installation & Usage

```bash
# Install dependencies
pip install -r requirements.txt

# Run the web application
streamlit run app.py

# Or use the streamlined motifs module directly
python -c "
import motifs_streamlined as motifs
results = motifs.all_motifs('GGGATTGGGATCGATGGGATTGGG')
print(f'Found {len(results)} motifs')
"
```

## 📈 Performance Improvements

- **Consolidated Detection**: Merged duplicate G4 modules, eliminating redundancy
- **Optimized Algorithms**: Reduced computational complexity in R-loop REZ detection
- **Memory Efficiency**: Streamlined data structures and reduced memory footprint
- **Modular Loading**: Lazy imports for optional features reduce startup time

## 📁 Project Structure

```
NBDFinder/
├── motifs/                    # Individual motif detection modules
├── motifs_streamlined.py      # New unified interface
├── archive/                   # Archived/backup files
├── MOTIF_CLASSIFICATION_TABLE.md  # Complete motif documentation
├── app.py                     # Main Streamlit application
└── requirements.txt           # Dependencies
```

## 🔗 Citation

If you use NBDFinder in your research, please cite:
```
NBDFinder: Advanced Non-B DNA Analysis Platform
Dr. Venkata Rajesh Yella, 2024
```

## 📄 License

Academic Use License - See LICENSE file for details.