import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re, io
from datetime import datetime
from motifs import all_motifs, find_hotspots
from utils import parse_fasta, wrap

# Example sequence
EXAMPLE_FASTA = """>Example
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG"""

# Configure page
st.set_page_config(
    page_title="Non-B DNA Motif Finder", 
    layout="wide",
    page_icon="🧬"
)

# Initialize session state
if 'seq' not in st.session_state:
    st.session_state['seq'] = ""
if 'df' not in st.session_state:
    st.session_state['df'] = pd.DataFrame()
if 'motif_results' not in st.session_state:
    st.session_state['motif_results'] = []
if 'analysis_status' not in st.session_state:
    st.session_state['analysis_status'] = ""
if 'stop_analysis' not in st.session_state:
    st.session_state['stop_analysis'] = False

# App pages
PAGES = {
    "Home": "Welcome and introduction",
    "Upload & Analyze": "Submit your DNA sequence for analysis",
    "Results": "View detected motifs and statistics",
    "Visualization": "Graphical representation of motifs",
    "Download": "Export your results",
    "Additional Information": "About the tool and methodology",
    "Motif Definitions Glossary": "Detailed motif descriptions"
}

# Sidebar navigation
st.sidebar.title("🧬 Navigation")
page = st.sidebar.radio("Go to", list(PAGES.keys()), 
                        help=PAGES[list(PAGES.keys())[0]])

st.sidebar.markdown("---")
st.sidebar.info(
    "Developed by Dr. Venkata Rajesh Yella & Chandrika Gummadi\n\n"
    "[GitHub Repository](https://github.com/VRYella/Non-B-DNA-Finder)"
)

# Main title
st.title("Non-B DNA Motif Finder")

def collect_all_motifs(seq, status_callback=None, stop_flag=None):
    """Wrapper function to collect all motifs with progress tracking"""
    if status_callback:
        status_callback("Scanning for non-B DNA motifs...")
    if stop_flag and stop_flag():
        return []
    return all_motifs(seq)

# Page: Home
if page == "Home":
    st.markdown("""
    <style>
        .home-header { 
            font-size: 2.5em; 
            font-weight: bold; 
            color: #1A5276; 
            text-align: center; 
            margin-bottom: 20px; 
        }
        .feature-card {
            background-color: #f8f9fa;
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 15px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
    </style>
    <div class='home-header'>Welcome to Non-B DNA Motif Finder</div>
    """, unsafe_allow_html=True)
    
    try:
        st.image("nbd.png", use_container_width=True)
    except FileNotFoundError:
        st.warning("Logo image not found. Using placeholder text.")
        st.markdown("### DNA Structure Analysis Tool")
    
    st.markdown("""
    <div class="feature-card">
    <h3>🔍 About This Tool</h3>
    <p>This application identifies and visualizes diverse non-B DNA structures in nucleotide sequences.</p>
    </div>
    """, unsafe_allow_html=True)
    
    cols = st.columns(2)
    with cols[0]:
        st.markdown("""
        <div class="feature-card">
        <h3>🧬 Detectable Structures</h3>
        <ul>
            <li>Curved DNA & A-Phased Repeats</li>
            <li>Z-DNA & Slipped DNA</li>
            <li>R-loops & Cruciforms</li>
            <li>Triplex DNA & H-DNA</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
    
    with cols[1]:
        st.markdown("""
        <div class="feature-card">
        <h3>🔬 Specialized Motifs</h3>
        <ul>
            <li>G-Quadruplexes (various types)</li>
            <li>i-Motifs & Hybrid Structures</li>
            <li>Non-B DNA Hotspots</li>
            <li>Multimeric structures</li>
        </ul>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("""
    <div style="text-align: center; margin-top: 20px;">
    <h4>Get started by uploading your sequence or using our example</h4>
    </div>
    """, unsafe_allow_html=True)

# Page: Upload & Analyze
elif page == "Upload & Analyze":
    st.markdown("<h2 style='color:#1A5276;'>Upload & Analyze</h2>", unsafe_allow_html=True)
    
    with st.expander("Sequence Input Options", expanded=True):
        tab1, tab2, tab3 = st.tabs(["File Upload", "Example Sequence", "Paste Sequence"])
        
        with tab1:
            fasta_file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"])
            if fasta_file:
                try:
                    seq = parse_fasta(fasta_file.read().decode("utf-8"))
                    st.session_state['seq'] = seq
                    st.success(f"FASTA loaded. Sequence length: {len(seq)} bp")
                except Exception as e:
                    st.error(f"Error reading file: {str(e)}")
        
        with tab2:
            if st.button("Load Example Sequence"):
                st.session_state['seq'] = parse_fasta(EXAMPLE_FASTA)
                st.success(f"Example loaded. Sequence length: {len(st.session_state['seq'])} bp")
                st.code(EXAMPLE_FASTA, language="fasta")
        
        with tab3:
            seq_input = st.text_area("Paste Sequence (FASTA or raw)", 
                                   value=st.session_state.get('seq', ''), 
                                   height=150)
            if seq_input:
                try:
                    st.session_state['seq'] = parse_fasta(seq_input)
                    st.success(f"Sequence parsed. Length: {len(st.session_state['seq'])} bp")
                except Exception as e:
                    st.error(f"Invalid sequence: {str(e)}")
    
    if st.session_state.get('seq'):
        st.markdown("---")
        st.subheader("Sequence Preview")
        st.text(wrap(st.session_state['seq'][:500]))  # Show first 500bp
        st.caption(f"Total length: {len(st.session_state['seq'])} bp | GC content: {gc_content(st.session_state['seq']):.1f}%")
        
        if st.button("Run Motif Analysis", type="primary"):
            seq = st.session_state['seq']
            if not seq or not re.match("^[ATGC]+$", seq.upper()):
                st.error("Invalid sequence. Only A/T/G/C characters allowed.")
            else:
                with st.spinner("Analyzing sequence for non-B DNA motifs..."):
                    try:
                        results = collect_all_motifs(seq)
                        if results:
                            st.session_state['motif_results'] = results
                            st.session_state['df'] = pd.DataFrame(results)
                            st.success(f"Analysis complete! Found {len(results)} motifs.")
                        else:
                            st.warning("Analysis complete but no motifs found.")
                            st.session_state['motif_results'] = []
                            st.session_state['df'] = pd.DataFrame()
                    except Exception as e:
                        st.error(f"Analysis failed: {str(e)}")

# Page: Results
elif page == "Results":
    st.markdown("<h2 style='color:#1A5276;'>Detected Motifs</h2>", unsafe_allow_html=True)
    
    if st.session_state.get('df', pd.DataFrame()).empty:
        st.info("No results available. Please run analysis first.")
    else:
        df = st.session_state['df']
        
        # Summary stats
        col1, col2, col3 = st.columns(3)
        col1.metric("Total Motifs", len(df))
        col2.metric("Unique Types", df['Subtype'].nunique())
        col3.metric("Sequence Coverage", 
                   f"{sum(df['Length']) / len(st.session_state['seq']) * 100:.1f}%")
        
        # Main dataframe
        st.dataframe(
            df[['Class', 'Subtype', 'Start', 'End', 'Length', 'Sequence', 'ScoreMethod', 'Score']],
            use_container_width=True,
            height=400
        )
        
        # Visualizations
        st.markdown("---")
        col1, col2 = st.columns(2)
        
        with col1:
            st.subheader("Motif Type Distribution")
            fig, ax = plt.subplots(figsize=(8, 4))
            sns.countplot(data=df, y='Subtype', ax=ax, 
                         order=df['Subtype'].value_counts().index)
            ax.set_xlabel("Count")
            ax.set_ylabel("")
            plt.tight_layout()
            st.pyplot(fig)
        
        with col2:
            st.subheader("Length Distribution")
            fig, ax = plt.subplots(figsize=(8, 4))
            sns.boxplot(data=df, x='Length', y='Class', ax=ax)
            ax.set_xlabel("Length (bp)")
            ax.set_ylabel("")
            plt.tight_layout()
            st.pyplot(fig)
        
        # Hotspots
        st.markdown("---")
        st.subheader("Hotspot Regions (≥3 motifs in 100 bp)")
        hotspots = find_hotspots(
            st.session_state['seq'], 
            st.session_state['motif_results'], 
            window=100, 
            min_count=3
        )
        
        if hotspots:
            hotspot_df = pd.DataFrame(hotspots)
            st.dataframe(hotspot_df)
            
            # Visualize hotspots on sequence
            fig, ax = plt.subplots(figsize=(10, 2))
            for _, row in hotspot_df.iterrows():
                ax.axvspan(row['RegionStart'], row['RegionEnd'], 
                           alpha=0.3, color='red')
            ax.set_xlim(0, len(st.session_state['seq']))
            ax.set_xlabel("Position (bp)")
            ax.set_yticks([])
            ax.set_title("Hotspot Regions Along Sequence")
            st.pyplot(fig)
        else:
            st.info("No hotspot regions found.")

# Page: Visualization
elif page == "Visualization":
    st.markdown("<h2 style='color:#1A5276;'>Motif Visualization</h2>", unsafe_allow_html=True)
    
    if st.session_state.get('df', pd.DataFrame()).empty:
        st.info("No results available. Please run analysis first.")
    else:
        df = st.session_state['df']
        seq = st.session_state['seq']
        
        # Create motif map
        motif_types = sorted(df['Subtype'].unique())
        color_palette = sns.color_palette('husl', n_colors=len(motif_types))
        color_map = {typ: color_palette[i] for i, typ in enumerate(motif_types)}
        y_map = {typ: i+1 for i, typ in enumerate(motif_types)}
        
        fig, ax = plt.subplots(figsize=(12, len(motif_types)*0.7 + 2))
        
        # Plot each motif
        for _, motif in df.iterrows():
            motif_type = motif['Subtype']
            y = y_map[motif_type]
            color = color_map[motif_type]
            ax.hlines(y, motif['Start'], motif['End'], 
                     color=color, linewidth=8, alpha=0.8)
        
        # Customize plot
        ax.set_yticks(list(y_map.values()))
        ax.set_yticklabels(list(y_map.keys()))
        ax.set_xlim(0, len(seq)+1)
        ax.set_xlabel('Position on Sequence (bp)')
        ax.set_title('Motif Map (Full Sequence)')
        sns.despine(left=False, bottom=False)
        plt.tight_layout()
        
        # Display in Streamlit
        st.pyplot(fig)
        
        # Add legend
        st.markdown("**Color Legend:**")
        cols = st.columns(4)
        for i, (motif, color) in enumerate(color_map.items()):
            with cols[i % 4]:
                st.markdown(
                    f"<span style='color:rgb({color[0]*255:.0f},{color[1]*255:.0f},{color[2]*255:.0f})'>"
                    f"■ {motif}</span>", 
                    unsafe_allow_html=True
                )

# Page: Download
elif page == "Download":
    st.markdown("<h2 style='color:#1A5276;'>Download Results</h2>", unsafe_allow_html=True)
    
    if st.session_state.get('df', pd.DataFrame()).empty:
        st.info("No results available to download. Please run analysis first.")
    else:
        df = st.session_state['df']
        
        # CSV Download
        csv = df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name="motif_results.csv",
            mime="text/csv",
            help="Download results as CSV file"
        )
        
        # Excel Download
        excel_buffer = io.BytesIO()
        with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
            df.to_excel(writer, index=False, sheet_name='Motif_Results')
            writer.close()
            
            st.download_button(
                label="Download Excel",
                data=excel_buffer.getvalue(),
                file_name="motif_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
                help="Download results as Excel file"
            )
        
        # FASTA Download
        if st.session_state.get('seq'):
            fasta_content = f">Analyzed_sequence\n{wrap(st.session_state['seq'])}"
            st.download_button(
                label="Download Sequence (FASTA)",
                data=fasta_content,
                file_name="analyzed_sequence.fasta",
                mime="text/plain",
                help="Download the analyzed sequence in FASTA format"
            )

# Page: Additional Information
elif page == "Additional Information":
    st.markdown("<h2 style='color:#1A5276;'>Additional Information</h2>", unsafe_allow_html=True)
    
    st.markdown("""
    ## About This Tool
    
    The Non-B DNA Motif Finder is designed to identify various non-canonical DNA structures 
    that play important roles in genomic regulation, stability, and disease.
    
    ### Key Features:
    - Comprehensive detection of 10+ non-B DNA structure types
    - Visualization of motif locations and hotspots
    - Export capabilities for further analysis
    - User-friendly interface with multiple input options
    
    ### Methodology:
    The tool uses regular expression patterns combined with specialized scoring algorithms 
    (like G4Hunter and Z-Seeker) to identify potential non-B DNA forming regions.
    
    ## References
    
    1. Yella VR, et al. (2018) *Nucleic Acids Research*  
    2. Puig Lombardi L, et al. (2019) *NAR Genomics and Bioinformatics*  
    3. Cer RZ, et al. (2013) *Nucleic Acids Research*
    
    ## Contact
    
    For questions or feedback, please contact:
    - Dr. Venkata Rajesh Yella: vryella@example.com
    - Chandrika Gummadi: cgummadi@example.com
    
    [GitHub Repository](https://github.com/VRYella/Non-B-DNA-Finder)
    """)

# Page: Motif Definitions Glossary
elif page == "Motif Definitions Glossary":
    st.markdown("<h2 style='color:#1A5276;'>Motif Definitions Glossary</h2>", unsafe_allow_html=True)
    
    with st.expander("Curved DNA / APRs", expanded=False):
        st.markdown("""
        **A-Phased Repeats (APRs):**  
        DNA sequences with periodically repeated AA/TT dinucleotides that introduce intrinsic curvature.
        Associated with nucleosome positioning and transcriptional regulation.
        """)
    
    with st.expander("Z-DNA", expanded=False):
        st.markdown("""
        **Z-DNA:**  
        Left-handed double helix formed by alternating purine-pyrimidine sequences (especially CG repeats).
        Important in transcriptional regulation and immune response.
        """)
    
    with st.expander("Slipped DNA", expanded=False):
        st.markdown("""
        **Slipped DNA:**  
        Structures formed when tandem repeats misalign during replication or repair.
        Associated with trinucleotide repeat expansion disorders.
        """)
    
    with st.expander("R-Loops", expanded=False):
        st.markdown("""
        **R-Loops:**  
        Three-stranded nucleic acid structures with an RNA-DNA hybrid and displaced single-stranded DNA.
        Play roles in transcription regulation and genome instability.
        """)
    
    with st.expander("Cruciforms & Hairpins", expanded=False):
        st.markdown("""
        **Cruciform DNA:**  
        Extruded secondary structures formed by inverted repeats.
        **Hairpins:**  
        Stem-loop structures formed by palindromic sequences.
        """)
    
    with st.expander("Triplex DNA & H-DNA", expanded=False):
        st.markdown("""
        **Triplex DNA:**  
        Triple-helical structures formed through Hoogsteen base pairing.
        **H-DNA:**  
        Intramolecular triplex formed by mirror repeat sequences.
        """)
    
    with st.expander("G-Quadruplexes", expanded=False):
        st.markdown("""
        **G-Quadruplexes (G4s):**  
        Four-stranded structures formed by stacked G-tetrads, stabilized by monovalent cations.
        Found in promoter regions, telomeres, and associated with cancer and neurodegenerative diseases.
        
        Variants include:
        - Canonical G4s
        - Relaxed G4s (longer loops)
        - Bulged G4s (interrupted G-runs)
        - Bipartite G4s (two separate G4 units)
        - Multimeric G4s (multiple contiguous G4 units)
        """)
    
    with st.expander("i-Motifs", expanded=False):
        st.markdown("""
        **i-Motifs:**  
        C-rich quadruplex structures formed in acidic conditions, complementary to G-quadruplexes.
        Potential roles in gene regulation and pH sensing.
        """)
    
    with st.expander("Hybrids & Hotspots", expanded=False):
        st.markdown("""
        **Hybrid Motifs:**  
        Genomic regions capable of forming multiple non-B DNA structures simultaneously.
        
        **Non-B Hotspots:**  
        Genomic regions with ≥3 different non-B motifs within 100 bp.
        Associated with genomic instability and disease.
        """)
