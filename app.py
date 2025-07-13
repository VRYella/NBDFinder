import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import io
from datetime import datetime
from typing import List, Dict
from PIL import Image

# Custom module imports with error handling
try:
    from motifs import all_motifs, find_hotspots
    from motifs import parse_fasta, wrap, gc_content, reverse_complement
except ImportError as e:
    st.error(f"Critical Import Error: {str(e)}")
    st.error("Please ensure motifs.py exists in the same directory with all required functions.")
    st.stop()

# Configure page
st.set_page_config(
    page_title="Non-B DNA Motif Finder",
    layout="wide",
    page_icon="🧬",
    menu_items={
        'About': "Non-B DNA Motif Finder | Developed by Dr. Venkata Rajesh Yella"
    }
)

# Example sequence with diverse motifs
EXAMPLE_FASTA = """>Example_Sequence
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA
CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
"""

# Initialize session state
if 'seq' not in st.session_state:
    st.session_state.seq = ""
if 'df' not in st.session_state:
    st.session_state.df = pd.DataFrame()
if 'motif_results' not in st.session_state:
    st.session_state.motif_results = []
if 'analysis_status' not in st.session_state:
    st.session_state.analysis_status = "Ready"
if 'hotspots' not in st.session_state:
    st.session_state.hotspots = []

# Define all 12 motif classes with colors
MOTIF_CLASSES = {
    "Curved_DNA": "#FF9AA2",
    "Z-DNA": "#FFB7B2",
    "Slipped_DNA": "#FFDAC1",
    "Cruciform": "#E2F0CB",
    "Triplex_DNA": "#B5EAD7",
    "G-Triplex": "#C7CEEA",
    "G4": "#A2D7D8",
    "i-Motif": "#B5EAD7",
    "R-Loop": "#FFD3B6",
    "Sticky_DNA": "#DCB8CB",
    "A-Phased_Repeat": "#A2C8CC",
    "Mirror_Repeat": "#D4A5A5"
}

# App pages
PAGES = {
    "Home": "Introduction and overview",
    "Upload & Analyze": "Submit DNA sequence for analysis",
    "Results": "View detected motifs and statistics", 
    "Visualization": "Graphical representation of motifs",
    "Download": "Export results for further analysis",
    "Documentation": "Scientific methods and references"
}

# Sidebar navigation
st.sidebar.title("🧬 Navigation")
page = st.sidebar.radio("Go to", list(PAGES.keys()) 

# Main title
st.title("Non-B DNA Motif Finder")
st.caption("Comprehensive detection of 12 non-canonical DNA structure types")

# Load and display the image on the home page
try:
    nbd_image = Image.open("nbd.PNG")
    st.image(nbd_image, caption="Non-B DNA Structures Overview", use_column_width=True)
except FileNotFoundError:
    st.warning("Image nbd.PNG not found. Please ensure it's in the same directory.")
except Exception as e:
    st.error(f"Error loading image: {str(e)}")

# Page: Home
if page == "Home":
    st.markdown("""
    ## Welcome to the Non-B DNA Motif Finder
    
    This tool identifies **12 classes** of non-canonical DNA structures using published algorithms:
    """)
    
    # Display motif classes with colors
    cols = st.columns(4)
    for i, (motif, color) in enumerate(MOTIF_CLASSES.items()):
        with cols[i % 4]:
            st.markdown(f"""
            <div style="background-color:{color}; padding:10px; border-radius:5px; margin-bottom:10px;">
                <b>{motif.replace('_', ' ')}</b>
            </div>
            """, unsafe_allow_html=True)
    
    st.markdown("""
    ### Key Features:
    - **Comprehensive Detection**: 12 non-B DNA structure types
    - **Scientific Validation**: Published algorithms and thresholds
    - **Interactive Visualization**: Genome browser-style display
    - **Export Capabilities**: CSV, Excel, and image exports
    
    ### How to Use:
    1. Upload or paste your DNA sequence
    2. Run the analysis
    3. Explore results through interactive visualizations
    4. Download data for further analysis
    """)

# Page: Upload & Analyze
elif page == "Upload & Analyze":
    st.header("Sequence Input")
    
    with st.expander("Input Options", expanded=True):
        input_method = st.radio("Select input method:", 
                              ["File Upload", "Example Sequence", "Paste Sequence"])
        
        if input_method == "File Upload":
            fasta_file = st.file_uploader("Upload FASTA file", 
                                         type=["fa", "fasta", "txt"])
            if fasta_file:
                try:
                    seq = parse_fasta(fasta_file.read().decode("utf-8"))
                    st.session_state.seq = seq
                    st.success(f"Loaded sequence: {len(seq):,} bp")
                except Exception as e:
                    st.error(f"Error reading file: {str(e)}")
        
        elif input_method == "Example Sequence":
            if st.button("Load Example"):
                st.session_state.seq = parse_fasta(EXAMPLE_FASTA)
                st.success(f"Example loaded: {len(st.session_state.seq):,} bp")
                st.code(EXAMPLE_FASTA, language="fasta")
        
        elif input_method == "Paste Sequence":
            seq_input = st.text_area("Paste DNA Sequence (FASTA or raw):", 
                                   height=150)
            if seq_input:
                try:
                    st.session_state.seq = parse_fasta(seq_input)
                    st.success(f"Sequence parsed: {len(st.session_state.seq):,} bp")
                except Exception as e:
                    st.error(f"Invalid sequence: {str(e)}")
    
    if st.session_state.seq:
        st.subheader("Sequence Preview")
        col1, col2 = st.columns(2)
        with col1:
            st.text(wrap(st.session_state.seq[:500]))  # First 500bp
        with col2:
            st.metric("GC Content", f"{gc_content(st.session_state.seq):.1f}%")
            st.metric("Sequence Length", f"{len(st.session_state.seq):,} bp")
        
        if st.button("Run Full Analysis", type="primary"):
            with st.spinner("Analyzing sequence for 12 motif types..."):
                try:
                    # Run motif detection
                    st.session_state.motif_results = all_motifs(st.session_state.seq)
                    st.session_state.df = pd.DataFrame(st.session_state.motif_results)
                    
                    # Find hotspots (corrected argument order)
                    st.session_state.hotspots = find_hotspots(
                        st.session_state.motif_results,
                        len(st.session_state.seq)
                    )
                    
                    if st.session_state.motif_results:
                        st.success(f"Found {len(st.session_state.motif_results)} motifs across {st.session_state.df['Class'].nunique()} classes")
                        st.session_state.analysis_status = "Complete"
                    else:
                        st.warning("No motifs detected")
                        st.session_state.analysis_status = "Complete"
                except Exception as e:
                    st.error(f"Analysis failed: {str(e)}")
                    st.session_state.analysis_status = "Error"

# [Rest of your pages (Results, Visualization, Download, Documentation)...]
# Continue with the same implementation as before for these pages
# Ensure all references to motif_results and hotspots use the session_state versions

# Page: Results
elif page == "Results":
    st.header("Analysis Results")
    
    if st.session_state.df.empty:
        st.info("No results available. Please run analysis first.")
    else:
        df = st.session_state.df.copy()
        
        # Convert Score column to numeric safely
        df['Score'] = pd.to_numeric(df['Score'], errors='coerce')
        df = df.dropna(subset=['Score'])
        
        # Summary statistics
        with st.expander("📊 Summary Statistics", expanded=True):
            cols = st.columns(4)
            cols[0].metric("Total Motifs", len(df))
            cols[1].metric("Unique Types", df['Subtype'].nunique())
            seq_coverage = sum(df['Length']) / len(st.session_state.seq) * 100
            cols[2].metric("Sequence Coverage", f"{seq_coverage:.1f}%")
            max_score = df['Score'].max()
            cols[3].metric("Top Score", f"{max_score:.2f}")
            
            st.progress(
                min(100, int(seq_coverage)),
                text=f"Sequence coverage: {seq_coverage:.1f}%"
            )

        # Main results table
        st.subheader("🧬 Detected Motifs")
        show_cols = ['Class', 'Subtype', 'Start', 'End', 'Length', 'Score', 'Sequence']
        st.dataframe(
            df[show_cols],
            use_container_width=True,
            height=400,
            column_config={
                "Score": st.column_config.ProgressColumn(
                    "Score",
                    format="%.2f",
                    min_value=0,
                    max_value=max(1, df['Score'].max())
                }
            }
        )

        # Visualization section
        # In your Results page section, replace the dataframe display code with:

# Main results table
st.subheader("🧬 Detected Motifs")
show_cols = ['Class', 'Subtype', 'Start', 'End', 'Length', 'Score', 'Sequence']
st.dataframe(
    df[show_cols],
    use_container_width=True,
    height=400,
    column_config={
        "Score": st.column_config.ProgressColumn(
            "Score",
            format="%.2f",
            min_value=0,
            max_value=max(1, df['Score'].max())
    }
)
        with tab2:
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.boxplot(
                data=df, 
                x='Length', 
                y='Class',
                palette=MOTIF_CLASSES.values()
            )
            ax.set_title("Motif Length Distribution")
            st.pyplot(fig)
            
        with tab3:
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.violinplot(
                data=df,
                x='Score',
                y='Class',
                palette=MOTIF_CLASSES.values()
            )
            ax.set_title("Score Distribution by Class")
            st.pyplot(fig)

        # Hotspot analysis
        st.subheader("🔥 Hotspot Regions")
        if st.session_state.hotspots:
            hotspot_df = pd.DataFrame(st.session_state.hotspots)
            st.dataframe(
                hotspot_df.sort_values('Score', ascending=False),
                use_container_width=True
            )
            
            fig, ax = plt.subplots(figsize=(12, 3))
            for _, row in hotspot_df.iterrows():
                ax.axvspan(
                    row['RegionStart'], 
                    row['RegionEnd'], 
                    alpha=0.3, 
                    color='red'
                )
            ax.set_xlim(0, len(st.session_state.seq))
            ax.set_xlabel("Sequence Position (bp)")
            ax.set_title("Hotspot Locations")
            st.pyplot(fig)
        else:
            st.info("No hotspot regions detected")

# [Continue with Visualization, Download, and Documentation pages...]
