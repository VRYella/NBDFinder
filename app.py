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

MOTIF_CLASSES = {
    "Curved_DNA": "#FF9AA2",
    "Z-DNA": "#FFB7B2",
    "Slipped_DNA": "#FFDAC1",
    "Cruciform": "#E2F0CB",
    "Triplex_DNA": "#B5EAD7",
    "Sticky_DNA": "#DCB8CB",
    "G-Triplex": "#C7CEEA",
    "G4": "#A2D7D8",
    "i-Motif": "#B0C4DE",
    "R-Loop": "#FFD3B6",
    "Hybrid": "#C1A192",         # New
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
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", list(PAGES.keys()))

# Main title
st.title("Non-B DNA Motif Finder")
st.caption("Comprehensive detection of 12 non-canonical DNA structure types")

# Load and display the image on the home page
try:
    nbd_image = Image.open("nbd.PNG")
    st.image(nbd_image, caption="Non-B DNA Structures Overview", use_container_width=True)
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
                )
            }
        )

        # Visualization section
        st.subheader("📈 Distribution Analysis")
        tab1, tab2, tab3 = st.tabs(["By Type", "By Length", "By Score"])
        
        with tab1:
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.countplot(
                data=df,
                y='Class',
                order=df['Class'].value_counts().index,
                palette=MOTIF_CLASSES.values()
            )
            ax.set_title("Motif Count by Class")
            st.pyplot(fig)
        
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

# Page: Visualization
elif page == "Visualization":
    st.header("Interactive Motif Visualization")
    
    if st.session_state.df.empty:
        st.info("No results to visualize. Please run analysis first.")
    else:
        df = st.session_state.df
        seq_len = len(st.session_state.seq)
        
        # Interactive controls
        st.sidebar.subheader("Visualization Settings")
        show_classes = st.sidebar.multiselect(
            "Select motif classes to display:",
            sorted(df['Class'].unique()),
            default=sorted(df['Class'].unique())
        )
        
        min_score = st.sidebar.slider(
            "Minimum confidence score:",
            0.0, 1.0, 0.5
        )
        
        position_range = st.sidebar.slider(
            "Sequence position range:",
            0, seq_len, (0, min(5000, seq_len))
        )
        
        # Filter data
        viz_df = df[
            (df['Class'].isin(show_classes)) & 
            (df['Score'] >= min_score) &
            (df['Start'] >= position_range[0]) & 
            (df['End'] <= position_range[1])
        ].copy()
        
        if viz_df.empty:
            st.warning("No motifs match the selected filters")
        else:
            # Create motif map
            st.subheader("Genome Browser View")
            
            # Calculate y-positions
            subtypes = sorted(viz_df['Subtype'].unique())
            y_pos = {subtype: i+1 for i, subtype in enumerate(subtypes)}
            
            # Create figure
            fig, ax = plt.subplots(figsize=(15, 8))
            
            # Plot each motif
            for _, row in viz_df.iterrows():
                ax.hlines(
                    y_pos[row['Subtype']],
                    row['Start'],
                    row['End'],
                    linewidth=10,
                    color=MOTIF_CLASSES[row['Class']],
                    alpha=0.8
                )
                # Add score as text
                ax.text(
                    (row['Start'] + row['End'])/2,
                    y_pos[row['Subtype']] + 0.1,
                    f"{row['Score']:.2f}",
                    ha='center',
                    fontsize=8
                )
            
            # Customize plot
            ax.set_yticks(list(y_pos.values()))
            ax.set_yticklabels(list(y_pos.keys()))
            ax.set_xlim(position_range[0], position_range[1])
            ax.set_xlabel("Sequence Position (bp)")
            ax.set_title(f"Non-B DNA Motifs ({position_range[0]:,}-{position_range[1]:,} bp)")
            plt.tight_layout()
            st.pyplot(fig)
            
            # Add interactive table below visualization
            st.dataframe(
                viz_df[['Class', 'Subtype', 'Start', 'End', 'Length', 'Score']],
                height=300
            )

# Page: Download
elif page == "Download":
    st.header("Download Results")
    
    if st.session_state.df.empty:
        st.info("No results available to download")
    else:
        # Prepare data
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        
        # CSV Download
        csv = st.session_state.df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="Download CSV (All Motifs)",
            data=csv,
            file_name=f"nonb_motifs_{timestamp}.csv",
            mime="text/csv"
        )
        
        # Excel Download with multiple sheets
        excel_buffer = io.BytesIO()
        with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
            st.session_state.df.to_excel(writer, index=False, sheet_name='Motifs')
            if st.session_state.hotspots:
                pd.DataFrame(st.session_state.hotspots).to_excel(
                    writer, index=False, sheet_name='Hotspots'
                )
            # Add sequence info
            pd.DataFrame({
                'Sequence Info': [
                    f"Length: {len(st.session_state.seq)} bp",
                    f"GC Content: {gc_content(st.session_state.seq):.1f}%",
                    f"Motifs Found: {len(st.session_state.df)}",
                    f"Analysis Date: {timestamp}"
                ]
            }).to_excel(writer, index=False, sheet_name='Summary')
            writer.close()
            
            st.download_button(
                label="Download Excel Workbook",
                data=excel_buffer.getvalue(),
                file_name=f"nonb_results_{timestamp}.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        
        # FASTA Download
        fasta_content = f">Analyzed_sequence\n{wrap(st.session_state.seq)}"
        st.download_button(
            label="Download Sequence (FASTA)",
            data=fasta_content,
            file_name=f"sequence_{timestamp}.fasta",
            mime="text/plain"
        )

# Page: Documentation
elif page == "Documentation":
    st.header("Scientific Documentation")
    
    with st.expander("All 12 Motif Types", expanded=True):
        st.markdown("""
        | Motif Type | Detection Method | Key References |
        |------------|------------------|----------------|
        | Curved DNA | A-tract phasing | Brukner et al. 1995 |
        | Z-DNA | Alternating Pu/Py | Ho et al. 2010 |
        | Slipped DNA | Direct repeats | Bacolla et al. 2006 |
        | Cruciform | Inverted repeats | Lilley 1985 |
        | Triplex DNA | Mirror repeats | Mirkin 1994 |
        | G-Triplex | Three G-runs | Karsisiotis 2011 |
        | G4 | G4Hunter | Bedrat et al. 2016 |
        | i-Motif | C-rich sequences | Zeraati et al. 2018 |
        | R-Loop | GC skew + G-clusters | Sanz et al. 2016 |
        | Sticky DNA | (GAA/TTC)n | Potaman et al. 2003 |
        | A-Phased Repeats | 10.5bp spacing | Trifonov 1980 |
        | Mirror Repeats | Self-complementary | Frank-Kamenetskii 1990 |
        """)
    
    with st.expander("Scoring Systems"):
        st.markdown("""
        ### G4Hunter Score
        ```
        score = mean(G-run contributions) - mean(C-run penalties)
        ```
        **Thresholds:**
        - Canonical: ≥1.2
        - Relaxed: ≥0.8
        - Bulged: ≥1.0
        
        ### Z-DNA Score
        ```
        score = (CG_pairs/5) + (total_alternating/15)
        ```
        
        ### i-Motif Score
        ```
        score = (sum(C_tracts)/16) + (C_content/2)
        ```
        """)
    
    st.markdown("""
    ## References
    1. Bedrat et al. (2016) Nucleic Acids Research  
    2. Ho et al. (2010) Nature Chemical Biology  
    3. Zeraati et al. (2018) Nature Chemistry  
    4. Bacolla et al. (2006) Nucleic Acids Research  
    5. Mirkin & Frank-Kamenetskii (1994) Annual Review of Biophysics
    """)
