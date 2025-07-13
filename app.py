import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
import io
from datetime import datetime
try:
    from .motifs import all_motifs, find_hotspots  # Relative import (if in package)
    from .utils import parse_fasta, wrap, gc_content
except ImportError:
    from motifs import all_motifs, find_hotspots  # Absolute import
    from utils import parse_fasta, wrap, gc_content
import os
import sys
from pathlib import Path
import streamlit as st



try:
    from motifs import all_motifs, find_hotspots
    from utils import parse_fasta, wrap, gc_content
except ImportError as e:
    st.error(f"Critical Import Error: {str(e)}")
    st.error("Please ensure motifs.py and utils.py exist in the same directory.")
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

# Example sequence
EXAMPLE_FASTA = """>Example
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG"""

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
page = st.sidebar.radio("Go to", list(PAGES.keys()), 
                       help=PAGES[list(PAGES.keys())[0]])

st.sidebar.markdown("---")
st.sidebar.info(
    "**Scientific Methods:**\n"
    "- G4Hunter algorithm (Bedrat et al. 2016)\n"
    "- Z-DNA predictor (Ho et al. 2010)\n"
    "- R-loop detection (Sanz et al. 2016)"
)

# Main title
st.title("Non-B DNA Motif Finder")
st.caption("Advanced detection of non-canonical DNA structures")

# Page: Home
if page == "Home":
    st.markdown("""
    <style>
        .feature-card {
            background-color: #f8f9fa;
            border-radius: 10px;
            padding: 20px;
            margin-bottom: 15px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
    </style>
    """, unsafe_allow_html=True)
    
    # Display logo or title
    try:
        st.image("nbd.PNG", use_container_width=True, 
                caption="Non-B DNA structural diversity")
    except:
        st.markdown("""
        <div style="text-align: center; margin-bottom: 30px;">
            <h1 style="color: #1A5276;">🧬 Non-B DNA Motif Finder</h1>
            <p>Advanced detection of non-canonical DNA structures</p>
        </div>
        """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="feature-card">
    <h3>🔬 Scientific Detection Methods</h3>
    <p>This tool identifies 12 classes of non-B DNA structures using published algorithms:</p>
    <ul>
        <li><b>G-Quadruplexes:</b> G4Hunter scoring (Bedrat et al. 2016)</li>
        <li><b>Z-DNA:</b> Alternating purine-pyrimidine detection</li>
        <li><b>R-loops:</b> GC-skew and G-cluster analysis</li>
        <li><b>i-Motifs:</b> C-rich sequence detection</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)
    
    st.markdown("""
    <div class="feature-card">
    <h3>📊 Analysis Pipeline</h3>
    <ol>
        <li>Sequence input (FASTA or raw sequence)</li>
        <li>Parallel motif scanning</li>
        <li>Structure validation</li>
        <li>Hotspot identification</li>
        <li>Interactive visualization</li>
    </ol>
    </div>
    """, unsafe_allow_html=True)

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
                    st.success(f"Loaded sequence: {len(seq)} bp")
                except Exception as e:
                    st.error(f"Error reading file: {str(e)}")
        
        elif input_method == "Example Sequence":
            if st.button("Load Example"):
                st.session_state.seq = parse_fasta(EXAMPLE_FASTA)
                st.success(f"Example loaded: {len(st.session_state.seq)} bp")
                st.code(EXAMPLE_FASTA, language="fasta")
        
        elif input_method == "Paste Sequence":
            seq_input = st.text_area("Paste DNA Sequence (FASTA or raw):", 
                                   height=150)
            if seq_input:
                try:
                    st.session_state.seq = parse_fasta(seq_input)
                    st.success(f"Sequence parsed: {len(st.session_state.seq)} bp")
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
            with st.spinner("Analyzing sequence..."):
                try:
                    # Run motif detection
                    st.session_state.motif_results = all_motifs(st.session_state.seq)
                    st.session_state.df = pd.DataFrame(st.session_state.motif_results)
                    
                    # Find hotspots
                    st.session_state.hotspots = find_hotspots(
                        st.session_state.seq, 
                        st.session_state.motif_results
                    )
                    
                    if st.session_state.motif_results:
                        st.success(f"Found {len(st.session_state.motif_results)} motifs")
                        st.session_state.analysis_status = "Complete"
                    else:
                        st.warning("No motifs detected")
                except Exception as e:
                    st.error(f"Analysis failed: {str(e)}")
                    st.session_state.analysis_status = "Error"

# Page: Results
elif page == "Results":
    st.header("Analysis Results")
    
    if st.session_state.df.empty or not hasattr(st.session_state, 'seq'):
        st.info("No results available. Please run analysis first.")
    else:
        df = st.session_state.df
        full_sequence = st.session_state.seq
        
        # Add sequence snippets to dataframe
        df['Sequence'] = df.apply(
            lambda row: full_sequence[row['Start']-1:row['End']], 
            axis=1
        )
        
        # Summary statistics with expandable details
        with st.expander("📊 Summary Statistics", expanded=True):
            col1, col2, col3, col4 = st.columns(4)
            col1.metric("Total Motifs", len(df))
            col2.metric("Unique Types", df['Subtype'].nunique())
            seq_coverage = sum(df['Length']) / len(full_sequence) * 100
            col3.metric("Sequence Coverage", f"{seq_coverage:.1f}%")
            col4.metric("Top Score", f"{df['Score'].max():.2f}")
            
            st.progress(min(100, int(seq_coverage)), 
                       text=f"Sequence covered by motifs: {seq_coverage:.1f}%")

        # Main results table with sequence display
        st.subheader("🧬 Detected Motifs")
        show_cols = st.multiselect(
            "Select columns to display:",
            df.columns.tolist(),
            default=['Class', 'Subtype', 'Start', 'End', 'Length', 'Score', 'Sequence']
        )
        
        st.dataframe(
            df[show_cols],
            use_container_width=True,
            height=400,
            column_config={
                "Score": st.column_config.ProgressColumn(
                    "Score",
                    help="Prediction confidence",
                    format="%.2f",
                    min_value=0,
                    max_value=1,
                ),
                "Sequence": st.column_config.TextColumn(
                    "Sequence",
                    help="Actual DNA sequence of motif",
                    width="medium"
                )
            }
        )

        # Interactive motif explorer
        st.subheader("🔍 Motif Explorer")
        selected_idx = st.selectbox(
            "Select a motif to inspect:",
            range(len(df)),
            format_func=lambda x: f"{df.iloc[x]['Class']} ({df.iloc[x]['Start']}-{df.iloc[x]['End']}) - Score: {df.iloc[x]['Score']:.2f}"
        )
        
        selected = df.iloc[selected_idx]
        context_size = st.slider("Context nucleotides to show", 0, 50, 10)
        
        # Display sequence with highlighted motif
        start_pos = max(0, selected['Start']-1 - context_size)
        end_pos = min(len(full_sequence), selected['End'] + context_size)
        context_seq = full_sequence[start_pos:end_pos]
        motif_start = selected['Start']-1 - start_pos
        motif_end = selected['End'] - start_pos
        
        st.markdown("**Sequence Context:**")
        st.code(
            f"{context_seq[:motif_start]}"
            f"\033[91m{context_seq[motif_start:motif_end]}\033[0m"
            f"{context_seq[motif_end:]}",
            language='text'
        )
        
        # Motif details in columns
        col1, col2, col3 = st.columns(3)
        with col1:
            st.markdown("**Position**")
            st.info(f"{selected['Start']}-{selected['End']}")
            st.markdown("**Length**")
            st.info(f"{selected['Length']} bp")
        with col2:
            st.markdown("**Type**")
            st.info(f"{selected['Class']} ({selected['Subtype']})")
            st.markdown("**Score**")
            st.info(f"{selected['Score']:.2f}")
        with col3:
            st.markdown("**Full Sequence**")
            st.text(wrap(selected['Sequence'], width=30))

        # Visualization section
        st.subheader("📈 Distribution Analysis")
        tab1, tab2, tab3 = st.tabs(["By Type", "By Length", "Score Distribution"])
        
        with tab1:
            fig = px.bar(
                df['Subtype'].value_counts().reset_index(),
                x='count',
                y='Subtype',
                orientation='h',
                title="Motif Count by Type"
            )
            st.plotly_chart(fig, use_container_width=True)
        
        with tab2:
            fig = px.box(
                df,
                x='Length',
                y='Class',
                title="Motif Length Distribution"
            )
            st.plotly_chart(fig, use_container_width=True)
            
        with tab3:
            fig = px.histogram(
                df,
                x='Score',
                color='Class',
                nbins=20,
                title="Score Distribution by Class"
            )
            st.plotly_chart(fig, use_container_width=True)

        # Hotspot analysis
        st.subheader("🔥 Hotspot Regions")
        if st.session_state.hotspots:
            hotspot_df = pd.DataFrame(st.session_state.hotspots)
            
            # Interactive hotspot explorer
            min_hotspot_score = st.slider(
                "Minimum hotspot score",
                0.0, 1.0, 0.5
            )
            filtered_hotspots = hotspot_df[hotspot_df['Score'].astype(float) >= min_hotspot_score]
            
            st.dataframe(
                filtered_hotspots.sort_values('Score', ascending=False),
                use_container_width=True
            )
            
            # Enhanced visualization
            fig = px.scatter(
                filtered_hotspots,
                x='RegionStart',
                y='Score',
                size='MotifCount',
                color='TypeDiversity',
                hover_data=['RegionStart', 'RegionEnd', 'MotifCount'],
                title="Hotspot Distribution"
            )
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.info("No hotspot regions detected")

        # Add sequence download option
        st.download_button(
            label="Download All Motif Sequences",
            data=df.to_csv(index=False),
            file_name="motif_sequences.csv",
            mime="text/csv"
        )
# Page: Visualization
elif page == "Visualization":
    st.header("Motif Visualization")
    
    if st.session_state.df.empty:
        st.info("No results to visualize. Please run analysis first.")
    else:
        df = st.session_state.df
        seq_len = len(st.session_state.seq)
        
        # Interactive controls
        st.sidebar.subheader("Visualization Settings")
        show_classes = st.sidebar.multiselect(
            "Select motif classes to display:",
            df['Class'].unique(),
            default=df['Class'].unique()
        )
        
        min_score = st.sidebar.slider(
            "Minimum confidence score:",
            0.0, 1.0, 0.7
        )
        
        # Filter data
        viz_df = df[
            (df['Class'].isin(show_classes)) & 
            (df['Score'].astype(float) >= min_score)
        ]
        
        if viz_df.empty:
            st.warning("No motifs match the selected filters")
        else:
            # Create motif map
            st.subheader("Motif Map")
            fig, ax = plt.subplots(figsize=(12, 6))
            
            # Assign y-positions for each subtype
            subtypes = sorted(viz_df['Subtype'].unique())
            y_pos = {subtype: i+1 for i, subtype in enumerate(subtypes)}
            colors = sns.color_palette("husl", len(subtypes))
            
            # Plot each motif
            for _, row in viz_df.iterrows():
                ax.hlines(
                    y_pos[row['Subtype']],
                    row['Start'],
                    row['End'],
                    linewidth=8,
                    color=colors[subtypes.index(row['Subtype'])],
                    alpha=0.7
                )
            
            # Customize plot
            ax.set_yticks(list(y_pos.values()))
            ax.set_yticklabels(list(y_pos.keys()))
            ax.set_xlim(0, seq_len)
            ax.set_xlabel("Sequence Position (bp)")
            ax.set_title("Non-B DNA Motif Distribution")
            plt.tight_layout()
            st.pyplot(fig)
            
            # Add legend
            st.markdown("**Color Legend:**")
            cols = st.columns(4)
            for i, subtype in enumerate(subtypes):
                with cols[i % 4]:
                    st.markdown(
                        f"<span style='color:{colors[i]}'>■ {subtype}</span>", 
                        unsafe_allow_html=True
                    )

# Page: Download
elif page == "Download":
    st.header("Download Results")
    
    if st.session_state.df.empty:
        st.info("No results available to download")
    else:
        # CSV Download
        csv = st.session_state.df.to_csv(index=False).encode('utf-8')
        st.download_button(
            label="Download CSV",
            data=csv,
            file_name="nonb_motifs.csv",
            mime="text/csv"
        )
        
        # Excel Download
        excel_buffer = io.BytesIO()
        with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
            st.session_state.df.to_excel(writer, index=False, sheet_name='Motifs')
            if st.session_state.hotspots:
                pd.DataFrame(st.session_state.hotspots).to_excel(
                    writer, index=False, sheet_name='Hotspots'
                )
            writer.close()
            st.download_button(
                label="Download Excel",
                data=excel_buffer.getvalue(),
                file_name="nonb_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        
        # FASTA Download
        if st.session_state.seq:
            fasta_content = f">Analyzed_sequence\n{wrap(st.session_state.seq)}"
            st.download_button(
                label="Download Sequence",
                data=fasta_content,
                file_name="sequence.fasta",
                mime="text/plain"
            )

# Page: Documentation
elif page == "Documentation":
    st.header("Scientific Documentation")
    
    with st.expander("Detection Methods", expanded=True):
        st.markdown("""
        ### G-Quadruplex Detection (G4Hunter)
        **Algorithm:**  
        Scores sequences based on G-richness and G-run arrangement:
        ```
        Score = mean([G-run scores] + [loop penalties])
        ```
        **Thresholds:**
        - Canonical: ≥1.2
        - Relaxed: ≥0.8
        - Bulged: ≥1.0 with penalty
        
        **Reference:**  
        Bedrat et al. (2016) *Nucleic Acids Research*
        """)
    
    with st.expander("Motif Definitions"):
        st.markdown("""
        | Motif Type | Pattern | Biological Significance |
        |------------|---------|-------------------------|
        | G4 | G≥3N1-7G≥3N1-7G≥3N1-7G≥3 | Gene regulation, telomere maintenance |
        | i-Motif | C≥3N0-7C≥3N0-7C≥3N0-7C≥3 | pH sensing, promoter regulation |
        | Z-DNA | (CG/GC/GT/TG/AC/CA)≥6 | Immune response, transcription |
        | R-loop | G-cluster + GC-rich region | Transcription termination, replication |
        """)
    
    with st.expander("Analysis Pipeline"):
        st.markdown("""
        ```mermaid
        graph TD
            A[Input Sequence] --> B(Preprocessing)
            B --> C[Parallel Motif Scanning]
            C --> D[Structure Validation]
            D --> E[Hotspot Identification]
            E --> F[Visualization]
            F --> G[Results Export]
        ```
        """)
    
    st.markdown("""
    ## References
    1. Bedrat et al. (2016) *Nucleic Acids Research*  
    2. Ho et al. (2010) *Nature Chemical Biology*  
    3. Zeraati et al. (2018) *Nature Chemistry*
    """)
