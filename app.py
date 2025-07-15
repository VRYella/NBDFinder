import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import io
from datetime import datetime
from PIL import Image

# Motif functions
try:
    from motifs import (
        all_motifs, 
        find_hotspots,
        parse_fasta, wrap, gc_content, reverse_complement,
        select_best_nonoverlapping_motifs
    )
except ImportError as e:
    st.error(f"Critical Import Error: {str(e)}")
    st.stop()

st.set_page_config(
    page_title="Non-B DNA Motif Finder",
    layout="wide",
    page_icon="🧬",
    menu_items={'About': "Non-B DNA Motif Finder | Developed by Dr. Venkata Rajesh Yella"}
)

EXAMPLE_FASTA = """>Example_Sequence
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA
CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
"""

# Session state initialization
for k, v in {
    'seq': "",
    'df': pd.DataFrame(),
    'motif_results': [],
    'motif_results_nonoverlap': [],
    'analysis_status': "Ready",
    'hotspots': []
}.items():
    if k not in st.session_state:
        st.session_state[k] = v

MOTIF_CLASSES = {
    "Curved_DNA": "#FF9AA2",
    "Z-DNA": "#FFB7B2",
    "Slipped_DNA": "#FFDAC1",
    "R-Loop": "#FFD3B6",
    "Cruciform": "#E2F0CB",
    "Triplex_DNA": "#B5EAD7",
    "Sticky_DNA": "#DCB8CB",
    "G-Triplex": "#C7CEEA",
    "G-quadruplex": "#A2D7D8",
    "i-Motif": "#B0C4DE",    
    "Hybrid": "#C1A192",
    "Non-B-DNA Clusters": "#A2C8CC"
    
}

PAGES = {
    "Home": "Introduction and overview",
    "Upload & Analyze": "Submit DNA sequence for analysis",
    "Results": "View detected motifs and statistics", 
    "Visualization": "Graphical representation of motifs",
    "Download": "Export results for further analysis",
    "Documentation": "Scientific methods and references"
}

# --- Sidebar navigation ---
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", list(PAGES.keys()))

# --- Developer info: fixed at sidebar bottom ---
st.markdown(
    """
    <style>
    .dev-footer {
        position: fixed;
        bottom: 0;
        left: 0;
        width: 22rem;
        padding: 14px 18px 14px 20px;
        background: #f7f7fc;
        border-top: 1px solid #e0e0e0;
        border-right: 1px solid #e0e0e0;
        border-bottom-left-radius: 12px;
        font-size: 15px;
        color: #1e293b;
        z-index: 100;
    }
    </style>
    <div class='dev-footer'>
        <b>Developed by</b><br>
        Dr. Venkata Rajesh Yella<br>
        <a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a><br>
        <a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
    </div>
    """,
    unsafe_allow_html=True
)

# --- Home page ---
if page == "Home":
    try:
        nbd_image = Image.open("nbd.PNG")
        st.image(nbd_image, use_container_width=True)
    except Exception:
        pass

    st.title("Non-B DNA Motif Finder")
    st.markdown("This tool identifies **12 classes** of non-canonical DNA structures using published algorithms:")

    # Motif cards
    cols = st.columns(4)
    for i, (motif, color) in enumerate(MOTIF_CLASSES.items()):
        with cols[i % 4]:
            st.markdown(
                f"<div style='background:{color};padding:10px 0 10px 0;border-radius:7px;margin-bottom:12px;text-align:center;font-weight:bold;font-size:1.1em;box-shadow:0 2px 10px #eee;'>"
                f"{motif.replace('_',' ')}</div>", unsafe_allow_html=True)

    # Features & usage paragraph
    st.markdown(
        """
        <div style='margin-top: 34px; font-size: 18px; line-height: 1.7; background: #f1f8fa; border-radius: 8px; padding: 22px 22px 18px 22px; box-shadow: 0px 2px 10px #e0e5ea;'>
        <b>The Non-B DNA Motif Finder</b> provides comprehensive detection of 12 distinct non-canonical DNA structure types, employing scientifically validated algorithms and established thresholds. Users benefit from interactive visualizations resembling a genome browser and versatile export options including CSV, Excel, and images.<br><br>
        <b>How to use:</b> Simply upload or paste your DNA sequence, execute the analysis, explore interactive visual representations, and download the detailed results for further examination.
        </div>
        """, unsafe_allow_html=True
    )

# --- Upload & Analyze ---
elif page == "Upload & Analyze":
    st.header("Sequence Input")
    with st.expander("Input Options", expanded=True):
        input_method = st.radio("Select input method:", 
            ["File Upload", "Example Sequence", "Paste Sequence"])
        if input_method == "File Upload":
            fasta_file = st.file_uploader("Upload FASTA file", type=["fa", "fasta", "txt"])
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
            seq_input = st.text_area("Paste DNA Sequence (FASTA or raw):", height=150)
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
            st.text(wrap(st.session_state.seq[:500]))
        with col2:
            st.metric("GC Content", f"{gc_content(st.session_state.seq):.1f}%")
            st.metric("Sequence Length", f"{len(st.session_state.seq):,} bp")
        if st.button("Run Full Analysis", type="primary"):
            with st.spinner("Analyzing sequence for 12 motif types..."):
                try:
                    st.session_state.motif_results = all_motifs(st.session_state.seq)
                    st.session_state.df = pd.DataFrame(st.session_state.motif_results)
                    st.session_state.motif_results_nonoverlap = select_best_nonoverlapping_motifs(
                        st.session_state.motif_results
                    )
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

# --- Results page with improved Hotspot visualization ---
elif page == "Results":
    st.header("Analysis Results")
    if not st.session_state.motif_results_nonoverlap:
        st.session_state.motif_results_nonoverlap = select_best_nonoverlapping_motifs(
            st.session_state.motif_results
        )
    df = pd.DataFrame(st.session_state.motif_results_nonoverlap)
    if df.empty:
        st.info("No results available. Please run analysis first.")
    else:
        df['Score'] = pd.to_numeric(df['Score'], errors='coerce')
        df = df.dropna(subset=['Score'])
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
        st.subheader("🧬 Detected Motifs")
        show_cols = ['Class', 'Subtype', 'Start', 'End', 'Length', 'Score', 'Sequence']
        st.dataframe(
            df[show_cols],
            use_container_width=True,
            height=400,
            column_config={
                "Score": st.column_config.ProgressColumn(
                    "Score", format="%.2f",
                    min_value=0,
                    max_value=max(1, df['Score'].max())
                )
            }
        )
        st.subheader("📈 Distribution Analysis")
        tab1, tab2, tab3, tab4 = st.tabs(["By Type", "By Length", "By Score", "Motif Density"])
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
        with tab4:
            fig, ax = plt.subplots(figsize=(12, 3))
            pos = []
            for _, row in df.iterrows():
                pos.extend(range(row['Start'], row['End'] + 1))
            if pos:
                sns.kdeplot(pos, fill=True, color='navy', ax=ax)
                ax.set_xlim(0, len(st.session_state.seq))
                ax.set_xlabel("Sequence Position (bp)")
                ax.set_ylabel("Motif Density")
                ax.set_title("Motif Density Along Sequence")
            st.pyplot(fig)
        st.subheader("🔥 Non-B DNA Clustered Regions (Hotspots)")
        if st.session_state.hotspots:
            hotspot_df = pd.DataFrame(st.session_state.hotspots)
            st.dataframe(
                hotspot_df.sort_values('Score', ascending=False),
                use_container_width=True
            )
            # Colorful, labeled blocks using tab20
            palette = sns.color_palette("tab20", len(hotspot_df))
            fig, ax = plt.subplots(figsize=(12, 3))
            for idx, (_, row) in enumerate(hotspot_df.iterrows()):
                ax.axvspan(
                    row['RegionStart'], 
                    row['RegionEnd'], 
                    alpha=0.5, 
                    color=palette[idx],
                    label=f"Cluster {idx+1}"
                )
                ax.text(
                    (row['RegionStart'] + row['RegionEnd'])/2,
                    0.7 + 0.1*(idx%2),
                    f"{row['RegionStart']}–{row['RegionEnd']}",
                    color='black',
                    fontsize=9,
                    ha='center',
                    va='bottom',
                    rotation=45
                )
            ax.set_xlim(0, len(st.session_state.seq))
            ax.set_xlabel("Sequence Position (bp)")
            ax.set_yticks([])
            ax.set_title("Non-B DNA Clustered Regions (Hotspots)")
            ax.legend(loc='upper right', fontsize=8)
            st.pyplot(fig)
            st.info("These regions represent clusters of non-B DNA motifs, i.e., hotspots where multiple motif types overlap or are concentrated.")
        else:
            st.info("No hotspot regions detected.")

# --- Visualization page: Colorful, score-free tracks ---
elif page == "Visualization":
    st.header("Genome Browser View (Interactive Motif Visualization)")
    if not st.session_state.motif_results_nonoverlap:
        st.session_state.motif_results_nonoverlap = select_best_nonoverlapping_motifs(
            st.session_state.motif_results
        )
    df = pd.DataFrame(st.session_state.motif_results_nonoverlap)
    if df.empty:
        st.info("No results to visualize. Please run analysis first.")
    else:
        df['Score'] = pd.to_numeric(df['Score'], errors='coerce')
        df = df.dropna(subset=['Score'])
        seq_len = len(st.session_state.seq)
        st.sidebar.subheader("Visualization Settings")
        show_classes = st.sidebar.multiselect(
            "Select motif classes to display:",
            sorted(df['Class'].unique()),
            default=sorted(df['Class'].unique())
        )
        min_score = st.sidebar.slider(
            "Minimum confidence score:",
            float(df['Score'].min()), float(df['Score'].max()), float(df['Score'].min())
        )
        position_range = st.sidebar.slider(
            "Sequence position range:",
            0, seq_len, (0, min(5000, seq_len))
        )
        viz_df = df[
            (df['Class'].isin(show_classes)) & 
            (df['Score'] >= min_score) &
            (df['Start'] >= position_range[0]) & 
            (df['End'] <= position_range[1])
        ].copy()
        if viz_df.empty:
            st.warning("No motifs match the selected filters")
        else:
            subtypes = sorted(viz_df['Subtype'].unique())
            y_pos = {subtype: i+1 for i, subtype in enumerate(subtypes)}
            color_map = dict(zip(subtypes, sns.color_palette("tab20", len(subtypes))))
            fig, ax = plt.subplots(figsize=(15, max(6, len(subtypes)//2+2)))
            for _, row in viz_df.iterrows():
                ax.hlines(
                    y_pos[row['Subtype']],
                    row['Start'],
                    row['End'],
                    linewidth=13,
                    color=color_map.get(row['Subtype'], "#888"),
                    alpha=0.88
                )
            ax.set_yticks(list(y_pos.values()))
            ax.set_yticklabels(list(y_pos.keys()))
            ax.set_xlim(position_range[0], position_range[1])
            ax.set_xlabel("Sequence Position (bp)")
            ax.set_title(f"Non-B DNA Motifs ({position_range[0]:,}-{position_range[1]:,} bp)")
            plt.tight_layout()
            st.pyplot(fig)
            st.dataframe(
                viz_df[['Class', 'Subtype', 'Start', 'End', 'Length', 'Score']],
                height=300
            )
            st.markdown(
                "<div style='display:flex;flex-wrap:wrap;gap:10px;'>"
                + "".join(
                    f"<div style='background:{color_map.get(s, '#ccc')};padding:7px 15px;border-radius:3px;margin:2px 5px;color:black;'>"
                    f"{s}</div>" for s in subtypes
                )
                + "</div>", unsafe_allow_html=True
            )

# --- Download page ---
elif page == "Download":
    st.header("Download Results")
    if not st.session_state.motif_results_nonoverlap:
        st.session_state.motif_results_nonoverlap = select_best_nonoverlapping_motifs(
            st.session_state.motif_results
        )
    df = pd.DataFrame(st.session_state.motif_results_nonoverlap)
    df['Score'] = pd.to_numeric(df['Score'], errors='coerce')
    df = df.dropna(subset=['Score'])
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    csv = df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Download CSV (Non-overlapping Motifs)",
        data=csv,
        file_name=f"nonb_motifs_nonoverlap_{timestamp}.csv",
        mime="text/csv"
    )
    st.markdown("For advanced analysis, you may download all motifs (overlapping):")
    all_df = pd.DataFrame(st.session_state.motif_results)
    all_df['Score'] = pd.to_numeric(all_df['Score'], errors='coerce')
    all_df = all_df.dropna(subset=['Score'])
    all_csv = all_df.to_csv(index=False).encode('utf-8')
    st.download_button(
        label="Download CSV (All Motifs - Overlapping)",
        data=all_csv,
        file_name=f"nonb_motifs_ALL_{timestamp}.csv",
        mime="text/csv"
    )
    # Excel
    excel_buffer = io.BytesIO()
    with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False, sheet_name='Motifs')
        if st.session_state.hotspots:
            pd.DataFrame(st.session_state.hotspots).to_excel(
                writer, index=False, sheet_name='Hotspots'
            )
        pd.DataFrame({
            'Sequence Info': [
                f"Length: {len(st.session_state.seq)} bp",
                f"GC Content: {gc_content(st.session_state.seq):.1f}%",
                f"Motifs Found: {len(df)}",
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
    # FASTA
    fasta_content = f">Analyzed_sequence\n{wrap(st.session_state.seq)}"
    st.download_button(
        label="Download Sequence (FASTA)",
        data=fasta_content,
        file_name=f"sequence_{timestamp}.fasta",
        mime="text/plain"
    )

# --- Documentation page ---
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
