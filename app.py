import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import io
from datetime import datetime
from PIL import Image

# ---- NCBI section ----
from Bio import Entrez, SeqIO

Entrez.email = "raazbiochem@gmail.com"
Entrez.api_key = None  # Optional: Add your API key here for higher rate limit

def fetch_ncbi_entry(db, query, query_type="accession", rettype="fasta", retmax=5):
    try:
        if query_type == "accession":
            handle = Entrez.efetch(db=db, id=query, rettype=rettype, retmode="text")
            records = list(SeqIO.parse(handle, rettype))
            handle.close()
            return records
        elif query_type == "gene":
            search = Entrez.esearch(db="gene", term=query, retmax=retmax)
            result = Entrez.read(search)
            if not result["IdList"]:
                return []
            gene_id = result["IdList"][0]
            link = Entrez.elink(dbfrom="gene", db=db, id=gene_id)
            link_result = Entrez.read(link)
            if not link_result[0]["LinkSetDb"]:
                return []
            linked_ids = [l['Id'] for l in link_result[0]['LinkSetDb'][0]['Link']]
            handle = Entrez.efetch(db=db, id=",".join(linked_ids), rettype=rettype, retmode="text")
            records = list(SeqIO.parse(handle, rettype))
            handle.close()
            return records
        elif query_type == "custom":
            search = Entrez.esearch(db=db, term=query, retmax=retmax)
            result = Entrez.read(search)
            ids = result["IdList"]
            if not ids:
                return []
            handle = Entrez.efetch(db=db, id=",".join(ids), rettype=rettype, retmode="text")
            records = list(SeqIO.parse(handle, rettype))
            handle.close()
            return records
        else:
            return []
    except Exception as e:
        st.error(f"NCBI fetch failed: {str(e)}")
        return []

def wrap(seq, width=70):
    return "\n".join([seq[i:i+width] for i in range(0, len(seq), width)])

# Import motif functions
try:
    from motifs import (
        all_motifs, 
        find_hotspots,
        parse_fasta, gc_content, reverse_complement,
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
    "Curved DNA": "#FF9AA2",
    "Z DNA": "#FFB7B2",
    "Slipped DNA": "#FFDAC1",
    "R Loop": "#FFD3B6",
    "Cruciform": "#E2F0CB",
    "Triplex DNA": "#B5EAD7",
    "Sticky DNA": "#DCB8CB",
    "G Triplex": "#C7CEEA",
    "G quadruplex": "#A2D7D8",
    "i Motif": "#B0C4DE",    
    "Hybrid": "#C1A192",
    "Non-B DNA Clusters": "#A2C8CC"
}

def motif_class_with_spaces(motif):
    return motif.replace("_", " ")

# ... (your existing imports)
from Bio import Entrez, SeqIO

# ... (your session state and motif class setup)

PAGES = {
    "Home": "Introduction and overview",
    "Upload & Analyze": "Submit DNA sequence for analysis",
    "Results": "View detected motifs and statistics", 
    "Visualization": "Graphical representation of motifs",
    "Download": "Export results for further analysis",
    "Documentation": "Scientific methods and references",
    "Advanced": "Advanced sequence/multifasta analysis"
}


# --- Sidebar navigation ---
st.sidebar.title("Navigation")
page = st.sidebar.radio("Go to", list(PAGES.keys()))

# --- Home page ---
if page == "Home":
    st.title("Non-B DNA Motif Finder")
    #st.markdown("This tool identifies **12 classes** of non-canonical DNA structures using published algorithms:")
    try:
        nbd_image = Image.open("nbd3.png")
        st.image(nbd_image, use_container_width=True)
    except Exception:
        pass



    st.markdown(
        """
        <div style='margin-top: 34px; font-size: 18px; line-height: 1.7; background: #f1f8fa; border-radius: 8px; padding: 22px 22px 18px 22px; box-shadow: 0px 2px 10px #e0e5ea;'>
        <b>The Non-B DNA Motif Finder</b> provides comprehensive detection of 12 distinct non-canonical DNA structure types, employing scientifically validated algorithms and established thresholds for exploratory genome-wide motif scanning (scientifically reasonable). Users benefit from interactive visualizations resembling a genome browser and versatile export options including CSV, Excel, and images.<br><br>
        <b>How to use:</b> Simply upload or paste your DNA sequence, execute the analysis, explore interactive visual representations, and download the detailed results for further examination.
        </div>
        """, unsafe_allow_html=True
    )

# --- Upload & Analyze ---
elif page == "Upload & Analyze":
    st.header("Sequence Input")
    with st.expander("Input Options", expanded=True):
        input_method = st.radio("Select input method:", 
            ["File Upload", "Example Sequence", "Paste Sequence", "NCBI Fetch"])
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
        elif input_method == "NCBI Fetch":
            st.markdown("Fetch a sequence from NCBI by accession, gene name, or custom query.")
            ncbi_db = st.selectbox("NCBI Database", ["nucleotide", "protein", "gene"])
            ncbi_query_type = st.radio("Query Type", ["Accession", "Gene Name", "Custom Query"])
            ncbi_query = st.text_input("Enter your query (accession, gene name, or Entrez term):")
            ncbi_rettype = st.selectbox("Return format", ["fasta", "gb"])
            ncbi_retmax = st.number_input("Max records (for gene/custom)", min_value=1, max_value=20, value=5)
            if st.button("Fetch from NCBI"):
                if ncbi_query:
                    qtype = {"Accession": "accession", "Gene Name": "gene", "Custom Query": "custom"}[ncbi_query_type]
                    with st.spinner("Contacting NCBI..."):
                        records = fetch_ncbi_entry(
                            db=ncbi_db,
                            query=ncbi_query,
                            query_type=qtype,
                            rettype=ncbi_rettype,
                            retmax=ncbi_retmax
                        )
                    if not records:
                        st.error("No sequence found for this query.")
                    else:
                        seq = str(records[0].seq)
                        st.session_state.seq = seq
                        st.success(f"Fetched sequence: {len(seq):,} bp")
                        st.code(f">{records[0].id}\n{wrap(seq)}", language="fasta")
                else:
                    st.warning("Enter a query before fetching.")

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
                    # Only keep non-overlapping results for all visual sections except Download
                    st.session_state.motif_results_nonoverlap = select_best_nonoverlapping_motifs(
                        st.session_state.motif_results
                    )
                    st.session_state.df = pd.DataFrame(st.session_state.motif_results_nonoverlap)
                    st.session_state.hotspots = find_hotspots(
                        st.session_state.motif_results_nonoverlap,
                        len(st.session_state.seq)
                    )
                    if st.session_state.motif_results_nonoverlap:
                        st.success(f"Found {len(st.session_state.motif_results_nonoverlap)} motifs across {st.session_state.df['Class'].nunique()} classes")
                        st.session_state.analysis_status = "Complete"
                    else:
                        st.warning("No motifs detected")
                        st.session_state.analysis_status = "Complete"
                except Exception as e:
                    st.error(f"Analysis failed: {str(e)}")
                    st.session_state.analysis_status = "Error"

# --- Results page (non-overlapping only) ---
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
        display_df = df.copy()
        display_df['Class'] = display_df['Class'].apply(lambda x: x.replace("_", " "))
        display_df['Subtype'] = display_df['Subtype'].apply(lambda x: x.replace("_", " "))
        st.dataframe(
            display_df[show_cols],
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
                data=display_df,
                y='Class',
                order=display_df['Class'].value_counts().index,
                palette=MOTIF_CLASSES.values()
            )
            ax.set_title("Motif Count by Class")
            st.pyplot(fig)
        with tab2:
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.boxplot(
                data=display_df, 
                x='Length', 
                y='Class',
                palette=MOTIF_CLASSES.values()
            )
            ax.set_title("Motif Length Distribution")
            st.pyplot(fig)
        with tab3:
            fig, ax = plt.subplots(figsize=(10, 6))
            sns.violinplot(
                data=display_df,
                x='Score',
                y='Class',
                palette=MOTIF_CLASSES.values()
            )
            ax.set_title("Score Distribution by Class")
            st.pyplot(fig)
        with tab4:
            fig, ax = plt.subplots(figsize=(12, 3))
            pos = []
            for _, row in display_df.iterrows():
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

# --- Visualization page: Non-overlapping motif tracks, no scores ---
elif page == "Visualization":
    st.header("Motif Map (Non-Overlapping Motif Tracks)")
    if not st.session_state.motif_results_nonoverlap:
        st.session_state.motif_results_nonoverlap = select_best_nonoverlapping_motifs(
            st.session_state.motif_results
        )
    df = pd.DataFrame(st.session_state.motif_results_nonoverlap)
    if df.empty:
        st.info("No results to visualize. Please run analysis first.")
    else:
        display_df = df.copy()
        display_df['Class'] = display_df['Class'].apply(lambda x: x.replace("_", " "))
        seq_len = len(st.session_state.seq)
        st.sidebar.subheader("Motif Map Settings")
        show_classes = st.sidebar.multiselect(
            "Select motif classes to display:",
            sorted(display_df['Class'].unique()),
            default=sorted(display_df['Class'].unique())
        )
        position_range = st.sidebar.slider(
            "Sequence position range:",
            0, seq_len, (0, min(5000, seq_len))
        )
        viz_df = display_df[
            (display_df['Class'].isin(show_classes)) &
            (display_df['Start'] >= position_range[0]) & 
            (display_df['End'] <= position_range[1])
        ].copy()
        if viz_df.empty:
            st.warning("No motifs match the selected filters")
        else:
            classes = sorted(viz_df['Class'].unique())
            y_pos = {motif:i+1 for i, motif in enumerate(classes)}
            color_map = dict(zip(classes, sns.color_palette("tab20", len(classes))))
            fig, ax = plt.subplots(figsize=(15, max(6, len(classes)//2+2)))
            for _, row in viz_df.iterrows():
                ax.hlines(
                    y_pos[row['Class']],
                    row['Start'],
                    row['End'],
                    linewidth=13,
                    color=color_map.get(row['Class'], "#888"),
                    alpha=0.88
                )
            ax.set_yticks(list(y_pos.values()))
            ax.set_yticklabels(list(y_pos.keys()))
            ax.set_xlim(position_range[0], position_range[1])
            ax.set_xlabel("Sequence Position (bp)")
            ax.set_title(f"Non-B DNA Motifs ({position_range[0]:,}-{position_range[1]:,} bp, non-overlapping)")
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
                    f"{s}</div>" for s in classes
                )
                + "</div>", unsafe_allow_html=True
            )

# --- Download page ---
# --- Advanced page ---
elif page == "Advanced":
    st.header("Advanced Analysis: All Non-B DNA Motifs from Multi-FASTA")
    st.markdown("""
    Upload a multi-FASTA file (all sequences must be of the same length).
    The app will automatically detect all 12 non-B DNA motif classes as defined in motifs.py and plot, for each class, a histogram of the motif positions across all sequences.
    """)

    fasta_file = st.file_uploader("Upload multi-FASTA file", type=["fa", "fasta", "txt"])

    if fasta_file:
        try:
            content = fasta_file.read().decode("utf-8")
            seqs, seq_names = [], []
            cur_seq, cur_name = "", ""
            for line in content.splitlines():
                if line.startswith(">"):
                    if cur_seq:
                        seqs.append(cur_seq.upper())
                        seq_names.append(cur_name)
                    cur_name = line.strip().lstrip(">")
                    cur_seq = ""
                else:
                    cur_seq += line.strip()
            if cur_seq:
                seqs.append(cur_seq.upper())
                seq_names.append(cur_name)
            if len(set(len(s) for s in seqs)) != 1:
                st.error("All FASTA sequences must be of the same length!")
            else:
                seq_len = len(seqs[0])
                st.success(f"Loaded {len(seqs)} sequences, each {seq_len} bp long.")
                # Collect motif occurrences for all classes
                motif_pos_dict = {}  # {class: [positions]}
                motif_count_dict = {}  # {class: hit count}
                motif_names_dict = {}  # {class: [subtype list]}
                for idx, seq in enumerate(seqs):
                    hits = all_motifs(seq)
                    for hit in hits:
                        motif_class = hit["Class"]
                        start = hit["Start"] - 1  # 0-based
                        if motif_class not in motif_pos_dict:
                            motif_pos_dict[motif_class] = []
                            motif_names_dict[motif_class] = set()
                            motif_count_dict[motif_class] = 0
                        motif_pos_dict[motif_class].append(start)
                        motif_names_dict[motif_class].add(hit.get("Subtype", ""))
                        motif_count_dict[motif_class] += 1
                if not motif_pos_dict:
                    st.warning("No motifs found in the provided sequences.")
                else:
                    motif_classes = sorted(motif_pos_dict.keys())
                    n = len(motif_classes)
                    fig, axes = plt.subplots(n, 1, figsize=(10, 3*n), sharex=True)
                    if n == 1:
                        axes = [axes]
                    for ax, motif_class in zip(axes, motif_classes):
                        positions = motif_pos_dict[motif_class]
                        ax.hist(positions, bins=np.arange(seq_len+1)-0.5, color='skyblue', edgecolor='black')
                        ax.set_title(
                            f"{motif_class} (Subtypes: {', '.join(sorted(motif_names_dict[motif_class]))})"
                        )
                        ax.set_ylabel("Count")
                        ax.set_xlim(-0.5, seq_len-0.5)
                        ax.grid(axis='y')
                    axes[-1].set_xlabel("Position in Sequence (0-based)")
                    plt.tight_layout()
                    st.pyplot(fig)
                    st.success("All motif histograms generated!")
                    # Show summary table
                    df_summary = pd.DataFrame({
                        "Motif Class": motif_classes,
                        "Total Hits": [motif_count_dict[m] for m in motif_classes],
                        "Unique Subtypes": [", ".join(sorted(motif_names_dict[m])) for m in motif_classes]
                    })
                    st.dataframe(df_summary)
        except Exception as e:
            st.error(f"Error processing file: {str(e)}")
    else:
        st.info("Please upload a multi-FASTA file.")
# ... (footer and remaining code)
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

# --- Footer: developer info, always at bottom ---
st.markdown("""
---
<div style='font-size: 15px; color: #1e293b; margin-top: 40px; text-align: left;'>
<b>Developed by</b><br>
Dr. Venkata Rajesh Yella<br>
<a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a><br>
<a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
</div>
""", unsafe_allow_html=True)
