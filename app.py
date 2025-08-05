import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import io
from collections import Counter
from PIL import Image
from Bio import Entrez, SeqIO

# Import motif functions (from your motifs.py)
from motifs import (
    all_motifs, 
    find_hotspots,
    parse_fasta, gc_content, reverse_complement,
    select_best_nonoverlapping_motifs, wrap
)

# --- Custom CSS for Big Tabs and Full Width ---
st.markdown("""
    <style>
    .stTabs [data-baseweb="tab"] {
        font-size: 1.65rem !important;
        font-weight: 700 !important;
        width: 18vw !important;
        min-width: 180px !important;
        padding-top: 18px !important;
        padding-bottom: 18px !important;
        justify-content: center !important;
    }
    .stTabs [data-baseweb="tab-list"] {
        justify-content: stretch !important;
        width: 100vw !important;
        margin-bottom: 8px !important;
    }
    /* Hide scroll bar for tabs */
    .stTabs [data-baseweb="tab-list"]::-webkit-scrollbar {display:none;}
    </style>
""", unsafe_allow_html=True)

# --- Page Config ---
st.set_page_config(
    page_title="Non-B DNA Motif Finder",
    layout="wide",
    page_icon="",
    menu_items={'About': "Non-B DNA Motif Finder | Developed by Dr. Venkata Rajesh Yella"}
)

# --- Constants ---
MOTIF_ORDER = [
    "Sticky DNA","Curved DNA","Z-DNA","eGZ (Extruded-G)","Slipped DNA","R-Loop",
    "Cruciform","Triplex DNA","G-Triplex","G4","Relaxed G4","Bulged G4","Bipartite G4",
    "Multimeric G4","i-Motif","AC-Motif","Hybrid","Non-B DNA Clusters"
]

MOTIF_COLORS = {
    "Curved DNA": "#FF9AA2","Z-DNA": "#FFB7B2","eGZ (Extruded-G)": "#6A4C93",
    "Slipped DNA": "#FFDAC1","R-Loop": "#FFD3B6","Cruciform": "#E2F0CB",
    "Triplex DNA": "#B5EAD7","Sticky DNA": "#DCB8CB","G-Triplex": "#C7CEEA",
    "G4": "#A2D7D8","Relaxed G4": "#A2D7B8","Bulged G4": "#A2A7D8",
    "Bipartite G4": "#A2D788","Multimeric G4": "#A2A7B8","i-Motif": "#B0C4DE",
    "Hybrid": "#C1A192","Non-B DNA Clusters": "#A2C8CC","AC-Motif": "#F5B041"
}

PAGES = {
    "Home": "Overview",
    "Upload & Analyze": "Sequence Upload and Motif Analysis",
    "Results": "Analysis Results and Visualization",
    "Download": "Export Data",
    "Documentation": "Scientific Documentation & References"
}

Entrez.email = "raazbiochem@gmail.com"
Entrez.api_key = None

EXAMPLE_FASTA = """>Example Sequence
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA
CTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCTCT
"""
EXAMPLE_MULTI_FASTA = """>Seq1
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
>Seq2
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
>Seq3
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC
GAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAAGAAA
"""

for k, v in {
    'seqs': [],
    'names': [],
    'results': [],
    'summary_df': pd.DataFrame(),
    'hotspots': [],
    'analysis_status': "Ready",
    'selected_motifs': MOTIF_ORDER,
}.items():
    if k not in st.session_state:
        st.session_state[k] = v

def get_basic_stats(seq, motifs=None):
    seq = seq.upper()
    length = len(seq)
    gc = gc_content(seq)
    at = (seq.count('A') + seq.count('T')) / length * 100 if length else 0
    stats = {
        "Length (bp)": length, "GC %": round(gc, 2), "AT %": round(at, 2),
        "A Count": seq.count('A'), "T Count": seq.count('T'), "G Count": seq.count('G'), "C Count": seq.count('C'),
    }
    if motifs is not None:
        covered = set()
        for m in motifs:
            covered.update(range(m['Start'], m['End']))
        coverage_pct = (len(covered) / length * 100) if length else 0
        stats["Motif Coverage (%)"] = round(coverage_pct, 2)
    return stats

# --- Horizontal Tabs Navigation ---
tabs = st.tabs(list(PAGES.keys()))
tab_pages = dict(zip(PAGES.keys(), tabs))

# --- HOME PAGE ---
with tab_pages["Home"]:
    col1, col2 = st.columns([1, 2])
    with col1:
        st.image("nbdcircle.JPG", width=200)
    with col2:
        st.markdown("""
        <div style='font-family:Montserrat, Arial; font-size:21px; color:#222; line-height:1.7;'>
        <b>Non-canonical DNA structures</b> play key roles in genome stability, regulation, and evolution. 
        This application detects and analyzes <b>18 distinct Non-B DNA motifs</b> in any DNA sequence or multi-FASTA file.<br>
        <br>
        Motifs are classified as follows:<br>
        <span style='color:#2c3e50;'><b>G-quadruplex-related</b> (G4, Relaxed G4, Bulged G4, Bipartite G4, Multimeric G4, G-Triplex, i-Motif, Hybrid), 
        <b>helix/curvature</b> (Z-DNA, eGZ (Extruded-G), Curved DNA, AC-Motif), <b>repeat/junction</b> (Slipped DNA, Cruciform, Sticky DNA, Triplex DNA), 
        <b>hybrid/cluster</b> (R-Loop, Non-B DNA Clusters).</span>
        <br><br>
        <b>Upload single or multi-FASTA files</b> and view interactive motif visualizations with downloadable results for further analysis.
        </div>
        """, unsafe_allow_html=True)

# --- UPLOAD & ANALYZE PAGE ---
with tab_pages["Upload & Analyze"]:
    st.markdown('<h2 style="font-family:Montserrat, Arial; font-weight:700; color:#002147; letter-spacing:1px; margin-bottom:18px;">Sequence Upload and Motif Analysis</h2>', unsafe_allow_html=True)
    st.markdown('<span style="font-family:Montserrat,Arial; font-size:16px;">Supports <b>multi-FASTA</b> and single FASTA. Paste, upload, select example, or fetch from NCBI.</span>', unsafe_allow_html=True)
    st.caption("Supported formats: .fa, .fasta, .txt | Limit: 200MB/file.")

    selected_motifs = st.multiselect(
        "Select Motif Classes for Analysis",
        MOTIF_ORDER,
        default=MOTIF_ORDER,
        help="Choose motif classes to analyze. Selecting 'Hybrid' or 'Non-B DNA Clusters' will run all motif modules."
    )
    st.session_state.selected_motifs = selected_motifs if selected_motifs else MOTIF_ORDER

    input_method = st.radio("Input Method:",
        ["Upload FASTA / Multi-FASTA File", "Paste Sequence(s)", "Example Sequence", "NCBI Fetch"],
        horizontal=True
    )
    
    seqs, names = [], []
    if input_method == "Upload FASTA / Multi-FASTA File":
        fasta_file = st.file_uploader("Drag and drop FASTA/multi-FASTA file here", type=["fa", "fasta", "txt"])
        if fasta_file:
            content = fasta_file.read().decode("utf-8")
            seqs, names = [], []
            cur_seq, cur_name = "", ""
            for line in content.splitlines():
                if line.startswith(">"):
                    if cur_seq:
                        seqs.append(parse_fasta(cur_seq))
                        names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                    cur_name = line.strip().lstrip(">")
                    cur_seq = ""
                else:
                    cur_seq += line.strip()
            if cur_seq:
                seqs.append(parse_fasta(cur_seq))
                names.append(cur_name if cur_name else f"Seq{len(seqs)}")
            if seqs:
                st.success(f"Loaded {len(seqs)} sequences.")
                for i, seq in enumerate(seqs[:3]):
                    stats = get_basic_stats(seq)
                    st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                    st.markdown(f"GC %: {stats['GC %']} | AT %: {stats['AT %']} | A: {stats['A Count']} | T: {stats['T Count']} | G: {stats['G Count']} | C: {stats['C Count']}")
                if len(seqs) > 3:
                    st.caption(f"...and {len(seqs)-3} more.")
            else:
                st.warning("No sequences found.")
    elif input_method == "Paste Sequence(s)":
        seq_input = st.text_area("Paste single or multi-FASTA here:", height=150)
        if seq_input:
            seqs, names = [], []
            cur_seq, cur_name = "", ""
            for line in seq_input.splitlines():
                if line.startswith(">"):
                    if cur_seq:
                        seqs.append(parse_fasta(cur_seq))
                        names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                    cur_name = line.strip().lstrip(">")
                    cur_seq = ""
                else:
                    cur_seq += line.strip()
            if cur_seq:
                seqs.append(parse_fasta(cur_seq))
                names.append(cur_name if cur_name else f"Seq{len(seqs)}")
            if seqs:
                st.success(f"Pasted {len(seqs)} sequences.")
                for i, seq in enumerate(seqs[:3]):
                    st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                    stats = get_basic_stats(seq)
                    st.markdown(f"GC %: {stats['GC %']} | AT %: {stats['AT %']} | A: {stats['A Count']} | T: {stats['T Count']} | G: {stats['G Count']} | C: {stats['C Count']}")
                if len(seqs) > 3:
                    st.caption(f"...and {len(seqs)-3} more.")
            else:
                st.warning("No sequences found.")
    elif input_method == "Example Sequence":
        ex_type = st.radio("Example Type:", ["Single Example", "Multi-FASTA Example"], horizontal=True)
        if ex_type == "Single Example":
            if st.button("Load Single Example"):
                seqs = [parse_fasta(EXAMPLE_FASTA)]
                names = ["Example Sequence"]
                st.success("Single example sequence loaded.")
                stats = get_basic_stats(seqs[0])
                st.code(EXAMPLE_FASTA, language="fasta")
                st.markdown(f"GC %: {stats['GC %']} | AT %: {stats['AT %']} | A: {stats['A Count']} | T: {stats['T Count']} | G: {stats['G Count']} | C: {stats['C Count']}")
        else:
            if st.button("Load Multi-FASTA Example"):
                seqs, names = [], []
                cur_seq, cur_name = "", ""
                for line in EXAMPLE_MULTI_FASTA.splitlines():
                    if line.startswith(">"):
                        if cur_seq:
                            seqs.append(parse_fasta(cur_seq))
                            names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                        cur_name = line.strip().lstrip(">")
                        cur_seq = ""
                    else:
                        cur_seq += line.strip()
                if cur_seq:
                    seqs.append(parse_fasta(cur_seq))
                    names.append(cur_name if cur_name else f"Seq{len(seqs)}")
                st.success(f"Multi-FASTA example loaded with {len(seqs)} sequences.")
                for i, seq in enumerate(seqs[:3]):
                    stats = get_basic_stats(seq)
                    st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                    st.markdown(f"GC %: {stats['GC %']} | AT %: {stats['AT %']} | A: {stats['A Count']} | T: {stats['T Count']} | G: {stats['G Count']} | C: {stats['C Count']}")
                st.code(EXAMPLE_MULTI_FASTA, language="fasta")
    elif input_method == "NCBI Fetch":
        db = st.selectbox("NCBI Database", ["nucleotide", "protein", "gene"])
        query_type = st.radio("Query Type", ["Accession", "Gene Name", "Custom Query"], horizontal=True)
        motif_examples = {
            "G-quadruplex": "NR_003287.2 (human telomerase RNA)",
            "Z-DNA": "NM_001126112.2 (human ADAR1 gene)",
            "R-loop": "NR_024540.1 (human SNRPN gene)",
            "eGZ-motif": "CGG repeat region",
            "AC-motif": "A-rich/C-rich consensus region"
        }
        with st.expander("Motif Example Queries"):
            for motif, example in motif_examples.items():
                st.write(f"**{motif}**: `{example}`")
        query = st.text_input("Enter query (accession, gene, etc.):")
        rettype = st.selectbox("Return Format", ["fasta", "gb"])
        retmax = st.number_input("Max Records", min_value=1, max_value=20, value=3)
        if st.button("Fetch from NCBI"):
            if query:
                qtype = {"Accession": "accession", "Gene Name": "gene", "Custom Query": "custom"}[query_type]
                with st.spinner("Contacting NCBI..."):
                    handle = Entrez.efetch(db=db, id=query, rettype=rettype, retmode="text")
                    records = list(SeqIO.parse(handle, rettype))
                    handle.close()
                    seqs = [str(rec.seq).upper().replace("U", "T") for rec in records]
                    names = [rec.id for rec in records]
                if seqs:
                    st.success(f"Fetched {len(seqs)} sequences.")
                    for i, seq in enumerate(seqs[:3]):
                        stats = get_basic_stats(seq)
                        st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC %: {stats['GC %']} | AT %: {stats['AT %']} | A: {stats['A Count']} | T: {stats['T Count']} | G: {stats['G Count']} | C: {stats['C Count']}")
            else:
                st.warning("Enter a query before fetching.")

    if seqs:
        st.session_state.seqs = seqs
        st.session_state.names = names
        st.session_state.results = []

    if st.session_state.seqs:
        st.subheader("Sequence Preview")
        for i, seq in enumerate(st.session_state.seqs[:2]):
            stats = get_basic_stats(seq)
            st.markdown(f"<b>{st.session_state.names[i]}</b> ({len(seq):,} bp) | GC %: {stats['GC %']} | AT %: {stats['AT %']} | A: {stats['A Count']} | T: {stats['T Count']} | G: {stats['G Count']} | C: {stats['C Count']}", unsafe_allow_html=True)
            st.code(wrap(seq[:400]), language="fasta")
        if len(st.session_state.seqs) > 2:
            st.caption(f"...and {len(st.session_state.seqs)-2} more.")

        # --- MOTIF ANALYSIS LOGIC ---
        run_all = any(m in st.session_state.selected_motifs for m in ["Hybrid", "Non-B DNA Clusters"])
        if st.button("Run Motif Analysis", type="primary"):
            st.session_state.analysis_status = "Running"
            motif_results = []
            for seq in st.session_state.seqs:
                if run_all:
                    motifs = all_motifs(seq)
                else:
                    motifs = [m for m in all_motifs(seq) if m['Class'] in st.session_state.selected_motifs]
                nonoverlapping = select_best_nonoverlapping_motifs(motifs)
                motif_results.append(nonoverlapping)
            st.session_state.results = motif_results

            # --- SUMMARY TABLE ---
            summary = []
            for i, motifs in enumerate(motif_results):
                stats = get_basic_stats(st.session_state.seqs[i], motifs)
                motif_types = Counter([m['Class'] if m['Class'] != "Z-DNA" or m.get("Subclass") != "eGZ (Extruded-G)" else "eGZ (Extruded-G)" for m in motifs])
                summary.append({
                    "Sequence Name": st.session_state.names[i],
                    "Length (bp)": stats['Length (bp)'],
                    "GC %": stats['GC %'],
                    "AT %": stats['AT %'],
                    "A Count": stats['A Count'],
                    "T Count": stats['T Count'],
                    "G Count": stats['G Count'],
                    "C Count": stats['C Count'],
                    "Motif Count": len(motifs),
                    "Motif Coverage (%)": stats["Motif Coverage (%)"],
                    "Motif Classes": ", ".join(f"{k} ({v})" for k, v in motif_types.items())
                })
            st.session_state.summary_df = pd.DataFrame(summary)
            st.success("Analysis complete! See 'Analysis Results and Visualization' tab for details.")
            st.session_state.analysis_status = "Complete"

# --- RESULTS PAGE ---
with tab_pages["Results"]:
    st.markdown('<h2 style="font-family:Montserrat, Arial; font-weight:700; color:#002147; letter-spacing:1px; margin-bottom:18px;">Analysis Results and Visualization</h2>', unsafe_allow_html=True)
    if not st.session_state.results:
        st.info("No analysis results. Please run motif analysis first.")
    else:
        st.dataframe(st.session_state.summary_df, use_container_width=True)
        seq_idx = 0
        if len(st.session_state.seqs) > 1:
            seq_idx = st.selectbox("Choose Sequence for Details:", range(len(st.session_state.seqs)), format_func=lambda i: st.session_state.names[i])
        motifs = st.session_state.results[seq_idx]
        if not motifs:
            st.warning("No motifs detected for this sequence.")
        else:
            df = pd.DataFrame(motifs)
            st.markdown(f"<h3 style='font-family:Montserrat, Arial; font-weight:600; color:#0A3D62;'>Motif Table for <b>{st.session_state.names[seq_idx]}</b></h3>", unsafe_allow_html=True)
            display_columns = [col for col in df.columns if col not in ['Sequence', 'SequenceName']]
            st.dataframe(df[display_columns], use_container_width=True, height=360)
            st.markdown('<span style="font-family:Montserrat,Arial;font-size:17px;"><b>Motif Type Distribution</b></span>', unsafe_allow_html=True)
            fig, ax = plt.subplots(figsize=(8,6))
            class_counts = df['Class'].value_counts().reindex(
                st.session_state.selected_motifs if not any(m in st.session_state.selected_motifs for m in ['Hybrid', 'Non-B DNA Clusters']) else MOTIF_ORDER, fill_value=0)
            if "Subclass" in df.columns:
                egz_count = (df["Subclass"] == "eGZ (Extruded-G)").sum()
                if "eGZ (Extruded-G)" in class_counts.index:
                    class_counts["eGZ (Extruded-G)"] = egz_count
            ax.barh(class_counts.index, class_counts.values, color=[MOTIF_COLORS.get(c, "#888") for c in class_counts.index])
            ax.set_xlabel("Motif Count")
            st.pyplot(fig)
            st.markdown('<span style="font-family:Montserrat,Arial;font-size:17px;"><b>Motif Map</b></span>', unsafe_allow_html=True)
            fig, ax = plt.subplots(figsize=(12,3))
            y = 1
            for _, row in df.iterrows():
                motif_class = row['Class']
                if motif_class == "Z-DNA" and row.get("Subclass", "") == "eGZ (Extruded-G)":
                    motif_class = "eGZ (Extruded-G)"
                color = MOTIF_COLORS.get(motif_class, "#888")
                ax.plot([row['Start'], row['End']], [y, y], lw=8, color=color, alpha=0.8)
                y += 0.12
            ax.set_yticks([])
            ax.set_xlabel("Sequence Position (bp)")
            ax.set_xlim(0, len(st.session_state.seqs[seq_idx]))
            ax.set_title(f"Motif Tracks: {st.session_state.names[seq_idx]}", fontweight='bold', fontsize=16)
            st.pyplot(fig)

# --- DOWNLOAD PAGE ---
with tab_pages["Download"]:
    st.header("Export Data")
    if not st.session_state.results:
        st.info("No results available to download.")
    else:
        df_all = []
        for i, motifs in enumerate(st.session_state.results):
            for m in motifs:
                m['Sequence Name'] = st.session_state.names[i]
                if m['Class'] == "Z-DNA" and m.get("Subclass", "") == "eGZ (Extruded-G)":
                    m['Class'] = "eGZ (Extruded-G)"
                df_all.append(m)
        df_all = pd.DataFrame(df_all)
        st.dataframe(df_all, use_container_width=True, height=350)
        csv_data = df_all.to_csv(index=False).encode("utf-8")
        st.download_button("Download CSV", data=csv_data, file_name="motif_results.csv", mime="text/csv")
        excel_data = io.BytesIO()
        with pd.ExcelWriter(excel_data, engine='xlsxwriter') as writer:
            df_all.to_excel(writer, index=False, sheet_name="Motifs")
        excel_data.seek(0)
        st.download_button("Download Excel", data=excel_data, file_name="motif_results.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

# --- DOCUMENTATION PAGE ---
with tab_pages["Documentation"]:
    st.header("Scientific Documentation & References")
    st.markdown("""
    <div style='background:#f4faff; border-radius:10px; padding:18px 18px 8px 18px; font-size:17px; font-family:Montserrat,Arial;'>
    <b>Motif Classes Detected:</b><br><br>
    <ul style="font-size:16px;">
        <li>
            <b>Curved DNA</b>: Identifies phased poly(A) or poly(T) tracts using regular expressions and spacing rules, reflecting regions of intrinsic curvature. Scoring is based on tract length and grouping of spaced tracts, indicating potential for DNA bending.
        </li>
        <li>
            <b>Z-DNA</b>: Detects alternating purine-pyrimidine patterns, especially GC-rich segments, using a windowed scoring algorithm (Z-Seeker). The regular expression strategy finds dinucleotide repeats; scoring reflects the weighted sum of favorable and unfavorable base steps.
        </li>
        <li>
            <b>eGZ-motif (Extruded-G Z-DNA)</b>: Searches for long (CGG)<sub>n</sub> runs via regex, representing extruded guanine Z-DNA variants. Scoring is normalized to repeat count, indicating stability of the motif.
        </li>
        <li>
            <b>Slipped DNA</b>: Recognizes direct repeats and short tandem repeats, using repeat-unit matching and regex for consecutive identical segments. Scoring is proportional to repeat length and unit copies, reflecting mutational susceptibility.
        </li>
        <li>
            <b>R-Loop</b>: Locates G-rich regions capable of forming stable RNA-DNA hybrids, using RLFS models defined by regular expressions. Thermodynamic scoring combines GC content and G-run frequency to estimate hybrid stability.
        </li>
        <li>
            <b>Cruciform</b>: Finds palindromic inverted repeats with short spacers using regex and reverse complement matching. Scoring factors in arm length and A/T content, indicating stem-loop formation propensity.
        </li>
        <li>
            <b>Triplex DNA / Mirror Repeat</b>: Detects purine- or pyrimidine-rich mirror repeats and triplex-forming motifs. Regular expressions identify repeated units separated by short spacers; scoring considers base composition and repeat purity.
        </li>
        <li>
            <b>Sticky DNA</b>: Searches for extended GAA or TTC repeats with regex. Scoring depends on repeat count, reflecting the likelihood of forming sticky junctions and associated instability.
        </li>
        <li>
            <b>G-Triplex</b>: Identifies three consecutive guanine runs separated by short loops. Regular expression and loop-length analysis determine motif stability; scoring incorporates G-run sum and loop penalty.
        </li>
        <li>
            <b>G4 (G-Quadruplex) and Variants</b>: Detects canonical, relaxed, bulged, bipartite, multimeric, and imperfect G4 motifs using specialized regex patterns for G-runs and loop sizes. The G4Hunter scoring system evaluates mean G/C content and structure favorability.
        </li>
        <li>
            <b>i-Motif</b>: Finds C-rich sequences capable of folding into i-motif structures under acidic conditions. Regex detects cytosine runs and loop lengths, while scoring incorporates run count, C-content, and loop properties.
        </li>
        <li>
            <b>AC-Motif</b>: Locates alternating A-rich and C-rich consensus regions using regex. Motif matches are scored by pattern presence, suggesting non-canonical secondary structure potential.
        </li>
        <li>
            <b>Hybrid Motif</b>: Represents DNA regions where two or more distinct motif classes overlap. Identified by interval intersection algorithms after all motif modules run, with scoring based on overlap diversity and region size.
        </li>
        <li>
            <b>Non-B DNA Clusters</b>: Represents hotspot regions where multiple non-B DNA motifs are found within a defined window using a sliding algorithm. The region is scored by motif count and type diversity, indicating local enrichment and possible functional significance.
        </li>
    </ul>
    <b>References:</b>
    <ul>
        <li>Bedrat et al., 2016 Nucleic Acids Research</li>
        <li>Ho et al., 2010 Nature Chemical Biology</li>
        <li>Kim et al., 2018 Nucleic Acids Research</li>
        <li>Zeraati et al., 2018 Nature Chemistry</li>
        <li>Bacolla et al., 2006 Nucleic Acids Research</li>
        <li>Mirkin & Frank-Kamenetskii, 1994 Annual Review of Biophysics</li>
        <li>New et al., 2020 Journal of DNA Structure</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)

st.markdown("""
---
<div style='font-size: 15px; color: #1e293b; margin-top: 40px; text-align: left; font-family:Montserrat,Arial;'>
<b>Developed by</b><br>
Dr. Venkata Rajesh Yella<br>
<a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a> |
<a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
</div>
""", unsafe_allow_html=True)
