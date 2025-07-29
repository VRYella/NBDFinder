import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import io
from collections import Counter
from datetime import datetime
from PIL import Image
from Bio import Entrez, SeqIO

# Import motif functions (from your motifs.py)
from motifs import (
    all_motifs, 
    find_hotspots,
    parse_fasta, gc_content, reverse_complement,
    select_best_nonoverlapping_motifs, wrap
)

Entrez.email = "raazbiochem@gmail.com"
Entrez.api_key = None

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

# --- Session state initialization ---
for k, v in {
    'seqs': [],   # List of sequences (single or multi)
    'names': [],  # Sequence names
    'results': [],  # List of motif results per sequence
    'summary_df': pd.DataFrame(),
    'hotspots': [],
    'analysis_status': "Ready",
}.items():
    if k not in st.session_state:
        st.session_state[k] = v

MOTIF_COLORS = {
    "Curved_DNA": "#FF9AA2", "Z-DNA": "#FFB7B2", "Slipped_DNA": "#FFDAC1", "R-Loop": "#FFD3B6",
    "Cruciform": "#E2F0CB", "Triplex_DNA": "#B5EAD7", "Sticky_DNA": "#DCB8CB", "G-Triplex": "#C7CEEA",
    "G4": "#A2D7D8", "i-Motif": "#B0C4DE", "Hybrid": "#C1A192", "Non-B DNA Clusters": "#A2C8CC"
}

PAGES = {
    "Home": "Overview",
    "Upload & Analyze": "Upload DNA (multi-FASTA supported!)",
    "Results": "Motif summary and visualization",
    "Download": "Export results",
    "Documentation": "Scientific methods & references"
}

def basic_stats(seq):
    seq = seq.upper()
    length = len(seq)
    gc = gc_content(seq)
    at = (seq.count('A') + seq.count('T')) / length * 100 if length else 0
    stats = {
        "Length": length,
        "GC%": round(gc, 2),
        "AT%": round(at, 2),
        "A": seq.count('A'),
        "T": seq.count('T'),
        "G": seq.count('G'),
        "C": seq.count('C'),
    }
    return stats

st.sidebar.title("Navigation")
page = st.sidebar.radio("", list(PAGES.keys()))

# --- Home page ---
if page == "Home":
    st.markdown("<h1 style='color:#0A3D62;'>Non-B DNA Motif Finder</h1>", unsafe_allow_html=True)
    try:
        st.image("nbd3.png", use_container_width=True)
    except Exception:
        pass
    st.markdown(
        """
        <div style='background: #eaf6fb; border-radius: 14px; padding: 22px 22px; box-shadow: 0px 4px 16px #e0e5ea; font-size: 18px;'>
        <b>Detect 12+ Non-Canonical DNA Motifs</b> in any DNA sequence or multi-FASTA file.<br>
        <ul>
        <li>Upload FASTA or multi-FASTA files (no need to specify single/multi format)</li>
        <li>Motifs detected include Z-DNA, G4, i-Motif, R-loop, Cruciform, Triplex, Hybrids, and more</li>
        <li>Interactive motif visualizations and downloadable results</li>
        </ul>
        <b>Developed by Dr. Venkata Rajesh Yella</b>
        </div>
        """, unsafe_allow_html=True
    )

# --- Upload & Analyze page ---
elif page == "Upload & Analyze":
    st.markdown("<h2 style='color:#0A3D62;'>Sequence Input</h2>", unsafe_allow_html=True)
    st.markdown("Supports <b>multi-FASTA</b> (multiple sequences) and single FASTA. Paste, upload, select example, or fetch from NCBI.", unsafe_allow_html=True)
    st.caption("Supported formats: .fa, .fasta, .txt | Limit: 200MB/file.")

    input_method = st.radio("Input method:", [
        "Upload FASTA / multi-FASTA file",
        "Paste Sequence(s)",
        "Example Sequence",
        "NCBI Fetch"
    ])
    
    seqs, names = [], []
    if input_method == "Upload FASTA / multi-FASTA file":
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
                    st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                    stats = basic_stats(seq)
                    st.markdown(f"GC%: {stats['GC%']} | AT%: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
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
                    stats = basic_stats(seq)
                    st.markdown(f"GC%: {stats['GC%']} | AT%: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                if len(seqs) > 3:
                    st.caption(f"...and {len(seqs)-3} more.")
            else:
                st.warning("No sequences found.")
    elif input_method == "Example Sequence":
        ex_type = st.radio("Example type:", ["Single Example", "Multi-FASTA Example"])
        if ex_type == "Single Example":
            if st.button("Load Single Example"):
                seqs = [parse_fasta(EXAMPLE_FASTA)]
                names = ["Example_Sequence"]
                st.success("Single example sequence loaded.")
                stats = basic_stats(seqs[0])
                st.code(EXAMPLE_FASTA, language="fasta")
                st.markdown(f"GC%: {stats['GC%']} | AT%: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
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
                    stats = basic_stats(seq)
                    st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                    st.markdown(f"GC%: {stats['GC%']} | AT%: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
                st.code(EXAMPLE_MULTI_FASTA, language="fasta")
    elif input_method == "NCBI Fetch":
        db = st.selectbox("NCBI Database", ["nucleotide", "protein", "gene"])
        query_type = st.radio("Query Type", ["Accession", "Gene Name", "Custom Query"])
        motif_examples = {
            "G-quadruplex": "NR_003287.2 (human telomerase RNA)",
            "Z-DNA": "NM_001126112.2 (human ADAR1 gene)",
            "R-loop": "NR_024540.1 (human SNRPN gene)"
        }
        with st.expander("Motif example queries"):
            for motif, example in motif_examples.items():
                st.write(f"**{motif}**: `{example}`")
        query = st.text_input("Enter query (accession, gene, etc.):")
        rettype = st.selectbox("Return format", ["fasta", "gb"])
        retmax = st.number_input("Max records", min_value=1, max_value=20, value=3)
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
                        stats = basic_stats(seq)
                        st.markdown(f"<b>{names[i]}</b>: <span style='color:#576574'>{len(seq):,} bp</span>", unsafe_allow_html=True)
                        st.markdown(f"GC%: {stats['GC%']} | AT%: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}")
            else:
                st.warning("Enter a query before fetching.")

    # Save to session state if sequences loaded
    if seqs:
        st.session_state.seqs = seqs
        st.session_state.names = names
        st.session_state.results = []  # reset previous results

    # Preview
    if st.session_state.seqs:
        st.subheader("Sequence Preview")
        for i, seq in enumerate(st.session_state.seqs[:2]):
            stats = basic_stats(seq)
            st.markdown(f"<b>{st.session_state.names[i]}</b> ({len(seq):,} bp) | GC%: {stats['GC%']} | AT%: {stats['AT%']} | A: {stats['A']} | T: {stats['T']} | G: {stats['G']} | C: {stats['C']}", unsafe_allow_html=True)
            st.code(wrap(seq[:400]), language="fasta")
        if len(st.session_state.seqs) > 2:
            st.caption(f"...and {len(st.session_state.seqs)-2} more.")

        if st.button("Run Motif Analysis", type="primary"):
            st.session_state.analysis_status = "Running"
            motif_results = []
            for seq in st.session_state.seqs:
                motifs = all_motifs(seq)
                nonoverlapping = select_best_nonoverlapping_motifs(motifs)
                motif_results.append(nonoverlapping)
            st.session_state.results = motif_results

            # Summary DataFrame
            summary = []
            for i, motifs in enumerate(motif_results):
                motif_types = Counter([m['Class'] for m in motifs])
                stats = basic_stats(st.session_state.seqs[i])
                summary.append({
                    "Sequence": st.session_state.names[i],
                    "Length": stats['Length'],
                    "GC%": stats['GC%'],
                    "AT%": stats['AT%'],
                    "A": stats['A'],
                    "T": stats['T'],
                    "G": stats['G'],
                    "C": stats['C'],
                    "Motif Count": len(motifs),
                    "Motif Types": ", ".join(f"{k}({v})" for k, v in motif_types.items())
                })
            st.session_state.summary_df = pd.DataFrame(summary)
            st.success("Analysis complete! See 'Results' tab for details.")
            st.session_state.analysis_status = "Complete"

# --- Results page ---
elif page == "Results":
    st.markdown("<h2 style='color:#0A3D62;'>Motif Summary & Visualization</h2>", unsafe_allow_html=True)
    if not st.session_state.results:
        st.info("No analysis results. Please run motif analysis first.")
    else:
        st.dataframe(st.session_state.summary_df, use_container_width=True)
        seq_idx = 0
        if len(st.session_state.seqs) > 1:
            seq_idx = st.selectbox("Choose sequence for details:", range(len(st.session_state.seqs)), format_func=lambda i: st.session_state.names[i])
        motifs = st.session_state.results[seq_idx]
        if not motifs:
            st.warning("No motifs detected for this sequence.")
        else:
            df = pd.DataFrame(motifs)
            st.markdown(f"### Motif Table for <b>{st.session_state.names[seq_idx]}</b>", unsafe_allow_html=True)
            st.dataframe(df, use_container_width=True, height=360)
            # Distribution plots
            motif_class_order = list(MOTIF_COLORS.keys())
            st.markdown("#### Motif Type Distribution")
            fig, ax = plt.subplots(figsize=(8,6))
            class_counts = df['Class'].value_counts().reindex(motif_class_order, fill_value=0)
            ax.barh(class_counts.index, class_counts.values, color=[MOTIF_COLORS.get(c, "#888") for c in class_counts.index])
            ax.set_xlabel("Motif Count")
            st.pyplot(fig)
            st.markdown("<br>", unsafe_allow_html=True)
            # Motif positions
            st.markdown("#### Motif Map")
            fig, ax = plt.subplots(figsize=(12,3))
            y = 1
            for _, row in df.iterrows():
                color = MOTIF_COLORS.get(row['Class'], "#888")
                ax.plot([row['Start'], row['End']], [y, y], lw=8, color=color, alpha=0.8)
                y += 0.12
            ax.set_yticks([])
            ax.set_xlabel("Sequence Position (bp)")
            ax.set_xlim(0, len(st.session_state.seqs[seq_idx]))
            ax.set_title(f"Motif tracks: {st.session_state.names[seq_idx]}")
            st.pyplot(fig)

# --- Download page ---
elif page == "Download":
    st.header("Download Results")
    if not st.session_state.results:
        st.info("No results available to download.")
    else:
        df_all = []
        for i, motifs in enumerate(st.session_state.results):
            for m in motifs:
                m['SequenceName'] = st.session_state.names[i]
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

# --- Documentation page ---
elif page == "Documentation":
    st.header("Scientific Documentation")
    st.markdown("""
    <div style='background:#f4faff; border-radius:10px; padding:18px 18px 8px 18px; font-size:17px;'>
    <b>Motif Classes Detected:</b><br>
    <ul>
        <li><b>Curved DNA:</b> Detects phasing of A/T tracts and global/local curvature. <i>Brukner et al., 1995</i>.</li>
        <li><b>Z-DNA:</b> Alternating purine/pyrimidine patterns via Z-Seeker scoring. <i>Ho et al., 2010</i>.</li>
        <li><b>Slipped DNA:</b> Direct repeats and short tandem repeats. <i>Bacolla et al., 2006</i>.</li>
        <li><b>R-Loop:</b> RLFS models for G-rich skew and thermodynamic stability. <i>Sanz et al., 2016</i>.</li>
        <li><b>Cruciform:</b> Inverted repeats with AT-rich arms. <i>Lilley, 1985</i>.</li>
        <li><b>Triplex DNA / Mirror Repeat:</b> Purine/pyrimidine-rich mirror repeats and triplex-forming motifs. <i>Mirkin, 1994</i>.</li>
        <li><b>Sticky DNA:</b> Extended GAA/TTC repeats. <i>Potaman et al., 2003</i>.</li>
        <li><b>G-Triplex & G4:</b> Canonical, relaxed, bulged, bipartite, multimeric, and imperfect G-quadruplexes. <i>Bedrat et al., 2016</i>.</li>
        <li><b>i-Motif:</b> C-rich, looped sequences. <i>Zeraati et al., 2018</i>.</li>
        <li><b>Hybrids:</b> Overlapping regions of two or more motif classes.</li>
        <li><b>Non-B DNA Clusters:</b> Hotspot regions with high motif density and diversity.</li>
    </ul>
    <b>Scoring:</b><br>
    - <b>G4Hunter</b> score for G4 and i-Motif, Z-Seeker for Z-DNA, tract length and A/T content for Curved DNA, repeat count for Sticky DNA, etc.<br>
    <b>References:</b>
    <ul>
        <li>Bedrat et al., 2016 Nucleic Acids Research</li>
        <li>Ho et al., 2010 Nature Chemical Biology</li>
        <li>Zeraati et al., 2018 Nature Chemistry</li>
        <li>Bacolla et al., 2006 Nucleic Acids Research</li>
        <li>Mirkin & Frank-Kamenetskii, 1994 Annual Review of Biophysics</li>
    </ul>
    </div>
    """, unsafe_allow_html=True)

# --- Footer ---
st.markdown("""
---
<div style='font-size: 15px; color: #1e293b; margin-top: 40px; text-align: left;'>
<b>Developed by</b><br>
Dr. Venkata Rajesh Yella<br>
<a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a> |
<a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
</div>
""", unsafe_allow_html=True)
