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

# Motif order as defined in motifs.py (all_motifs function)
MOTIF_ORDER = [
    "Sticky_DNA",
    "Curved_DNA",
    "Z-DNA",
    "eGZ (extruded-G)",
    "Slipped_DNA",
    "R-Loop",
    "Cruciform",
    "Triplex_DNA",
    "G-Triplex",
    "G4",
    "Relaxed_G4",
    "Bulged_G4",
    "Bipartite_G4",
    "Multimeric_G4",
    "i-Motif",
    "AC-Motif",
    "Hybrid",
    "Non-B DNA Clusters",
]

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
    'seqs': [],
    'names': [],
    'results': [],
    'summary_df': pd.DataFrame(),
    'hotspots': [],
    'analysis_status': "Ready",
}.items():
    if k not in st.session_state:
        st.session_state[k] = v

MOTIF_COLORS = {
    "Curved_DNA": "#FF9AA2",
    "Z-DNA": "#FFB7B2",
    "Slipped_DNA": "#FFDAC1",
    "R-Loop": "#FFD3B6",
    "Cruciform": "#E2F0CB",
    "Triplex_DNA": "#B5EAD7",
    "Sticky_DNA": "#DCB8CB",
    "G-Triplex": "#C7CEEA",
    "G4": "#A2D7D8",
    "Relaxed_G4": "#A2D7B8",
    "Bulged_G4": "#A2A7D8",
    "Bipartite_G4": "#A2D788",
    "Multimeric_G4": "#A2A7B8",
    "i-Motif": "#B0C4DE",
    "Hybrid": "#C1A192",
    "Non-B DNA Clusters": "#A2C8CC",
    "eGZ (extruded-G)": "#6A4C93",
    "AC-Motif": "#F5B041"
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
        <b>Detect 14+ Non-Canonical DNA Motifs</b> in any DNA sequence or multi-FASTA file.<br>
        <ul>
        <li>Upload FASTA or multi-FASTA files (no need to specify single/multi format)</li>
        <li>Motifs detected include Z-DNA, eGZ-motif, AC-motif, G4, i-Motif, R-loop, Cruciform, Triplex, Hybrids, and more</li>
        <li>Interactive motif visualizations and downloadable results</li>
        </ul>
        <b>Developed by Dr. Venkata Rajesh Yella</b>
        </div>
        """, unsafe_allow_html=True
    )

# --- Upload & Analyze page ---
elif page == "Upload & Analyze":
    # ... (same as your original code for upload/input/preview) ...
    pass

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
            motif_class_order = MOTIF_ORDER
            st.markdown("#### Motif Type Distribution")
            fig, ax = plt.subplots(figsize=(8,6))
            class_counts = df['Class'].value_counts().reindex(motif_class_order, fill_value=0)
            if "Subclass" in df.columns:
                egz_count = (df["Subclass"] == "eGZ (extruded-G)").sum()
                class_counts["eGZ (extruded-G)"] = egz_count
            ax.barh(class_counts.index, class_counts.values, color=[MOTIF_COLORS.get(c, "#888") for c in class_counts.index])
            ax.set_xlabel("Motif Count")
            st.pyplot(fig)
            st.markdown("<br>", unsafe_allow_html=True)
            # Motif positions
            st.markdown("#### Motif Map")
            fig, ax = plt.subplots(figsize=(12,3))
            y = 1
            for _, row in df.iterrows():
                motif_class = row['Class']
                if motif_class == "Z-DNA" and row.get("Subclass", "") == "eGZ (extruded-G)":
                    motif_class = "eGZ (extruded-G)"
                color = MOTIF_COLORS.get(motif_class, "#888")
                ax.plot([row['Start'], row['End']], [y, y], lw=8, color=color, alpha=0.8)
                y += 0.12
            ax.set_yticks([])
            ax.set_xlabel("Sequence Position (bp)")
            ax.set_xlim(0, len(st.session_state.seqs[seq_idx]))
            ax.set_title(f"Motif tracks: {st.session_state.names[seq_idx]}")
            st.pyplot(fig)

# --- Download/Documentation/Footer pages ---
# ... (same as your original code) ...

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
                if m['Class'] == "Z-DNA" and m.get("Subclass", "") == "eGZ (extruded-G)":
                    m['Class'] = "eGZ (extruded-G)"
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
        <li><b>eGZ-motif (extruded-G Z-DNA):</b> Long (CGG)<sub>n</sub> runs, a special Z-DNA variant. <i>Kim et al., 2018</i>.</li>
        <li><b>Slipped DNA:</b> Direct repeats and short tandem repeats. <i>Bacolla et al., 2006</i>.</li>
        <li><b>R-Loop:</b> RLFS models for G-rich skew and thermodynamic stability. <i>Sanz et al., 2016</i>.</li>
        <li><b>Cruciform:</b> Inverted repeats with AT-rich arms. <i>Lilley, 1985</i>.</li>
        <li><b>Triplex DNA / Mirror Repeat:</b> Purine/pyrimidine-rich mirror repeats and triplex-forming motifs. <i>Mirkin, 1994</i>.</li>
        <li><b>Sticky DNA:</b> Extended GAA/TTC repeats. <i>Potaman et al., 2003</i>.</li>
        <li><b>G-Triplex & G4:</b> Canonical, relaxed, bulged, bipartite, multimeric, and imperfect G-quadruplexes. <i>Bedrat et al., 2016</i>.</li>
        <li><b>i-Motif:</b> C-rich, looped sequences. <i>Zeraati et al., 2018</i>.</li>
        <li><b>AC-motif:</b> Consensus A-rich/C-rich motif. <i>New et al., 2020</i>.</li>
        <li><b>Hybrids:</b> Overlapping regions of two or more motif classes.</li>
        <li><b>Non-B DNA Clusters:</b> Hotspot regions with high motif density and diversity.</li>
    </ul>
    <b>Scoring:</b><br>
    - <b>G4Hunter</b> score for G4 and i-Motif, Z-Seeker for Z-DNA, tract length and A/T content for Curved DNA, repeat count for Sticky DNA, normalized repeats for eGZ-motif, pattern match for AC-motif, etc.<br>
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
