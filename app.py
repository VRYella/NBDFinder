import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import io
import time
from collections import Counter
from Bio import Entrez, SeqIO

from motifs import (
    all_motifs, 
    find_hotspots,
    parse_fasta, gc_content, reverse_complement,
    select_best_nonoverlapping_motifs, wrap
)

# ---------- PATCH: Ensure every motif has Subtype ----------
def ensure_subtype(motif):
    """Guarantee every motif has a string 'Subtype'"""
    if isinstance(motif, dict):
        if 'Subtype' not in motif or motif['Subtype'] is None:
            motif['Subtype'] = 'Other'
        return motif
    else:
        # Handle non-dict motifs gracefully (could log/warn here)
        return {'Subtype': 'Other', 'Motif': motif}


# ---------- PROFESSIONAL CSS FOR BALANCED DESIGN ----------
st.markdown("""
    <style>
    body, [data-testid="stAppViewContainer"], .main {
        background: #f7fafd !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
    }
    /* Tabs: medium-large, bold, clean */
    .stTabs [data-baseweb="tab-list"] {
        width: 100vw !important;
        justify-content: stretch !important;
        border-bottom: 2px solid #1565c0;
        background: linear-gradient(90deg,#eaf3fa 0%,#f7fafd 100%) !important;
        box-shadow: 0 2px 8px #dae5f2;
        margin-bottom: 0;
    }
    .stTabs [data-baseweb="tab"] {
        font-size: 1.45rem !important;
        font-weight: 700 !important;
        flex: 1 1 0%;
        min-width: 0 !important;
        padding: 15px 0 15px 0 !important;
        text-align: center;
        color: #1565c0 !important;
        background: #eaf3fa !important;
        border-right: 1px solid #eee !important;
        letter-spacing: 0.03em;
    }
    .stTabs [aria-selected="true"] {
        color: #002147 !important;
        border-bottom: 5px solid #1565c0 !important;
        background: #f7fafd !important;
        box-shadow: 0 4px 8px #e0e5ea;
    }
    .stTabs [data-baseweb="tab"]:last-child {
        border-right: none !important;
    }
    /* Headings: harmonized medium size */
    h1, h2, h3, h4 {
        font-family: 'Montserrat', Arial, sans-serif !important;
        color: #1565c0 !important;
        font-weight: 800 !important;
        margin-top: 0.8em;
        margin-bottom: 0.8em;
    }
    h1 { font-size:2.05rem !important; }
    h2 { font-size:1.55rem !important; }
    h3 { font-size:1.19rem !important; }
    h4 { font-size:1.09rem !important; }
    /* Markdown/text: medium size, easy reading */
    .stMarkdown, .markdown-text-container, .stText, p, span, label, input, .stTextInput>div>div>input, .stSelectbox>div>div>div, .stMultiSelect>div>div>div, .stRadio>div>div>label>div {
        font-size: 1.08rem !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
    }
    /* Buttons: modern, medium */
    .stButton>button {
        font-size: 1.08rem !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
        padding: 0.45em 1.2em !important;
        background: linear-gradient(90deg,#1565c0 0%,#2e8bda 100%) !important;
        color: #fff !important;
        border-radius: 7px !important;
        border: none !important;
        font-weight: 600 !important;
        box-shadow: 0 2px 8px #b5cbe6;
        transition: background 0.2s;
    }
    .stButton>button:hover {
        background: linear-gradient(90deg,#2e8bda 0%,#1565c0 100%) !important;
    }
    /* DataFrame font */
    .stDataFrame, .stTable {
        font-size: 1.05rem !important;
        font-family: 'Montserrat', Arial, sans-serif !important;
    }
    </style>
""", unsafe_allow_html=True)

# ---------- PAGE CONFIG ----------
st.set_page_config(
    page_title="Non-B DNA Motif Finder",
    layout="wide",
    page_icon="🧬",
    menu_items={'About': "Non-B DNA Motif Finder | Developed by Dr. Venkata Rajesh Yella"}
)

# ---------- CONSTANTS ----------
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
    'analysis_start_time': None,
    'analysis_progress': 0,
    'is_analyzing': False,
    'stop_analysis': False,
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

def create_progress_tracker():
    """Create a progress tracking container"""
    if st.session_state.is_analyzing:
        progress_container = st.container()
        with progress_container:
            col1, col2, col3 = st.columns([2, 1, 1])
            
            with col1:
                progress_bar = st.progress(st.session_state.analysis_progress)
                if st.session_state.analysis_start_time:
                    elapsed = time.time() - st.session_state.analysis_start_time
                    st.write(f"⏱️ Analysis time: {elapsed:.1f}s | Status: {st.session_state.analysis_status}")
            
            with col2:
                if st.button("⏸️ Stop Analysis", key="stop_btn"):
                    st.session_state.stop_analysis = True
                    st.session_state.is_analyzing = False
                    st.warning("Analysis stopped by user")
                    
            with col3:
                st.write("🔄 Processing...")
        
        return progress_container
    return None

def analyze_sequence_with_progress(seq, seq_name, selected_motifs):
    """Analyze sequence with progress updates"""
    total_steps = len(selected_motifs) if selected_motifs else 18
    
    # Initialize progress
    st.session_state.analysis_progress = 0
    st.session_state.analysis_status = "Starting analysis..."
    
    # Check if we need to run all motifs
    run_all = any(m in selected_motifs for m in ["Hybrid", "Non-B DNA Clusters"])
    
    if run_all:
        st.session_state.analysis_status = "Running comprehensive motif analysis..."
        motifs = all_motifs(seq)
        st.session_state.analysis_progress = 0.8
    else:
        motifs = all_motifs(seq)
        # Filter by selected motifs
        motifs = [m for m in motifs if m['Class'] in selected_motifs]
        st.session_state.analysis_progress = 0.7
    
    st.session_state.analysis_status = "Processing results..."
    # PATCH: Ensure every motif has a 'Subtype'
    motifs = [ensure_subtype(m) for m in motifs]
    nonoverlapping = select_best_nonoverlapping_motifs(motifs)
    
    st.session_state.analysis_progress = 1.0
    st.session_state.analysis_status = "Analysis complete!"
    
    return nonoverlapping

def create_enhanced_motif_visualization(motifs, seq_name, seq_length):
    """Create enhanced interactive motif visualization using Plotly"""
    if not motifs:
        st.warning("No motifs to visualize")
        return
    
    # Create interactive plot
    fig = go.Figure()
    
    # Color mapping for motifs
    y_positions = {}
    y_counter = 0
    
    for i, motif in enumerate(motifs):
        motif_class = motif['Class']
        if motif_class == "Z-DNA" and motif.get("Subclass", "") == "eGZ (Extruded-G)":
            motif_class = "eGZ (Extruded-G)"
        
        if motif_class not in y_positions:
            y_positions[motif_class] = y_counter
            y_counter += 1
        
        color = MOTIF_COLORS.get(motif_class, "#888888")
        
        # Add motif as a horizontal bar
        fig.add_trace(go.Scatter(
            x=[motif['Start'], motif['End']],
            y=[y_positions[motif_class], y_positions[motif_class]],
            mode='lines',
            line=dict(color=color, width=8),
            name=motif_class,
            showlegend=motif_class not in [trace.name for trace in fig.data],
            hovertemplate=f"<b>{motif_class}</b><br>" +
                         f"Position: {motif['Start']}-{motif['End']}<br>" +
                         f"Length: {motif['Length']} bp<br>" +
                         f"Score: {motif.get('Score', 'N/A')}<br>" +
                         "<extra></extra>"
        ))
    
    # Update layout
    fig.update_layout(
        title=f"Interactive Motif Map: {seq_name}",
        xaxis_title="Sequence Position (bp)",
        yaxis_title="Motif Classes",
        yaxis=dict(
            tickmode='array',
            tickvals=list(y_positions.values()),
            ticktext=list(y_positions.keys()),
            showgrid=True
        ),
        height=max(400, len(y_positions) * 50),
        hovermode='closest',
        showlegend=True
    )
    
    return fig

# ---- TABS ----
tabs = st.tabs(list(PAGES.keys()))
tab_pages = dict(zip(PAGES.keys(), tabs))

# ---------- HOME ----------
with tab_pages["Home"]:
    st.markdown("<h1>Non-B DNA Motif Finder</h1>", unsafe_allow_html=True)
    left, right = st.columns([1,1])
    with left:
        st.image("nbdcircle.JPG", use_container_width=True)
    with right:
        st.markdown("""
        <div style='font-family:Montserrat, Arial; font-size:1.14rem; color:#222; line-height:1.7; padding:18px; background:#f8f9fa; border-radius:14px; box-shadow:0 2px 8px #eee;'>
        <b>Non-canonical DNA structures</b> play key roles in genome stability, regulation, and evolution.<br>
        This application detects and analyzes <b>18 distinct Non-B DNA motifs</b> in any DNA sequence or multi-FASTA file.<br>
        <b>Motif Classes:</b><br>
        <span style='color:#1565c0;'>
            <b>G-quadruplex-related</b> (G4, Relaxed G4, Bulged G4, Bipartite G4, Multimeric G4, G-Triplex, i-Motif, Hybrid),<br>
            <b>helix/curvature</b> (Z-DNA, eGZ (Extruded-G), Curved DNA, AC-Motif),<br>
            <b>repeat/junction</b> (Slipped DNA, Cruciform, Sticky DNA, Triplex DNA),<br>
            <b>hybrid/cluster</b> (R-Loop, Non-B DNA Clusters).
        </span>
        <br>
        <b>Upload single or multi-FASTA files...</b>
        </div>
        """, unsafe_allow_html=True)

# ---------- UPLOAD & ANALYZE ----------
with tab_pages["Upload & Analyze"]:
    st.markdown("<h2>Sequence Upload and Motif Analysis</h2>", unsafe_allow_html=True)
    st.markdown('<span style="font-family:Montserrat,Arial; font-size:1.12rem;">Supports multi-FASTA and single FASTA. Paste, upload, select example, or fetch from NCBI.</span>', unsafe_allow_html=True)
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

        # Enhanced analysis section with progress tracking
        st.markdown("---")
        col1, col2 = st.columns([3, 1])
        
        with col1:
            run_all = any(m in st.session_state.selected_motifs for m in ["Hybrid", "Non-B DNA Clusters"])
            
            if not st.session_state.is_analyzing:
                if st.button("🚀 Start Motif Analysis", type="primary", key="start_analysis"):
                    st.session_state.is_analyzing = True
                    st.session_state.analysis_start_time = time.time()
                    st.session_state.analysis_progress = 0
                    st.session_state.stop_analysis = False
                    st.rerun()
            else:
                st.button("🚀 Analysis Running...", disabled=True, key="analysis_running")
        
        with col2:
            if st.session_state.is_analyzing:
                st.metric("Status", "🔄 Running")
            else:
                st.metric("Status", "⏸️ Ready")
        
        # Progress tracking container
        progress_container = create_progress_tracker()
        
        # Run analysis if triggered
        if st.session_state.is_analyzing and not st.session_state.stop_analysis:
            with st.spinner("Analyzing sequences..."):
                motif_results = []
                total_seqs = len(st.session_state.seqs)
                
                for idx, seq in enumerate(st.session_state.seqs):
                    if st.session_state.stop_analysis:
                        break
                        
                    st.session_state.analysis_status = f"Processing sequence {idx+1}/{total_seqs}..."
                    st.session_state.analysis_progress = idx / total_seqs
                    
                    try:
                        motifs = analyze_sequence_with_progress(seq, st.session_state.names[idx], st.session_state.selected_motifs)
                        motif_results.append(motifs)
                    except Exception as e:
                        st.error(f"Error analyzing sequence {idx+1}: {str(e)}")
                        motif_results.append([])
                
                if not st.session_state.stop_analysis:
                    st.session_state.results = motif_results
                    
                    # Generate summary statistics
                    summary = []
                    for i, motifs in enumerate(motif_results):
                        stats = get_basic_stats(st.session_state.seqs[i], motifs)
                        motif_types = Counter([m['Class'] if m['Class'] != "Z-DNA" or m.get("Subclass") != "eGZ (Extruded-G)" else "eGZ (Extruded-G)" for m in motifs])
                        summary.append({
                            "Sequence Name": st.session_state.names[i],
                            "Length (bp)": stats['Length (bp)'],
                            "GC %": stats['GC %'],
                            "Motif Count": len(motifs),
                            "Motif Coverage (%)": stats["Motif Coverage (%)"],
                            "Top Motifs": ", ".join([f"{k}({v})" for k, v in sorted(motif_types.items(), key=lambda x: x[1], reverse=True)[:3]])
                        })
                    
                    st.session_state.summary_df = pd.DataFrame(summary)
                    st.success("✅ Analysis complete! View results in the 'Results' tab.")
                
                # Reset analysis state
                st.session_state.is_analyzing = False
                st.session_state.analysis_status = "Ready"
                st.rerun()

# ---------- RESULTS ----------
with tab_pages["Results"]:
    st.markdown('<h2>Analysis Results and Visualization</h2>', unsafe_allow_html=True)
    if not st.session_state.results:
        st.info("No analysis results. Please run motif analysis first.")
    else:
        # Enhanced summary table with key metrics only
        st.subheader("📊 Analysis Summary")
        summary_cols = ["Sequence Name", "Length (bp)", "GC %", "Motif Count", "Motif Coverage (%)", "Top Motifs"]
        st.dataframe(st.session_state.summary_df[summary_cols], use_container_width=True)
        
        seq_idx = 0
        if len(st.session_state.seqs) > 1:
            seq_idx = st.selectbox("Choose Sequence for Detailed Analysis:", range(len(st.session_state.seqs)), format_func=lambda i: st.session_state.names[i])
        
        motifs = st.session_state.results[seq_idx]
        if not motifs:
            st.warning("No motifs detected for this sequence.")
        else:
            df = pd.DataFrame(motifs)
            
            # Enhanced motif table with essential columns only
            st.markdown(f"<h3>🧬 Detailed Motifs for <b>{st.session_state.names[seq_idx]}</b></h3>", unsafe_allow_html=True)
            
            # Key columns for display (removing arm length, adding predicted sequence preview)
            essential_columns = ['Class', 'Subtype', 'Start', 'End', 'Length', 'Score']
            
            # Add predicted sequence preview (first 50 chars)
            if 'Sequence' in df.columns:
                df['Sequence Preview'] = df['Sequence'].apply(lambda x: x.replace('\n', '')[:50] + '...' if len(x.replace('\n', '')) > 50 else x.replace('\n', ''))
                essential_columns.append('Sequence Preview')
            
            display_df = df[essential_columns].copy()
            
            # Format numeric columns
            if 'Score' in display_df.columns:
                display_df['Score'] = display_df['Score'].apply(lambda x: f"{x:.2f}" if isinstance(x, (int, float)) else str(x))
            
            st.dataframe(display_df, use_container_width=True, height=400)
            
            # Enhanced visualizations
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown('<h4>📈 Motif Type Distribution</h4>', unsafe_allow_html=True)
                
                # Create interactive bar chart
                class_counts = df['Class'].value_counts()
                if "Subclass" in df.columns:
                    egz_count = (df["Subclass"] == "eGZ (Extruded-G)").sum()
                    if egz_count > 0:
                        class_counts["eGZ (Extruded-G)"] = egz_count
                
                fig_bar = px.bar(
                    x=class_counts.values, 
                    y=class_counts.index,
                    orientation='h',
                    color=class_counts.index,
                    color_discrete_map=MOTIF_COLORS,
                    title="Motif Counts by Type"
                )
                fig_bar.update_layout(height=400, showlegend=False)
                st.plotly_chart(fig_bar, use_container_width=True)
            
            with col2:
                st.markdown('<h4>📊 Motif Length Distribution</h4>', unsafe_allow_html=True)
                
                # Create motif length histogram
                if 'Length' in df.columns:
                    fig_hist = px.histogram(
                        df, 
                        x='Length', 
                        color='Class',
                        color_discrete_map=MOTIF_COLORS,
                        title="Distribution of Motif Lengths"
                    )
                    fig_hist.update_layout(height=400)
                    st.plotly_chart(fig_hist, use_container_width=True)
            
            # Enhanced interactive motif map
            st.markdown('<h4>🗺️ Interactive Motif Map</h4>', unsafe_allow_html=True)
            
            # Create the enhanced visualization
            interactive_fig = create_enhanced_motif_visualization(motifs, st.session_state.names[seq_idx], len(st.session_state.seqs[seq_idx]))
            if interactive_fig:
                st.plotly_chart(interactive_fig, use_container_width=True)
            
            # Additional sequence composition analysis
            st.markdown('<h4>🔬 Sequence Composition Analysis</h4>', unsafe_allow_html=True)
            seq = st.session_state.seqs[seq_idx]
            stats = get_basic_stats(seq, motifs)
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("GC Content", f"{stats['GC %']}%")
            with col2:
                st.metric("AT Content", f"{stats['AT %']}%")
            with col3:
                st.metric("Motif Coverage", f"{stats['Motif Coverage (%)']}%")
            with col4:
                st.metric("Total Motifs", len(motifs))
                
            # Motif density heatmap
            if len(motifs) > 0:
                st.markdown('<h4>🌡️ Motif Density Heatmap</h4>', unsafe_allow_html=True)
                
                # Create density array
                seq_len = len(seq)
                window_size = max(50, seq_len // 20)  # Adaptive window size
                density = []
                positions = []
                
                for i in range(0, seq_len, window_size):
                    window_end = min(i + window_size, seq_len)
                    window_motifs = [m for m in motifs if m['Start'] >= i and m['End'] <= window_end]
                    density.append(len(window_motifs))
                    positions.append(f"{i}-{window_end}")
                
                if density:
                    fig_heatmap = px.bar(
                        x=positions,
                        y=density,
                        title=f"Motif Density (window size: {window_size} bp)",
                        labels={'x': 'Sequence Position', 'y': 'Motif Count'}
                    )
                    fig_heatmap.update_layout(height=300)
                    st.plotly_chart(fig_heatmap, use_container_width=True)

# ---------- DOWNLOAD ----------
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

# ---------- DOCUMENTATION ----------
with tab_pages["Documentation"]:
    st.header("Scientific Documentation & References")
    st.markdown("""
    <div style='background:#f4faff; border-radius:12px; padding:18px; font-size:1.08rem; font-family:Montserrat,Arial;'>
    <b>Motif Classes Detected:</b><br><br>
    <ul>
        <li><b>Curved DNA</b>: Identifies phased poly(A) or poly(T) tracts using regex and spacing rules, reflecting intrinsic curvature. Scoring is based on tract length/grouping.</li>
        <li><b>Z-DNA</b>: Detects alternating purine-pyrimidine patterns, GC-rich segments. Uses windowed scoring; regex finds dinucleotide repeats.</li>
        <li><b>eGZ-motif (Extruded-G Z-DNA)</b>: Searches for long (CGG)<sub>n</sub> runs via regex. Scored by repeat count.</li>
        <li><b>Slipped DNA</b>: Recognizes direct/tandem repeats by repeat-unit matching and regex. Scoring by length and unit copies.</li>
        <li><b>R-Loop</b>: Finds G-rich regions for stable RNA-DNA hybrids; RLFS model and regex. Thermodynamic scoring for hybrid stability.</li>
        <li><b>Cruciform</b>: Finds palindromic inverted repeats with spacers, regex and reverse complement. Scoring by arm length and A/T content.</li>
        <li><b>Triplex DNA / Mirror Repeat</b>: Detects purine/pyrimidine mirror repeats/triplex motifs. Regex identifies units; scoring by composition/purity.</li>
        <li><b>Sticky DNA</b>: Searches extended GAA/TTC repeats. Scoring by repeat count.</li>
        <li><b>G-Triplex</b>: Finds three consecutive guanine runs by regex and loop length. Scoring by G-run sum and loop penalty.</li>
        <li><b>G4 (G-Quadruplex) and Variants</b>: Detects canonical/variant G4 motifs by G-run/loop regex. G4Hunter scoring for content/structure.</li>
        <li><b>i-Motif</b>: C-rich sequences for i-motif under acid. Regex for C runs/loops; scoring by run count and content.</li>
        <li><b>AC-Motif</b>: Alternating A-rich/C-rich consensus regions by regex. Scoring by pattern presence.</li>
        <li><b>Hybrid Motif</b>: Regions where motif classes overlap; found by interval intersection, scored on diversity/size.</li>
        <li><b>Non-B DNA Clusters</b>: Hotspots with multiple motifs in a window; sliding algorithm, scored by motif count/diversity.</li>
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
<div style='font-size: 1.05rem; color: #1e293b; margin-top: 30px; text-align: left; font-family:Montserrat,Arial;'>
<b>Developed by</b><br>
Dr. Venkata Rajesh Yella<br>
<a href='mailto:yvrajesh_bt@kluniversity.in'>yvrajesh_bt@kluniversity.in</a> |
<a href='https://github.com/VRYella' target='_blank'>GitHub: VRYella</a>
</div>
""", unsafe_allow_html=True)
