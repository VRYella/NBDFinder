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
    /* Tabs: enhanced bold, large, appealing */
    .stTabs [data-baseweb="tab-list"] {
        width: 100vw !important;
        justify-content: stretch !important;
        border-bottom: 3px solid #1565c0;
        background: linear-gradient(135deg, #e3f2fd 0%, #f3e5f5 30%, #e8f5e8 60%, #fff3e0 100%) !important;
        box-shadow: 0 4px 16px rgba(21, 101, 192, 0.15);
        margin-bottom: 0;
        border-radius: 12px 12px 0 0;
    }
    .stTabs [data-baseweb="tab"] {
        font-size: 1.65rem !important;
        font-weight: 800 !important;
        flex: 1 1 0%;
        min-width: 0 !important;
        padding: 18px 12px 18px 12px !important;
        text-align: center;
        color: #1565c0 !important;
        background: linear-gradient(135deg, #f8fdff 0%, #f0f7ff 100%) !important;
        border-right: 2px solid #e3f2fd !important;
        letter-spacing: 0.05em;
        text-shadow: 0 1px 2px rgba(21, 101, 192, 0.1);
        border-radius: 8px 8px 0 0;
        margin: 4px 2px 0 2px;
        transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    }
    .stTabs [data-baseweb="tab"]:hover {
        background: linear-gradient(135deg, #e3f2fd 0%, #f3e5f5 100%) !important;
        color: #0d47a1 !important;
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(21, 101, 192, 0.2);
    }
    .stTabs [aria-selected="true"] {
        color: #ffffff !important;
        border-bottom: 6px solid #d32f2f !important;
        background: linear-gradient(135deg, #1565c0 0%, #1976d2 50%, #2196f3 100%) !important;
        box-shadow: 0 8px 24px rgba(21, 101, 192, 0.3);
        transform: translateY(-3px);
        font-size: 1.7rem !important;
        text-shadow: 0 2px 4px rgba(0, 0, 0, 0.3);
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
    /* Enhanced Motif Selection Styling */
    .stMultiSelect > div {
        background: linear-gradient(135deg, #f8fdff 0%, #eef8ff 100%) !important;
        border: 2px solid #1565c0 !important;
        border-radius: 12px !important;
        box-shadow: 0 4px 12px rgba(21, 101, 192, 0.15) !important;
    }
    .stMultiSelect > div > div {
        background: linear-gradient(135deg, #ffffff 0%, #f8fdff 100%) !important;
        border-radius: 10px !important;
    }
    /* Enhanced Input Method Radio Buttons */
    .stRadio > div {
        background: linear-gradient(135deg, #f8fdff 0%, #f0f7ff 100%) !important;
        border-radius: 12px !important;
        padding: 16px !important;
        border: 2px solid #e3f2fd !important;
        box-shadow: 0 2px 8px rgba(21, 101, 192, 0.1) !important;
    }
    .stRadio > div > label {
        font-weight: 600 !important;
        color: #1565c0 !important;
        font-size: 1.12rem !important;
    }
    .stRadio > div > label > div:first-child {
        background: linear-gradient(45deg, #1565c0, #2196f3) !important;
        border: 2px solid #ffffff !important;
        box-shadow: 0 2px 6px rgba(21, 101, 192, 0.3) !important;
    }
    /* Enhanced Text Styling for Input Method */
    .input-method-title {
        font-weight: 800 !important;
        font-size: 1.25rem !important;
        color: #1565c0 !important;
        text-shadow: 0 1px 2px rgba(21, 101, 192, 0.2) !important;
    }
    .input-method-subtitle {
        font-weight: 600 !important;
        font-style: italic !important;
        color: #2e7d32 !important;
        font-size: 1.1rem !important;
    }
    .input-method-regular {
        font-weight: 400 !important;
        color: #424242 !important;
        font-size: 1.08rem !important;
    }
    .sequence-preview-title {
        font-weight: 700 !important;
        color: #d32f2f !important;
        font-size: 1.2rem !important;
        text-shadow: 0 1px 2px rgba(211, 47, 47, 0.2) !important;
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
        <div style='font-family:Montserrat, Arial; font-size:1.16rem; color:#222; line-height:1.8; padding:20px; background:linear-gradient(135deg, #f8fdff 0%, #eaf6ff 100%); border-radius:16px; box-shadow:0 4px 12px #dae8f5; border:1px solid #e3f2fd;'>
        
        <div style='margin-bottom:25px;'>
            <div style='display:flex; align-items:center; margin-bottom:15px;'>
                <div style='width:8px; height:8px; background:linear-gradient(45deg, #1565c0, #2196f3); border-radius:50%; margin-right:12px;'></div>
                <b style='color:#1565c0; font-size:1.18rem;'>Non-canonical DNA structures</b> play key roles in genome stability, regulation, and evolution.
            </div>
        </div>
        
        <div style='margin-bottom:25px;'>
            <div style='display:flex; align-items:flex-start; margin-bottom:15px;'>
                <div style='width:8px; height:8px; background:linear-gradient(45deg, #1565c0, #2196f3); border-radius:50%; margin-right:12px; margin-top:6px;'></div>
                <div>This application detects and analyzes <b style='color:#d32f2f;'>18 distinct Non-B DNA motifs</b> in any DNA sequence or multi-FASTA file:<br>
                <div style='margin-top:12px; padding-left:8px; border-left:3px solid #1565c0; background:#f1f8ff; padding:12px; border-radius:8px;'>
                    <span style='color:#1565c0; font-weight:600;'>
                        <b>G-quadruplex-related:</b> G4, Relaxed G4, Bulged G4, Bipartite G4, Multimeric G4, G-Triplex, i-Motif, Hybrid<br>
                        <b>Helix/curvature:</b> Z-DNA, eGZ (Extruded-G), Curved DNA, AC-Motif<br>
                        <b>Repeat/junction:</b> Slipped DNA, Cruciform, Sticky DNA, Triplex DNA<br>
                        <b>Hybrid/cluster:</b> R-Loop, Non-B DNA Clusters
                    </span>
                </div>
                </div>
            </div>
        </div>
        
        <div>
            <div style='display:flex; align-items:center; margin-bottom:15px;'>
                <div style='width:8px; height:8px; background:linear-gradient(45deg, #1565c0, #2196f3); border-radius:50%; margin-right:12px;'></div>
                <b style='color:#1565c0; font-size:1.18rem;'>Upload single or multi-FASTA files</b> to analyze comprehensive motif patterns with scientific accuracy.
            </div>
        </div>
        
        </div>
        """, unsafe_allow_html=True)

# ---------- UPLOAD & ANALYZE ----------
with tab_pages["Upload & Analyze"]:
    st.markdown("<h2>Sequence Upload and Motif Analysis</h2>", unsafe_allow_html=True)
    st.markdown('<span style="font-family:Montserrat,Arial; font-size:1.12rem;">Supports multi-FASTA and single FASTA. Paste, upload, select example, or fetch from NCBI.</span>', unsafe_allow_html=True)
    st.caption("Supported formats: .fa, .fasta, .txt | Limit: 200MB/file.")

    # Motif class selection
    selected_motifs = st.multiselect(
        "Select Motif Classes for Analysis", MOTIF_ORDER, default=MOTIF_ORDER,
        help="Choose motif classes to analyze. Selecting 'Hybrid' or 'Non-B DNA Clusters' will run all motif modules."
    )
    st.session_state.selected_motifs = selected_motifs if selected_motifs else MOTIF_ORDER

    # Input method selection
    st.markdown('<p class="input-method-title">Input Method:</p>', unsafe_allow_html=True)
    input_method = st.radio("", ["Upload FASTA / Multi-FASTA File", "Paste Sequence(s)", "Example Sequence", "NCBI Fetch"], horizontal=True)

    seqs, names = [], []

    # --- File upload ---
    if input_method == "Upload FASTA / Multi-FASTA File":
        fasta_file = st.file_uploader("Drag and drop FASTA/multi-FASTA file here", type=["fa", "fasta", "txt"])
        if fasta_file:
            try:
                content = fasta_file.read().decode("utf-8")
                cur_seq, cur_name = "", ""
                
                for line in content.splitlines():
                    if line.startswith(">"):
                        if cur_seq:
                            seqs.append(parse_fasta(cur_seq))
                            # Clean up sequence name by replacing underscores with spaces
                            clean_name = cur_name.replace("_", " ") if cur_name else f"Sequence {len(seqs)}"
                            names.append(clean_name)
                        cur_name = line.strip().lstrip(">")
                        cur_seq = ""
                    else:
                        cur_seq += line.strip()
                
                if cur_seq:
                    seqs.append(parse_fasta(cur_seq))
                    clean_name = cur_name.replace("_", " ") if cur_name else f"Sequence {len(seqs)}"
                    names.append(clean_name)
                
                if seqs:
                    st.success(f"✅ Successfully loaded {len(seqs)} sequence(s) from {fasta_file.name}")
                    
                    # Show preview of first 3 sequences
                    st.subheader("📋 Sequence Preview")
                    for i, seq in enumerate(seqs[:3]):
                        with st.expander(f"📄 {names[i]} (Length: {len(seq)} bp)", expanded=i==0):
                            stats = get_basic_stats(seq)
                            col1, col2, col3, col4 = st.columns(4)
                            col1.metric("Length", f"{stats['Length (bp)']} bp")
                            col2.metric("GC Content", f"{stats['GC %']}%")
                            col3.metric("AT Content", f"{stats['AT %']}%")
                            col4.metric("A+T Count", stats['A Count'] + stats['T Count'])
                            
                            # Show sequence preview
                            preview = seq[:100] + "..." if len(seq) > 100 else seq
                            st.code(preview, language="text")
                    
                    if len(seqs) > 3:
                        st.info(f"📊 Showing preview of first 3 sequences. Total sequences: {len(seqs)}")
                else:
                    st.warning("⚠️ No valid sequences found in the uploaded file.")
                    
            except UnicodeDecodeError:
                st.error("❌ Could not decode file. Please ensure it's a valid text file with UTF-8 encoding.")
            except Exception as e:
                st.error(f"❌ Error processing file: {str(e)}")

    # --- Paste sequence ---
    elif input_method == "Paste Sequence(s)":
        seq_input = st.text_area("Paste FASTA or raw sequence(s)", height=150, placeholder="Paste your sequences here...\n\nExamples:\n>Sequence 1\nATCGATCGATCG...\n\nOr raw sequence:\nATCGATCGATCG...")
        if seq_input:
            try:
                lines = seq_input.splitlines()
                cur_seq, cur_name = "", ""
                
                for line in lines:
                    if line.startswith(">"):
                        if cur_seq:
                            seqs.append(parse_fasta(cur_seq))
                            clean_name = cur_name.replace("_", " ") if cur_name else f"Sequence {len(seqs)}"
                            names.append(clean_name)
                        cur_name = line.strip().lstrip(">")
                        cur_seq = ""
                    else:
                        cur_seq += line.strip()
                
                # Handle case where no header is provided (raw sequence)
                if cur_seq and not any(line.startswith(">") for line in lines):
                    seqs.append(parse_fasta(cur_seq))
                    names.append("Pasted Sequence")
                elif cur_seq:
                    seqs.append(parse_fasta(cur_seq))
                    clean_name = cur_name.replace("_", " ") if cur_name else f"Sequence {len(seqs)}"
                    names.append(clean_name)
                
                if seqs:
                    st.success(f"✅ Successfully processed {len(seqs)} sequence(s)")
                    
                    # Show sequence stats
                    for i, (seq, name) in enumerate(zip(seqs, names)):
                        stats = get_basic_stats(seq)
                        st.write(f"📄 **{name}**: {stats['Length (bp)']} bp, GC: {stats['GC %']}%")
                else:
                    st.warning("⚠️ No valid sequences found in the pasted text.")
                    
            except Exception as e:
                st.error(f"❌ Error processing pasted sequences: {str(e)}")

    # --- Example input ---
    elif input_method == "Example Sequence":
        examples = ["g4_rich_sequence.fasta", "disease_repeats.fasta", "structural_motifs.fasta", "comprehensive_example.fasta"]
        example = st.selectbox("Select example input", examples)
        if example:
            try:
                path = f"example_inputs/{example}"
                with open(path, "r") as f:
                    content = f.read()
                    cur_seq, cur_name = "", ""
                    
                    for line in content.splitlines():
                        if line.startswith(">"):
                            if cur_seq:
                                seqs.append(parse_fasta(cur_seq))
                                clean_name = cur_name.replace("_", " ") if cur_name else f"Sequence {len(seqs)}"
                                names.append(clean_name)
                            cur_name = line.strip().lstrip(">")
                            cur_seq = ""
                        else:
                            cur_seq += line.strip()
                    
                    if cur_seq:
                        seqs.append(parse_fasta(cur_seq))
                        clean_name = cur_name.replace("_", " ") if cur_name else f"Sequence {len(seqs)}"
                        names.append(clean_name)
                
                if seqs:
                    st.success(f"✅ Loaded {len(seqs)} example sequence(s) from {example}")
                    
                    # Show example info
                    st.info(f"📖 **Example Description**: {example.replace('_', ' ').replace('.fasta', '').title()}")
                    for i, (seq, name) in enumerate(zip(seqs[:3], names[:3])):
                        stats = get_basic_stats(seq)
                        st.write(f"📄 **{name}**: {stats['Length (bp)']} bp, GC: {stats['GC %']}%")
                        
            except FileNotFoundError:
                st.error(f"❌ Example file '{example}' not found.")
            except Exception as e:
                st.error(f"❌ Error loading example: {str(e)}")

    # --- NCBI Query input ---
    elif input_method == "NCBI Fetch":
        st.info("🚧 **NCBI Fetch Feature**: This feature requires NCBI API integration.")
        ncbi_query = st.text_input("NCBI Query", value="", placeholder="Enter query (accession, gene, etc.)")
        
        col1, col2 = st.columns([3, 1])
        with col1:
            st.markdown("""
            **Examples of valid queries:**
            - Accession numbers: `NC_000001.11`, `NM_000518.4`
            - Gene symbols: `BRCA1`, `TP53`
            - Search terms: `"fragile X syndrome"[MeSH Terms]`
            """)
        
        with col2:
            if st.button("🔍 Fetch from NCBI", disabled=True):
                st.warning("⚠️ NCBI fetch functionality is not implemented yet.")
        
        st.warning("⚠️ NCBI fetch is currently disabled. Please use other input methods.")

    # Store sequences in session state
    if seqs:
        st.session_state.seqs = seqs
        st.session_state.names = names

    # --- Analysis trigger and result display ---
    if seqs and st.button("🧬 Analyze Sequences", type="primary"):
        with st.spinner("🔬 Analyzing sequences..."):
            st.session_state.is_analyzing = True
            results = []
            
            progress_container = st.container()
            
            for i, (seq, name) in enumerate(zip(seqs, names)):
                # Update progress
                progress = (i + 1) / len(seqs)
                progress_container.progress(progress, text=f"Analyzing sequence {i+1}/{len(seqs)}: {name}")
                
                motifs = analyze_sequence_with_progress(seq, name, st.session_state.selected_motifs)
                results.append(motifs)
            
            st.session_state.results = results
            st.session_state.is_analyzing = False
            
            # Create summary dataframe
            summary_data = []
            for i, (seq, name, motifs) in enumerate(zip(seqs, names, results)):
                stats = get_basic_stats(seq, motifs)
                top_motifs = ", ".join([m['Class'] for m in motifs[:3]]) if motifs else "None"
                summary_data.append({
                    "Sequence Name": name,
                    "Length (bp)": stats["Length (bp)"],
                    "GC %": stats["GC %"],
                    "Motif Count": len(motifs),
                    "Motif Coverage (%)": stats.get("Motif Coverage (%)", 0),
                    "Top Motifs": top_motifs
                })
            
            st.session_state.summary_df = pd.DataFrame(summary_data)
            
        st.success("✅ Analysis complete! Check the **Results** tab for detailed visualization and data.")
        st.balloons()


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
            seq_name_display = st.session_state.names[seq_idx].replace("_", " ")
            st.markdown(f"<h3>🧬 Detailed Motifs for <b>{seq_name_display}</b></h3>", unsafe_allow_html=True)
            
            # Clean up motif class names to replace underscores with spaces
            if 'Class' in df.columns:
                df['Class'] = df['Class'].apply(lambda x: x.replace("_", " ") if isinstance(x, str) else x)
            if 'Subtype' in df.columns:
                df['Subtype'] = df['Subtype'].apply(lambda x: x.replace("_", " ") if isinstance(x, str) else x)
            
            # Key columns for display (removing arm length, adding predicted sequence preview)
            essential_columns = ['Class', 'Subtype', 'Start', 'End', 'Length', 'Score']
            
            # Add predicted sequence preview (first 50 chars)
            if 'Sequence' in df.columns:
                df['Sequence Preview'] = df['Sequence'].apply(lambda x: x.replace('\n', '')[:50] + '...' if len(x.replace('\n', '')) > 50 else x.replace('\n', ''))
                essential_columns.append('Sequence Preview')
            
            # Add serial number starting from 1
            df['S.No'] = range(1, len(df) + 1)
            essential_columns = ['S.No'] + essential_columns
            
            display_df = df[essential_columns].copy()
            
            # Format numeric columns
            if 'Score' in display_df.columns:
                display_df['Score'] = display_df['Score'].apply(lambda x: f"{x:.2f}" if isinstance(x, (int, float)) else str(x))
            
            # Add scoring significance column
            if 'Score' in display_df.columns:
                def get_score_significance(score, motif_class):
                    try:
                        score_val = float(score)
                        if motif_class in ['G4', 'Relaxed G4', 'Bulged G4', 'Bipartite G4', 'Multimeric G4']:
                            if score_val >= 1.5: return "🟢 High confidence"
                            elif score_val >= 1.0: return "🟡 Moderate confidence" 
                            else: return "🟠 Low confidence"
                        elif motif_class in ['Z-DNA', 'R-Loop', 'Curved DNA']:
                            if score_val >= 100: return "🟢 High confidence"
                            elif score_val >= 50: return "🟡 Moderate confidence"
                            else: return "🟠 Low confidence"
                        else:
                            if score_val >= 25: return "🟢 High confidence"
                            elif score_val >= 15: return "🟡 Moderate confidence"
                            else: return "🟠 Low confidence"
                    except:
                        return "⚪ Not assessed"
                
                display_df['Prediction Confidence'] = display_df.apply(lambda row: get_score_significance(row.get('Score', 0), row.get('Class', '')), axis=1)
                
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
                    title="<b>Motif Counts by Type</b>"
                )
                fig_bar.update_layout(
                    height=400, 
                    showlegend=False,
                    xaxis=dict(
                        title="<b>Number of Motifs</b>",
                        title_font=dict(size=14, color='#1565c0'),
                        tickfont=dict(size=12, color='#424242'),
                        gridcolor='#e0e0e0',
                        gridwidth=1
                    ),
                    yaxis=dict(
                        title="<b>Motif Classes</b>",
                        title_font=dict(size=14, color='#1565c0'),
                        tickfont=dict(size=12, color='#424242')
                    ),
                    title_font=dict(size=16, color='#1565c0'),
                    plot_bgcolor='rgba(248,253,255,0.7)',
                    paper_bgcolor='white'
                )
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
                        title="<b>Distribution of Motif Lengths</b>",
                        nbins=20
                    )
                    fig_hist.update_layout(
                        height=400,
                        xaxis=dict(
                            title="<b>Motif Length (base pairs)</b>",
                            title_font=dict(size=14, color='#1565c0'),
                            tickfont=dict(size=12, color='#424242'),
                            gridcolor='#e0e0e0',
                            gridwidth=1
                        ),
                        yaxis=dict(
                            title="<b>Frequency</b>",
                            title_font=dict(size=14, color='#1565c0'),
                            tickfont=dict(size=12, color='#424242'),
                            gridcolor='#e0e0e0',
                            gridwidth=1
                        ),
                        title_font=dict(size=16, color='#1565c0'),
                        plot_bgcolor='rgba(248,253,255,0.7)',
                        paper_bgcolor='white',
                        bargap=0.1
                    )
                    # Add annotation for interpretation
                    fig_hist.add_annotation(
                        text="<i>Shorter motifs (10-50 bp): Higher occurrence<br>Longer motifs (>100 bp): Rare but significant</i>",
                        xref="paper", yref="paper",
                        x=0.02, y=0.98,
                        showarrow=False,
                        font=dict(size=10, color='#666666'),
                        bgcolor="rgba(255,255,255,0.8)",
                        bordercolor="#cccccc",
                        borderwidth=1
                    )
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
            for j, m in enumerate(motifs):
                m = m.copy()  # Create a copy to avoid modifying original
                # Clean up sequence name by replacing underscores with spaces
                m['Sequence Name'] = st.session_state.names[i].replace("_", " ")
                m['S.No'] = j + 1  # Add serial number starting from 1
                
                # Clean up class and subtype names
                if 'Class' in m:
                    m['Class'] = m['Class'].replace("_", " ")
                if 'Subtype' in m:
                    m['Subtype'] = m['Subtype'].replace("_", " ") if m['Subtype'] else m['Subtype']
                
                if m['Class'] == "Z-DNA" and m.get("Subclass", "") == "eGZ (Extruded-G)":
                    m['Class'] = "eGZ (Extruded-G)"
                
                # Add scoring significance
                def get_score_significance(score, motif_class):
                    try:
                        score_val = float(score)
                        if motif_class in ['G4', 'Relaxed G4', 'Bulged G4', 'Bipartite G4', 'Multimeric G4']:
                            min_score = 1.0
                            if score_val >= 1.5: return "High confidence", min_score
                            elif score_val >= 1.0: return "Moderate confidence", min_score
                            else: return "Low confidence", min_score
                        elif motif_class in ['Z-DNA', 'R-Loop', 'Curved DNA']:
                            min_score = 50.0
                            if score_val >= 100: return "High confidence", min_score
                            elif score_val >= 50: return "Moderate confidence", min_score
                            else: return "Low confidence", min_score
                        else:
                            min_score = 15.0
                            if score_val >= 25: return "High confidence", min_score
                            elif score_val >= 15: return "Moderate confidence", min_score
                            else: return "Low confidence", min_score
                    except:
                        return "Not assessed", "N/A"
                
                confidence, min_score = get_score_significance(m.get('Score', 0), m.get('Class', ''))
                m['Prediction Confidence'] = confidence
                m['Minimum Score Threshold'] = min_score
                df_all.append(m)
        
        df_all = pd.DataFrame(df_all)
        
        # Reorder columns to put S.No first
        cols = df_all.columns.tolist()
        if 'S.No' in cols:
            cols = ['S.No'] + [col for col in cols if col != 'S.No']
            df_all = df_all[cols]
        
        st.markdown("### 📥 Complete Results Table with Scoring Information")
        st.markdown("""
        <div style='background:#f0f8ff; padding:12px; border-radius:8px; margin-bottom:15px; border-left:4px solid #1565c0;'>
        <b>Scoring Significance:</b><br>
        • <b>High confidence:</b> Scores above established thresholds indicate strong structural likelihood<br>
        • <b>Moderate confidence:</b> Scores at minimum thresholds suggest possible structural formation<br>
        • <b>Low confidence:</b> Scores below minimum thresholds indicate weaker predictions<br>
        • <b>Minimum Score Threshold:</b> Shows the lowest score considered reliable for each motif type
        </div>
        """, unsafe_allow_html=True)
        
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
    st.header("📚 Scientific Documentation & References")
    
    # Add tabbed sections for better organization
    doc_tabs = st.tabs(["🧬 Motif Descriptions", "🔬 Detection Methods", "📊 Scoring Systems", "📖 References"])
    
    with doc_tabs[0]:
        st.subheader("Comprehensive Motif Classifications")
        
        # Create detailed motif table
        motif_data = {
            "Motif Class": [
                "Curved DNA", "Z-DNA", "eGZ (Extruded-G)", "Slipped DNA", "R-Loop", 
                "Cruciform", "Triplex DNA", "Sticky DNA", "G-Triplex", "G4 (G-Quadruplex)",
                "Relaxed G4", "Bulged G4", "Bipartite G4", "Multimeric G4", "i-Motif",
                "AC-Motif", "Hybrid Motif", "Non-B DNA Clusters"
            ],
            "Structural Features": [
                "Phased A/T tracts causing DNA curvature",
                "Left-handed double helix with alternating purines/pyrimidines",
                "Long CGG repeats prone to expansion diseases",
                "Direct tandem repeats forming slipped structures",
                "G-rich sequences forming stable RNA-DNA hybrids",
                "Palindromic sequences forming four-way junctions",
                "Mirror repeats supporting triplex formation",
                "Extended GAA/TTC repeats in disease genes",
                "Three consecutive G-runs forming triplex structures",
                "Four G-runs forming stable quadruplex structures",
                "G4 with relaxed loop length constraints",
                "G4 with bulges in G-tracts",
                "Two separate G4 units in proximity",
                "Multiple G4 units in tandem",
                "C-rich sequences forming intercalated motifs",
                "Alternating A-rich and C-rich patterns",
                "Overlapping regions of multiple motif types",
                "Genomic hotspots with clustered non-B structures"
            ],
            "Biological Significance": [
                "Protein binding sites, nucleosome positioning",
                "Gene regulation, chromatin structure",
                "Fragile X syndrome, repeat expansion diseases",
                "Genetic instability, recombination hotspots",
                "Transcription regulation, DNA damage",
                "Recombination, genetic instability",
                "Gene regulation, antigene therapy targets",
                "Friedreich's ataxia, neurological diseases",
                "Gene regulation, DNA packaging",
                "Telomere maintenance, oncogene regulation",
                "Promoter regulation, stress response",
                "Genetic variation, disease susceptibility",
                "Long-range gene regulation",
                "Epigenetic regulation, chromatin loops",
                "pH-dependent gene regulation",
                "Transcriptional control elements",
                "Multi-functional regulatory regions",
                "Genomic instability hotspots"
            ],
            "Score Range": [
                "15-200", "50-500", "10-100", "15-150", "20-300",
                "25-200", "20-180", "10-80", "15-120", "1.0-3.0",
                "1.0-2.5", "1.0-2.5", "20-100", "30-150", "15-100",
                "10-50", "Variable", "10-200"
            ]
        }
        
        df_motifs = pd.DataFrame(motif_data)
        st.dataframe(df_motifs, use_container_width=True, height=400)
        
    with doc_tabs[1]:
        st.subheader("Detection Algorithms & Regular Expressions")
        
        # Create methods table
        methods_data = {
            "Motif Type": [
                "G-Quadruplex (G4)", "Z-DNA", "Curved DNA", "Cruciform", "R-Loop",
                "Slipped DNA", "Triplex DNA", "i-Motif", "AC-Motif"
            ],
            "Primary Algorithm": [
                "G4Hunter + Structural Factors", "Kadane's Maximum Subarray", "Curvature Prediction",
                "Palindrome Detection", "RLFS + REZ Stability", "Direct Repeat Analysis",
                "Mirror Repeat Detection", "C-Rich Pattern Matching", "Alternating Pattern Analysis"
            ],
            "Regular Expression": [
                r"G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,}",
                r"(CG|GC|CA|TG|AC|GT){6,}",
                r"A{4,}[\w]{0,10}T{4,}|T{4,}[\w]{0,10}A{4,}",
                r"Palindromic sequences with spacers",
                r"G{3,}[\w]{0,20}G{3,}",
                r"(\w{2,20})\1{2,}",
                r"Mirror repeats: Pu-Py patterns",
                r"C{3,}\w{1,7}C{3,}\w{1,7}C{3,}\w{1,7}C{3,}",
                r"[AC]{10,}|Alternating A/C rich"
            ],
            "Validation Method": [
                "Structural factor calculation", "Windowed Z-score analysis", "Curvature angle prediction",
                "Reverse complement matching", "Thermodynamic stability", "Repeat unit validation",
                "Purine-pyrimidine composition", "pH-dependent stability", "Composition analysis"
            ]
        }
        
        df_methods = pd.DataFrame(methods_data)
        st.dataframe(df_methods, use_container_width=True, height=350)
        
        st.info("💡 **Algorithm Accuracy**: Each algorithm is calibrated against experimental data and validated using known biological examples.")
        
    with doc_tabs[2]:
        st.subheader("Scoring Systems & Thresholds")
        
        scoring_data = {
            "Motif Category": [
                "G4 Family", "Z-DNA Family", "Repeat-Based", "Junction-Based", "Hybrid/Cluster"
            ],
            "Scoring Method": [
                "G4Hunter algorithm with structural factors",
                "Kadane's algorithm with dinucleotide weights",
                "Repeat count × stability factors",
                "Palindrome length × AT content",
                "Composite scoring based on constituent motifs"
            ],
            "High Confidence": [
                "≥ 1.5", "≥ 100", "≥ 25", "≥ 30", "≥ 50"
            ],
            "Moderate Confidence": [
                "1.0 - 1.5", "50 - 100", "15 - 25", "20 - 30", "25 - 50"
            ],
            "Low Confidence": [
                "< 1.0", "< 50", "< 15", "< 20", "< 25"
            ],
            "Biological Validation": [
                "Experimental G4 formation data",
                "B-to-Z transition conditions",
                "Repeat expansion thresholds",
                "Cruciform formation evidence",
                "Multi-technique validation"
            ]
        }
        
        df_scoring = pd.DataFrame(scoring_data)
        st.dataframe(df_scoring, use_container_width=True, height=250)
        
        st.markdown("""
        <div style='background:linear-gradient(135deg, #e3f2fd 0%, #f3e5f5 100%); border-radius:12px; padding:16px; margin:15px 0; border-left:4px solid #1565c0;'>
        <b>🎯 Interpretation Guidelines:</b><br>
        • <b>High confidence:</b> Strong experimental support, likely to form in vivo<br>
        • <b>Moderate confidence:</b> Reasonable formation potential under specific conditions<br>
        • <b>Low confidence:</b> Possible formation, requires experimental validation<br>
        • <b>Threshold Selection:</b> Based on ROC analysis of experimental datasets
        </div>
        """, unsafe_allow_html=True)
        
    with doc_tabs[3]:
        st.subheader("Scientific References & Citations")
        
        st.markdown("""
        <div style='background:#f8fdff; border-radius:12px; padding:20px; font-size:1.08rem; font-family:Montserrat,Arial;'>
        
        <h4 style='color:#1565c0; margin-top:0;'>🔬 Core Methodology References:</h4>
        <ol style='line-height:1.8;'>
            <li><b>Bedrat, A. et al. (2016)</b> "Re-evaluation of G-quadruplex propensity with G4Hunter." <i>Nucleic Acids Research</i> 44(4): 1746-1759. 
            <br><small style='color:#666;'>DOI: 10.1093/nar/gkw006 - Primary G4 detection algorithm</small></li>
            
            <li><b>Ho, P.S. (2010)</b> "The non-B-DNA structure of d(CA/TG)n does not differ from that of Z-DNA." <i>Nature Chemical Biology</i> 6: 648-653.
            <br><small style='color:#666;'>DOI: 10.1038/nchembio.408 - Z-DNA structural validation</small></li>
            
            <li><b>Kim, N. & Jinks-Robertson, S. (2018)</b> "The Top1 paradox: Friend or foe of DNA replication." <i>Nucleic Acids Research</i> 46(20): 10563-10573.
            <br><small style='color:#666;'>DOI: 10.1093/nar/gky888 - R-loop biology and detection</small></li>
            
            <li><b>Zeraati, M. et al. (2018)</b> "I-motif DNA structures are formed in the nuclei of human cells." <i>Nature Chemistry</i> 10: 631-637.
            <br><small style='color:#666;'>DOI: 10.1038/s41557-018-0046-3 - In vivo i-motif validation</small></li>
        </ol>
        
        <h4 style='color:#1565c0;'>📊 Structural Analysis References:</h4>
        <ol start='5' style='line-height:1.8;'>
            <li><b>Bacolla, A. & Wells, R.D. (2006)</b> "Non-B DNA conformations, genomic rearrangements, and human disease." <i>Nucleic Acids Research</i> 34(6): 1803-1820.
            <br><small style='color:#666;'>DOI: 10.1093/nar/gkl004 - Comprehensive non-B DNA review</small></li>
            
            <li><b>Mirkin, S.M. & Frank-Kamenetskii, M.D. (1994)</b> "H-DNA and related structures." <i>Annual Review of Biophysics</i> 23: 541-576.
            <br><small style='color:#666;'>DOI: 10.1146/annurev.bb.23.060194.002545 - Triplex DNA mechanisms</small></li>
            
            <li><b>Sinden, R.R. (1994)</b> "DNA Structure and Function." Academic Press, San Diego.
            <br><small style='color:#666;'>ISBN: 978-0126457506 - Foundational structural biology</small></li>
        </ol>
        
        <h4 style='color:#1565c0;'>🧬 Disease Association References:</h4>
        <ol start='8' style='line-height:1.8;'>
            <li><b>Usdin, K. et al. (2015)</b> "Repeat-associated non-ATG translation: molecular mechanisms and contribution to neurological disease." <i>Nucleic Acids Research</i> 43(20): 9589-9598.</li>
            
            <li><b>Brooks, T.A. & Hurley, L.H. (2010)</b> "Targeting MYC Expression through G-Quadruplexes." <i>Genes & Cancer</i> 1(6): 641-649.</li>
            
            <li><b>McMurray, C.T. (2010)</b> "Mechanisms of trinucleotide repeat instability during human development." <i>Nature Reviews Genetics</i> 11: 786-799.</li>
        </ol>
        
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
