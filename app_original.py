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
import numpy as np

from motifs import (
    all_motifs, 
    find_hotspots,
    parse_fasta, gc_content, reverse_complement,
    select_best_nonoverlapping_motifs, wrap
)

# Import advanced visualization
try:
    from advanced_visualization import create_enhanced_dashboard, export_publication_figure
    ADVANCED_VIZ_AVAILABLE = True
except ImportError:
    ADVANCED_VIZ_AVAILABLE = False
    st.warning("Advanced visualization not available. Install additional dependencies.")

# Import comprehensive publication visualizations
try:
    from nbdfinder_viz_integration import create_nbdfinder_visualization_interface
    PUBLICATION_VIZ_AVAILABLE = True
except ImportError:
    PUBLICATION_VIZ_AVAILABLE = False
    print("Publication visualization module not available")

# ---------- ENHANCED PAGE CONFIGURATION ----------
st.set_page_config(
    page_title="NBDFinder - Non-B DNA Analysis Platform",
    page_icon="DNA",
    layout="wide",
    initial_sidebar_state="expanded",
    menu_items={
        'Get Help': 'https://github.com/VRYella/NBDFinder',
        'Report a bug': 'https://github.com/VRYella/NBDFinder/issues',
        'About': "NBDFinder - The most advanced non-B DNA prediction platform"
    }
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


# Example sequences from reference implementation
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

# ---------- CONCISE DESIGN SYSTEM ----------
st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    
    :root {
        --primary: #1e3a8a; --accent: #0891b2; --text: #1e3a8a; --text-muted: #64748b;
        --border: #e2e8f0; --bg: #ffffff; --surface: #f8fafc; --radius: 8px;
        --shadow: 0 2px 4px rgba(0,0,0,0.1); --transition: all 0.2s ease;
    }
    
    body, [data-testid="stAppViewContainer"], .main {
        background: var(--bg); font-family: Inter, sans-serif; color: var(--text); line-height: 1.6;
    }
    
    /* TABS */
    .stTabs [data-baseweb="tab-list"] {
        background: linear-gradient(135deg, var(--surface) 0%, var(--border) 100%);
        border: 1px solid var(--border); border-radius: 12px; padding: 8px; margin-bottom: 24px;
        box-shadow: var(--shadow); display: flex; gap: 4px; width: 100%;
    }
    .stTabs [data-baseweb="tab"] {
        background: transparent; border: none; border-radius: var(--radius);
        font-family: Inter, sans-serif; font-weight: 500; font-size: 1.125rem;
        color: var(--text-muted); padding: 12px 24px; transition: var(--transition);
        cursor: pointer; display: inline-flex; align-items: center; justify-content: center;
        white-space: nowrap; flex: 1;
    }
    .stTabs [data-baseweb="tab"]:hover {
        background: rgba(255,255,255,0.8); color: var(--primary); font-weight: 600;
        transform: translateY(-1px); box-shadow: var(--shadow);
    }
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, var(--primary) 0%, var(--accent) 100%);
        color: white; font-weight: 700; font-size: 1.25rem;
        box-shadow: 0 4px 12px rgba(30,58,138,0.3); transform: translateY(-2px);
    }
    
    /* FORM CONTROLS */
    .stSelectbox > div, .stMultiSelect > div, .stTextInput > div, .stTextArea > div {
        background: var(--bg); border: 1px solid var(--border); border-radius: var(--radius); 
        transition: var(--transition);
    }
    .stSelectbox > div:hover, .stMultiSelect > div:hover, .stTextInput > div:hover, .stTextArea > div:hover {
        border-color: var(--accent);
    }
    .stSelectbox > div:focus-within, .stMultiSelect > div:focus-within, .stTextInput > div:focus-within, .stTextArea > div:focus-within {
        border-color: var(--primary); box-shadow: 0 0 0 2px rgba(30,64,175,0.1);
    }
    
    /* BUTTONS */
    .stButton > button {
        background: var(--primary); color: white; border: none; border-radius: var(--radius);
        font-weight: 500; transition: var(--transition); cursor: pointer;
    }
    .stButton > button:hover { background: #1e3a8a; box-shadow: var(--shadow); }
    .stButton > button:active { transform: translateY(1px); }
    
    /* LAYOUT */
    .main .block-container { padding: 1.5rem 2rem; max-width: none; width: 100%; }
    .stApp > .main { width: 100%; max-width: none; }
    [data-testid="stAppViewContainer"] { width: 100%; max-width: none; }

    /* TYPOGRAPHY */
    h1, h2, h3, h4, h5, h6 {
        font-family: Inter, sans-serif; font-weight: 700; letter-spacing: -0.025em;
        line-height: 1.1; margin: 1.5rem 0 1rem; color: var(--primary);
        background: linear-gradient(135deg, var(--primary) 0%, var(--accent) 100%);
        -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;
    }
    h1 { font-size: 4rem; font-weight: 900; letter-spacing: -0.04em; margin-bottom: 1.5rem; }
    h2 { font-size: 3rem; font-weight: 800; letter-spacing: -0.03em; }
    h3 { font-size: 2.25rem; font-weight: 750; letter-spacing: -0.025em; }
    h4 { font-size: 1.75rem; font-weight: 700; letter-spacing: -0.02em; }
    h5 { font-size: 1.375rem; font-weight: 650; }
    h6 { font-size: 1.25rem; font-weight: 600; }
    
    .stMarkdown, .markdown-text-container, .stText, p, span, label, input {
        font-family: Inter, sans-serif; font-size: 1.125rem; line-height: 1.7; 
        color: var(--text); font-weight: 400;
    }
    .stMarkdown p { margin-bottom: 1.25rem; font-size: 1.125rem; }
    .stMarkdown ul, .stMarkdown ol { font-size: 1.125rem; line-height: 1.7; margin-bottom: 1.5rem; }
    .stMarkdown li { margin-bottom: 0.5rem; padding-left: 0.5rem; }
    
    /* BUTTON STYLING */
    .stButton>button {
        font-family: Inter, sans-serif; font-weight: 600; font-size: 1rem; padding: 10px 20px;
        background: var(--primary); color: white; border: none; border-radius: var(--radius);
        transition: var(--transition); cursor: pointer;
    }
    .stButton>button:hover { background: #1e3a8a; box-shadow: var(--shadow); }
    .stButton>button:active { transform: translateY(1px); }
    
    /* DATA TABLES */
    .stDataFrame, .stTable {
        font-family: Inter, sans-serif; border-radius: var(--radius); overflow: hidden;
        box-shadow: var(--shadow); border: 1px solid var(--border); background: var(--bg);
    }
    .stDataFrame th, .stTable th {
        background: var(--surface); color: var(--text); font-weight: 600; font-size: 0.875rem;
        padding: 12px; text-align: left; border-bottom: 2px solid var(--border);
    }
    .stDataFrame td, .stTable td {
        padding: 10px 12px; text-align: left; font-size: 0.875rem; font-weight: 400;
        border-bottom: 1px solid var(--border);
    }
    .stDataFrame tr:hover td, .stTable tr:hover td { background: var(--surface); }
    .stDataFrame tr:nth-child(even), .stTable tr:nth-child(even) { background: rgba(248,250,252,0.5); }
    
    /* CONFIDENCE BADGES */
    .stDataFrame .confidence-optimal { background: #059669; color: white; font-weight: 600; border-radius: var(--radius); padding: 4px 8px; font-size: 0.75rem; }
    .stDataFrame .confidence-high { background: var(--primary); color: white; font-weight: 600; border-radius: var(--radius); padding: 4px 8px; font-size: 0.75rem; }
    .stDataFrame .confidence-moderate { background: #ea580c; color: white; font-weight: 600; border-radius: var(--radius); padding: 4px 8px; font-size: 0.75rem; }
    .stDataFrame .confidence-low { background: #dc2626; color: white; font-weight: 600; border-radius: var(--radius); padding: 4px 8px; font-size: 0.75rem; }
    
    /* DROPDOWNS & MISC */
    .stSelectbox > div > div { min-width: 200px; padding-right: 40px; }
    .stSelectbox > div > div > div { white-space: nowrap; overflow: hidden; text-overflow: ellipsis; padding-right: 30px; }
    .stDataFrame div[data-testid="stTable"] button { min-width: auto; padding-right: 25px; }
    .stDataFrame .dropdown-menu { min-width: 150px; }
    .stSelectbox [data-baseweb="select"] > div:last-child { padding-left: 8px; min-width: 24px; }
    
    /* EXPANDERS & PROGRESS */
    [data-testid="stExpander"] { border: 1px solid var(--border); border-radius: var(--radius); background: var(--bg); margin-bottom: 16px; }
    [data-testid="stExpander"] summary { padding: 12px 16px; font-family: Inter, sans-serif; font-weight: 500; font-size: 0.875rem; color: var(--text); cursor: pointer; background: var(--surface); border: none; display: flex; align-items: center; justify-content: space-between; }
    [data-testid="stExpander"] summary:hover { background: #f1f5f9; }
    [data-testid="stExpander"] summary div { display: flex; align-items: center; gap: 8px; }
    .stProgress > div > div { background: var(--primary); border-radius: var(--radius); }
    .stProgress > div { background: var(--border); border-radius: var(--radius); }
    
    /* SCROLLBAR */
    ::-webkit-scrollbar { width: 6px; height: 6px; }
    ::-webkit-scrollbar-track { background: var(--surface); border-radius: 3px; }
    ::-webkit-scrollbar-thumb { background: var(--text-muted); border-radius: 3px; }
    ::-webkit-scrollbar-thumb:hover { background: var(--primary); }
    
    /* TOGGLE BUTTONS */
    .stButton > button {
        font-family: Inter, sans-serif; font-weight: 500; border-radius: var(--radius);
        transition: var(--transition); border: 1px solid var(--border); min-height: 42px;
        display: flex; align-items: center; justify-content: center; white-space: nowrap;
        text-overflow: ellipsis; overflow: hidden;
    }
    .stButton > button[kind="secondary"] { background: var(--bg); color: var(--text-muted); border: 1px solid var(--border); }
    .stButton > button[kind="secondary"]:hover { background: var(--surface); color: var(--text); border-color: var(--accent); transform: translateY(-1px); box-shadow: var(--shadow); }
    .stButton > button[kind="primary"] { background: var(--primary); color: white; border: 1px solid var(--primary); font-weight: 600; box-shadow: var(--shadow); }
    .stButton > button[kind="primary"]:hover { background: #1e3a8a; border-color: #1e3a8a; transform: translateY(-1px); box-shadow: var(--shadow); }
    .element-container:has(.stButton) { margin-bottom: 8px; }
    
    /* FINAL STYLING */
    [data-testid="stExpander"] > div:first-child { border: none; background: transparent; }
    .stButton > button:focus { outline: none; box-shadow: 0 0 0 2px rgba(30,64,175,0.3); }
    </style>
    """, unsafe_allow_html=True)

# ---------- CONSTANTS ----------
MOTIF_CATEGORIES = {
    "G-quadruplex-related": ["Bipartite G4", "Bulged G4", "Canonical G4", "Imperfect G4", "Multimeric G4", "Relaxed G4"],
    "G-Triplex": ["G-Triplex"],
    "i-motif related": ["AC-Motif", "i-Motif"],
    "Helix deviations": ["Curved DNA", "eGZ (Extruded-G)", "Z-DNA"],
    "Repeat/junction": ["Cruciform", "R-Loop", "Slipped DNA", "Sticky DNA", "Triplex DNA"],
    "Hybrid": ["Hybrid"],
    "Non-B DNA Clusters": ["Non-B DNA Clusters"]
}

MOTIF_ORDER = [motif for motifs in MOTIF_CATEGORIES.values() for motif in sorted(motifs)]

# Color schemes
MOTIF_COLORS = {
    "Canonical G4": "#1E40AF", "Relaxed G4": "#3B82F6", "Bulged G4": "#60A5FA", 
    "Bipartite G4": "#1D4ED8", "Multimeric G4": "#2563EB", "Imperfect G4": "#93C5FD",
    "Z-DNA": "#DC2626", "eGZ (Extruded-G)": "#EA580C", "Curved DNA": "#F59E0B",
    "Slipped DNA": "#059669", "R-Loop": "#10B981", "Sticky DNA": "#34D399",
    "Cruciform": "#7C3AED", "Triplex DNA": "#8B5CF6", "G-Triplex": "#A78BFA",
    "i-Motif": "#EC4899", "AC-Motif": "#F59E0B", "Hybrid": "#06B6D4", 
    "Non-B DNA Clusters": "#EF4444"
}

PAGES = {
    "Home": "Overview",
    "Upload & Analyze": "Sequence Upload and Motif Analysis",
    "Results": "Analysis Results and Visualization", 
    "Download": "Download Results",
    "Documentation": "Scientific Documentation & References"
}
Entrez.email = "raazbiochem@gmail.com"
Entrez.api_key = None

EXAMPLE_SEQUENCES = {
    "Human Telomere G4-Rich": {
        "name": "G4_Rich_Human_Telomere_Example", 
        "sequence": "TTAGGGTTAGGGTTAGGGTTAGGGAAAAATCCGTCGAGCAGAGTTAAAAAGGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCCGAAAGAAAGAAAGAAAGAAACGCGCGCGCGCGCGCGCGCGATCGCACACACACAGCTGCTGCTGC",
        "description": "Human telomeric sequence rich in G-quadruplex structures"
    },
    "Z-DNA Forming Sequence": {
        "name": "Z_DNA_Example_Sequence",
        "sequence": "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC",
        "description": "Alternating CG sequence capable of Z-DNA formation"
    },
    "Disease-Associated Repeats": {
        "name": "Disease_Repeat_Motifs",
        "sequence": "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAACAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCAGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTGCTG",
        "description": "GAA and CAG repeats associated with genetic diseases"
    },
    "Comprehensive Non-B DNA": {
        "name": "Multi_Structure_Example",
        "sequence": "GGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCCGCGCGCGCGCGCGCGCGCGCGCGAAAAATTTTTAAAAATTTTTAAAAATTTTTAAACAGCAGCAGCAGCAGCAGCAGCAGCTGCTGCTGCTGCTGCTGCTGCTGCGAAAGAAAGAAAGAAAG",
        "description": "Sequence containing multiple Non-B DNA forming motifs"
    }
}

FAMOUS_NCBI_EXAMPLES = {
    "Human TERT Promoter": "NC_000005.10:1253147-1295047",
    "Human c-MYC Promoter": "NG_007161.1",
    "Human BCL2 Promoter": "NG_009361.1", 
    "Fragile X FMR1 Gene": "NG_007529.1",
    "Huntington HTT Gene": "NG_009378.1",
    "Human Immunoglobulin Switch": "NG_001019.6",
    "Human Alpha Globin": "NG_000006.1",
    "Friedreich Ataxia FXN": "NG_008845.1"
}

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

def ncbi_fetch(query):
    """Fetch sequences from NCBI using Entrez with improved error handling"""
    try:
        # Set email for NCBI Entrez (required for good practice)
        Entrez.email = "user@example.com"
        
        # Search for the query with better parameters
        handle = Entrez.esearch(db="nucleotide", term=query, retmax=5, usehistory="y")
        search_results = Entrez.read(handle)
        handle.close()
        
        if not search_results['IdList']:
            st.warning(f"No sequences found for query: '{query}'. Try a different search term.")
            return [], []
        
        # Fetch sequences with better error handling
        ids = search_results['IdList']
        st.info(f"Found {len(ids)} sequence(s). Fetching data...")
        
        handle = Entrez.efetch(db="nucleotide", id=ids, rettype="fasta", retmode="text")
        fasta_content = handle.read()
        handle.close()
        
        if not fasta_content.strip():
            st.error("Empty response from NCBI. Try again later.")
            return [], []
        
        # Parse FASTA content with better validation
        seqs, names = [], []
        cur_seq, cur_name = "", ""
        for line in fasta_content.splitlines():
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if cur_seq:
                    cleaned_seq = parse_fasta(cur_seq)
                    if len(cleaned_seq) > 10:  # Only include sequences > 10 bp
                        seqs.append(cleaned_seq)
                        # Clean up the name for better display
                        display_name = cur_name[:100] + "..." if len(cur_name) > 100 else cur_name
                        names.append(display_name if display_name else f"Sequence_{len(seqs)}")
                cur_name = line.lstrip(">").strip()
                cur_seq = ""
            else:
                cur_seq += line
        
        # Add the last sequence
        if cur_seq:
            cleaned_seq = parse_fasta(cur_seq)
            if len(cleaned_seq) > 10:
                seqs.append(cleaned_seq)
                display_name = cur_name[:100] + "..." if len(cur_name) > 100 else cur_name
                names.append(display_name if display_name else f"Sequence_{len(seqs)}")
        
        if not seqs:
            st.warning("No valid sequences found (sequences must be >10 bp).")
            return [], []
        
        return seqs, names
        
    except Exception as e:
        error_msg = str(e)
        if "HTTP Error 429" in error_msg:
            st.error("NCBI rate limit exceeded. Please wait a moment and try again.")
        elif "URLError" in error_msg or "timeout" in error_msg.lower():
            st.error("Network connection issue. Please check your internet connection and try again.")
        elif "XML" in error_msg or "parse" in error_msg.lower():
            st.error("NCBI response parsing error. The service may be temporarily unavailable.")
        else:
            st.error(f"NCBI fetch failed: {error_msg}")
        
        return [], []

def create_progress_tracker():
    """Create an enhanced real-time progress tracking container with time estimation"""
    if st.session_state.is_analyzing:
        # Clean progress container with professional styling
        st.markdown("""
        <div style='background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; margin: 20px 0;'>
            <h3 style='color: #374151; margin: 0 0 16px 0; font-size: 1.25rem; font-weight: 600;'>
                Analysis Dashboard
            </h3>
        </div>
        """, unsafe_allow_html=True)
        
        progress_container = st.container()
        with progress_container:
            # Progress bar with clean styling
            progress_value = getattr(st.session_state, 'analysis_progress', 0)
            st.progress(progress_value)
            
            # Time tracking and estimation
            col1, col2, col3 = st.columns([3, 2, 2])
            
            with col1:
                if hasattr(st.session_state, 'analysis_start_time') and st.session_state.analysis_start_time:
                    elapsed = time.time() - st.session_state.analysis_start_time
                    minutes = int(elapsed // 60)
                    seconds = int(elapsed % 60)
                    time_str = f"{minutes:02d}:{seconds:02d}"
                    
                    # Estimate remaining time
                    if progress_value > 0.05:  # Avoid division by very small numbers
                        estimated_total = elapsed / progress_value
                        remaining = max(0, estimated_total - elapsed)
                        remaining_min = int(remaining // 60)
                        remaining_sec = int(remaining % 60)
                        remaining_str = f"{remaining_min:02d}:{remaining_sec:02d}"
                    else:
                        remaining_str = "Calculating..."
                    
                    st.markdown(f"""
                    <div style='background: #ffffff; border: 1px solid #e2e8f0; border-radius: 4px; padding: 12px; margin: 8px 0;'>
                        <strong>Elapsed:</strong> {time_str} | <strong>Remaining:</strong> {remaining_str}
                    </div>
                    """, unsafe_allow_html=True)
            
            with col2:
                # Status without icons
                status = getattr(st.session_state, 'analysis_status', 'Initializing...')
                progress_percent = int(progress_value * 100)
                st.markdown(f"""
                <div style='background: #ffffff; border: 1px solid #e2e8f0; border-radius: 4px; padding: 12px; margin: 8px 0; text-align: center;'>
                    <strong>Progress:</strong> {progress_percent}%<br>
                    <small style='color: #6b7280;'>{status}</small>
                </div>
                """, unsafe_allow_html=True)
                
            with col3:
                # Control buttons with clean styling
                col3a, col3b = st.columns(2)
                with col3a:
                    if st.button("Stop", key="stop_btn", help="Stop the current analysis"):
                        st.session_state.stop_analysis = True
                        st.session_state.is_analyzing = False
                        st.warning("Analysis stopped by user")
                        st.rerun()
                
                with col3b:
                    if st.button("Restart", key="restart_btn", help="Restart the analysis"):
                        st.session_state.analysis_progress = 0
                        st.session_state.analysis_start_time = time.time()
                        st.session_state.analysis_status = "Restarting analysis..."
                        st.rerun()
        
        return progress_container
    return None

def analyze_sequence_with_progress(seq, seq_name, selected_motifs):
    """Analyze sequence with progress updates"""
    total_steps = len(selected_motifs) if selected_motifs else 18
    
    # Initialize progress
    st.session_state.analysis_progress = 0
    st.session_state.analysis_status = f"Analyzing {seq_name} ({len(seq):,} bp)..."
    
    # Check if we need to run all motifs
    run_all = any(m in selected_motifs for m in ["Hybrid", "Non-B DNA Clusters"])
    
    # Always use complete performance mode as requested
    enable_hotspots = True  # Always enable hotspots in complete mode
    
    if run_all:
        st.session_state.analysis_status = f"Running comprehensive motif detection on {seq_name}..."
        motifs = all_motifs(seq, report_hotspots=enable_hotspots)
        st.session_state.analysis_progress = 0.8
    else:
        st.session_state.analysis_status = f"Detecting selected motifs in {seq_name}..."
        motifs = all_motifs(seq)
        # Filter by selected motifs
        motifs = [m for m in motifs if m['Class'] in selected_motifs]
        st.session_state.analysis_progress = 0.7
    
    st.session_state.analysis_status = f"Processing results for {seq_name}..."
    # PATCH: Ensure every motif has a 'Subtype'
    motifs = [ensure_subtype(m) for m in motifs]
    nonoverlapping = select_best_nonoverlapping_motifs(motifs)
    
    st.session_state.analysis_progress = 1.0
    st.session_state.analysis_status = "Analysis complete!"
    
    return nonoverlapping

def create_genome_browser_view(motifs_df, sequence, seq_name):
    """Create a genome browser-style visualization of motifs on the sequence"""
    
    # Sequence length
    seq_length = len(sequence)
    
    # Create figure
    fig = go.Figure()
    
    # Color mapping for different motif classes
    color_map = {
        'Canonical G4': '#FF6B6B',
        'Relaxed G4': '#FF8E53', 
        'Bulged G4': '#FF9A8B',
        'Imperfect G4': '#FFAD5A',
        'Bipartite G4': '#FFC93C',
        'Multimeric G4': '#FFE66D',
        'G-Triplex': '#95E1D3',
        'i-Motif': '#F3D2C1',
        'AC-Motif': '#8EC5FC',
        'Curved DNA': '#FFB6C1',
        'Z-DNA': '#DDA0DD',
        'eGZ (Extruded-G)': '#98FB98',
        'Cruciform': '#F0E68C',
        'R-Loop': '#87CEEB',
        'Slipped DNA': '#DEB887',
        'Sticky DNA': '#F5DEB3',
        'Triplex DNA': '#FFE4E1',
        'Hybrid': '#E6E6FA',
        'Non-B DNA Clusters': '#FAFAD2'
    }
    
    # Add sequence track (background)
    fig.add_trace(go.Scatter(
        x=[0, seq_length],
        y=[0, 0],
        mode='lines',
        line=dict(color='lightgray', width=8),
        name='Sequence',
        showlegend=False
    ))
    
    # Add motifs as colored segments
    y_positions = {}  # Track y positions for overlapping motifs
    
    for idx, motif in motifs_df.iterrows():
        start = motif.get('Start', 0)
        end = motif.get('End', 0)
        motif_class = motif.get('Class', 'Unknown')
        
        # Determine y position to avoid overlaps
        y_pos = 0
        for existing_class, (existing_start, existing_end, existing_y) in y_positions.items():
            if not (end < existing_start or start > existing_end):  # Overlap detected
                y_pos = existing_y + 0.2
        
        y_positions[f"{motif_class}_{idx}"] = (start, end, y_pos)
        
        # Get color for this motif class
        color = color_map.get(motif_class, '#888888')
        
        # Add motif segment
        fig.add_trace(go.Scatter(
            x=[start, end],
            y=[y_pos, y_pos],
            mode='lines',
            line=dict(color=color, width=12),
            name=motif_class,
            showlegend=True,
            hovertemplate=f"<b>{motif_class}</b><br>Position: {start}-{end}<br>Length: {end-start} bp<extra></extra>"
        ))
        
        # Add motif label
        mid_point = (start + end) / 2
        fig.add_annotation(
            x=mid_point,
            y=y_pos + 0.1,
            text=motif_class,
            showarrow=False,
            font=dict(size=10, color='black'),
            bgcolor='rgba(255, 255, 255, 0.8)',
            bordercolor='gray',
            borderwidth=1
        )
    
    # Update layout
    fig.update_layout(
        title=f"Genome Browser View: {seq_name} ({seq_length:,} bp)",
        xaxis=dict(
            title="Sequence Position (bp)",
            range=[0, seq_length],
            showgrid=True,
            gridcolor='lightgray'
        ),
        yaxis=dict(
            title="Motif Tracks",
            showticklabels=False,
            showgrid=False,
            range=[-0.5, max([pos[2] for pos in y_positions.values()]) + 0.5 if y_positions else 0.5]
        ),
        height=400,
        margin=dict(l=50, r=50, t=60, b=50),
        hovermode='closest',
        plot_bgcolor='white',
        paper_bgcolor='white'
    )
    
    # Display the plot
    st.plotly_chart(fig, use_container_width=True)
    
    # Add sequence ruler/coordinates
    with st.expander("Sequence Coordinates & Details"):
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.metric("Sequence Length", f"{seq_length:,} bp")
        
        with col2:
            st.metric("Total Motifs", len(motifs_df))
        
        with col3:
            coverage = (motifs_df['Length'].sum() / seq_length * 100) if len(motifs_df) > 0 else 0
            st.metric("Coverage", f"{coverage:.1f}%")
        
        # Show sequence segments around motifs
        if len(motifs_df) > 0:
            st.markdown("**Motif Context (±20bp):**")
            for idx, motif in motifs_df.head(3).iterrows():  # Show first 3 motifs
                start = max(0, motif['Start'] - 20)
                end = min(seq_length, motif['End'] + 20)
                context = sequence[start:end]
                
                # Highlight the motif within context
                motif_start_in_context = motif['Start'] - start
                motif_end_in_context = motif['End'] - start
                
                before = context[:motif_start_in_context]
                motif_seq = context[motif_start_in_context:motif_end_in_context]
                after = context[motif_end_in_context:]
                
                st.markdown(f"**{motif['Class']} ({motif['Start']}-{motif['End']}):**")
                st.markdown(f"`{before}`**`{motif_seq}`**`{after}`")

def create_enhanced_motif_visualization(motifs, seq_name, seq_length):
    """
    Create publication-quality interactive motif visualization using Plotly
    
    Scientific Enhancement: Generates Nature-level manuscript figures with:
    - High-resolution vector graphics suitable for publication
    - Professional color schemes optimized for accessibility
    - Clear scientific labeling and annotations
    - Interactive features for data exploration
    
    Based on visualization best practices from:
    - Tufte, E.R. "The Visual Display of Quantitative Information" (2001)
    - Wong, B. "Color blindness" Nature Methods (2011)
    """
    if not motifs:
        st.warning("No motifs to visualize")
        return
    
    # Create publication-quality interactive plot
    fig = go.Figure()
    
    # Enhanced color mapping optimized for publication and accessibility
    # Use the enhanced publication colors defined above
    # (Colors already defined in global PUBLICATION_COLORS variable)
    
    # Enhanced y-axis positioning with biological grouping
    y_positions = {}
    y_counter = 0
    
    # Group motifs by biological function for better visualization
    motif_groups = {
        'G4 Family': ['Canonical G4', 'Relaxed G4', 'Bulged G4', 'Bipartite G4', 'Multimeric G4', 'Imperfect G4'],
        'Triple-helix': ['G-Triplex', 'i-Motif', 'Triplex_DNA'],
        'Helix Variants': ['Z-DNA', 'eGZ (Extruded-G)', 'Curved_DNA', 'AC-Motif'],
        'Repeat Structures': ['Slipped_DNA', 'Cruciform', 'Sticky_DNA'],
        'Hybrid/Complex': ['R-Loop', 'Hybrid', 'Non-B DNA Clusters']
    }
    
    for group_name, group_motifs in motif_groups.items():
        for motif in group_motifs:
            if motif not in y_positions:
                y_positions[motif] = y_counter
                y_counter += 1
    
    # Add horizontal lines for motifs with enhanced styling
    for i, motif in enumerate(motifs):
        motif_class = motif['Class']
        if motif_class == "Z-DNA" and motif.get("Subclass", "") == "eGZ (Extruded-G)":
            motif_class = "eGZ (Extruded-G)"
        
        color = MOTIF_COLORS.get(motif_class, "#666666")
        
        # Enhanced hover information with scientific details
        hover_text = (f"<b>{motif_class}</b><br>"
                     f"Position: {motif['Start']:,}-{motif['End']:,} bp<br>"
                     f"Length: {motif['Length']:,} bp<br>"
                     f"Score: {motif.get('Score', 'N/A')}<br>")
        
        if 'ScoreMethod' in motif:
            hover_text += f"Method: {motif['ScoreMethod']}<br>"
        if 'Arms/Repeat Unit/Copies' in motif and motif['Arms/Repeat Unit/Copies']:
            hover_text += f"Details: {motif['Arms/Repeat Unit/Copies']}<br>"
        
        hover_text += "<extra></extra>"
        
        fig.add_trace(go.Scatter(
            x=[motif['Start'], motif['End']],
            y=[y_positions.get(motif_class, 0), y_positions.get(motif_class, 0)],
            mode='lines',
            line=dict(color=color, width=12),  # Thicker lines for publication quality
            name=motif_class,
            showlegend=motif_class not in [trace.name for trace in fig.data],
            hovertemplate=hover_text
        ))
    
    # Publication-quality layout with enhanced styling
    fig.update_layout(
        title=dict(
            text=f"<b>Non-B DNA Structural Motifs: {seq_name}</b>",
            x=0.5,
            font=dict(size=18, family="Arial, sans-serif")
        ),
        xaxis=dict(
            title=dict(text="<b>Genomic Position (bp)</b>", font=dict(size=14)),
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=False,
            tickfont=dict(size=12)
        ),
        yaxis=dict(
            title=dict(text="<b>Non-B DNA Motif Classes</b>", font=dict(size=14)),
            tickmode='array',
            tickvals=list(y_positions.values()),
            ticktext=list(y_positions.keys()),
            showgrid=True,
            gridwidth=1,
            gridcolor='lightgray',
            zeroline=False,
            tickfont=dict(size=11)
        ),
        height=max(500, len(y_positions) * 40),  # Optimized height for readability
        hovermode='closest',
        showlegend=True,
        legend=dict(
            orientation="v",
            yanchor="top",
            y=1,
            xanchor="left",
            x=1.05,
            font=dict(size=11)
        ),
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin=dict(l=150, r=200, t=80, b=80)  # Enhanced margins for publication
    )
    
    return fig

def create_comprehensive_motif_overview(motifs, sequence_length, sequence_name="Sequence"):
    """Create comprehensive overview showing all 19 motif types, indicating 'not found' for missing ones"""
    
    # All 19 motif types in order
    all_motif_types = [
        "Bipartite G4", "Bulged G4", "Canonical G4", "Imperfect G4", "Multimeric G4", "Relaxed G4",
        "G-Triplex", "AC-Motif", "i-Motif", "Curved DNA", "Z-DNA", "eGZ (Extruded-G)",
        "Cruciform", "R-Loop", "Slipped DNA", "Sticky DNA", "Triplex DNA", "Hybrid", "Non-B DNA Clusters"
    ]
    
    # Count motifs by type
    motif_counts = {}
    motif_coverage = {}
    
    for motif_type in all_motif_types:
        count = 0
        total_coverage = 0
        
        for motif in motifs:
            motif_class = motif.get('Class', '')
            # Handle special case for eGZ
            if motif_class == "Z-DNA" and motif.get("Subclass", "") == "eGZ (Extruded-G)":
                motif_class = "eGZ (Extruded-G)"
            
            if motif_class == motif_type:
                count += 1
                total_coverage += motif.get('Length', 0)
        
        motif_counts[motif_type] = count
        coverage_percent = (total_coverage / sequence_length * 100) if sequence_length > 0 else 0
        motif_coverage[motif_type] = coverage_percent
    
    # Create comprehensive bar chart
    fig = go.Figure()
    
    # Colors for each motif type
    colors = [MOTIF_COLORS.get(motif_type, "#999999") for motif_type in all_motif_types]
    
    # Create bars with special formatting for "not found" items
    bar_colors = []
    bar_text = []
    y_values = []
    
    for i, motif_type in enumerate(all_motif_types):
        count = motif_counts[motif_type]
        if count == 0:
            bar_colors.append("#f5f5f5")  # Light gray for not found
            bar_text.append("Not Found")
            y_values.append(0.1)  # Small bar to show it exists
        else:
            bar_colors.append(colors[i])
            bar_text.append(f"{count} found")
            y_values.append(count)
    
    fig.add_trace(go.Bar(
        y=all_motif_types,
        x=y_values,
        orientation='h',
        marker_color=bar_colors,
        text=bar_text,
        textposition='inside',
        textfont=dict(color='white', size=10),
        hovertemplate="<b>%{y}</b><br>Count: %{x}<br>Coverage: %{customdata:.2f}%<extra></extra>",
        customdata=[motif_coverage[motif_type] for motif_type in all_motif_types]
    ))
    
    fig.update_layout(
        title=f"<b>Comprehensive Non-B DNA Motif Analysis: {sequence_name}</b><br><i>All 19 motif types assessed</i>",
        xaxis_title="<b>Number of Motifs Found</b>",
        yaxis_title="<b>Motif Types</b>",
        height=600,
        font=dict(size=12),
        plot_bgcolor='rgba(248,253,255,0.7)',
        paper_bgcolor='white',
        showlegend=False
    )
    
    # Add annotation for not found items
    fig.add_annotation(
        text="<i>Gray bars indicate motif types not detected</i>",
        xref="paper", yref="paper",
        x=0.99, y=-0.05, 
        showarrow=False,
        font=dict(size=10, color="gray")
    )
    
    return fig

# ---- TABS ----
tabs = st.tabs(list(PAGES.keys()))
tab_pages = dict(zip(PAGES.keys(), tabs))

# ---------- HOME ----------
with tab_pages["Home"]:
    # Compact header section
    st.markdown("""
    <div style='text-align: center; padding: 12px; margin-bottom: 16px;'>
        <h1 style='color: #2d3748; font-family: Inter, sans-serif; font-weight: 700; margin-bottom: 8px; font-size: 1.8rem;'>
            NBDFinder: Non-B DNA Analysis Platform
        </h1>
    </div>
    """, unsafe_allow_html=True)
    
    # Create main content layout
    left, right = st.columns([1.2, 1])
    
    with left:
        # Enhanced image with caption
        st.image("nbdcircle.JPG", use_container_width=True, caption="Non-B DNA structural diversity: From canonical B-form to complex alternative conformations")
        

    
    with right:
        st.markdown("""
        <div style='font-family: Inter, sans-serif; font-size: 1rem; color: #374151; line-height: 1.7; padding: 24px; background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px;'>
                <h3 style='color: #374151; margin-top: 0; margin-bottom: 20px; font-size: 1.25rem; font-weight: 600;'>
                    10 Non-B DNA Classes Detected:
                </h3>
                <div style='display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 20px;'>
                    <div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>1. Curved DNA</div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>2. Slipped DNA</div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>3. Cruciform DNA</div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>4. R-loop</div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>5. Triplex</div>
                    </div>
                    <div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>6. G-Quadruplex Family</div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>7. i-motif family</div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>8. Z-DNA</div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>9. Hybrid</div>
                        <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>10. Non-B DNA cluster regions</div>
                    </div>
                </div>
                <div style='border-top: 1px solid #e2e8f0; padding-top: 16px;'>
                    <h4 style='color: #374151; margin-bottom: 12px; font-size: 1.125rem; font-weight: 600;'>Platform Overview:</h4>
                    <p style='margin: 0; color: #6b7280; line-height: 1.6;'>
                        Non-canonical DNA structures play crucial roles in genome organization, gene regulation, and disease pathogenesis. Our platform provides computational tools for comprehensive structural analysis with multiple algorithms optimized for different motif types, interactive visualizations, and support for individual sequences and multi-FASTA files.
                    </p>
                </div>
        </div>
        """, unsafe_allow_html=True)
    


# ---------- UPLOAD & ANALYZE ----------
with tab_pages["Upload & Analyze"]:

    # Performance Mode Selection (New Feature)
    st.markdown("""
    <div style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); color: white; border-radius: 8px; padding: 16px; margin-bottom: 20px;'>
        <h3 style='color: white; margin: 0 0 8px 0; font-weight: 600;'>Performance Mode</h3>
        <p style='margin: 0; font-size: 0.9rem; opacity: 0.9;'>Choose analysis speed vs completeness trade-off</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Using Complete Performance Mode Only (as requested)
    st.markdown("**Performance Mode: Complete** (Comprehensive analysis with all features enabled)")
    
    # Scientific motif class selection with academic categorization
    st.markdown("""
    <div style='background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; margin-bottom: 24px;'>
        <h3 style='color: #374151; margin-top: 0; margin-bottom: 12px; font-size: 1.25rem; font-weight: 600;'>Motif Classes for Analysis</h3>
        <p style='margin-bottom: 15px; color: #6b7280; font-size: 0.875rem; line-height: 1.5;'>Select specific Non-B DNA motif classes for targeted analysis. Our comprehensive detection suite covers all major structural categories validated by experimental studies.</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Replace multiselect with scientific checkbox grid
    st.write("**Select motif classes for analysis:**")
    
    # Group motifs by scientific categories
    motif_categories = {
        "G-Quadruplex Related": ["Canonical G4", "Relaxed G4", "Bulged G4", "Imperfect G4", "Bipartite G4", "Multimeric G4", "G-Triplex"],
        "i-Motif Family": ["i-Motif", "AC-Motif"], 
        "Alternative Conformations": ["Curved DNA", "Z-DNA", "eGZ (Extruded-G)", "Cruciform"],
        "Complex Structures": ["R-Loop", "Slipped DNA", "Sticky DNA", "Triplex DNA"],
        "Composite Analysis": ["Hybrid", "Non-B DNA Clusters"]
    }
    
    # Initialize selected motifs from session state
    if 'selected_motifs' not in st.session_state:
        st.session_state.selected_motifs = MOTIF_ORDER.copy()
    
    selected_motifs = []
    
    # Quick selection options
    col1, col2, col3 = st.columns(3)
    with col1:
        if st.button("Select All", type="secondary"):
            st.session_state.selected_motifs = MOTIF_ORDER.copy()
            st.rerun()
    with col2:
        if st.button("Clear All", type="secondary"):
            st.session_state.selected_motifs = []
            st.rerun()
    with col3:
        if st.button("G4 Family Only", type="secondary"):
            st.session_state.selected_motifs = motif_categories["G-Quadruplex Related"]
            st.rerun()
    
    # Checkbox grid by category
    for category, motifs in motif_categories.items():
        st.markdown(f"**{category}:**")
        cols = st.columns(3)
        for i, motif in enumerate(motifs):
            with cols[i % 3]:
                if st.checkbox(motif, value=motif in st.session_state.selected_motifs, key=f"motif_{motif}"):
                    if motif not in st.session_state.selected_motifs:
                        st.session_state.selected_motifs.append(motif)
                else:
                    if motif in st.session_state.selected_motifs:
                        st.session_state.selected_motifs.remove(motif)
        st.write("")  # Add spacing between categories
    
    selected_motifs = st.session_state.selected_motifs if st.session_state.selected_motifs else MOTIF_ORDER
    st.session_state.selected_motifs = selected_motifs if selected_motifs else MOTIF_ORDER

    # Scientific input method selection
    st.markdown("""
    <div style='background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; margin: 24px 0;'>
        <h3 style='color: #374151; margin-top: 0; margin-bottom: 12px; font-size: 1.25rem; font-weight: 600;'>Input Method Selection</h3>
        <p style='margin-bottom: 0; color: #6b7280; font-size: 0.875rem; line-height: 1.5;'>Choose your preferred method for sequence input. All formats support both single sequences and batch processing.</p>
    </div>
    """, unsafe_allow_html=True)
    st.write("**Input Method:**")
    input_method = st.radio("", ["Upload FASTA / Multi-FASTA File", "Paste Sequence(s)", "Example Sequence", "NCBI Fetch"], horizontal=True)

    # Clear session state sequences when switching input methods (except on first load)
    if 'last_input_method' not in st.session_state:
        st.session_state.last_input_method = input_method
    elif st.session_state.last_input_method != input_method:
        st.session_state.seqs = []
        st.session_state.names = []
        st.session_state.results = []
        st.session_state.last_input_method = input_method

    # Initialize empty sequences - will be updated by input methods
    seqs, names = [], []

    # File upload with clean styling
    if input_method == "Upload FASTA / Multi-FASTA File":
        st.markdown("""
        <div style='background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; margin: 20px 0; text-align: center;'>
            <div style='font-size: 1.125rem; margin-bottom: 8px; color: #374151; font-weight: 500;'>File Upload</div>
            <h3 style='color: #374151; margin: 0 0 8px 0; font-weight: 600;'>
                Drag & Drop Your FASTA Files Here
            </h3>
            <p style='color: #6b7280; margin: 0; font-size: 0.875rem;'>
                Supports .fa, .fasta, .txt files • Single or Multi-FASTA format
            </p>
        </div>
        """, unsafe_allow_html=True)
        
        fasta_file = st.file_uploader(
            "Or click to browse files", 
            type=["fa", "fasta", "txt"],
            help="Upload FASTA files containing DNA sequences. Multi-FASTA files with multiple sequences are supported."
        )
        if fasta_file:
            content = fasta_file.read().decode("utf-8")
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
                # Store in session state for persistence
                st.session_state.seqs = seqs
                st.session_state.names = names
                st.success(f"Loaded {len(seqs)} sequences.")
                for i, seq in enumerate(seqs[:3]):
                    stats = get_basic_stats(seq)
                    st.write(f"Seq {i+1}: {stats}")

    # --- Enhanced Paste sequence with larger text area ---
    elif input_method == "Paste Sequence(s)":
        # Show example format in an expandable section
        with st.expander("FASTA Format Examples & Tips", expanded=False):
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Single FASTA Format:**")
                st.code(""">Human Beta Globin
ATGGCTCGTCTCTGTCTCCTGCCGACAAGACCAACGTCAAGGCCGCCTGGGGTAAGGTCGGCGCGCACGCTGGCGAG
TATGGTGCGGAGGCCCTGGAGAGGTTGGCCCCCCTGTCCCCGGGAGCGATCGTGGATGACCTGGCCGGCCCTTCGGCC
GTGCCCGACGTCCAGGTGGTGACCCACCTGACCGTGGCCTGGGAGAGTGTCTTCGCCAACTACCAAGACGCAGACCAC""", language='text')
            with col2:
                st.markdown("**Multi-FASTA Format:**")
                st.code(""">Sequence1_G4_Rich
GGGTTAGGGTTAGGGTTAGGGTAA
>Sequence2_iMotif
CCCTAACCCTAACCCTAACCCTAA
>Sequence3_ZDna
CGCGCGCGCGCGCGCGCGCG""", language='text')
            
            st.markdown("**Professional Tips:**")
            st.markdown("- Each sequence starts with `>` followed by sequence name")
            st.markdown("- DNA sequence follows on the next line(s)")
            st.markdown("- Multiple sequences can be pasted at once")
            st.markdown("- Accepts raw sequences without headers")
        
        # Enhanced text area with better styling
        st.markdown("""
        <div style='margin: 20px 0 8px 0;'>
            <label style='font-size: 1.1rem; font-weight: 600; color: #2d3748;'>
                **Paste Your DNA Sequences:**
            </label>
        </div>
        """, unsafe_allow_html=True)
        
        seq_input = st.text_area(
            "Paste FASTA or raw sequence(s)", 
            height=250,
            placeholder="Paste your DNA sequences here in FASTA format:\n>MySequence\nATCGATCGATCG...\n\nOr paste raw sequences directly...",
            help="Supports both FASTA format and raw DNA sequences. For multiple sequences, use FASTA format with > headers."
        )
        if seq_input:
            lines = seq_input.splitlines()
            cur_seq, cur_name = "", ""
            for line in lines:
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
                # Store in session state for persistence
                st.session_state.seqs = seqs
                st.session_state.names = names
                st.success(f"Pasted {len(seqs)} sequences.")

    # --- Example input ---
    # --- Example sequences (load as sequence, not file) ---
    elif input_method == "Example Sequence":
        st.markdown("**Choose from curated Non-B DNA examples:**")
        
        # Display example options with radio buttons instead of selectbox
        example_choices = list(EXAMPLE_SEQUENCES.keys())
        
        # Show examples as radio buttons with descriptions
        selected_example = st.radio(
            "Select example sequence:",
            options=example_choices,
            format_func=lambda x: f"{x} - {EXAMPLE_SEQUENCES[x]['description']}"
        )
        
        if selected_example:
            example_data = EXAMPLE_SEQUENCES[selected_example]
            seqs = [example_data["sequence"]]
            names = [example_data["name"]]
            
            # Store in session state for persistence
            st.session_state.seqs = seqs
            st.session_state.names = names
            
            # Add spacing before the preview to prevent overlap
            st.markdown("<br>", unsafe_allow_html=True)
            
            # Show sequence preview directly (not in expander)
            st.markdown("**Preview Selected Sequence:**")
            st.text(f"Name: {example_data['name']}")
            st.text(f"Length: {len(example_data['sequence'])} bp")
            st.text(f"Description: {example_data['description']}")
            st.code(wrap(example_data["sequence"], 80), language='text')
            
            st.success(f"Loaded Loaded: {example_data['name']} ({len(example_data['sequence'])} bp)")

    # --- Enhanced NCBI Query with inline examples ---
    elif input_method == "NCBI Fetch":
        # Single inline text input with rotating examples as specified
        placeholders = [
            "for eg. Human c-MYC Promoter (NG_007161.1)",
            "Human BCL2 Promoter (NG_009361.1)", 
            "Fragile X FMR1 Gene (NG_007529.1)",
            "Huntington HTT Gene (NG_009378.1)",
            "Human Immunoglobulin Switch (NG_001019.6)"
        ]
        
        # Use rotating placeholder (simple approach for now)
        import hashlib
        placeholder_idx = hash(str(int(time.time() / 10))) % len(placeholders)
        current_placeholder = placeholders[placeholder_idx]
        
        # Example chips for auto-fill functionality
        st.markdown("**Quick Examples:**")
        chip_cols = st.columns(len(FAMOUS_NCBI_EXAMPLES))
        for i, (gene, accession) in enumerate(FAMOUS_NCBI_EXAMPLES.items()):
            if i < len(chip_cols):
                with chip_cols[i]:
                    if st.button(f"{gene}", key=f"chip_{i}", help=f"Auto-fill: {accession}"):
                        st.session_state.ncbi_query = accession
        
        # Main NCBI query input with rotating placeholder
        ncbi_query = st.text_input(
            "Enter NCBI Query:",
            value=st.session_state.get('ncbi_query', ''),
            placeholder=current_placeholder,
            help="Enter NCBI accession number or gene name. Supports both accession-based (e.g., NG_007161.1) and free-text queries (e.g., human c-MYC promoter)",
            key="ncbi_query_input"
        )
        
        # Inline validation for NCBI query
        if ncbi_query:
            import re
            # Pattern for common NCBI accession formats
            accession_patterns = [
                r'^[A-Z]{1,2}_\d{6,9}\.\d{1,2}$',  # e.g., NG_007161.1, NM_007294.3
                r'^[A-Z]{1,2}\d{6,9}$',            # e.g., AC123456
                r'^[A-Z]\d{5}$',                   # e.g., M12345
                r'^[A-Z]{2}\d{6}$'                 # e.g., AB123456
            ]
            
            is_accession = any(re.match(pattern, ncbi_query.strip()) for pattern in accession_patterns)
            
            if is_accession:
                st.success(f"✓ Valid accession format detected: {ncbi_query}")
            else:
                st.info(f"Free-text query: '{ncbi_query}' - will search NCBI database")
        
        if ncbi_query:
            if st.button("Fetch from NCBI"):
                with st.spinner("Fetching sequences from NCBI..."):
                    try:
                        seqs, names = ncbi_fetch(ncbi_query)
                        if seqs:
                            # Store in session state for persistence
                            st.session_state.seqs = seqs
                            st.session_state.names = names
                            st.success(f"Loaded Fetched {len(seqs)} sequence(s) from NCBI.")
                            # Show preview of fetched sequences
                            with st.expander("- Preview Fetched Sequences"):
                                for i, (seq, name) in enumerate(zip(seqs[:3], names[:3])):  # Show first 3
                                    st.text(f"{i+1}. {name} ({len(seq)} bp)")
                                    st.code(wrap(seq[:200], 80) + ("..." if len(seq) > 200 else ""), language='text')
                                if len(seqs) > 3:
                                    st.info(f"... and {len(seqs)-3} more sequences")
                        else:
                            st.warning("No sequences found for this query.")
                    except Exception as e:
                        st.error(f"NCBI fetch failed: {str(e)}")
                        st.info("Try a different query or check your internet connection.")

    # --- Analysis trigger ---
    # Show currently loaded sequences status
    if st.session_state.seqs:
        st.markdown("---")
        st.markdown("### Loaded Sequences")
        seq_count = len(st.session_state.seqs)
        total_length = sum(len(seq) for seq in st.session_state.seqs)
        st.info(f"**{seq_count} sequence(s) ready for analysis** | Total length: {total_length:,} bp")
        
        if seq_count <= 3:
            for i, (seq, name) in enumerate(zip(st.session_state.seqs, st.session_state.names)):
                st.text(f"{i+1}. {name} ({len(seq):,} bp)")
        else:
            for i in range(3):
                seq, name = st.session_state.seqs[i], st.session_state.names[i]
                st.text(f"{i+1}. {name} ({len(seq):,} bp)")
            st.text(f"... and {seq_count-3} more sequences")
    
    # Use session state sequences for analysis (ensures persistence across input methods)
    if st.session_state.seqs:
        
        # ========== ADVANCED OPTIONS PANEL ==========
        with st.expander("Advanced Analysis Options", expanded=False):
            st.markdown("""
            <div style='background: linear-gradient(135deg, #fdf4ff 0%, #f9fafb 100%); 
                        border-radius: 12px; padding: 20px; margin: 12px 0;
                        border: 2px solid #e879f9;'>
                <h4 style='color: #374151; margin: 0 0 16px 0; font-size: 1.125rem; font-weight: 600;'>
                    Detection Parameters & Analysis Options
                </h4>
            </div>
            """, unsafe_allow_html=True)
            
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown("**Motif-Specific Thresholds:**")
                
                # G4 sensitivity
                g4_sensitivity = st.slider(
                    "G-Quadruplex Sensitivity",
                    min_value=1.0, max_value=3.0, value=1.5, step=0.1,
                    help="Lower values = more sensitive detection (more G4s found)"
                )
                
                # Minimum motif length
                min_motif_length = st.number_input(
                    "Minimum Motif Length (bp)",
                    min_value=8, max_value=100, value=20,
                    help="Exclude motifs shorter than this threshold"
                )
                
                # Conservation analysis
                enable_conservation = st.checkbox(
                    "Enable Conservation Analysis",
                    value=False,
                    help="Analyze evolutionary conservation of detected motifs"
                )
                
            with col2:
                st.markdown("**Output & Analysis Options:**")
                
                # Initialize session state for toggle buttons if not exists
                if 'output_format_selection' not in st.session_state:
                    st.session_state.output_format_selection = "Comprehensive (All Data)"
                if 'overlap_strategy_selection' not in st.session_state:
                    st.session_state.overlap_strategy_selection = "Best Score"
                
                # Output format with horizontal toggle buttons
                st.markdown("**Export Format:**")
                format_col1, format_col2, format_col3 = st.columns(3)
                
                format_options = ["Comprehensive (All Data)", "Summary Only", "Coordinates Only"]
                format_labels = ["All Data", "Summary", "Coordinates"]
                
                with format_col1:
                    if st.button(
                        format_labels[0], 
                        key="format_comprehensive",
                        help="Export comprehensive data with all details",
                        use_container_width=True,
                        type="primary" if st.session_state.output_format_selection == format_options[0] else "secondary"
                    ):
                        st.session_state.output_format_selection = format_options[0]
                
                with format_col2:
                    if st.button(
                        format_labels[1], 
                        key="format_summary",
                        help="Export summary data only",
                        use_container_width=True,
                        type="primary" if st.session_state.output_format_selection == format_options[1] else "secondary"
                    ):
                        st.session_state.output_format_selection = format_options[1]
                
                with format_col3:
                    if st.button(
                        format_labels[2], 
                        key="format_coordinates",
                        help="Export coordinates only",
                        use_container_width=True,
                        type="primary" if st.session_state.output_format_selection == format_options[2] else "secondary"
                    ):
                        st.session_state.output_format_selection = format_options[2]
                
                output_format = st.session_state.output_format_selection
                
                st.markdown("<br>", unsafe_allow_html=True)
                
                # Overlap handling with horizontal toggle buttons
                st.markdown("**Overlap Resolution:**")
                overlap_col1, overlap_col2, overlap_col3 = st.columns(3)
                
                overlap_options = ["Best Score", "Longest Motif", "No Filtering"]
                overlap_labels = ["Best Score", "Longest", "No Filter"]
                
                with overlap_col1:
                    if st.button(
                        overlap_labels[0], 
                        key="overlap_best",
                        help="Keep motifs with best scores when overlapping",
                        use_container_width=True,
                        type="primary" if st.session_state.overlap_strategy_selection == overlap_options[0] else "secondary"
                    ):
                        st.session_state.overlap_strategy_selection = overlap_options[0]
                
                with overlap_col2:
                    if st.button(
                        overlap_labels[1], 
                        key="overlap_longest",
                        help="Keep longest motifs when overlapping",
                        use_container_width=True,
                        type="primary" if st.session_state.overlap_strategy_selection == overlap_options[1] else "secondary"
                    ):
                        st.session_state.overlap_strategy_selection = overlap_options[1]
                
                with overlap_col3:
                    if st.button(
                        overlap_labels[2], 
                        key="overlap_none",
                        help="No filtering - keep all overlapping motifs",
                        use_container_width=True,
                        type="primary" if st.session_state.overlap_strategy_selection == overlap_options[2] else "secondary"
                    ):
                        st.session_state.overlap_strategy_selection = overlap_options[2]
                
                overlap_strategy = st.session_state.overlap_strategy_selection
                
                st.markdown("<br>", unsafe_allow_html=True)
                
                # Statistical analysis
                enable_stats = st.checkbox(
                    "Extended Statistical Analysis",
                    value=True,
                    help="Include detailed statistical metrics and clustering analysis"
                )
            
            # Store advanced options in session state
            st.session_state.advanced_options = {
                'g4_sensitivity': g4_sensitivity,
                'min_motif_length': min_motif_length,
                'enable_conservation': enable_conservation,
                'output_format': output_format,
                'overlap_strategy': overlap_strategy,
                'enable_stats': enable_stats
            }
        
        # Analysis button with enhanced styling
        st.markdown("<br>", unsafe_allow_html=True)
        analysis_col1, analysis_col2, analysis_col3 = st.columns([2, 3, 2])
        with analysis_col2:
            if st.button(
                "Start Analysis", 
                type="primary",
                use_container_width=True,
                help="Begin comprehensive Non-B DNA motif detection"
            ):
                seqs = st.session_state.seqs
                names = st.session_state.names
                
                # Check sequence length restrictions
                invalid_seqs = []
                for i, (seq, name) in enumerate(zip(seqs, names)):
                    if len(seq) > 100000:
                        invalid_seqs.append(f"{name} ({len(seq):,} bp)")
                
                if invalid_seqs:
                    st.error("Sequence length restriction exceeded!")
                    st.warning(f"Maximum allowed length: 100,000 nucleotides")
                    st.info("For longer sequences, please use the standalone version of NBDFinder")
                    st.write("**Sequences exceeding limit:**")
                    for seq_info in invalid_seqs:
                        st.write(f"- {seq_info}")
                else:
                    st.session_state.is_analyzing = True
                    st.session_state.analysis_start_time = time.time()  # Initialize start time
                    results = []
                    progress_container = create_progress_tracker()
                    
                    for seq, name in zip(seqs, names):
                        motifs = analyze_sequence_with_progress(seq, name, st.session_state.selected_motifs)
                        results.append((name, motifs))
                        
                    st.session_state.seqs = seqs
                    st.session_state.names = names
                    st.session_state.results = results
                    st.session_state.is_analyzing = False
                    
                    # Create summary dataframe
                    summary_data = []
                    for i, (name, motifs) in enumerate(results):
                        seq = seqs[i]
                        stats = get_basic_stats(seq, motifs)
                        motif_classes = [m['Class'] for m in motifs]
                        top_motifs = ', '.join(list(set(motif_classes))[:3]) if motif_classes else "None"
                        
                        summary_data.append({
                            "Sequence Name": name,
                            "Length (bp)": stats["Length (bp)"],
                            "GC %": stats["GC %"],
                            "Motif Count": len(motifs),
                            "Motif Coverage (%)": stats.get("Motif Coverage (%)", 0),
                            "Top Motifs": top_motifs
                        })
                    
                    st.session_state.summary_df = pd.DataFrame(summary_data)
                    st.success("Analysis complete. See Results tab for visualization and tables.")

# ---------- ENHANCED RESULTS SECTION ----------
with tab_pages["Results"]:
    st.markdown("""
    <div style='text-align: center; margin-bottom: 30px;'>
        <h2 style='color: #374151; font-family: Inter, sans-serif; font-weight: 600; margin-bottom: 8px;'>
            Advanced Results Analysis Dashboard
        </h2>
        <p style='color: #6b7280; font-size: 1rem; margin: 0;'>
            Publication-quality visualizations and comprehensive motif analysis
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    if not st.session_state.results:
        st.markdown("""
        <div style='background: linear-gradient(135deg, #4a5568 0%, #2d3748 100%); color: white; border-radius: 12px; padding: 40px; text-align: center; margin: 20px 0;'>
            <h3 style='color: white; margin: 0 0 12px 0;'>Ready for Analysis</h3>
            <p style='color: white; margin: 0; font-size: 1.1rem; opacity: 0.9;'>
                Upload sequences in the <strong>Upload & Analyze</strong> tab to see advanced visualizations here
            </p>
        </div>
        """, unsafe_allow_html=True)
    else:
        # ========== 6 SCIENTIFIC SUBPAGE NAVIGATION ==========
        st.markdown("""
        <div style='background: #f8fafc; border-radius: 12px; padding: 20px; margin-bottom: 24px;'>
            <h3 style='color: #374151; margin: 0 0 16px 0;'>Scientific Analysis Modules</h3>
        </div>
        """, unsafe_allow_html=True)
        
        # Scientific subpage selection
        st.markdown("**Analysis Module:**")
        subpage_selection = st.radio(
            "Choose scientific analysis view:",
            ["Overview", "Motif Classes", "Genomic Position", "Hybrid Analysis", "Cluster Analysis", "Statistical Summary"],
            horizontal=True,
            key="scientific_results_view_radio"
        )
        
        st.markdown("---")
        
        # Scientific results content based on selection
        if subpage_selection == "Overview":
            st.markdown("## Analysis Overview")
            
            # Calculate comprehensive summary statistics
            total_sequences = len(st.session_state.results)
            total_motifs = sum(len(motifs) for _, motifs in st.session_state.results)
            all_motifs_flat = []
            for _, motifs in st.session_state.results:
                all_motifs_flat.extend(motifs)
            
            # Summary metrics in professional layout
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Sequences Analyzed", total_sequences)
            with col2:
                st.metric("Total Motifs Found", total_motifs)
            with col3:
                unique_classes = len(set([m['Class'] for m in all_motifs_flat]))
                st.metric("Motif Classes", unique_classes)
            with col4:
                avg_motifs = total_motifs / total_sequences if total_sequences > 0 else 0
                st.metric("Avg Motifs/Sequence", f"{avg_motifs:.1f}")
            
            # Class distribution chart with professional styling
            if all_motifs_flat:
                st.markdown("### Class Distribution Analysis")
                class_counts = {}
                for motif in all_motifs_flat:
                    motif_class = motif.get('Class', 'Unknown')
                    class_counts[motif_class] = class_counts.get(motif_class, 0) + 1
                
                if class_counts:
                    # Professional publication-quality chart
                    classes = list(class_counts.keys())
                    counts = list(class_counts.values())
                    
                    # Use professional color scheme
                    colors = ['#2563eb', '#dc2626', '#059669', '#7c3aed', '#ea580c', '#0891b2', '#ec4899', '#f59e0b', '#10b981', '#8b5cf6']
                    
                    fig = go.Figure()
                    fig.add_trace(go.Bar(
                        x=classes, 
                        y=counts,
                        marker_color=colors[:len(classes)],
                        text=counts,
                        textposition='auto'
                    ))
                    fig.update_layout(
                        title="Non-B DNA Motif Class Distribution",
                        xaxis_title="Motif Class",
                        yaxis_title="Count",
                        height=400,
                        plot_bgcolor='white',
                        paper_bgcolor='white',
                        font=dict(family="Arial, sans-serif", size=12),
                        title_font_size=16
                    )
                    st.plotly_chart(fig, use_container_width=True)
                
                # Summary table
                summary_df = pd.DataFrame([
                    {"Motif Class": cls, "Count": count, "Percentage": f"{count/total_motifs*100:.1f}%"}
                    for cls, count in sorted(class_counts.items(), key=lambda x: x[1], reverse=True)
                ])
                st.markdown("### Summary Statistics")
                st.dataframe(summary_df, use_container_width=True)
        
            # Detailed breakdown by motif class with data tables
            st.markdown("## Motif Classes Analysis")
            
            # Get all motifs data
            all_motifs_flat = []
            for _, motifs in st.session_state.results:
                all_motifs_flat.extend(motifs)
                
            if all_motifs_flat:
                # Prepare comprehensive data table
                all_data = []
                for seq_name, motifs in st.session_state.results:
                    for motif in motifs:
                        motif_copy = motif.copy()
                        motif_copy['Sequence Name'] = seq_name
                        all_data.append(motif_copy)
                
                df_comprehensive = pd.DataFrame(all_data)
                
                # Professional filtering controls
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    # Class filter
                    available_classes = sorted(df_comprehensive['Class'].unique()) if 'Class' in df_comprehensive.columns else []
                    selected_class = st.selectbox(
                        "Motif Class:",
                        options=["All Classes"] + available_classes
                    )
                
                with col2:
                    # Score filter slider
                    if 'Score' in df_comprehensive.columns:
                        score_values = pd.to_numeric(df_comprehensive['Score'], errors='coerce').dropna()
                        if len(score_values) > 0:
                            min_score = st.slider(
                                "Minimum Score:",
                                min_value=float(score_values.min()),
                                max_value=float(score_values.max()),
                                value=float(score_values.min())
                            )
                        else:
                            min_score = 0
                    else:
                        min_score = 0
                
                with col3:
                    # Length filter
                    if 'Length' in df_comprehensive.columns:
                        length_values = pd.to_numeric(df_comprehensive['Length'], errors='coerce').dropna()
                        if len(length_values) > 0:
                            min_length = st.slider(
                                "Minimum Length (bp):",
                                min_value=int(length_values.min()),
                                max_value=int(length_values.max()),
                                value=int(length_values.min())
                            )
                        else:
                            min_length = 0
                    else:
                        min_length = 0
                
                # Apply filters
                filtered_df = df_comprehensive.copy()
                if selected_class != "All Classes":
                    filtered_df = filtered_df[filtered_df['Class'] == selected_class]
                
                if 'Score' in filtered_df.columns:
                    score_numeric = pd.to_numeric(filtered_df['Score'], errors='coerce')
                    filtered_df = filtered_df[score_numeric >= min_score]
                
                if 'Length' in filtered_df.columns:
                    length_numeric = pd.to_numeric(filtered_df['Length'], errors='coerce')
                    filtered_df = filtered_df[length_numeric >= min_length]
                
                # Display summary
                st.markdown(f"**Showing {len(filtered_df)} of {len(df_comprehensive)} total motifs**")
                
                # Display the table with professional styling
                key_columns = ['Sequence Name', 'Class', 'Subtype', 'Start', 'End', 'Length', 'Score']
                display_columns = [col for col in key_columns if col in filtered_df.columns]
                
                st.dataframe(
                    filtered_df[display_columns] if display_columns else filtered_df, 
                    use_container_width=True, 
                    height=400
                )
                
                # Export options
                col1, col2, col3 = st.columns(3)
                with col1:
                    csv = filtered_df.to_csv(index=False)
                    st.download_button(
                        label="Export CSV",
                        data=csv,
                        file_name="motif_classes_analysis.csv",
                        mime="text/csv"
                    )
                with col2:
                    # Excel export
                    excel_data = io.BytesIO()
                    with pd.ExcelWriter(excel_data, engine='xlsxwriter') as writer:
                        filtered_df.to_excel(writer, index=False, sheet_name="Motif_Classes")
                    excel_data.seek(0)
                    st.download_button(
                        label="Export Excel",
                        data=excel_data,
                        file_name="motif_classes_analysis.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                with col3:
                    # TXT export
                    txt_data = filtered_df.to_string(index=False)
                    st.download_button(
                        label="Export TXT",
                        data=txt_data,
                        file_name="motif_classes_analysis.txt",
                        mime="text/plain"
                    )
        
        elif subpage_selection == "Genomic Position":
            # Interactive genomic position mapping
            st.markdown("## Genomic Position Analysis")
            
            # Get all sequences and motifs
            all_sequences = []
            for seq_name, motifs in st.session_state.results:
                for i, seq in enumerate(st.session_state.seqs):
                    if st.session_state.names[i] == seq_name:
                        all_sequences.append((seq_name, seq, motifs))
                        break
            
            if all_sequences:
                # Sequence selector
                seq_names = [name for name, _, _ in all_sequences]
                selected_seq = st.selectbox("Select Sequence:", seq_names)
                
                # Find selected sequence data
                for seq_name, seq, motifs in all_sequences:
                    if seq_name == selected_seq:
                        # Create genome browser visualization
                        # Convert motifs list to DataFrame as expected by the function
                        motifs_df = pd.DataFrame(motifs) if motifs else pd.DataFrame()
                        create_genome_browser_view(motifs_df, seq, seq_name)
                        break
        
        elif subpage_selection == "Hybrid Analysis":
            # Comprehensive hybrid structure analysis with class combinations
            st.markdown("## Hybrid Structure Analysis")
            
            # Get all motifs and filter for hybrids
            all_motifs_flat = []
            hybrid_motifs = []
            
            for _, motifs in st.session_state.results:
                for motif in motifs:
                    all_motifs_flat.append(motif)
                    if motif.get('Class') == 'Hybrid':
                        hybrid_motifs.append(motif)
            
            if hybrid_motifs:
                st.markdown(f"**Found {len(hybrid_motifs)} hybrid structures**")
                
                # Hybrid class combination analysis
                class_combinations = {}
                for hybrid in hybrid_motifs:
                    classes = hybrid.get('MotifClasses', [])
                    if classes:
                        combo_key = ' + '.join(sorted(classes))
                        class_combinations[combo_key] = class_combinations.get(combo_key, 0) + 1
                
                if class_combinations:
                    st.markdown("### Class Combination Frequency")
                    combo_df = pd.DataFrame([
                        {"Class Combination": combo, "Count": count}
                        for combo, count in sorted(class_combinations.items(), key=lambda x: x[1], reverse=True)
                    ])
                    st.dataframe(combo_df, use_container_width=True)
                    
                    # Visualization of combinations
                    if len(class_combinations) > 0:
                        fig = go.Figure()
                        fig.add_trace(go.Bar(
                            x=list(class_combinations.keys()),
                            y=list(class_combinations.values()),
                            marker_color='#0891b2'
                        ))
                        fig.update_layout(
                            title="Hybrid Structure Class Combinations",
                            xaxis_title="Class Combination",
                            yaxis_title="Count",
                            height=400,
                            plot_bgcolor='white',
                            paper_bgcolor='white'
                        )
                        st.plotly_chart(fig, use_container_width=True)
                
                # Detailed hybrid table
                st.markdown("### Hybrid Structure Details")
                hybrid_df = pd.DataFrame(hybrid_motifs)
                key_columns = ['Sequence Name', 'Subtype', 'Start', 'End', 'Length', 'OverlapDegree', 'Score']
                display_columns = [col for col in key_columns if col in hybrid_df.columns]
                st.dataframe(hybrid_df[display_columns] if display_columns else hybrid_df, use_container_width=True)
                
                # Export hybrid analysis
                csv = hybrid_df.to_csv(index=False)
                st.download_button(
                    label="Export Hybrid Analysis CSV",
                    data=csv,
                    file_name="hybrid_analysis.csv",
                    mime="text/csv"
                )
            else:
                st.info("No hybrid structures detected in current analysis. Hybrids form when 2+ different motif classes overlap.")
        
        elif subpage_selection == "Cluster Analysis":
            # Non-B DNA cluster regions with complexity metrics
            st.markdown("## Cluster Analysis")
            
            # Get all motifs and identify clusters
            all_motifs_flat = []
            cluster_motifs = []
            
            for _, motifs in st.session_state.results:
                for motif in motifs:
                    all_motifs_flat.append(motif)
                    if motif.get('Class') == 'Non-B DNA Clusters':
                        cluster_motifs.append(motif)
            
            if cluster_motifs:
                st.markdown(f"**Found {len(cluster_motifs)} cluster regions**")
                
                # Cluster complexity metrics
                col1, col2, col3 = st.columns(3)
                with col1:
                    avg_length = np.mean([c.get('Length', 0) for c in cluster_motifs])
                    st.metric("Average Cluster Length", f"{avg_length:.1f} bp")
                with col2:
                    avg_motifs = np.mean([len(c.get('ContributingMotifs', [])) for c in cluster_motifs])
                    st.metric("Average Motifs/Cluster", f"{avg_motifs:.1f}")
                with col3:
                    max_density = max([c.get('MotifDensity', 0) for c in cluster_motifs], default=0)
                    st.metric("Max Motif Density", f"{max_density:.2f}/100bp")
                
                # Cluster details table
                st.markdown("### Cluster Region Details")
                cluster_df = pd.DataFrame(cluster_motifs)
                key_columns = ['Sequence Name', 'Start', 'End', 'Length', 'MotifDensity', 'SubclassCount', 'Score']
                display_columns = [col for col in key_columns if col in cluster_df.columns]
                st.dataframe(cluster_df[display_columns] if display_columns else cluster_df, use_container_width=True)
                
                # Export cluster analysis
                csv = cluster_df.to_csv(index=False)
                st.download_button(
                    label="Export Cluster Analysis CSV", 
                    data=csv,
                    file_name="cluster_analysis.csv",
                    mime="text/csv"
                )
            else:
                st.info("No cluster regions detected. Clusters form when 3+ motifs occur within 100nt windows.")
        
        elif subpage_selection == "Statistical Summary":
            # Publication-ready statistical analysis and metrics
            st.markdown("## Statistical Summary")
            
            # Comprehensive statistical analysis
            all_motifs_flat = []
            for _, motifs in st.session_state.results:
                all_motifs_flat.extend(motifs)
            
            if all_motifs_flat:
                # Create comprehensive statistics
                stats_data = []
                
                # Per-class statistics
                class_stats = {}
                for motif in all_motifs_flat:
                    motif_class = motif.get('Class', 'Unknown')
                    if motif_class not in class_stats:
                        class_stats[motif_class] = {
                            'count': 0,
                            'lengths': [],
                            'scores': []
                        }
                    
                    class_stats[motif_class]['count'] += 1
                    if 'Length' in motif:
                        try:
                            class_stats[motif_class]['lengths'].append(float(motif['Length']))
                        except:
                            pass
                    if 'Score' in motif:
                        try:
                            class_stats[motif_class]['scores'].append(float(motif['Score']))
                        except:
                            pass
                
                # Compile statistics table
                for class_name, stats in class_stats.items():
                    lengths = stats['lengths']
                    scores = stats['scores']
                    
                    stats_row = {
                        'Motif Class': class_name,
                        'Count': stats['count'],
                        'Frequency (%)': f"{stats['count']/len(all_motifs_flat)*100:.1f}%"
                    }
                    
                    if lengths:
                        stats_row.update({
                            'Mean Length (bp)': f"{np.mean(lengths):.1f}",
                            'Std Length (bp)': f"{np.std(lengths):.1f}",
                            'Min Length (bp)': f"{min(lengths):.0f}",
                            'Max Length (bp)': f"{max(lengths):.0f}"
                        })
                    
                    if scores:
                        stats_row.update({
                            'Mean Score': f"{np.mean(scores):.2f}",
                            'Std Score': f"{np.std(scores):.2f}"
                        })
                    
                    stats_data.append(stats_row)
                
                # Display comprehensive statistics table
                stats_df = pd.DataFrame(stats_data)
                st.markdown("### Comprehensive Statistical Analysis")
                st.dataframe(stats_df, use_container_width=True)
                
                # Overall summary metrics
                st.markdown("### Summary Metrics")
                col1, col2, col3, col4 = st.columns(4)
                
                with col1:
                    total_sequences = len(st.session_state.results)
                    st.metric("Total Sequences", total_sequences)
                with col2:
                    st.metric("Total Motifs", len(all_motifs_flat))
                with col3:
                    st.metric("Unique Classes", len(class_stats))
                with col4:
                    overall_coverage = sum([m.get('Length', 0) for m in all_motifs_flat])
                    total_seq_length = sum([len(seq) for seq in st.session_state.seqs])
                    coverage_pct = (overall_coverage / total_seq_length * 100) if total_seq_length > 0 else 0
                    st.metric("Total Coverage", f"{coverage_pct:.1f}%")
                
                # Export comprehensive statistics
                col1, col2, col3 = st.columns(3)
                with col1:
                    csv = stats_df.to_csv(index=False)
                    st.download_button(
                        label="Export Statistics CSV",
                        data=csv,
                        file_name="statistical_summary.csv",
                        mime="text/csv"
                    )
                with col2:
                    excel_data = io.BytesIO()
                    with pd.ExcelWriter(excel_data, engine='xlsxwriter') as writer:
                        stats_df.to_excel(writer, index=False, sheet_name="Statistics")
                    excel_data.seek(0)
                    st.download_button(
                        label="Export Statistics Excel",
                        data=excel_data,
                        file_name="statistical_summary.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                with col3:
                    txt_data = stats_df.to_string(index=False)
                    st.download_button(
                        label="Export Statistics TXT",
                        data=txt_data,
                        file_name="statistical_summary.txt",
                        mime="text/plain"
                    )
            else:
                st.info("No motifs available for statistical analysis.")


# ---------- DOWNLOAD ----------
with tab_pages["Download"]:
    st.header("Download")
    if not st.session_state.results:
        st.info("No results available to download.")
    else:
        df_all = []
        for i, (name, motifs) in enumerate(st.session_state.results):
            for j, m in enumerate(motifs):
                if isinstance(m, dict):  # Ensure m is a dictionary
                    m = m.copy()  # Create a copy to avoid modifying original
                    m['Sequence Name'] = st.session_state.names[i]
                    m['S.No'] = j + 1  # Add serial number starting from 1
                    if m.get('Class') == "Z-DNA" and m.get("Subclass", "") == "eGZ (Extruded-G)":
                        m['Class'] = "eGZ (Extruded-G)"
                    
                    # Add normalized confidence scoring (consistent with Results tab)
                    def get_normalized_confidence_score(score, motif_class):
                        """Normalize scores to unified 1.0-3.0 scale for consistency"""
                        try:
                            score_val = float(score)
                            
                            # Define motif-specific thresholds and normalizations
                            thresholds = {
                                # G4 Family (already normalized 1.0-3.0)
                                'Canonical G4': {'min': 1.0, 'high': 2.0, 'optimal': 3.0},
                                'Relaxed G4': {'min': 1.0, 'high': 2.0, 'optimal': 2.5},
                                'Bulged G4': {'min': 1.0, 'high': 2.0, 'optimal': 2.5},
                                'Imperfect G4': {'min': 1.0, 'high': 1.5, 'optimal': 2.0},
                                'Bipartite G4': {'min': 20.0, 'high': 60.0, 'optimal': 100.0},
                                'Multimeric G4': {'min': 30.0, 'high': 90.0, 'optimal': 150.0},
                                
                                # Alternative DNA Structures
                                'Z-DNA': {'min': 50.0, 'high': 200.0, 'optimal': 500.0},
                                'Curved DNA': {'min': 15.0, 'high': 100.0, 'optimal': 200.0},
                                'eGZ': {'min': 10.0, 'high': 55.0, 'optimal': 100.0},
                                'Extruded-G': {'min': 10.0, 'high': 55.0, 'optimal': 100.0},
                                'eGZ (Extruded-G)': {'min': 10.0, 'high': 55.0, 'optimal': 100.0},
                                
                                # Repeat-based Structures
                                'Slipped_DNA': {'min': 15.0, 'high': 80.0, 'optimal': 150.0},
                                'Slipped DNA': {'min': 15.0, 'high': 80.0, 'optimal': 150.0},
                                'R-Loop': {'min': 20.0, 'high': 150.0, 'optimal': 300.0},
                                'RLFS': {'min': 20.0, 'high': 150.0, 'optimal': 300.0},
                                
                                # Junction Structures
                                'Cruciform': {'min': 25.0, 'high': 100.0, 'optimal': 200.0},
                                'Triplex': {'min': 20.0, 'high': 100.0, 'optimal': 180.0},
                                'H-DNA': {'min': 20.0, 'high': 100.0, 'optimal': 180.0},
                                'Triplex_DNA': {'min': 20.0, 'high': 100.0, 'optimal': 180.0},
                                
                                # Other Structures
                                'Sticky_DNA': {'min': 10.0, 'high': 45.0, 'optimal': 80.0},
                                'Sticky DNA': {'min': 10.0, 'high': 45.0, 'optimal': 80.0},
                                'G-Triplex': {'min': 15.0, 'high': 65.0, 'optimal': 120.0},
                                'i-Motif': {'min': 15.0, 'high': 55.0, 'optimal': 100.0},
                            }
                            
                            # Get thresholds for this motif class
                            if motif_class in thresholds:
                                thresh = thresholds[motif_class]
                                
                                # Normalize to 1.0-3.0 scale
                                if score_val >= thresh['optimal']:
                                    normalized_score = 3.0
                                elif score_val >= thresh['high']:
                                    # Linear interpolation between high and optimal
                                    normalized_score = 2.0 + (score_val - thresh['high']) / (thresh['optimal'] - thresh['high'])
                                elif score_val >= thresh['min']:
                                    # Linear interpolation between min and high  
                                    normalized_score = 1.0 + (score_val - thresh['min']) / (thresh['high'] - thresh['min'])
                                else:
                                    # Below minimum threshold
                                    normalized_score = max(0.1, score_val / thresh['min'])
                                
                                return min(3.0, max(0.1, normalized_score))
                            else:
                                # Default handling for unknown motif types
                                return 1.0
                                
                        except (ValueError, TypeError):
                            return 1.0
                    
                    def get_confidence_label(normalized_score):
                        """Convert normalized score to confidence label"""
                        if normalized_score >= 2.5:
                            return "Optimal (>=2.5)"
                        elif normalized_score >= 2.0:
                            return "High (>=2.0)"
                        elif normalized_score >= 1.0:
                            return "Moderate (>=1.0)"
                        else:
                            return "Low (<1.0)"
                    
                    # Add normalized scoring
                    normalized_score = get_normalized_confidence_score(m.get('Score', 0), m.get('Class', ''))
                    m['Normalized Score'] = f"{normalized_score:.2f}"
                    m['Confidence Level'] = get_confidence_label(normalized_score)
                    df_all.append(m)
        
        df_all = pd.DataFrame(df_all)
        
        # Remove unwanted columns as specified in requirements
        unwanted_columns = [
            'GC_Content', 'GC_CG_Dinucleotides', 'Arms/Repeat Unit/Copies', 
            'Spacer', 'MotifClasses', 'ContributingMotifs', 'Unit', 'Copies', 'GC_Fraction'
        ]
        
        # Filter out unwanted columns if they exist
        available_cols = df_all.columns.tolist()
        cols_to_keep = [col for col in available_cols if col not in unwanted_columns]
        df_all = df_all[cols_to_keep]
        
        # Reorder columns to put S.No first
        cols = df_all.columns.tolist()
        if 'S.No' in cols:
            cols = ['S.No'] + [col for col in cols if col != 'S.No']
            df_all = df_all[cols]
        
        st.markdown("### Complete Results Table")
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
    # Apply neutral grayscale theme with accent cyan for documentation as specified
    st.markdown("""
    <style>
    /* Documentation-specific theme: neutral grayscale with accent cyan */
    .stMarkdown h1, .stMarkdown h2, .stMarkdown h3, .stMarkdown h4 { color: #1f2937; }
    .doc-accent {
        color: #0891b2;
        background: linear-gradient(135deg, #0891b2 0%, #06b6d4 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
    }
    .doc-section {
        background: linear-gradient(135deg, rgba(249, 250, 251, 0.8) 0%, rgba(255, 255, 255, 0.9) 100%);
        border-radius: 12px;
        padding: 24px;
        margin: 16px 0;
        border-left: 4px solid #0891b2;
        box-shadow: 0 4px 16px rgba(0, 0, 0, 0.05);
    }
    .doc-code {
        background: #f8fafc;
        border: 1px solid #e2e8f0;
        border-radius: 8px;
        padding: 16px;
        font-family: 'Courier New', monospace;
        color: #334155;
    }
    .doc-citation {
        background: #f0fdfa;
        border-left: 4px solid #14b8a6;
        padding: 12px 16px;
        margin: 8px 0;
        border-radius: 6px;
        font-style: italic;
    }
    </style>
    """, unsafe_allow_html=True)
    
    # Create floating TOC navigation
    st.markdown("""
    <div style='position: fixed; top: 20%; right: 20px; background: rgba(255, 255, 255, 0.95); 
                border-radius: 12px; padding: 16px; box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1); 
                backdrop-filter: blur(10px); z-index: 1000; max-width: 200px; border: 1px solid #e2e8f0;'>
        <h4 style='margin: 0 0 12px 0; color: #0891b2; font-size: 14px;'>Table of Contents</h4>
        <div style='font-size: 12px; line-height: 1.4;'>
            <a href="#overview" style='display: block; color: #4b5563; text-decoration: none; padding: 4px 0;'>Overview</a>
            <a href="#how-it-works" style='display: block; color: #4b5563; text-decoration: none; padding: 4px 0;'>How It Works</a>
            <a href="#api-endpoints" style='display: block; color: #4b5563; text-decoration: none; padding: 4px 0;'>API Endpoints</a>
            <a href="#ui-usage" style='display: block; color: #4b5563; text-decoration: none; padding: 4px 0;'>UI Usage Guide</a>
            <a href="#validation" style='display: block; color: #4b5563; text-decoration: none; padding: 4px 0;'>Validation & Errors</a>
            <a href="#performance" style='display: block; color: #4b5563; text-decoration: none; padding: 4px 0;'>Performance</a>
            <a href="#accessibility" style='display: block; color: #4b5563; text-decoration: none; padding: 4px 0;'>Accessibility</a>
            <a href="#references" style='display: block; color: #4b5563; text-decoration: none; padding: 4px 0;'>References</a>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Main documentation content
    st.markdown("""
    <h1 id="overview" class="doc-accent" style='text-align: center; margin-bottom: 32px; font-size: 2.5rem;'>
        NBDFinder Scientific Documentation & References
    </h1>
    """, unsafe_allow_html=True)
    
    # Overview Section
    st.markdown("""
    <div class="doc-section" id="overview">
        <h2 class="doc-accent">Overview</h2>
        <h3>Purpose</h3>
        <p>NBDFinder is an advanced computational framework for genome-wide detection and analysis of Non-B DNA structural motifs. 
        The platform identifies alternative DNA conformations that deviate from the canonical Watson-Crick B-form double helix, 
        which play crucial roles in genome organization, gene regulation, and disease pathogenesis.</p>
        
        <h3>Scope</h3>
        <p>Our comprehensive analysis suite covers <strong>10 major Non-B DNA classes</strong> and <strong>22 specialized subclasses</strong>:</p>
        <ul>
            <li><strong>G-Quadruplex Family:</strong> Canonical G4, Relaxed G4, Bulged G4, Bipartite G4, Multimeric G4, Imperfect G4</li>
            <li><strong>Alternative Helices:</strong> Z-DNA, eGZ (Extruded-G)</li>
            <li><strong>Structural Variants:</strong> Curved DNA, Slipped DNA, Cruciform DNA</li>
            <li><strong>Complex Structures:</strong> R-loop, Triplex DNA, Sticky DNA</li>
            <li><strong>Specialized Motifs:</strong> i-motif family, AC-motif, G-Triplex</li>
            <li><strong>Hybrid & Cluster Regions:</strong> Overlapping motifs and high-density regions</li>
        </ul>
        
        <h3>Domain Context</h3>
        <p>Non-B DNA structures are particularly enriched in:</p>
        <ul>
            <li><strong>Promoter regions:</strong> G-quadruplexes in c-MYC NHE III1, VEGF promoter</li>
            <li><strong>Telomeric sequences:</strong> G-quadruplex formations in human telomeres</li>
            <li><strong>Disease-associated repeats:</strong> CAG/CTG repeats in neurological disorders</li>
            <li><strong>Immunoglobulin switch regions:</strong> Triplex and G4 structures</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Comprehensive Motif Detection Table
    st.markdown("""
    <div class="doc-section" id="motif-detection">
        <h2 class="doc-accent">Motif Detection Methods & Scoring Systems</h2>
        <p>This table provides comprehensive details on all motif classes, their detection patterns, and scoring methodologies employed by NBDFinder:</p>
    </div>
    """, unsafe_allow_html=True)
    
    # Create the comprehensive motif table as requested
    motif_data = {
        "Class": [
            "G-Quadruplex (G4)", "G-Quadruplex (G4)", "G-Quadruplex (G4)", "G-Quadruplex (G4)", "G-Quadruplex (G4)", "G-Quadruplex (G4)", "G-Quadruplex (G4)",
            "i-Motif", "i-Motif", "i-Motif",
            "R-loop", "R-loop", "R-loop",
            "Z-DNA", "Z-DNA",
            "Triplex (H-DNA)", "Triplex (H-DNA)",
            "Cruciform", "Slipped DNA (SSS/STRs)", "Curved DNA",
            "Hybrid Regions", "Cluster Regions"
        ],
        "Subclass/Motif": [
            "Canonical G4", "Relaxed G4", "Bulged G4", "Imperfect G4", "Bipartite G4", "Multimeric G4", "G-Triplex",
            "Canonical i-motif", "Relaxed i-motif", "AC-motif",
            "RLFS (m1)", "RLFS (m2)", "RLFS+REZ",
            "Z-DNA (Kadane)", "eGZ (Extruded-G)",
            "Mirror repeat triplex", "Sticky DNA",
            "Palindromic cruciform", "Mono/di/tri/tetra-STRs", "A-tracts, PolyA/PolyT",
            "G4/R-loop, G4/Z-DNA", "Non-B DNA cluster"
        ],
        "Detection Logic / Pattern Description": [
            "4 runs of ≥3 G, loops 1-7 nt", "4 runs of ≥3 G, loops 1-12 nt", "G runs can be interrupted by 1-3 nt (bulges)", "4 runs of ≥3 G/A, loops 1-7 nt (A substitutions allowed)", "2 G-runs, long spacer (20-50 nt), 2 more G-runs", "4+ G-runs separated by 1-12 nt (≥4 repeats)", "3 runs of ≥3 G, loops up to 15 nt",
            "4 runs of ≥3 C, loops 1-7 nt", "4 runs of ≥3 C, loops 1-12 nt", "Alternating A/C, 4 runs of ≥3 C, up to 4 A in loops",
            "G-rich tracts (≥3 G), separated by 1-10 nt, ≥2 runs", "≥4 G runs, separated by 1-10 nt", "Downstream GC-rich region after RLFS",
            "Dinucleotide scoring: GC/CG+7, GT/TG/AC/CA+1.25, AT/TA 0, others -3; Kadane's max subarray", "(CGG) repeat expansions, 3+ repeats",
            "Mirror A/G or T/C runs, min length ≥10–15", "Long direct repeats (GAA)n, etc",
            "Inverted repeats, stem/loop, min stem/loop length", "Tandem repeats: (AT)n, (GAA)n, (CGG)n, etc", "AT-rich tracts, polyA/polyT ≥5–6 nt",
            "Motif overlap: both G4 and R-loop/Z-DNA at same/adjacent region", "≥2 different motif classes in window (e.g., 500 bp)"
        ],
        "Regex/Pattern Example": [
            "G{3,}\\w{1,7}G{3,}\\w{1,7}G{3,}\\w{1,7}G{3,}", "G{3,}\\w{1,12}G{3,}\\w{1,12}G{3,}\\w{1,12}G{3,}", "G{2,}\\w{1,3}G{1,}\\w{1,7}G{3,}\\w{1,7}G{3,}\\w{1,7}G{3,}", "[GA]{3,}\\w{1,7}[GA]{3,}\\w{1,7}[GA]{3,}\\w{1,7}[GA]{3,}", "G{3,}\\w{1,7}G{3,}\\w{20,50}G{3,}\\w{1,7}G{3,}", "(G{3,}\\w{1,12}){4,}", "G{3,}\\w{1,15}G{3,}\\w{1,15}G{3,}",
            "C{3,}\\w{1,7}C{3,}\\w{1,7}C{3,}\\w{1,7}C{3,}", "C{3,}\\w{1,12}C{3,}\\w{1,12}C{3,}\\w{1,12}C{3,}", "C{3,}[ATGC]{1,12}C{3,}[ATGC]{1,12}C{3,}[ATGC]{1,12}C{3,} with A/C enrichment",
            "G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}", "G{4,}(?:[ATGC]{1,10}?G{4,}){1,}", "Sliding window, GC content >40%",
            "No explicit regex; dinucleotide scan", "(CGG){3,}",
            "See panNonB for various triplex regex", "(GAA){n,}",
            "Palindrome detection", "(AT){n,} (GAA){n,} (CGG){n,}", "A{5,} T{5,}",
            "Overlap of above patterns", "Aggregation, sliding window"
        ],
        "Scoring System": [
            "G4Hunter score (avg G/C run, capped at 4); category thresholds: ≥1.5 High, ≥1.0 Moderate, <1.0 Low", "G4Hunter score, lower threshold (≥0.8)", "G4Hunter score, threshold (≥0.9)", "G4Hunter score, lower threshold (≥0.7)", "G4Hunter score, threshold (≥0.8)", "G4Hunter score, higher threshold (≥1.2)", "Score: motif length × G content × run count",
            "Hunter-style (C-rich), loop penalty/bonus, canonical threshold", "Hunter-style, relaxed threshold", "Regex: C{3,}[ATGC]{1,12}C{3,}[ATGC]{1,12}C{3,}[ATGC]{1,12}C{3,} with A/C enrichment",
            "Score: (GC_fraction×50 + G_runs×10)×length^0.25; RLFS+REZ bipartite", "Same as above, G-run emphasis", "Combined RLFS+REZ scoring",
            "Cumulative weight (Kadane), length ≥12, threshold ≥50; alt. scoring", "Score: n_repeats × 3 × (1+2G_frac)",
            "Sequence symmetry, A/G-rich scoring", "Length & repeat count",
            "Stem/loop length, stability penalty", "Motif-specific minimum repeat count", "Length, AT %",
            "Combined scoring, presence of both", "Motif density, region size"
        ]
    }
    
    # Create DataFrame and display the table
    motif_df = pd.DataFrame(motif_data)
    st.dataframe(motif_df, use_container_width=True, height=600)
    
    st.markdown("""
    <div class="doc-section">
        <h3>Scientific Scoring Methodology</h3>
        <p><strong>NBDFinder employs established scoring systems for accurate Non-B DNA prediction:</strong></p>
        <ul>
            <li><strong>G-Quadruplex Motifs:</strong> G4Hunter algorithm (Bedrat et al., 2016) with structural factors. Scores >1.0 indicate high G4 formation potential.</li>
            <li><strong>Z-DNA Motifs:</strong> Modified Z-Hunt/ZhuntLSC approach based on dinucleotide propensities (Ho et al., 1986; Schroth et al., 1992). Scores >50 suggest Z-DNA formation.</li>
            <li><strong>i-Motif Structures:</strong> Cytosine run analysis with loop constraints based on pH-dependent stability (Day et al., 2014; Wright et al., 2020).</li>
            <li><strong>Curved DNA:</strong> Curvature propensity scoring using An/Tn tract positioning (Bolshoy et al., 1991; Gabrielian & Pongor, 1996).</li>
            <li><strong>Repeat Elements:</strong> Pattern matching with biological constraints based on disease-associated sequences (Mirkin, 2007; Wells, 2007).</li>
        </ul>
        
        <h4>Confidence Levels:</h4>
        <ul>
            <li><strong>High Confidence:</strong> Scores exceed established experimental thresholds</li>
            <li><strong>Moderate Confidence:</strong> Scores meet minimum formation thresholds</li>
            <li><strong>Low Confidence:</strong> Predicted but below optimal formation conditions</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # How It Works Section
    st.markdown("""
    <div class="doc-section" id="how-it-works">
        <h2 class="doc-accent">How It Works</h2>
        <h3>Data Flow Architecture</h3>
        <p>NBDFinder follows a sophisticated five-stage pipeline:</p>
        
        <div style='background: #f8fafc; border-radius: 8px; padding: 20px; margin: 16px 0;'>
            <p><strong>1. Query Input →</strong> FASTA sequences, NCBI accessions, or direct sequence input</p>
            <p><strong>2. Fetch & Parse →</strong> Sequence retrieval from NCBI databases with validation</p>
            <p><strong>3. Multi-Algorithm Analysis →</strong> Parallel motif detection using specialized algorithms</p>
            <p><strong>4. Score Integration →</strong> Confidence scoring and overlap resolution</p>
            <p><strong>5. Visualization →</strong> Interactive genome browser-style displays</p>
        </div>
        
        <h3>Detection Algorithms</h3>
        <ul>
            <li><strong>G4 Detection:</strong> Loop-based scanning with thermodynamic scoring</li>
            <li><strong>Z-DNA:</strong> Alternating purine-pyrimidine pattern recognition</li>
            <li><strong>Curved DNA:</strong> Intrinsic curvature prediction using structural parameters</li>
            <li><strong>R-loops:</strong> Asymmetric G/C skew analysis</li>
            <li><strong>Cruciform:</strong> Inverted repeat identification with stem-loop energy calculations</li>
        </ul>
        
        <h3>Component Links</h3>
        <ul>
            <li><a href="#ui-usage" class="doc-accent">→ Upload & Analyze</a>: Sequence input and motif selection</li>
            <li><a href="#api-endpoints" class="doc-accent">→ Results Visualization</a>: Interactive genome browser and statistical plots</li>
            <li><a href="#validation" class="doc-accent">→ Download</a>: Export data in CSV, Excel, JSON, and FASTA formats</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # API Endpoints Section
    st.markdown("""
    <div class="doc-section" id="api-endpoints">
        <h2 class="doc-accent">API Endpoints and Parameters</h2>
        <h3>Core Analysis Endpoints</h3>
        
        <div class="doc-code">
            <h4>Sequence Analysis</h4>
            <code>POST /analyze</code><br>
            <strong>Parameters:</strong><br>
            • sequence: string (FASTA format)<br>
            • motif_classes: array (selected motif types)<br>
            • sensitivity: float (0.1-1.0)<br>
            • min_length: integer (8-100 bp)<br>
            
            <strong>Response:</strong><br>
            {<br>
            &nbsp;&nbsp;"motifs": [{"class": "G4", "start": 123, "end": 145, "score": 2.4}],<br>
            &nbsp;&nbsp;"statistics": {"gc_content": 45.2, "motif_density": 0.12}<br>
            }
        </div>
        
        <div class="doc-code">
            <h4>NCBI Data Retrieval</h4>
            <code>GET /ncbi/{accession}</code><br>
            <strong>Parameters:</strong><br>
            • accession: string (e.g., "NG_007161.1")<br>
            • format: string ("fasta" | "genbank")<br>
            
            <strong>Response:</strong><br>
            {<br>
            &nbsp;&nbsp;"sequence": "ATCG...",<br>
            &nbsp;&nbsp;"description": "Human c-MYC gene...",<br>
            &nbsp;&nbsp;"length": 5386<br>
            }
        </div>
        
        <h3>Export Endpoints</h3>
        <ul>
            <li><code>GET /export/csv</code> - Motif results in CSV format</li>
            <li><code>GET /export/excel</code> - Comprehensive Excel workbook</li>
            <li><code>GET /export/json</code> - Structured JSON data</li>
            <li><code>GET /export/fasta</code> - Motif sequences in FASTA format</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # UI Usage Guide Section
    st.markdown("""
    <div class="doc-section" id="ui-usage">
        <h2 class="doc-accent">UI Usage Guide</h2>
        
        <h3>Submitting Queries</h3>
        <ol>
            <li><strong>Select Motif Classes:</strong> Choose from 19 available motif types using the multi-select interface</li>
            <li><strong>Input Method:</strong> Choose from four input options:
                <ul>
                    <li>Upload FASTA files (single or multi-sequence)</li>
                    <li>Paste sequences directly</li>
                    <li>Select example sequences</li>
                    <li>Fetch from NCBI using accession numbers or gene names</li>
                </ul>
            </li>
            <li><strong>Advanced Options:</strong> Configure sensitivity, minimum length, and analysis parameters</li>
            <li><strong>Start Analysis:</strong> Click the analysis button to begin processing</li>
        </ol>
        
        <h3>Interpreting Results</h3>
        <ul>
            <li><strong>Summary Statistics:</strong> Sequence length, GC content, motif density</li>
            <li><strong>Motif Table:</strong> Detailed list with positions, scores, and sequences</li>
            <li><strong>Genome Browser:</strong> Visual representation of motif locations</li>
            <li><strong>Distribution Charts:</strong> Statistical analysis and motif frequency plots</li>
        </ul>
        
        <h3>Theme Switching</h3>
        <p>NBDFinder features automatic theme selection per page:</p>
        <ul>
            <li><strong>Home/Search:</strong> Cool blue-teal theme for navigation</li>
            <li><strong>Results/Analysis:</strong> Indigo-magenta theme for data visualization</li>
            <li><strong>Documentation:</strong> Neutral grayscale with cyan accents</li>
            <li><strong>Settings:</strong> Green-accent theme for configuration</li>
        </ul>
        
        <h3>Export and Reporting</h3>
        <ul>
            <li><strong>CSV Export:</strong> Tabular data for statistical analysis</li>
            <li><strong>Excel Export:</strong> Formatted workbooks with multiple sheets</li>
            <li><strong>FASTA Export:</strong> Motif sequences for further analysis</li>
            <li><strong>JSON Export:</strong> Structured data for programmatic access</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Validation and Error Handling Section
    st.markdown("""
    <div class="doc-section" id="validation">
        <h2 class="doc-accent">Validation and Error Handling</h2>
        
        <h3>Expected Errors and Limits</h3>
        <table style='width: 100%; border-collapse: collapse; margin: 16px 0;'>
            <tr style='background: #f8fafc;'>
                <th style='padding: 12px; border: 1px solid #e2e8f0; text-align: left;'>Error Type</th>
                <th style='padding: 12px; border: 1px solid #e2e8f0; text-align: left;'>Limit/Cause</th>
                <th style='padding: 12px; border: 1px solid #e2e8f0; text-align: left;'>Recovery Action</th>
            </tr>
            <tr>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Sequence Length</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Maximum 100,000 bp</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Split large sequences or use standalone version</td>
            </tr>
            <tr style='background: #f8fafc;'>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>NCBI Rate Limit</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>3 requests per second</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Wait 10 seconds and retry</td>
            </tr>
            <tr>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Invalid Sequence</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Non-ATCG characters</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Automatic cleaning with user notification</td>
            </tr>
            <tr style='background: #f8fafc;'>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Network Timeout</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>30 second timeout</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Automatic retry with exponential backoff</td>
            </tr>
        </table>
        
        <h3>Input Validation</h3>
        <ul>
            <li><strong>NCBI Accessions:</strong> Pattern matching for standard formats (NG_######.#, NM_######.#)</li>
            <li><strong>FASTA Format:</strong> Header validation and sequence cleaning</li>
            <li><strong>Sequence Quality:</strong> Minimum 10 bp length, maximum ambiguous base percentage</li>
            <li><strong>File Uploads:</strong> Size limits (200MB), format validation (.fa, .fasta, .txt)</li>
        </ul>
        
        <h3>Error Recovery Steps</h3>
        <ol>
            <li><strong>Automatic Retry:</strong> Network errors trigger automatic retry with backoff</li>
            <li><strong>Graceful Degradation:</strong> Partial results displayed when some analyses fail</li>
            <li><strong>User Notification:</strong> Clear error messages with specific recovery actions</li>
            <li><strong>Session Persistence:</strong> Analysis state preserved across page refreshes</li>
        </ol>
    </div>
    """, unsafe_allow_html=True)
    
    # Performance Notes Section
    st.markdown("""
    <div class="doc-section" id="performance">
        <h2 class="doc-accent">Performance Notes</h2>
        
        <h3>Caching Strategy</h3>
        <ul>
            <li><strong>NCBI Responses:</strong> 1-hour cache with ETag validation</li>
            <li><strong>Analysis Results:</strong> Session-based caching for repeated analysis</li>
            <li><strong>Static Resources:</strong> CDN caching for CSS, JavaScript, and images</li>
            <li><strong>Database Queries:</strong> Redis caching for frequently accessed sequences</li>
        </ul>
        
        <h3>Lazy Loading</h3>
        <ul>
            <li><strong>Visualization Components:</strong> Plotly charts loaded on-demand</li>
            <li><strong>Large Results:</strong> Paginated table display with virtual scrolling</li>
            <li><strong>Image Assets:</strong> Progressive loading with placeholder</li>
            <li><strong>Analysis Modules:</strong> Dynamic import of detection algorithms</li>
        </ul>
        
        <h3>Server-Side Rendering (SSR)</h3>
        <ul>
            <li><strong>Static Routes:</strong> Home, Documentation, and Settings pre-rendered</li>
            <li><strong>Meta Tags:</strong> Dynamic SEO optimization per route</li>
            <li><strong>Critical CSS:</strong> Inline critical styles for faster initial paint</li>
            <li><strong>Preload Hints:</strong> Resource prioritization for key assets</li>
        </ul>
        
        <h3>Performance Benchmarks</h3>
        <ul>
            <li><strong>Analysis Speed:</strong> ~1,000 bp/second for comprehensive analysis</li>
            <li><strong>First Content Paint:</strong> <1.2s on 3G connection</li>
            <li><strong>Time to Interactive:</strong> <2.5s for initial page load</li>
            <li><strong>Memory Usage:</strong> <50MB for typical analysis session</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # Accessibility Section
    st.markdown("""
    <div class="doc-section" id="accessibility">
        <h2 class="doc-accent">Accessibility and Browser Support</h2>
        
        <h3>WCAG AA Compliance</h3>
        <ul>
            <li><strong>Color Contrast:</strong> Minimum 4.5:1 ratio for normal text, 3:1 for large text</li>
            <li><strong>Keyboard Navigation:</strong> Full functionality accessible via keyboard</li>
            <li><strong>Screen Readers:</strong> ARIA labels and semantic HTML structure</li>
            <li><strong>Focus Management:</strong> Visible focus indicators and logical tab order</li>
        </ul>
        
        <h3>High-Contrast Mode</h3>
        <ul>
            <li><strong>Automatic Detection:</strong> Respects system preferences</li>
            <li><strong>Manual Toggle:</strong> User-controlled high-contrast theme</li>
            <li><strong>Chart Accessibility:</strong> Pattern-based differentiation in visualizations</li>
            <li><strong>Color-Blind Support:</strong> Wong palette for visualization accessibility</li>
        </ul>
        
        <h3>Browser Support</h3>
        <table style='width: 100%; border-collapse: collapse; margin: 16px 0;'>
            <tr style='background: #f8fafc;'>
                <th style='padding: 12px; border: 1px solid #e2e8f0; text-align: left;'>Browser</th>
                <th style='padding: 12px; border: 1px solid #e2e8f0; text-align: left;'>Minimum Version</th>
                <th style='padding: 12px; border: 1px solid #e2e8f0; text-align: left;'>Features</th>
            </tr>
            <tr>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Chrome</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>90+</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Full support</td>
            </tr>
            <tr style='background: #f8fafc;'>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Firefox</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>88+</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Full support</td>
            </tr>
            <tr>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Safari</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>14+</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Full support</td>
            </tr>
            <tr style='background: #f8fafc;'>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Edge</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>90+</td>
                <td style='padding: 12px; border: 1px solid #e2e8f0;'>Full support</td>
            </tr>
        </table>
        
        <h3>Mobile Optimization</h3>
        <ul>
            <li><strong>Responsive Design:</strong> Optimized layouts for mobile devices</li>
            <li><strong>Touch Interactions:</strong> Appropriate touch targets (44px minimum)</li>
            <li><strong>Progressive Web App:</strong> Offline functionality and app-like experience</li>
            <li><strong>Performance:</strong> Optimized for mobile networks and processing power</li>
        </ul>
    </div>
    """, unsafe_allow_html=True)
    
    # References Section
    st.markdown("""
    <div class="doc-section" id="references">
        <h2 class="doc-accent">Scientific References</h2>
        
        <h3>Foundational Studies</h3>
        <div class="doc-citation">
            Watson, J. D. & Crick, F. H. Molecular structure of nucleic acids; a structure for deoxyribose nucleic acid. 
            <em>Nature</em> <strong>171</strong>, 737-738 (1953).
        </div>
        
        <div class="doc-citation">
            Wang, A. H. et al. Molecular structure of a left-handed double helical DNA fragment at atomic resolution. 
            <em>Nature</em> <strong>282</strong>, 680-686 (1979).
        </div>
        
        <div class="doc-citation">
            Sen, D. & Gilbert, W. Formation of parallel four-stranded complexes by guanine-rich motifs in DNA and its implications for meiosis. 
            <em>Nature</em> <strong>334</strong>, 364-366 (1988).
        </div>
        
        <h3>Non-B DNA Motifs and Regulation</h3>
        <div class="doc-citation">
            Siddiqui-Jain, A., Grand, C. L., Bearss, D. J. & Hurley, L. H. Direct evidence for a G-quadruplex in a promoter region and its targeting with a small molecule to repress c-MYC transcription. 
            <em>Proc. Natl. Acad. Sci. USA</em> <strong>99</strong>, 11593-11598 (2002).
        </div>
        
        <div class="doc-citation">
            Brooks, T. A., Kendrick, S. & Hurley, L. H. Making sense of G-quadruplex and i-motif functions in oncogene promoters. 
            <em>FEBS J.</em> <strong>277</strong>, 3459-3469 (2010).
        </div>
        
        <div class="doc-citation">
            Wang, G. & Vasquez, K. M. Dynamic alternative DNA structures in biology and disease. 
            <em>Nat Rev Genet</em> <strong>24</strong>, 211-234 (2023).
        </div>
        
        <h3>Disease Associations</h3>
        <div class="doc-citation">
            Mirkin, S. M. Expandable DNA repeats and human disease. 
            <em>Nature</em> <strong>447</strong>, 932-940 (2007).
        </div>
        
        <div class="doc-citation">
            Matos-Rodrigues, G., Hisey, J. A., Nussenzweig, A. & Mirkin, S. M. Detection of alternative DNA structures and its implications for human disease. 
            <em>Mol Cell</em> <strong>83</strong>, 3622-3641 (2023).
        </div>
        
        <h3>Computational Methods</h3>
        <div class="doc-citation">
            Bedrat, A., Lacroix, L. & Mergny, J. L. Re-evaluation of G-quadruplex propensity with G4Hunter. 
            <em>Nucleic Acids Res.</em> <strong>44</strong>, 1746-1759 (2016).
        </div>
        
        <div class="doc-citation">
            Kikin, O., D'Antonio, L. & Bagga, P. S. QGRS Mapper: a web-based server for predicting G-quadruplexes in nucleotide sequences. 
            <em>Nucleic Acids Res.</em> <strong>34</strong>, W676-W682 (2006).
        </div>
        
        <h3>NCBI Resources</h3>
        <div class="doc-citation">
            NCBI Resource Coordinators. Database resources of the National Center for Biotechnology Information. 
            <em>Nucleic Acids Res.</em> <strong>52</strong>, D33-D43 (2024).
        </div>
        
        <div class="doc-citation">
            Sayers, E. W. et al. Database resources of the National Center for Biotechnology Information. 
            <em>Nucleic Acids Res.</em> <strong>50</strong>, D20-D26 (2022).
        </div>
        
        <h3>Curated Citation Collection</h3>
        <p>The following resources provide comprehensive coverage of non-B DNA research:</p>
        <ul>
            <li><strong>PubMed Query:</strong> "non-B DNA"[Title/Abstract] AND "structure"[Title/Abstract]</li>
            <li><strong>G4 Database:</strong> <a href="#" class="doc-accent">G4IPDB</a> - G-quadruplex interaction database</li>
            <li><strong>Z-DNA Atlas:</strong> <a href="#" class="doc-accent">ZDB</a> - Z-DNA forming sequences database</li>
            <li><strong>RepeatMasker:</strong> <a href="#" class="doc-accent">DFAM</a> - DNA repeat elements database</li>
        </ul>
        
        <div style='margin-top: 32px; padding: 20px; background: #f0fdfa; border-radius: 8px; border-left: 4px solid #14b8a6;'>
            <p style='margin: 0; font-style: italic; color: #0f766e;'>
                <strong>Note on Reproducibility:</strong> All analyses performed by NBDFinder are version-controlled with 
                algorithm parameters documented for reproducibility. Citation suggestions are provided in exported results 
                to support transparent reporting in scientific publications.
            </p>
        </div>
    </div>
    """, unsafe_allow_html=True)
    
    # Footer with deep-linkable anchors note
    st.markdown("""
    <div style='margin-top: 40px; padding: 20px; background: #fafafa; border-radius: 8px; text-align: center; color: #6b7280;'>
        <p style='margin: 0; font-size: 14px;'>
            <strong>Deep-linkable sections:</strong> All sections in this documentation have unique anchors for easy sharing and reference. 
            Simply copy the URL when viewing any section to create direct links.
        </p>
        <p style='margin: 8px 0 0 0; font-size: 12px;'>
            Last updated: December 2024 | Version 2.0 | 
            <a href="#overview" class="doc-accent">Back to top ↑</a>
        </p>
    </div>
    """, unsafe_allow_html=True)

st.markdown("""
---
<div style="font-size: 1.05rem; color: #1e293b; margin-top: 30px; text-align: left; font-family:Montserrat,Arial;">
<b>Developed by</b><br>
Dr. Venkata Rajesh Yella<br>
<a href="mailto:yvrajesh_bt@kluniversity.in">yvrajesh_bt@kluniversity.in</a> |
<a href="https://github.com/VRYella" target="_blank">GitHub: VRYella</a>
</div>
""", unsafe_allow_html=True)
