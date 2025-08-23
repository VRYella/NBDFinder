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
    print("Advanced visualization not available. Install additional dependencies.")

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
    page_icon="🧬",
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
    """Guarantee every motif has a string 'Subtype' using proper classification mapping"""
    from motifs.classification_config import get_official_classification, LEGACY_TO_OFFICIAL_MAPPING
    
    if isinstance(motif, dict):
        if 'Subtype' not in motif or motif['Subtype'] is None or motif['Subtype'] == '':
            # Try to get proper classification based on Class
            motif_class = motif.get('Class', '')
            if motif_class in LEGACY_TO_OFFICIAL_MAPPING:
                # Use proper subtype from mapping
                official_class, official_subtype = get_official_classification(motif_class, '')
                motif['Subtype'] = official_subtype if official_subtype else 'Canonical'
            else:
                # Generic fallback that's more meaningful than 'Other'
                motif['Subtype'] = 'Canonical'
        return motif
    else:
        # Handle non-dict motifs gracefully
        return {'Subtype': 'Canonical', 'Motif': motif}

# ---------- FLOATING COMPONENTS ----------
def add_floating_back_to_top():
    """Add floating back to top button"""
    st.markdown("""
    <div id="back-to-top" style="
        position: fixed; 
        bottom: 20px; 
        right: 20px; 
        background: linear-gradient(135deg, #1e3a8a 0%, #0891b2 100%);
        color: white; 
        padding: 12px 16px; 
        border-radius: 50px; 
        cursor: pointer; 
        box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        z-index: 1000;
        font-weight: 600;
        font-size: 14px;
        transition: all 0.3s ease;
        border: none;
        backdrop-filter: blur(10px);
    " onclick="window.scrollTo({top: 0, behavior: 'smooth'})">
        ↑ Top
    </div>
    
    <style>
    #back-to-top:hover {
        transform: translateY(-2px);
        box-shadow: 0 6px 20px rgba(0,0,0,0.2);
    }
    </style>
    """, unsafe_allow_html=True)



# ---------- DESIGN SYSTEM ----------
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
        -webkit-background-clip: text; -webkit-text-fill-color: transparent;
        background-clip: text;
    }
    h1 { font-size: 3.5rem; font-weight: 900; letter-spacing: -0.04em; margin-bottom: 1.5rem; }
    h2 { font-size: 2.75rem; font-weight: 800; letter-spacing: -0.03em; }
    h3 { font-size: 2.25rem; font-weight: 750; letter-spacing: -0.025em; }
    h4 { font-size: 1.875rem; font-weight: 700; letter-spacing: -0.02em; }
    h5 { font-size: 1.5rem; font-weight: 650; }
    h6 { font-size: 1.375rem; font-weight: 600; }
    
    .stMarkdown, .markdown-text-container, .stText, p, span, label, input {
        font-family: Inter, sans-serif; font-size: 1.25rem; line-height: 1.7; 
        color: var(--text); font-weight: 400;
    }
    .stMarkdown p { margin-bottom: 1.25rem; font-size: 1.25rem; }
    
    /* Better readability for all text elements */
    .stSelectbox label, .stMultiSelect label, .stTextInput label, .stTextArea label,
    .stSlider label, .stRadio label, .stCheckbox label {
        font-size: 1.125rem !important; font-weight: 500;
    }
    
    /* Form inputs with larger text */
    .stSelectbox div[data-baseweb="select"] > div, 
    .stTextInput input, .stTextArea textarea {
        font-size: 1.125rem !important;
    }
    
    /* EXPANDERS */
    .streamlit-expanderHeader {
        background: var(--surface); border: 1px solid var(--border); border-radius: var(--radius);
        font-weight: 500; color: var(--text);
    }
    .streamlit-expanderContent {
        background: var(--bg); border: 1px solid var(--border); border-top: none;
        border-radius: 0 0 var(--radius) var(--radius);
    }
    
    /* CUSTOM CARDS */
    .feature-card {
        background: var(--bg); border: 1px solid var(--border); border-radius: var(--radius);
        padding: 24px; margin: 16px 0; box-shadow: var(--shadow); transition: var(--transition);
    }
    .feature-card:hover {
        transform: translateY(-2px); box-shadow: 0 4px 12px rgba(0,0,0,0.1);
    }
    </style>
    """, unsafe_allow_html=True)

# ---------- CONSTANTS ----------
MOTIF_COLORS = {
    "Curved DNA": "#3B82F6",
    "Slipped DNA": "#10B981", 
    "Cruciform": "#F59E0B",
    "R-loop": "#EF4444",
    "Triplex": "#8B5CF6",
    "G-Quadruplex": "#06B6D4",
    "i-motif": "#EC4899",
    "Z-DNA": "#84CC16",
    "Hybrid": "#F97316",
    "Non-B DNA Clusters": "#EF4444"
}

# Main application pages structure
MAIN_PAGES = {
    "Home": "Home / Welcome",
    "Upload & Analyze": "Sequence Upload and Analysis",
    "Results": "Analysis Results and Visualization", 
    "Visualization": "Publication-Ready Visualizations",
    "Clinical/Disease": "Disease Annotation and Clinical Data",
    "Download & Export": "Download Results and Export Data",
    "Documentation": "Scientific Documentation & References",
    "Settings": "Application Settings and Preferences"
}

# Initialize session state
if 'current_page' not in st.session_state:
    st.session_state.current_page = "Home"
if 'seqs' not in st.session_state:
    st.session_state.seqs = []
if 'names' not in st.session_state:
    st.session_state.names = []
if 'results' not in st.session_state:
    st.session_state.results = []
if 'analysis_settings' not in st.session_state:
    st.session_state.analysis_settings = {}

# Set up Entrez for NCBI access
Entrez.email = "raazbiochem@gmail.com"
Entrez.api_key = None

# Example sequences
EXAMPLE_SEQUENCES = {
    "Human Telomere G4-Rich": {
        "name": "G4_Rich_Human_Telomere_Example", 
        "sequence": "TTAGGGTTAGGGTTAGGGTTAGGGAAAAATCCGTCGAGCAGAGTTAAAAAGGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCCGAAAGAAAGAAAGAAAGAAACGCGCGCGCGCGCGCGCGCGATCGCACACACACAGCTGCTGCTGC",
        "description": "Human telomeric sequence rich in G-quadruplex structures"
    },
    "Z-DNA Forming Sequence": {
        "name": "Z_DNA_Example_Sequence",
        "sequence": "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
        "description": "Alternating CG sequence prone to Z-DNA formation"
    },
    "Disease-Associated GAA Repeats": {
        "name": "Disease_GAA_Repeats_Example",
        "sequence": "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA",
        "description": "Frataxin GAA repeat expansion associated with Friedreich's ataxia"
    }
}

# ---- MAIN APPLICATION ----

# Add floating components
add_floating_back_to_top()

# Header
st.markdown("""
<div style='text-align: center; padding: 20px; margin-bottom: 30px; background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 100%); border-radius: 12px; border: 1px solid #e2e8f0;'>
    <h1 style='color: #1e3a8a; font-family: Inter, sans-serif; font-weight: 700; margin-bottom: 8px; font-size: 2.5rem;'>
        NBDFinder: Non-B DNA Analysis Platform
    </h1>
    <p style='color: #64748b; font-size: 1.125rem; margin: 0;'>Publication-ready computational framework for genome-wide detection and analysis of non-B DNA structural motifs</p>
</div>
""", unsafe_allow_html=True)

# Main navigation tabs
main_tabs = st.tabs(list(MAIN_PAGES.keys()))
tab_dict = dict(zip(MAIN_PAGES.keys(), main_tabs))

# =====================================================
# HOME / WELCOME PAGE
# =====================================================
with tab_dict["Home"]:
    col1, col2 = st.columns([1.2, 1])
    
    with col1:
        # Display main image
        try:
            st.image("nbdcircle.JPG", use_container_width=True, 
                    caption="Non-B DNA structural diversity: From canonical B-form to complex alternative conformations")
        except:
            st.info("Main image not found. Please ensure nbdcircle.JPG is in the working directory.")
    
    with col2:
        st.markdown("""
        <div class="feature-card">
            <h3 style='color: #374151; margin-top: 0; margin-bottom: 20px; font-size: 1.25rem; font-weight: 600;'>
                10 Non-B DNA Classes Detected
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
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>10. Non-B DNA clusters</div>
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # Version and citation info
    st.markdown("---")
    col1, col2, col3 = st.columns(3)
    
    with col1:
        st.markdown("""
        <div class="feature-card">
            <h4>Version Information</h4>
            <p><strong>Version:</strong> 2.0</p>
            <p><strong>Last Updated:</strong> December 2024</p>
            <p><strong>License:</strong> Academic Use</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="feature-card">
            <h4>Citation</h4>
            <p style="font-size: 0.9rem; font-style: italic;">
            Yella, V.R. (2024). NBDFinder: A Comprehensive Computational Framework 
            for Genome-Wide Detection and Analysis of Non-B DNA Structural Motifs. 
            GitHub Repository.
            </p>
        </div>
        """, unsafe_allow_html=True)
    
    with col3:
        st.markdown("""
        <div class="feature-card">
            <h4>Platform Overview</h4>
            <p>High-performance analysis platform for non-canonical DNA structures with:</p>
            <ul style="font-size: 0.9rem;">
                <li>Real-time detection algorithms</li>
                <li>Publication-quality visualizations</li>
                <li>Clinical disease annotations</li>
                <li>Multi-format data export</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)

# =====================================================
# UPLOAD & ANALYZE PAGE
# =====================================================
with tab_dict["Upload & Analyze"]:
    # Subtabs for Upload & Analyze
    upload_subtabs = st.tabs(["Upload Sequence", "Parameter Settings", "Run Analysis"])
    
    # ---- UPLOAD SEQUENCE SUBTAB ----
    with upload_subtabs[0]:
        st.markdown("""
        <div style='background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; margin-bottom: 24px;'>
            <h3 style='color: #374151; margin-top: 0; margin-bottom: 12px;'>Sequence Input</h3>
            <p style='margin-bottom: 0; color: #6b7280; font-size: 0.875rem;'>Choose your preferred method for sequence input. All formats support both single sequences and batch processing.</p>
        </div>
        """, unsafe_allow_html=True)
        
        input_method = st.radio("**Input Method:**", 
                               ["Upload FASTA / Multi-FASTA File", "Paste Sequence(s)", "Example Sequence", "NCBI Fetch"], 
                               horizontal=True)
        
        # Clear session state when switching methods
        if 'last_input_method' not in st.session_state:
            st.session_state.last_input_method = input_method
        elif st.session_state.last_input_method != input_method:
            st.session_state.seqs = []
            st.session_state.names = []
            st.session_state.results = []
            st.session_state.last_input_method = input_method
        
        seqs, names = [], []
        
        if input_method == "Upload FASTA / Multi-FASTA File":
            st.markdown("""
            <div class="feature-card">
                <h4>File Upload</h4>
                <p>Drag & drop your FASTA files here</p>
                <p style="font-size: 0.875rem; color: #6b7280;">Supports .fa, .fasta, .txt files • Single or Multi-FASTA format</p>
            </div>
            """, unsafe_allow_html=True)
            
            fasta_file = st.file_uploader(
                "Choose file", 
                type=['fa', 'fasta', 'txt'],
                help="Upload single or multi-FASTA files for analysis"
            )
            
            if fasta_file:
                try:
                    content = fasta_file.read().decode('utf-8')
                    seqs, names = parse_fasta(content)
                    if seqs:
                        st.success(f"Successfully loaded {len(seqs)} sequence(s)")
                        for i, (name, seq) in enumerate(zip(names, seqs)):
                            with st.expander(f"Sequence {i+1}: {name[:50]}..."):
                                st.text(f"Length: {len(seq)} bp")
                                st.text(f"GC Content: {gc_content(seq):.1f}%")
                                st.text_area("Sequence:", seq[:200] + ("..." if len(seq) > 200 else ""), height=100)
                    else:
                        st.error("❌ No valid sequences found in file")
                except Exception as e:
                    st.error(f"❌ Error reading file: {e}")
        
        elif input_method == "Paste Sequence(s)":
            st.markdown("""
            <div class="feature-card">
                <h4>📝 Manual Sequence Input</h4>
                <p>Paste your sequence(s) in FASTA format or raw sequence</p>
            </div>
            """, unsafe_allow_html=True)
            
            sequence_input = st.text_area(
                "Enter sequence(s):",
                height=200,
                help="Paste FASTA format sequences or raw DNA sequences",
                placeholder=">My_Sequence_1\nATCGATCGATCG...\n>My_Sequence_2\nGCGCGCGCGCGC..."
            )
            
            if sequence_input:
                if sequence_input.startswith('>'):
                    seqs, names = parse_fasta(sequence_input)
                else:
                    # Treat as raw sequence
                    cleaned_seq = ''.join(sequence_input.upper().split())
                    if all(c in 'ATCGRYSWKMBDHVN' for c in cleaned_seq):
                        seqs = [cleaned_seq]
                        names = ["User_Input_Sequence"]
                    else:
                        st.error("❌ Invalid sequence characters detected")
                
                if seqs:
                    st.success(f"Parsed {len(seqs)} sequence(s)")
        
        elif input_method == "Example Sequence":
            st.markdown("""
            <div class="feature-card">
                <h4>Example Datasets</h4>
                <p>Choose from pre-loaded example sequences for testing</p>
            </div>
            """, unsafe_allow_html=True)
            
            selected_example = st.selectbox(
                "Select example sequence:",
                list(EXAMPLE_SEQUENCES.keys()),
                help="Choose an example sequence to demonstrate NBDFinder capabilities"
            )
            
            if selected_example:
                example_data = EXAMPLE_SEQUENCES[selected_example]
                st.info(f"**{selected_example}**: {example_data['description']}")
                seqs = [example_data['sequence']]
                names = [example_data['name']]
                
                with st.expander("View sequence details"):
                    st.text(f"Length: {len(seqs[0])} bp")
                    st.text(f"GC Content: {gc_content(seqs[0]):.1f}%")
                    st.text_area("Sequence:", seqs[0], height=100)
        
        elif input_method == "NCBI Fetch":
            st.markdown("""
            <div class="feature-card">
                <h4>NCBI Database Fetch</h4>
                <p>Retrieve sequences directly from NCBI databases</p>
            </div>
            """, unsafe_allow_html=True)
            
            accession = st.text_input(
                "Enter NCBI Accession Number:",
                help="Format: NG_123456.1, NM_123456.1, etc.",
                placeholder="NG_123456.1"
            )
            
            if accession and st.button("Fetch Sequence"):
                try:
                    with st.spinner("Fetching sequence from NCBI..."):
                        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
                        fasta_content = handle.read()
                        handle.close()
                        
                        seqs, names = parse_fasta(fasta_content)
                        if seqs:
                            st.success(f"Successfully fetched: {names[0]}")
                        else:
                            st.error("❌ No sequence data found")
                except Exception as e:
                    st.error(f"❌ Error fetching sequence: {e}")
        
        # Store sequences in session state
        if seqs:
            st.session_state.seqs = seqs
            st.session_state.names = names
    
    # ---- PARAMETER SETTINGS SUBTAB ----
    with upload_subtabs[1]:
        st.markdown("""
        <div style='background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; margin-bottom: 24px;'>
            <h3 style='color: #374151; margin-top: 0; margin-bottom: 12px;'>Analysis Parameters</h3>
            <p style='margin-bottom: 0; color: #6b7280; font-size: 0.875rem;'>Configure motif detection settings and sensitivity thresholds</p>
        </div>
        """, unsafe_allow_html=True)
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Motif Classes Selection**")
            motif_classes = st.multiselect(
                "Select motif classes to analyze:",
                list(MOTIF_COLORS.keys()),
                default=list(MOTIF_COLORS.keys()),
                help="Choose which non-B DNA motif types to detect"
            )
            
            st.markdown("**Sensitivity Settings**")
            sensitivity = st.select_slider(
                "Detection sensitivity:",
                options=["Low", "Medium", "High", "Maximum"],
                value="High",
                help="Higher sensitivity detects more motifs but may include false positives"
            )
        
        with col2:
            st.markdown("**Advanced Options**")
            with st.expander("Advanced Detection Parameters"):
                min_motif_length = st.number_input("Minimum motif length (bp):", 
                                                 min_value=3, max_value=100, value=10)
                max_motif_length = st.number_input("Maximum motif length (bp):", 
                                                 min_value=10, max_value=1000, value=200)
                overlap_threshold = st.slider("Overlap threshold (%):", 
                                             min_value=0, max_value=100, value=50)
                quality_filter = st.checkbox("Apply quality filtering", value=True)
                
            with st.expander("Performance Settings"):
                parallel_processing = st.checkbox("Enable parallel processing", value=True)
                memory_optimization = st.checkbox("Memory optimization", value=True)
        
        # Store settings in session state
        st.session_state.analysis_settings = {
            'motif_classes': motif_classes,
            'sensitivity': sensitivity,
            'min_motif_length': min_motif_length,
            'max_motif_length': max_motif_length,
            'overlap_threshold': overlap_threshold,
            'quality_filter': quality_filter,
            'parallel_processing': parallel_processing,
            'memory_optimization': memory_optimization
        }
    
    # ---- RUN ANALYSIS SUBTAB ----
    with upload_subtabs[2]:
        st.markdown("""
        <div style='background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; margin-bottom: 24px;'>
            <h3 style='color: #374151; margin-top: 0; margin-bottom: 12px;'>Analysis Execution</h3>
            <p style='margin-bottom: 0; color: #6b7280; font-size: 0.875rem;'>Review settings and launch motif detection analysis</p>
        </div>
        """, unsafe_allow_html=True)
        
        # Settings summary
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Analysis Summary**")
            if st.session_state.seqs:
                st.success(f"{len(st.session_state.seqs)} sequence(s) loaded")
                for i, name in enumerate(st.session_state.names):
                    st.text(f"• {name[:50]}... ({len(st.session_state.seqs[i])} bp)")
            else:
                st.warning("No sequences loaded. Please upload sequences first.")
        
        with col2:
            st.markdown("**Current Settings**")
            if 'analysis_settings' in st.session_state:
                settings = st.session_state.analysis_settings
                st.text(f"• Motif classes: {len(settings.get('motif_classes', []))} selected")
                st.text(f"• Sensitivity: {settings.get('sensitivity', 'High')}")
                st.text(f"• Length range: {settings.get('min_motif_length', 10)}-{settings.get('max_motif_length', 200)} bp")
                st.text(f"• Quality filter: {'✓' if settings.get('quality_filter', True) else '✗'}")
        
        st.markdown("---")
        
        # Analysis button
        if st.session_state.seqs and not st.session_state.get('analysis_running', False):
            if st.button("Start Analysis", type="primary", use_container_width=True):
                st.session_state.analysis_running = True
                st.rerun()
        
        # Analysis progress
        if st.session_state.get('analysis_running', False):
            st.markdown("**Analysis in Progress**")
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            # Simulate analysis process
            for i in range(101):
                time.sleep(0.02)  # Simulate processing time
                progress_bar.progress(i)
                if i < 20:
                    status_text.text("Initializing analysis...")
                elif i < 40:
                    status_text.text("Processing sequences...")
                elif i < 70:
                    status_text.text("Detecting motifs...")
                elif i < 90:
                    status_text.text("Generating results...")
                else:
                    status_text.text("Finalizing analysis...")
            
            # Run actual analysis (placeholder - integrate with existing motif detection)
            try:
                results = []
                for seq, name in zip(st.session_state.seqs, st.session_state.names):
                    # This would call the actual motif detection functions
                    motifs = all_motifs(seq)  # Using existing function
                    results.append({
                        'sequence_name': name,
                        'sequence': seq,
                        'motifs': motifs,
                        'total_motifs': len(motifs),
                        'sequence_length': len(seq)
                    })
                
                st.session_state.results = results
                st.session_state.analysis_running = False
                st.success("Analysis completed successfully!")
                st.balloons()
                
            except Exception as e:
                st.session_state.analysis_running = False
                st.error(f"❌ Analysis failed: {e}")

# =====================================================
# RESULTS PAGE
# =====================================================
with tab_dict["Results"]:
    if not st.session_state.results:
        st.info("No analysis results available. Please run an analysis first.")
    else:
        # Results subtabs
        results_subtabs = st.tabs(["Overview", "Motif Classes", "Genomic Position", "Hybrid & Cluster Analysis", "Statistical Summary"])
        
        # ---- OVERVIEW SUBTAB ----
        with results_subtabs[0]:
            st.markdown("### Analysis Overview")
            
            # Summary metrics
            total_sequences = len(st.session_state.results)
            total_motifs = sum(r['total_motifs'] for r in st.session_state.results)
            
            col1, col2, col3, col4 = st.columns(4)
            with col1:
                st.metric("Total Sequences", total_sequences)
            with col2:
                st.metric("Total Motifs Found", total_motifs)
            with col3:
                st.metric("Average Motifs/Sequence", f"{total_motifs/total_sequences:.1f}" if total_sequences > 0 else "0")
            with col4:
                # Calculate motif types found
                unique_types = set()
                for result in st.session_state.results:
                    for motif in result['motifs']:
                        unique_types.add(motif.get('Motif', 'Unknown'))
                st.metric("Motif Types Found", len(unique_types))
            
            # Results table
            st.markdown("### Results Summary")
            results_data = []
            for result in st.session_state.results:
                results_data.append({
                    'Sequence Name': result['sequence_name'],
                    'Length (bp)': result['sequence_length'],
                    'Total Motifs': result['total_motifs'],
                    'Motifs/kb': f"{(result['total_motifs'] / result['sequence_length'] * 1000):.2f}"
                })
            
            st.dataframe(pd.DataFrame(results_data), use_container_width=True)
        
        # ---- MOTIF CLASSES SUBTAB ----
        with results_subtabs[1]:
            st.markdown("### Motif Classes Analysis")
            
            # Motif type distribution
            motif_counts = {}
            for result in st.session_state.results:
                for motif in result['motifs']:
                    motif_type = motif.get('Motif', 'Unknown')
                    motif_counts[motif_type] = motif_counts.get(motif_type, 0) + 1
            
            if motif_counts:
                # Create bar chart
                fig = px.bar(
                    x=list(motif_counts.keys()),
                    y=list(motif_counts.values()),
                    title="Motif Type Distribution",
                    labels={'x': 'Motif Type', 'y': 'Count'}
                )
                fig.update_layout(height=400)
                st.plotly_chart(fig, use_container_width=True)
                
                # Detailed table
                st.markdown("### Detailed Motif Breakdown")
                motif_df = pd.DataFrame(list(motif_counts.items()), columns=['Motif Type', 'Count'])
                motif_df['Percentage'] = (motif_df['Count'] / motif_df['Count'].sum() * 100).round(2)
                st.dataframe(motif_df, use_container_width=True)
        
        # ---- GENOMIC POSITION SUBTAB ----
        with results_subtabs[2]:
            st.markdown("### Genomic Position Analysis")
            
            # Sequence selection for detailed view
            if len(st.session_state.results) > 1:
                selected_seq = st.selectbox(
                    "Select sequence for position analysis:",
                    range(len(st.session_state.results)),
                    format_func=lambda x: st.session_state.results[x]['sequence_name']
                )
            else:
                selected_seq = 0
            
            result = st.session_state.results[selected_seq]
            
            # Create position plot
            if result['motifs']:
                positions = []
                types = []
                for motif in result['motifs']:
                    positions.append(motif.get('Start', 0))
                    types.append(motif.get('Motif', 'Unknown'))
                
                # Scatter plot of motif positions
                fig = px.scatter(
                    x=positions,
                    y=types,
                    title=f"Motif Positions in {result['sequence_name']}",
                    labels={'x': 'Position (bp)', 'y': 'Motif Type'},
                    height=400
                )
                fig.update_layout(showlegend=False)
                st.plotly_chart(fig, use_container_width=True)
                
                # Position table
                position_df = pd.DataFrame({
                    'Motif Type': types,
                    'Position': positions,
                    'End': [motif.get('End', motif.get('Start', 0)) for motif in result['motifs']]
                })
                st.dataframe(position_df, use_container_width=True)
        
        # ---- HYBRID & CLUSTER ANALYSIS SUBTAB ----
        with results_subtabs[3]:
            st.markdown("### Hybrid & Cluster Analysis")
            st.info("Hybrid motif and cluster analysis functionality will be implemented here.")
            
            # Placeholder for hybrid analysis
            st.markdown("""
            **Coming Soon:**
            - Overlapping motif detection
            - Cluster region identification
            - Hybrid structure analysis
            - Co-occurrence patterns
            """)
        
        # ---- STATISTICAL SUMMARY SUBTAB ----
        with results_subtabs[4]:
            st.markdown("### Statistical Summary")
            
            # Advanced statistics
            if st.session_state.results:
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**Distribution Statistics**")
                    motif_counts_per_seq = [r['total_motifs'] for r in st.session_state.results]
                    st.text(f"Mean motifs per sequence: {np.mean(motif_counts_per_seq):.2f}")
                    st.text(f"Median motifs per sequence: {np.median(motif_counts_per_seq):.2f}")
                    st.text(f"Standard deviation: {np.std(motif_counts_per_seq):.2f}")
                    st.text(f"Min motifs: {np.min(motif_counts_per_seq)}")
                    st.text(f"Max motifs: {np.max(motif_counts_per_seq)}")
                
                with col2:
                    st.markdown("**Sequence Statistics**")
                    seq_lengths = [r['sequence_length'] for r in st.session_state.results]
                    st.text(f"Mean sequence length: {np.mean(seq_lengths):.0f} bp")
                    st.text(f"Median sequence length: {np.median(seq_lengths):.0f} bp")
                    st.text(f"Total analyzed: {np.sum(seq_lengths):,} bp")

# =====================================================
# VISUALIZATION PAGE
# =====================================================
with tab_dict["Visualization"]:
    if not st.session_state.results:
        st.info("No analysis results available for visualization. Please run an analysis first.")
    else:
        # Visualization subtabs
        viz_subtabs = st.tabs(["Motif Distributions", "Publication Figures", "3D/Advanced Views"])
        
        # ---- MOTIF DISTRIBUTIONS SUBTAB ----
        with viz_subtabs[0]:
            st.markdown("### Interactive Motif Distribution Plots")
            
            # Use the existing visualization integration if available
            if PUBLICATION_VIZ_AVAILABLE:
                # Prepare motif data for visualization
                all_motifs_for_viz = []
                for result in st.session_state.results:
                    all_motifs_for_viz.extend(result['motifs'])
                
                # Call the integrated visualization function
                create_nbdfinder_visualization_interface(all_motifs_for_viz, 
                                                       sum(r['sequence_length'] for r in st.session_state.results))
            else:
                st.info("Advanced publication visualizations not available. Basic charts shown below.")
                
                # Basic visualization fallback
                motif_counts = {}
                for result in st.session_state.results:
                    for motif in result['motifs']:
                        motif_type = motif.get('Motif', 'Unknown')
                        motif_counts[motif_type] = motif_counts.get(motif_type, 0) + 1
                
                if motif_counts:
                    fig = px.pie(
                        values=list(motif_counts.values()),
                        names=list(motif_counts.keys()),
                        title="Motif Type Distribution"
                    )
                    st.plotly_chart(fig, use_container_width=True)
        
        # ---- PUBLICATION FIGURES SUBTAB ----
        with viz_subtabs[1]:
            st.markdown("### Publication-Ready Figures")
            
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Export Options:**")
                export_format = st.selectbox("Format:", ["PNG (300 DPI)", "PDF", "SVG"])
                figure_size = st.selectbox("Size:", ["Small (8x6)", "Medium (12x8)", "Large (16x12)"])
                
            with col2:
                st.markdown("**Style Options:**")
                color_scheme = st.selectbox("Color scheme:", ["Default", "Colorblind-friendly", "Grayscale"])
                include_title = st.checkbox("Include title", value=True)
            
            if st.button("Generate Publication Figure"):
                st.info("Publication figure generation functionality will be implemented here.")
        
        # ---- 3D/ADVANCED VIEWS SUBTAB ----
        with viz_subtabs[2]:
            st.markdown("### 3D and Advanced Visualizations")
            
            if ADVANCED_VIZ_AVAILABLE:
                st.info("Advanced 3D visualization functionality will be integrated here.")
            else:
                st.info("3D visualization requires additional dependencies. Please install the advanced visualization package.")

# =====================================================
# CLINICAL/DISEASE ANNOTATION PAGE
# =====================================================
with tab_dict["Clinical/Disease"]:
    # Clinical subtabs
    clinical_subtabs = st.tabs(["Disease Motif Results", "Clinical Summary"])
    
    # ---- DISEASE MOTIF RESULTS SUBTAB ----
    with clinical_subtabs[0]:
        st.markdown("### 🩺 Disease-Associated Motif Analysis")
        
        if not st.session_state.results:
            st.info("No analysis results available. Please run an analysis first.")
        else:
            st.markdown("""
            <div class="feature-card">
                <h4>Disease Annotation Status</h4>
                <p>Clinical disease annotation functionality is being integrated.</p>
                <p><strong>Planned Features:</strong></p>
                <ul>
                    <li>Disease-associated repeat expansions</li>
                    <li>Pathogenic motif classifications</li>
                    <li>Clinical significance scoring</li>
                    <li>Literature references</li>
                </ul>
            </div>
            """, unsafe_allow_html=True)
    
    # ---- CLINICAL SUMMARY SUBTAB ----
    with clinical_subtabs[1]:
        st.markdown("### Clinical Interpretation Summary")
        
        st.markdown("""
        <div class="feature-card">
            <h4>Clinical Data Integration</h4>
            <p>This section will provide:</p>
            <ul>
                <li>Risk assessment for detected motifs</li>
                <li>Links to clinical databases (ClinVar, OMIM)</li>
                <li>Therapeutic implications</li>
                <li>Patient counseling information</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)

# =====================================================
# DOWNLOAD & EXPORT PAGE
# =====================================================
with tab_dict["Download & Export"]:
    st.markdown("### Download Analysis Results")
    
    if not st.session_state.results:
        st.info("No analysis results available for download. Please run an analysis first.")
    else:
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Data Export Options**")
            
            # Results table download
            if st.button("Download Results Table (CSV)"):
                results_data = []
                for result in st.session_state.results:
                    for motif in result['motifs']:
                        results_data.append({
                            'Sequence_Name': result['sequence_name'],
                            'Motif_Type': motif.get('Motif', 'Unknown'),
                            'Start_Position': motif.get('Start', 0),
                            'End_Position': motif.get('End', 0),
                            'Score': motif.get('Score', 0),
                            'Strand': motif.get('Strand', '+')
                        })
                
                df = pd.DataFrame(results_data)
                csv = df.to_csv(index=False)
                st.download_button(
                    label="Download CSV",
                    data=csv,
                    file_name="nbdfinder_results.csv",
                    mime="text/csv"
                )
            
            # Excel export
            if st.button("Download Results (Excel)"):
                st.info("Excel export functionality will be implemented.")
        
        with col2:
            st.markdown("**Visualization Export**")
            
            if st.button("Download Plots (PNG)"):
                st.info("Plot export functionality will be implemented.")
            
            if st.button("📄 Generate Analysis Report"):
                st.info("PDF report generation will be implemented.")
        
        st.markdown("---")
        st.markdown("**Session Data**")
        
        # Session summary
        col1, col2, col3 = st.columns(3)
        with col1:
            st.metric("Sequences Analyzed", len(st.session_state.results))
        with col2:
            st.metric("Total Motifs", sum(r['total_motifs'] for r in st.session_state.results))
        with col3:
            analysis_time = "< 1 min"  # Placeholder
            st.metric("Analysis Time", analysis_time)

# =====================================================
# DOCUMENTATION PAGE
# =====================================================
with tab_dict["Documentation"]:
    # Documentation subtabs
    doc_subtabs = st.tabs(["Motif Definitions", "Detection Algorithms", "Scoring Systems", "User Guide", "References"])
    
    # ---- MOTIF DEFINITIONS SUBTAB ----
    with doc_subtabs[0]:
        st.markdown("### Non-B DNA Motif Definitions")
        
        motif_definitions = {
            "Curved DNA": "DNA sequences that adopt intrinsic curvature due to specific base arrangements",
            "Slipped DNA": "Structures formed by repetitive sequences that can form hairpin loops",
            "Cruciform": "Four-way junctions formed by inverted repeat sequences",
            "R-loop": "Three-stranded nucleic acid structures with RNA-DNA hybrid",
            "Triplex": "Three-stranded DNA structures formed by Hoogsteen base pairing",
            "G-Quadruplex": "Four-stranded structures formed by guanine-rich sequences",
            "i-motif": "Four-stranded structures formed by cytosine-rich sequences",
            "Z-DNA": "Left-handed double helix formed by alternating purine-pyrimidine sequences",
            "Hybrid": "Complex structures with mixed motif characteristics",
            "Clusters": "Regions with high density of multiple non-B DNA motifs"
        }
        
        for motif, definition in motif_definitions.items():
            with st.expander(f"{motif}"):
                st.write(definition)
                st.info(f"Detailed information about {motif} detection algorithms and biological significance.")
    
    # ---- DETECTION ALGORITHMS SUBTAB ----
    with doc_subtabs[1]:
        st.markdown("### Detection Algorithms")
        
        st.markdown("""
        <div class="feature-card">
            <h4>Algorithm Overview</h4>
            <p>NBDFinder employs state-of-the-art computational algorithms for each motif type:</p>
            <ul>
                <li><strong>Sequence-based detection:</strong> Pattern matching with biological constraints</li>
                <li><strong>Energy-based models:</strong> Thermodynamic stability calculations</li>
                <li><strong>Machine learning:</strong> AI-powered classification for complex motifs</li>
                <li><strong>Hybrid approaches:</strong> Combined methods for maximum accuracy</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # ---- SCORING SYSTEMS SUBTAB ----
    with doc_subtabs[2]:
        st.markdown("### Scoring Systems")
        
        st.markdown("""
        <div class="feature-card">
            <h4>Motif Scoring Methodology</h4>
            <p>Each detected motif is assigned confidence scores based on:</p>
            <ul>
                <li><strong>Sequence specificity:</strong> How well the sequence matches known patterns</li>
                <li><strong>Structural stability:</strong> Predicted thermodynamic stability</li>
                <li><strong>Conservation:</strong> Evolutionary conservation across species</li>
                <li><strong>Experimental validation:</strong> Support from literature and databases</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)
    
    # ---- USER GUIDE SUBTAB ----
    with doc_subtabs[3]:
        st.markdown("### User Guide")
        
        st.markdown("""
        <div class="feature-card">
            <h4>Quick Start Guide</h4>
            <ol>
                <li><strong>Upload Sequences:</strong> Use the Upload & Analyze tab to input your sequences</li>
                <li><strong>Configure Parameters:</strong> Set motif classes and sensitivity in Parameter Settings</li>
                <li><strong>Run Analysis:</strong> Execute the analysis and monitor progress</li>
                <li><strong>View Results:</strong> Explore comprehensive results in the Results tab</li>
                <li><strong>Generate Visualizations:</strong> Create publication-ready figures</li>
                <li><strong>Download Data:</strong> Export results in multiple formats</li>
            </ol>
        </div>
        """, unsafe_allow_html=True)
    
    # ---- REFERENCES SUBTAB ----
    with doc_subtabs[4]:
        st.markdown("### Scientific References")
        
        st.markdown("""
        <div class="feature-card">
            <h4>📄 Key Publications</h4>
            <p><strong>Primary Citation:</strong></p>
            <blockquote>
            Yella, V.R. (2024). NBDFinder: A Comprehensive Computational Framework for Genome-Wide Detection and Analysis of Non-B DNA Structural Motifs. GitHub Repository.
            </blockquote>
            
            <p><strong>Related Literature:</strong></p>
            <ul>
                <li>Brascher et al. (2023). Non-B DNA structures and genome instability</li>
                <li>Guiblet et al. (2021). Non-B DNA secondary structures and disease</li>
                <li>Wells et al. (2020). Non-canonical DNA structures in biology</li>
            </ul>
        </div>
        """, unsafe_allow_html=True)

# =====================================================
# SETTINGS PAGE
# =====================================================
with tab_dict["Settings"]:
    st.markdown("### Application Settings & Customization")
    
    # Create subtabs for different settings categories
    settings_subtabs = st.tabs(["Appearance", "Typography", "Analysis Defaults", "Data Export"])
    
    # ---- APPEARANCE SUBTAB ----
    with settings_subtabs[0]:
        st.markdown("#### Theme & Color Settings")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Color Theme**")
            
            theme_mode = st.selectbox(
                "Color theme:",
                ["Auto", "Light", "Dark", "High Contrast"],
                help="Choose your preferred color theme"
            )
            
            background_color = st.color_picker(
                "Background color:",
                value="#FFFFFF",
                help="Set custom background color"
            )
            
            primary_color = st.color_picker(
                "Primary accent color:",
                value="#1E3A8A",
                help="Main accent color for buttons and highlights"
            )
        
        with col2:
            st.markdown("**Text Colors**")
            
            text_color = st.color_picker(
                "Main text color:",
                value="#1E3A8A",
                help="Primary text color"
            )
            
            muted_text_color = st.color_picker(
                "Secondary text color:",
                value="#64748B",
                help="Color for secondary/muted text"
            )
            
            border_color = st.color_picker(
                "Border color:",
                value="#E2E8F0",
                help="Color for borders and dividers"
            )
    
    # ---- TYPOGRAPHY SUBTAB ----
    with settings_subtabs[1]:
        st.markdown("#### Font & Typography Settings")
        
        col1, col2, col3 = st.columns(3)
        
        with col1:
            st.markdown("**Font Family**")
            
            font_family = st.selectbox(
                "Primary font:",
                ["Inter", "Arial", "Helvetica", "Times New Roman", "Georgia", "Verdana", "Roboto", "Open Sans"],
                help="Choose the main font for the application"
            )
            
            fallback_fonts = st.multiselect(
                "Fallback fonts:",
                ["Arial", "Helvetica", "Times New Roman", "Georgia", "Verdana", "sans-serif", "serif"],
                default=["Arial", "sans-serif"],
                help="Backup fonts if primary font fails to load"
            )
        
        with col2:
            st.markdown("**Text Sizes**")
            
            base_font_size = st.slider(
                "Base text size (rem):",
                min_value=0.8,
                max_value=2.0,
                value=1.25,
                step=0.05,
                help="Base font size for regular text"
            )
            
            heading_scale = st.slider(
                "Heading size multiplier:",
                min_value=1.2,
                max_value=3.0,
                value=2.0,
                step=0.1,
                help="How much larger headings should be"
            )
            
            label_font_size = st.slider(
                "Label text size (rem):",
                min_value=0.9,
                max_value=1.5,
                value=1.125,
                step=0.025,
                help="Font size for form labels"
            )
        
        with col3:
            st.markdown("**Text Styling**")
            
            line_height = st.slider(
                "Line height:",
                min_value=1.0,
                max_value=2.5,
                value=1.7,
                step=0.1,
                help="Spacing between lines of text"
            )
            
            font_weight = st.selectbox(
                "Text weight:",
                ["300 (Light)", "400 (Normal)", "500 (Medium)", "600 (Semi-bold)", "700 (Bold)"],
                index=1,
                help="Default font weight"
            )
            
            letter_spacing = st.slider(
                "Letter spacing (em):",
                min_value=-0.05,
                max_value=0.1,
                value=0.0,
                step=0.005,
                help="Space between letters"
            )
    
    # ---- ANALYSIS DEFAULTS SUBTAB ----
    with settings_subtabs[2]:
        st.markdown("#### Default Analysis Settings")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**Detection Parameters**")
            
            default_sensitivity = st.selectbox(
                "Default sensitivity:",
                ["Low", "Medium", "High", "Maximum"],
                index=2,
                help="Default sensitivity for motif detection"
            )
            
            default_motif_classes = st.multiselect(
                "Default motif classes:",
                ["Curved DNA", "Slipped DNA", "Cruciform", "R-loop", "Triplex", 
                 "G-Quadruplex", "i-motif", "Z-DNA", "Hybrid", "Disease-Associated"],
                default=["G-Quadruplex", "i-motif", "Z-DNA"],
                help="Motif classes to detect by default"
            )
            
        with col2:
            st.markdown("**Output Preferences**")
            
            default_export_format = st.selectbox(
                "Default export format:",
                ["Excel (.xlsx)", "CSV", "TSV", "JSON"],
                help="Default format for data export"
            )
            
            show_sequence_details = st.checkbox(
                "Show sequence details by default",
                value=True,
                help="Display sequence information in results"
            )
            
            auto_generate_visualizations = st.checkbox(
                "Auto-generate visualizations",
                value=False,
                help="Automatically create plots after analysis"
            )
    
    # ---- DATA EXPORT SUBTAB ----
    with settings_subtabs[3]:
        st.markdown("#### Data Export Settings")
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**File Naming**")
            
            filename_prefix = st.text_input(
                "Default filename prefix:",
                value="NBDFinder_Results",
                help="Prefix for exported files"
            )
            
            include_timestamp = st.checkbox(
                "Include timestamp in filenames",
                value=True,
                help="Add timestamp to exported filenames"
            )
            
        with col2:
            st.markdown("**Export Options**")
            
            include_metadata = st.checkbox(
                "Include analysis metadata",
                value=True,
                help="Include analysis parameters in export"
            )
            
            compress_exports = st.checkbox(
                "Compress large exports",
                value=False,
                help="Compress files larger than 10MB"
            )
    
    st.markdown("---")
    
    # Action buttons
    col1, col2, col3 = st.columns([1, 1, 1])
    
    with col1:
        if st.button("💾 Save Settings", use_container_width=True):
            # Store settings in session state
            st.session_state.user_settings = {
                'theme_mode': theme_mode,
                'background_color': background_color,
                'primary_color': primary_color,
                'text_color': text_color,
                'muted_text_color': muted_text_color,
                'border_color': border_color,
                'font_family': font_family,
                'fallback_fonts': fallback_fonts,
                'base_font_size': base_font_size,
                'heading_scale': heading_scale,
                'label_font_size': label_font_size,
                'line_height': line_height,
                'font_weight': font_weight,
                'letter_spacing': letter_spacing,
                'default_sensitivity': default_sensitivity,
                'default_motif_classes': default_motif_classes,
                'default_export_format': default_export_format,
                'show_sequence_details': show_sequence_details,
                'auto_generate_visualizations': auto_generate_visualizations,
                'filename_prefix': filename_prefix,
                'include_timestamp': include_timestamp,
                'include_metadata': include_metadata,
                'compress_exports': compress_exports
            }
            st.success("✅ Settings saved successfully!")
    
    with col2:
        if st.button("🔄 Reset to Defaults", use_container_width=True):
            if 'user_settings' in st.session_state:
                del st.session_state.user_settings
            st.info("🔄 Settings reset to default values")
    
    with col3:
        if st.button("📥 Export Settings", use_container_width=True):
            if 'user_settings' in st.session_state:
                import json
                settings_json = json.dumps(st.session_state.user_settings, indent=2)
                st.download_button(
                    label="⬇️ Download settings.json",
                    data=settings_json,
                    file_name="nbdfinder_settings.json",
                    mime="application/json"
                )
            else:
                st.warning("No custom settings to export")
    
    # Live preview section
    st.markdown("---")
    st.markdown("#### 🎨 Live Preview")
    
    # Apply custom CSS based on current settings
    custom_css = f"""
    <style>
    :root {{
        --custom-bg: {background_color};
        --custom-primary: {primary_color};
        --custom-text: {text_color};
        --custom-muted: {muted_text_color};
        --custom-border: {border_color};
        --custom-font: {font_family};
        --custom-base-size: {base_font_size}rem;
        --custom-label-size: {label_font_size}rem;
        --custom-line-height: {line_height};
        --custom-letter-spacing: {letter_spacing}em;
    }}
    
    .preview-card {{
        background: var(--custom-bg);
        border: 2px solid var(--custom-border);
        border-radius: 8px;
        padding: 20px;
        font-family: var(--custom-font), {', '.join(fallback_fonts)};
        color: var(--custom-text);
        font-size: var(--custom-base-size);
        line-height: var(--custom-line-height);
        letter-spacing: var(--custom-letter-spacing);
    }}
    
    .preview-card h3 {{
        color: var(--custom-primary);
        font-size: calc(var(--custom-base-size) * {heading_scale});
        margin-top: 0;
    }}
    
    .preview-card .label {{
        color: var(--custom-muted);
        font-size: var(--custom-label-size);
        font-weight: {font_weight.split()[0]};
    }}
    </style>
    
    <div class="preview-card">
        <h3>Sample Text Preview</h3>
        <div class="label">Label text example</div>
        <p>This is how regular paragraph text will appear with your current settings. The font family is {font_family} with a base size of {base_font_size}rem and line height of {line_height}.</p>
        <p style="color: var(--custom-muted);">This is secondary text using the muted color you selected.</p>
    </div>
    """
    
    st.markdown(custom_css, unsafe_allow_html=True)

# Footer
st.markdown("---")
st.markdown("""
<div style="text-align: center; padding: 20px; color: #6b7280; font-size: 0.875rem;">
    <p><strong>NBDFinder v2.0</strong> | Developed by Dr. Venkata Rajesh Yella | 
    <a href="https://github.com/VRYella/NBDFinder" style="color: #0891b2;">GitHub</a> | 
    <a href="mailto:yvrajesh_bt@kluniversity.in" style="color: #0891b2;">Contact</a></p>
    <p>© 2024 NBDFinder. Academic use license.</p>
</div>
""", unsafe_allow_html=True)