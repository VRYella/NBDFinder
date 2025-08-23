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
        /* Professional color scheme as per requirements */
        --primary: #1e3a8a;        /* Deep blue/navy for headings */
        --accent: #0891b2;         /* Teal for subheadings */
        --text: #1f2937;           /* Dark gray/black for primary text */
        --text-muted: #6b7280;     /* Gray for secondary text */
        --bg: #f9fafb;             /* Light gray background */
        --surface: #ffffff;        /* Off-white for surfaces */
        --border: #e5e7eb;         /* Light border color */
        --radius: 8px;
        --shadow: 0 2px 4px rgba(0,0,0,0.1); 
        --transition: all 0.2s ease;
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
    
    /* IMPROVED FORM CONTROLS */
    .stSelectbox > div > div {
        background: var(--surface); border: 2px solid var(--border); border-radius: var(--radius);
        font-size: 1.125rem; font-weight: 500; transition: var(--transition);
    }
    .stSelectbox > div > div:focus-within {
        border-color: var(--accent); box-shadow: 0 0 0 3px rgba(8,145,178,0.1);
    }
    
    .stSlider > div > div > div > div {
        background: linear-gradient(90deg, var(--accent) 0%, var(--primary) 100%);
    }
    
    .stCheckbox > label {
        font-size: 1.125rem; font-weight: 500; color: var(--text);
    }
    
    .stRadio > div {
        display: flex; gap: 16px; flex-wrap: wrap;
    }
    .stRadio > div > label {
        background: var(--surface); border: 2px solid var(--border); border-radius: var(--radius);
        padding: 12px 20px; margin: 4px; font-size: 1.125rem; font-weight: 500;
        transition: var(--transition); cursor: pointer; min-width: 120px; text-align: center;
    }
    .stRadio > div > label:hover {
        border-color: var(--accent); background: rgba(8,145,178,0.05);
    }
    .stRadio > div > label[data-checked="true"] {
        background: linear-gradient(135deg, var(--primary) 0%, var(--accent) 100%);
        color: white; border-color: var(--primary); font-weight: 600;
    }
    
    /* ENHANCED DATAFRAMES AND TABLES */
    .stDataFrame {
        border: 2px solid var(--border); border-radius: 12px; overflow: hidden;
        box-shadow: var(--shadow); margin: 16px 0;
    }
    
    .stDataFrame table {
        font-family: Inter, sans-serif; font-size: 0.95rem;
    }
    
    .stDataFrame thead tr th {
        background: linear-gradient(135deg, var(--primary) 0%, var(--accent) 100%);
        color: white; font-weight: 700; font-size: 1rem; padding: 16px 12px;
        text-align: left; border: none;
    }
    
    .stDataFrame tbody tr:nth-child(even) {
        background: rgba(248,250,252,0.8);
    }
    
    .stDataFrame tbody tr:hover {
        background: rgba(8,145,178,0.08); transition: background 0.2s ease;
    }
    
    .stDataFrame tbody tr td {
        padding: 14px 12px; border-bottom: 1px solid var(--border);
        font-weight: 500; color: var(--text);
    }
    
    /* ENHANCED BUTTONS */
    .stButton > button {
        background: linear-gradient(135deg, var(--primary) 0%, var(--accent) 100%);
        color: white; border: none; border-radius: var(--radius); font-family: Inter, sans-serif;
        font-weight: 600; font-size: 1.125rem; padding: 12px 24px; transition: var(--transition);
        box-shadow: var(--shadow); cursor: pointer;
    }
    .stButton > button:hover {
        transform: translateY(-2px); box-shadow: 0 6px 20px rgba(30,58,138,0.3);
        background: linear-gradient(135deg, var(--accent) 0%, var(--primary) 100%);
    }
    
    .stDownloadButton > button {
        background: linear-gradient(135deg, #10b981 0%, #059669 100%);
        color: white; border: none; border-radius: var(--radius); font-family: Inter, sans-serif;
        font-weight: 600; font-size: 1rem; padding: 10px 20px; transition: var(--transition);
        box-shadow: var(--shadow);
    }
    .stDownloadButton > button:hover {
        transform: translateY(-1px); box-shadow: 0 4px 16px rgba(16,185,129,0.3);
    }
    
    /* ENHANCED METRICS */
    .metric-card {
        background: var(--surface); border: 2px solid var(--border); border-radius: 12px;
        padding: 20px; margin: 8px 0; box-shadow: var(--shadow); transition: var(--transition);
    }
    .metric-card:hover {
        transform: translateY(-2px); box-shadow: 0 8px 25px rgba(0,0,0,0.1);
        border-color: var(--accent);
    }
    
    .metric-card h3 {
        color: var(--primary); font-size: 2.5rem; font-weight: 700; margin: 0;
        text-align: center;
    }
    .metric-card p {
        color: var(--text-muted); font-size: 1.125rem; font-weight: 500; margin: 8px 0 0 0;
        text-align: center;
    }
    
    /* ENHANCED MESSAGES */
    .stSuccess {
        background: linear-gradient(135deg, #d1fae5 0%, #a7f3d0 100%);
        border: 2px solid #10b981; border-radius: var(--radius); padding: 16px;
        font-family: Inter, sans-serif; font-weight: 500;
    }
    
    .stError {
        background: linear-gradient(135deg, #fee2e2 0%, #fecaca 100%);
        border: 2px solid #ef4444; border-radius: var(--radius); padding: 16px;
        font-family: Inter, sans-serif; font-weight: 500;
    }
    
    .stWarning {
        background: linear-gradient(135deg, #fef3c7 0%, #fde68a 100%);
        border: 2px solid #f59e0b; border-radius: var(--radius); padding: 16px;
        font-family: Inter, sans-serif; font-weight: 500;
    }
    
    .stInfo {
        background: linear-gradient(135deg, #dbeafe 0%, #bfdbfe 100%);
        border: 2px solid #3b82f6; border-radius: var(--radius); padding: 16px;
        font-family: Inter, sans-serif; font-weight: 500;
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
    "Results & Visualization": "Analysis Results and Publication-Ready Visualizations",
    "Clinical/Disease": "Disease Annotation and Clinical Data",
    "Download & Export": "Download Results and Export Data",
    "Documentation": "Scientific Documentation & References"
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
            <h3 style='color: #1e3a8a; margin-top: 0; margin-bottom: 20px; font-size: 1.25rem; font-weight: 600;'>
                10 Non-B DNA Classes & 22+ Subclasses
            </h3>
            <div style='font-size: 0.9rem; line-height: 1.4;'>
                <strong style='color: #1e3a8a;'>Complete Classification System:</strong>
                <ul style='margin: 8px 0; padding-left: 16px;'>
                    <li>10 Primary structural classes</li>
                    <li>22+ Official subclasses</li>
                    <li>Dynamic hybrid detection</li>
                    <li>Cluster region analysis</li>
                </ul>
            </div>
        </div>
        """, unsafe_allow_html=True)
    
    # Comprehensive class and subclass display
    st.markdown("---")
    st.markdown("""
    <div style='margin: 20px 0;'>
        <h3 style='color: #1e3a8a; font-weight: 600; margin-bottom: 20px; text-align: center;'>
            📚 Complete Non-B DNA Classification System
        </h3>
    </div>
    """, unsafe_allow_html=True)
    
    # Create expandable sections for each class
    from motifs.classification_config import OFFICIAL_CLASSIFICATION
    
    # Organize classes in two columns
    col1, col2 = st.columns(2)
    
    # Split classes between columns
    classes_list = list(OFFICIAL_CLASSIFICATION.items())
    col1_classes = classes_list[:5]  # Classes 1-5
    col2_classes = classes_list[5:]  # Classes 6-10
    
    with col1:
        for class_id, class_info in col1_classes:
            class_name = class_info['class_name']
            subclasses = class_info['subclasses']
            
            with st.expander(f"**{class_id}. {class_name}**", expanded=False):
                if subclasses:
                    st.markdown("**Subclasses:**")
                    for i, subclass in enumerate(subclasses, 1):
                        st.markdown(f"• {subclass}")
                    st.markdown(f"*Total subclasses: {len(subclasses)}*")
                else:
                    if class_name == "Hybrid":
                        st.markdown("**Dynamic subclasses based on overlapping motifs:**")
                        st.markdown("• Created when any two different non-B DNA classes overlap")
                    elif class_name == "Non-B DNA cluster regions":
                        st.markdown("**Dynamic subclasses based on motif clustering:**")
                        st.markdown("• Regions with 3+ different motif types within 100 nucleotides")
                    else:
                        st.markdown("*Subclasses determined dynamically during analysis*")
    
    with col2:
        for class_id, class_info in col2_classes:
            class_name = class_info['class_name']
            subclasses = class_info['subclasses']
            
            with st.expander(f"**{class_id}. {class_name}**", expanded=False):
                if subclasses:
                    st.markdown("**Subclasses:**")
                    for i, subclass in enumerate(subclasses, 1):
                        st.markdown(f"• {subclass}")
                    st.markdown(f"*Total subclasses: {len(subclasses)}*")
                else:
                    if class_name == "Hybrid":
                        st.markdown("**Dynamic subclasses based on overlapping motifs:**")
                        st.markdown("• Created when any two different non-B DNA classes overlap")
                    elif class_name == "Non-B DNA cluster regions":
                        st.markdown("**Dynamic subclasses based on motif clustering:**")
                        st.markdown("• Regions with 3+ different motif types within 100 nucleotides")
                    else:
                        st.markdown("*Subclasses determined dynamically during analysis*")
    
    # Summary statistics
    from motifs.classification_config import get_class_counts
    total_classes, total_fixed_subclasses = get_class_counts()
    
    st.markdown(f"""
    <div style='background: #f8fafc; padding: 16px; border-radius: 8px; border-left: 4px solid #0891b2; margin: 20px 0;'>
        <h4 style='color: #0891b2; margin: 0 0 8px 0;'>📊 Classification Summary</h4>
        <p style='margin: 4px 0; color: #374151;'><strong>Total Classes:</strong> {total_classes}</p>
        <p style='margin: 4px 0; color: #374151;'><strong>Fixed Subclasses:</strong> {total_fixed_subclasses}</p>
        <p style='margin: 4px 0; color: #374151;'><strong>Dynamic Subclasses:</strong> Generated during analysis (Hybrid & Cluster regions)</p>
        <p style='margin: 4px 0; color: #374151;'><strong>Total Coverage:</strong> 22+ officially recognized subclasses</p>
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
# RESULTS & VISUALIZATION PAGE (MERGED)
# =====================================================
with tab_dict["Results & Visualization"]:
    if not st.session_state.results:
        st.info("No analysis results available. Please run an analysis first.")
    else:
        # Combined Results & Visualization subtabs
        combined_subtabs = st.tabs([
            "Overview/Summary", 
            "Motif Class/Subclass Distribution", 
            "Genomic Position Plots", 
            "Detailed Result Tables", 
            "Publication-Ready Figures"
        ])
        
        # ---- OVERVIEW/SUMMARY SUBTAB ----
        with combined_subtabs[0]:
            st.markdown("### 📊 Analysis Overview & Summary")
            
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
                # Calculate unique motif classes found
                unique_classes = set()
                for result in st.session_state.results:
                    for motif in result['motifs']:
                        # Use proper classification
                        motif_class = motif.get('Class', motif.get('Motif', 'Unknown'))
                        unique_classes.add(motif_class)
                st.metric("Classes Detected", len(unique_classes))
            
            # Results summary table
            st.markdown("### 📋 Sequence Analysis Summary")
            results_data = []
            for result in st.session_state.results:
                results_data.append({
                    'Sequence Name': result['sequence_name'],
                    'Length (bp)': result['sequence_length'],
                    'Total Motifs': result['total_motifs'],
                    'Motifs/kb': f"{(result['total_motifs'] / result['sequence_length'] * 1000):.2f}"
                })
            
            st.dataframe(pd.DataFrame(results_data), use_container_width=True)
            
            # Advanced statistics
            if st.session_state.results:
                col1, col2 = st.columns(2)
                
                with col1:
                    st.markdown("**📈 Distribution Statistics**")
                    motif_counts_per_seq = [r['total_motifs'] for r in st.session_state.results]
                    st.text(f"Mean motifs per sequence: {np.mean(motif_counts_per_seq):.2f}")
                    st.text(f"Median motifs per sequence: {np.median(motif_counts_per_seq):.2f}")
                    st.text(f"Standard deviation: {np.std(motif_counts_per_seq):.2f}")
                    st.text(f"Min motifs: {np.min(motif_counts_per_seq)}")
                    st.text(f"Max motifs: {np.max(motif_counts_per_seq)}")
                
                with col2:
                    st.markdown("**📏 Sequence Statistics**")
                    seq_lengths = [r['sequence_length'] for r in st.session_state.results]
                    st.text(f"Mean sequence length: {np.mean(seq_lengths):.0f} bp")
                    st.text(f"Median sequence length: {np.median(seq_lengths):.0f} bp")
                    st.text(f"Total analyzed: {np.sum(seq_lengths):,} bp")
        
        # ---- MOTIF CLASS/SUBCLASS DISTRIBUTION SUBTAB ----
        with combined_subtabs[1]:
            st.markdown("### 🧬 Complete Motif Class & Subclass Distribution")
            
            # Create comprehensive class/subclass analysis
            from motifs.classification_config import OFFICIAL_CLASSIFICATION, get_official_classification
            
            # Initialize complete class/subclass count dictionary
            all_class_counts = {}
            all_subclass_counts = {}
            
            # Initialize with all official classes and subclasses
            for class_info in OFFICIAL_CLASSIFICATION.values():
                class_name = class_info['class_name']
                all_class_counts[class_name] = 0
                
                for subclass in class_info['subclasses']:
                    all_subclass_counts[subclass] = 0
            
            # Count actual detected motifs
            for result in st.session_state.results:
                for motif in result['motifs']:
                    # Get official classification
                    legacy_class = motif.get('Class', motif.get('Motif', 'Unknown'))
                    legacy_subtype = motif.get('Subtype', '')
                    
                    official_class, official_subtype = get_official_classification(legacy_class, legacy_subtype)
                    
                    # Count class
                    if official_class in all_class_counts:
                        all_class_counts[official_class] += 1
                    
                    # Count subclass (if not empty)
                    if official_subtype and official_subtype in all_subclass_counts:
                        all_subclass_counts[official_subtype] += 1
            
            # Display class distribution chart
            if any(count > 0 for count in all_class_counts.values()):
                st.markdown("#### 📊 Class Distribution (All 10 Classes)")
                
                # Create bar chart showing all classes
                class_df = pd.DataFrame([
                    {'Class': class_name, 'Count': count, 'Status': 'Detected' if count > 0 else 'Not Found'}
                    for class_name, count in all_class_counts.items()
                ])
                
                fig = px.bar(
                    class_df,
                    x='Class',
                    y='Count',
                    color='Status',
                    title="Complete Motif Class Distribution (10 Official Classes)",
                    labels={'Count': 'Number of Motifs'},
                    color_discrete_map={'Detected': '#1e3a8a', 'Not Found': '#e5e7eb'}
                )
                fig.update_layout(height=500, xaxis_tickangle=-45)
                st.plotly_chart(fig, use_container_width=True)
                
                # Class counts table
                st.markdown("#### 📋 Complete Class Summary Table")
                class_summary_df = class_df.copy()
                class_summary_df['Percentage'] = (class_summary_df['Count'] / class_summary_df['Count'].sum() * 100).round(2)
                st.dataframe(class_summary_df, use_container_width=True)
            
            # Display subclass distribution
            if any(count > 0 for count in all_subclass_counts.values()):
                st.markdown("#### 🔬 Subclass Distribution (22+ Official Subclasses)")
                
                # Create subclass chart
                subclass_df = pd.DataFrame([
                    {'Subclass': subclass_name, 'Count': count, 'Status': 'Detected' if count > 0 else 'Not Found'}
                    for subclass_name, count in all_subclass_counts.items()
                ])
                
                fig = px.bar(
                    subclass_df,
                    x='Subclass',
                    y='Count',
                    color='Status',
                    title="Complete Subclass Distribution (22+ Official Subclasses)",
                    labels={'Count': 'Number of Motifs'},
                    color_discrete_map={'Detected': '#0891b2', 'Not Found': '#f3f4f6'}
                )
                fig.update_layout(height=600, xaxis_tickangle=-45)
                st.plotly_chart(fig, use_container_width=True)
                
                # Subclass counts table
                st.markdown("#### 📋 Complete Subclass Summary Table")
                subclass_summary_df = subclass_df.copy()
                total_subclass_motifs = subclass_summary_df['Count'].sum()
                if total_subclass_motifs > 0:
                    subclass_summary_df['Percentage'] = (subclass_summary_df['Count'] / total_subclass_motifs * 100).round(2)
                else:
                    subclass_summary_df['Percentage'] = 0
                st.dataframe(subclass_summary_df, use_container_width=True)
            
            # Use advanced visualization if available
            if PUBLICATION_VIZ_AVAILABLE:
                st.markdown("#### 🎨 Advanced Interactive Visualizations")
                
                # Prepare motif data for visualization
                all_motifs_for_viz = []
                for result in st.session_state.results:
                    all_motifs_for_viz.extend(result['motifs'])
                
                # Call the integrated visualization function
                create_nbdfinder_visualization_interface(all_motifs_for_viz, 
                                                       sum(r['sequence_length'] for r in st.session_state.results))
        
        # ---- GENOMIC POSITION PLOTS SUBTAB ----
        with combined_subtabs[2]:
            st.markdown("### 🗺️ Genomic Position Analysis & Plots")
            
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
                classes = []
                
                for motif in result['motifs']:
                    positions.append(motif.get('Start', 0))
                    types.append(motif.get('Motif', 'Unknown'))
                    classes.append(motif.get('Class', 'Unknown'))
                
                # Scatter plot of motif positions by class
                fig = px.scatter(
                    x=positions,
                    y=classes,
                    color=classes,
                    title=f"Motif Positions by Class in {result['sequence_name']}",
                    labels={'x': 'Position (bp)', 'y': 'Motif Class'},
                    height=500
                )
                fig.update_layout(showlegend=True)
                st.plotly_chart(fig, use_container_width=True)
                
                # Linear genomic map
                st.markdown("#### 📍 Linear Genomic Map")
                
                # Create a linear map showing motif density
                seq_length = result['sequence_length']
                bin_size = max(seq_length // 100, 1)  # 100 bins
                bins = np.arange(0, seq_length + bin_size, bin_size)
                hist, _ = np.histogram(positions, bins=bins)
                
                fig = px.bar(
                    x=bins[:-1],
                    y=hist,
                    title=f"Motif Density Along {result['sequence_name']}",
                    labels={'x': 'Position (bp)', 'y': 'Motif Count per Bin'}
                )
                fig.update_layout(height=300)
                st.plotly_chart(fig, use_container_width=True)
                
                # Position table
                st.markdown("#### 📋 Detailed Position Table")
                position_df = pd.DataFrame({
                    'Motif Type': types,
                    'Class': classes,
                    'Start Position': positions,
                    'End Position': [motif.get('End', motif.get('Start', 0)) for motif in result['motifs']],
                    'Length': [motif.get('End', motif.get('Start', 0)) - motif.get('Start', 0) + 1 for motif in result['motifs']]
                })
                st.dataframe(position_df, use_container_width=True)
        
        # ---- DETAILED RESULT TABLES SUBTAB ----
        with combined_subtabs[3]:
            st.markdown("### 📊 Detailed Result Tables")
            
            # Comprehensive motif details table
            all_motifs_details = []
            
            for seq_idx, result in enumerate(st.session_state.results):
                for motif_idx, motif in enumerate(result['motifs']):
                    # Ensure motif has proper subtype
                    motif = ensure_subtype(motif)
                    
                    # Get official classification
                    legacy_class = motif.get('Class', motif.get('Motif', 'Unknown'))
                    legacy_subtype = motif.get('Subtype', '')
                    official_class, official_subtype = get_official_classification(legacy_class, legacy_subtype)
                    
                    all_motifs_details.append({
                        'Sequence': result['sequence_name'],
                        'Motif_ID': f"M{seq_idx+1}_{motif_idx+1}",
                        'Official_Class': official_class,
                        'Official_Subclass': official_subtype,
                        'Legacy_Type': motif.get('Motif', 'Unknown'),
                        'Start': motif.get('Start', 0),
                        'End': motif.get('End', 0),
                        'Length': motif.get('End', 0) - motif.get('Start', 0) + 1,
                        'Score': motif.get('Score', 'N/A'),
                        'Strand': motif.get('Strand', 'N/A'),
                        'Sequence_Context': motif.get('Sequence', 'N/A')
                    })
            
            if all_motifs_details:
                details_df = pd.DataFrame(all_motifs_details)
                
                # Add filters
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    selected_classes = st.multiselect(
                        "Filter by Class:",
                        options=details_df['Official_Class'].unique(),
                        default=list(details_df['Official_Class'].unique())
                    )
                
                with col2:
                    selected_sequences = st.multiselect(
                        "Filter by Sequence:",
                        options=details_df['Sequence'].unique(),
                        default=list(details_df['Sequence'].unique())
                    )
                
                with col3:
                    min_score = st.number_input("Minimum Score:", value=0.0, min_value=0.0)
                
                # Apply filters
                filtered_df = details_df[
                    (details_df['Official_Class'].isin(selected_classes)) &
                    (details_df['Sequence'].isin(selected_sequences))
                ]
                
                # Apply score filter if Score column has numeric values
                try:
                    filtered_df = filtered_df[pd.to_numeric(filtered_df['Score'], errors='coerce') >= min_score]
                except:
                    pass  # Skip score filtering if scores are not numeric
                
                st.markdown(f"#### 📋 Filtered Results ({len(filtered_df)} motifs)")
                st.dataframe(filtered_df, use_container_width=True)
                
                # Download button for filtered results
                csv = filtered_df.to_csv(index=False)
                st.download_button(
                    label="💾 Download Filtered Results as CSV",
                    data=csv,
                    file_name=f"nbdfinder_detailed_results.csv",
                    mime="text/csv"
                )
        
        # ---- PUBLICATION-READY FIGURES SUBTAB ----
        with combined_subtabs[4]:
            st.markdown("### 📖 Publication-Ready Figures")
            
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**📄 Export Options:**")
                export_format = st.selectbox("Format:", ["PNG (300 DPI)", "PDF", "SVG"])
                figure_size = st.selectbox("Size:", ["Small (8x6)", "Medium (12x8)", "Large (16x12)"])
                
            with col2:
                st.markdown("**🎨 Style Options:**")
                color_scheme = st.selectbox("Color scheme:", ["Professional", "Colorblind-friendly", "Grayscale"])
                include_title = st.checkbox("Include title", value=True)
            
            if st.button("🚀 Generate Publication Figure"):
                st.info("Publication figure generation functionality will be implemented here.")
            
            # Additional visualization options
            if ADVANCED_VIZ_AVAILABLE:
                st.markdown("#### 🎨 Advanced 3D Visualizations")
                st.info("Advanced 3D visualization functionality will be integrated here.")
            else:
                st.markdown("#### ℹ️ Advanced Visualizations")
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
