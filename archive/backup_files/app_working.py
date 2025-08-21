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

# ---------- ENHANCED PAGE CONFIGURATION ----------
st.set_page_config(
    page_title="NBDFinder - Advanced Non-B DNA Analysis Platform",
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


# ---------- PREMIUM TOP-1-PERCENTILE DESIGN SYSTEM ----------
st.markdown("""
    <style>
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@300;400;500;600;700;800;900&family=JetBrains+Mono:wght@400;500;600&display=swap');
    
    :root {
        --primary-gradient: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
        --secondary-gradient: linear-gradient(135deg, #f093fb 0%, #f5576c 100%);
        --success-gradient: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%);
        --warning-gradient: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%);
        --glass-bg: rgba(255, 255, 255, 0.08);
        --glass-border: rgba(255, 255, 255, 0.2);
        --shadow-soft: 0 8px 32px rgba(31, 38, 135, 0.15);
        --shadow-strong: 0 16px 64px rgba(31, 38, 135, 0.25);
        --text-primary: #2d3748;
        --text-secondary: #4a5568;
        --border-radius: 16px;
        --transition: all 0.3s cubic-bezier(0.4, 0, 0.2, 1);
    }
    
    body, [data-testid="stAppViewContainer"], .main {
        background: linear-gradient(135deg, #f7fafc 0%, #edf2f7 25%, #e2e8f0 50%, #cbd5e0 75%, #a0aec0 100%) !important;
        font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif !important;
        color: var(--text-primary) !important;
    }
    /* PREMIUM TAB SYSTEM - GLASSMORPHISM DESIGN */
    .stTabs [data-baseweb="tab-list"] {
        background: linear-gradient(135deg, rgba(255, 255, 255, 0.25) 0%, rgba(255, 255, 255, 0.1) 100%) !important;
        backdrop-filter: blur(20px) !important;
        border: 1px solid var(--glass-border) !important;
        border-radius: 24px !important;
        padding: 8px !important;
        margin-bottom: 32px !important;
        box-shadow: var(--shadow-soft) !important;
        justify-content: space-evenly !important;
    }
    .stTabs [data-baseweb="tab"] {
        background: transparent !important;
        border: none !important;
        border-radius: 16px !important;
        font-family: 'Inter', sans-serif !important;
        font-weight: 600 !important;
        font-size: 0.95rem !important;
        color: var(--text-secondary) !important;
        padding: 12px 24px !important;
        margin: 4px !important;
        transition: var(--transition) !important;
        position: relative !important;
        overflow: hidden !important;
    }
    .stTabs [data-baseweb="tab"]:hover {
        background: rgba(255, 255, 255, 0.2) !important;
        color: var(--text-primary) !important;
        transform: translateY(-2px) !important;
        box-shadow: 0 8px 25px rgba(0, 0, 0, 0.1) !important;
    }
    .stTabs [aria-selected="true"] {
        background: linear-gradient(135deg, #667eea 0%, #764ba2 100%) !important;
        color: white !important;
        transform: translateY(-3px) !important;
        box-shadow: 0 12px 30px rgba(102, 126, 234, 0.3) !important;
        font-weight: 700 !important;
    }
    .stTabs [aria-selected="true"]::before {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background: linear-gradient(135deg, rgba(255, 255, 255, 0.2) 0%, transparent 50%) !important;
        pointer-events: none;
    }
    /* PREMIUM CONTENT CONTAINERS */
    .stContainer, .main .block-container {
        background: linear-gradient(135deg, rgba(255, 255, 255, 0.95) 0%, rgba(255, 255, 255, 0.8) 100%) !important;
        backdrop-filter: blur(20px) !important;
        border-radius: var(--border-radius) !important;
        border: 1px solid var(--glass-border) !important;
        box-shadow: var(--shadow-soft) !important;
        padding: 32px !important;
        margin: 24px 0 !important;
        position: relative !important;
        overflow: hidden !important;
    }
    .stContainer::before, .main .block-container::before {
        content: '';
        position: absolute;
        top: 0;
        left: 0;
        right: 0;
        bottom: 0;
        background: linear-gradient(135deg, rgba(255, 255, 255, 0.1) 0%, transparent 50%) !important;
        pointer-events: none;
    }
        background: linear-gradient(135deg, #e3f2fd 0%, #f3e5f5 25%, #e8f5e8 50%, #fff3e0 75%, #fce4ec 100%) !important;
        box-shadow: 0 6px 24px rgba(13, 71, 161, 0.25);
        margin-bottom: 0;
        border-radius: 16px 16px 0 0;
    }
    .stTabs [data-baseweb="tab"] {
        font-size: 1.85rem !important;
        font-weight: 900 !important;
        flex: 1 1 0%;
        min-width: 0 !important;
        padding: 22px 15px 22px 15px !important;
        text-align: center;
        color: #0d47a1 !important;
        background: linear-gradient(135deg, #ffffff 0%, #f8fdff 50%, #f0f9ff 100%) !important;
        border-right: 3px solid #e3f2fd !important;
        letter-spacing: 0.08em;
        text-shadow: 0 2px 4px rgba(13, 71, 161, 0.15);
        border-radius: 12px 12px 0 0;
        margin: 6px 3px 0 3px;
        transition: all 0.4s cubic-bezier(0.4, 0, 0.2, 1);
    }
    .stTabs [data-baseweb="tab"]:hover {
        background: linear-gradient(135deg, #e3f2fd 0%, #f3e5f5 50%, #e8f5e8 100%) !important;
        color: #0d47a1 !important;
        transform: translateY(-4px);
        box-shadow: 0 8px 32px rgba(13, 71, 161, 0.3);
        font-size: 1.9rem !important;
    }
    .stTabs [aria-selected="true"] {
        color: #ffffff !important;
        border-bottom: 8px solid #ff6f00 !important;
        background: linear-gradient(135deg, #0d47a1 0%, #1565c0 30%, #1976d2 60%, #42a5f5 100%) !important;
        box-shadow: 0 12px 40px rgba(13, 71, 161, 0.4);
        transform: translateY(-6px);
        font-size: 2.0rem !important;
        text-shadow: 0 3px 6px rgba(0, 0, 0, 0.4);
        border-left: 4px solid #ff6f00 !important;
        border-right: 4px solid #ff6f00 !important;
    }
    .stTabs [data-baseweb="tab"]:last-child {
        border-right: none !important;
    }
    /* PREMIUM TYPOGRAPHY SYSTEM */
    h1, h2, h3, h4, h5, h6 {
        font-family: 'Inter', sans-serif !important;
        font-weight: 800 !important;
        letter-spacing: -0.025em !important;
        line-height: 1.2 !important;
        margin-top: 2rem !important;
        margin-bottom: 1rem !important;
        background: var(--primary-gradient) !important;
        -webkit-background-clip: text !important;
        -webkit-text-fill-color: transparent !important;
        background-clip: text !important;
        position: relative !important;
    }
    h1 { 
        font-size: 3.5rem !important; 
        font-weight: 900 !important;
        letter-spacing: -0.05em !important;
    }
    h2 { 
        font-size: 2.75rem !important; 
        font-weight: 800 !important;
    }
    h3 { 
        font-size: 2.25rem !important; 
        font-weight: 700 !important;
    }
    h4 { 
        font-size: 1.875rem !important; 
        font-weight: 600 !important;
    }
    h5 { 
        font-size: 1.5rem !important; 
        font-weight: 600 !important;
    }
    h6 { 
        font-size: 1.25rem !important; 
        font-weight: 600 !important;
    }
    /* PREMIUM TEXT AND INTERACTIVE ELEMENTS */
    .stMarkdown, .markdown-text-container, .stText, p, span, label, input {
        font-family: 'Inter', sans-serif !important;
        font-size: 1rem !important;
        line-height: 1.7 !important;
        color: var(--text-primary) !important;
        font-weight: 400 !important;
    }
    
    /* PREMIUM BUTTON SYSTEM */
    .stButton>button {
        font-family: 'Inter', sans-serif !important;
        font-weight: 600 !important;
        font-size: 1rem !important;
        padding: 12px 32px !important;
        background: var(--primary-gradient) !important;
        color: white !important;
        border: none !important;
        border-radius: var(--border-radius) !important;
        box-shadow: var(--shadow-soft) !important;
        transition: var(--transition) !important;
        position: relative !important;
        overflow: hidden !important;
        cursor: pointer !important;
    }
    .stButton>button:hover {
        transform: translateY(-2px) !important;
        box-shadow: var(--shadow-strong) !important;
        background: linear-gradient(135deg, #764ba2 0%, #667eea 100%) !important;
    }
    .stButton>button:active {
        transform: translateY(0) !important;
        box-shadow: var(--shadow-soft) !important;
    }
    .stButton>button::before {
        content: '';
        position: absolute;
        top: 0;
        left: -100%;
        width: 100%;
        height: 100%;
        background: linear-gradient(90deg, transparent, rgba(255, 255, 255, 0.2), transparent);
        transition: var(--transition);
    }
    .stButton>button:hover::before {
        left: 100%;
    }
    /* PREMIUM DATA TABLES */
    .stDataFrame, .stTable {
        font-family: 'Inter', sans-serif !important;
        border-radius: var(--border-radius) !important;
        overflow: hidden !important;
        box-shadow: var(--shadow-soft) !important;
        border: 1px solid var(--glass-border) !important;
        background: rgba(255, 255, 255, 0.95) !important;
        backdrop-filter: blur(10px) !important;
    }
    .stDataFrame th, .stTable th {
        background: var(--primary-gradient) !important;
        color: white !important;
        font-weight: 700 !important;
        font-size: 0.875rem !important;
        padding: 16px 12px !important;
        text-align: center !important;
        border: none !important;
        letter-spacing: 0.05em !important;
        text-transform: uppercase !important;
    }
    .stDataFrame td, .stTable td {
        padding: 12px !important;
        text-align: center !important;
        font-size: 0.875rem !important;
        font-weight: 500 !important;
        border-bottom: 1px solid rgba(0, 0, 0, 0.05) !important;
        transition: var(--transition) !important;
    }
    .stDataFrame tr:hover td, .stTable tr:hover td {
        background: rgba(102, 126, 234, 0.05) !important;
    }
    .stDataFrame tr:nth-child(even), .stTable tr:nth-child(even) {
        background: rgba(248, 250, 252, 0.8) !important;
    }
    .stDataFrame tr:nth-child(odd), .stTable tr:nth-child(odd) {
        background: rgba(255, 255, 255, 0.9) !important;
    }
    
    /* Enhanced confidence level styling with vibrant color coding */
    .stDataFrame .confidence-optimal {
        background: linear-gradient(135deg, #2e7d32, #43a047, #66bb6a) !important;
        color: white !important;
        font-weight: 800 !important;
        border-radius: 8px !important;
        padding: 6px 12px !important;
        box-shadow: 0 3px 8px rgba(46, 125, 50, 0.4);
        text-shadow: 0 1px 2px rgba(0, 0, 0, 0.3);
    }
    
    .stDataFrame .confidence-high {
        background: linear-gradient(135deg, #1565c0, #1976d2, #42a5f5) !important;
        color: white !important;
        font-weight: 800 !important;
        border-radius: 8px !important;
        padding: 6px 12px !important;
        box-shadow: 0 3px 8px rgba(21, 101, 192, 0.4);
        text-shadow: 0 1px 2px rgba(0, 0, 0, 0.3);
    }
    
    .stDataFrame .confidence-moderate {
        background: linear-gradient(135deg, #f57f17, #ffa726, #ffcc02) !important;
        color: white !important;
        font-weight: 800 !important;
        border-radius: 8px !important;
        padding: 6px 12px !important;
        box-shadow: 0 3px 8px rgba(245, 127, 23, 0.4);
        text-shadow: 0 1px 2px rgba(0, 0, 0, 0.3);
    }
    
    .stDataFrame .confidence-low {
        background: linear-gradient(135deg, #d32f2f, #f44336, #ef5350) !important;
        color: white !important;
        font-weight: 800 !important;
        border-radius: 8px !important;
        padding: 6px 12px !important;
        box-shadow: 0 3px 8px rgba(211, 47, 47, 0.4);
        text-shadow: 0 1px 2px rgba(0, 0, 0, 0.3);
    }
    /* Fix dropdown text overlap issues */
    .stSelectbox > div > div {
        min-width: 200px !important;
        padding-right: 40px !important;
    }
    .stSelectbox > div > div > div {
        white-space: nowrap !important;
        overflow: hidden !important;
        text-overflow: ellipsis !important;
        padding-right: 30px !important;
    }
    /* Fix data table dropdown menus */
    .stDataFrame div[data-testid="stTable"] button {
        min-width: auto !important;
        padding-right: 25px !important;
    }
    .stDataFrame .dropdown-menu {
        min-width: 150px !important;
    }
    /* Fix selectbox arrow positioning */
    .stSelectbox [data-baseweb="select"] > div:last-child {
        padding-left: 8px !important;
        min-width: 24px !important;
    }
    /* PREMIUM FORM CONTROLS */
    .stSelectbox > div, .stMultiSelect > div, .stTextInput > div, .stTextArea > div {
        background: rgba(255, 255, 255, 0.9) !important;
        backdrop-filter: blur(10px) !important;
        border: 1px solid var(--glass-border) !important;
        border-radius: var(--border-radius) !important;
        box-shadow: var(--shadow-soft) !important;
        transition: var(--transition) !important;
    }
    .stSelectbox > div:hover, .stMultiSelect > div:hover, .stTextInput > div:hover, .stTextArea > div:hover {
        border-color: rgba(102, 126, 234, 0.5) !important;
        box-shadow: 0 8px 25px rgba(102, 126, 234, 0.15) !important;
    }
    .stSelectbox > div > div, .stMultiSelect > div > div, .stTextInput > div > div > input, .stTextArea > div > div > textarea {
        background: transparent !important;
        border: none !important;
        font-family: 'Inter', sans-serif !important;
        font-size: 0.875rem !important;
        color: var(--text-primary) !important;
        padding: 12px 16px !important;
    }
    
    /* ENHANCED VISUAL ELEMENTS */
    .stProgress > div > div {
        background: var(--primary-gradient) !important;
        border-radius: 10px !important;
        box-shadow: 0 2px 10px rgba(102, 126, 234, 0.3) !important;
    }
    .stProgress > div {
        background: rgba(255, 255, 255, 0.2) !important;
        border-radius: 10px !important;
        backdrop-filter: blur(5px) !important;
    }
    
    /* PREMIUM SIDEBAR */
    .css-1d391kg {
        background: linear-gradient(135deg, rgba(255, 255, 255, 0.95) 0%, rgba(255, 255, 255, 0.8) 100%) !important;
        backdrop-filter: blur(20px) !important;
        border-right: 1px solid var(--glass-border) !important;
    }
    
    /* CUSTOM ANIMATIONS */
    @keyframes fadeInUp {
        from {
            opacity: 0;
            transform: translateY(30px);
        }
        to {
            opacity: 1;
            transform: translateY(0);
        }
    }
    .stContainer {
        animation: fadeInUp 0.6s ease-out;
    }
    
    /* SCROLLBAR STYLING */
    ::-webkit-scrollbar {
        width: 8px;
        height: 8px;
    }
    ::-webkit-scrollbar-track {
        background: rgba(255, 255, 255, 0.1);
        border-radius: 10px;
    }
    ::-webkit-scrollbar-thumb {
        background: var(--primary-gradient);
        border-radius: 10px;
        box-shadow: inset 0 0 5px rgba(0, 0, 0, 0.2);
    }
    ::-webkit-scrollbar-thumb:hover {
        background: linear-gradient(135deg, #764ba2 0%, #667eea 100%);
    }
    </style>
""", unsafe_allow_html=True)


# ---------- PAGE CONFIG ----------
st.set_page_config(
    page_title="Non-B DNA Motif Finder",
    layout="wide",
    page_icon="DNA",
    menu_items={'About': "Non-B DNA Motif Finder | Developed by Dr. Venkata Rajesh Yella"}
)

# ---------- CONSTANTS ----------
# Motifs organized by category and alphabetically within each category as per requirements
MOTIF_CATEGORIES = {
    "G-quadruplex-related": [
        "Bipartite G4", "Bulged G4", "Canonical G4", "Imperfect G4", "Multimeric G4", "Relaxed G4"
    ],
    "G-Triplex": [
        "G-Triplex"
    ],
    "i-motif related": [
        "AC-Motif", "i-Motif"
    ],
    "Helix deviations": [
        "Curved DNA", "eGZ (Extruded-G)", "Z-DNA"
    ],
    "Repeat/junction": [
        "Cruciform", "R-Loop", "Slipped DNA", "Sticky DNA", "Triplex DNA"
    ],
    "Hybrid": [
        "Hybrid"
    ],
    "Non-B DNA Clusters": [
        "Non-B DNA Clusters"
    ]
}

# Flatten categories to create alphabetical order within categories
MOTIF_ORDER = []
for category, motifs in MOTIF_CATEGORIES.items():
    MOTIF_ORDER.extend(sorted(motifs))

# Enhanced color scheme from the reference implementation
# Enhanced premium color schemes for top-tier visualization and accessibility
MOTIF_COLORS = {
    # G4 Family - Sophisticated blues with perfect contrast
    "Canonical G4": "#1E40AF", "Relaxed G4": "#3B82F6", "Bulged G4": "#60A5FA", 
    "Bipartite G4": "#1D4ED8", "Multimeric G4": "#2563EB", "Imperfect G4": "#93C5FD",
    
    # Alternative Structures - Elegant warm tones  
    "Z-DNA": "#DC2626", "eGZ (Extruded-G)": "#EA580C", "Curved DNA": "#F59E0B",
    
    # Repeats - Professional greens
    "Slipped DNA": "#059669", "R-Loop": "#10B981", "Sticky DNA": "#34D399",
    
    # Junctions - Rich purples and magentas
    "Cruciform": "#7C3AED", "Triplex DNA": "#8B5CF6", "G-Triplex": "#A78BFA",
    
    # Special Motifs - Distinctive premium colors
    "i-Motif": "#EC4899", "AC-Motif": "#F59E0B", "Hybrid": "#06B6D4", 
    "Non-B DNA Clusters": "#EF4444"
}

# Premium publication-quality color scheme for scientific figures
PUBLICATION_COLORS = {
    # G4 Family - Deep professional blues
    'Canonical G4': '#1E3A8A',      
    'Relaxed G4': '#3B82F6',        
    'Bulged G4': '#60A5FA',         
    'Bipartite G4': '#1E40AF',      
    'Multimeric G4': '#1D4ED8',     
    'Imperfect G4': '#93C5FD',      
    
    # G-related structures - Distinguished colors
    'G-Triplex': '#7C2D92',         
    'i-Motif': '#EC4899',           
    
    # Z-DNA Family - Rich purples (left-handed helix)
    'Z-DNA': '#7C2D92',             
    'eGZ (Extruded-G)': '#A855F7',  
    
    # Structural variants - Sophisticated warm colors
    'Curved_DNA': '#DC2626',        
    'Slipped_DNA': '#EA580C',       
    'Cruciform': '#F59E0B',         
    
    # Complex structures - Premium earth tones
    'Triplex_DNA': '#059669',       
    'Sticky_DNA': '#0891B2',        
    'R-Loop': '#7C3AED',            
    
    # Specialized motifs - Distinctive premium colors
    'AC-Motif': '#F59E0B',          
    'Hybrid': '#06B6D4',            
    'Non-B DNA Clusters': '#EF4444' 
}
PAGES = {
    "Home": "Overview",
    "Upload & Analyze": "Sequence Upload and Motif Analysis",
    "Results": "Analysis Results and Visualization", 
    "Advanced Disease Results": "Disease-Associated Non-B DNA Motifs Analysis",
    "Download": "Export Data",
    "Documentation": "Scientific Documentation & References"
}
Entrez.email = "raazbiochem@gmail.com"
Entrez.api_key = None

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

# Enhanced example sequences with specific Non-B DNA structures
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

# Famous genes/sequences known for Non-B DNA motifs  
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
    """Create a progress tracking container"""
    if st.session_state.is_analyzing:
        progress_container = st.container()
        with progress_container:
            col1, col2, col3 = st.columns([2, 1, 1])
            
            with col1:
                progress_bar = st.progress(st.session_state.analysis_progress)
                if st.session_state.analysis_start_time:
                    elapsed = time.time() - st.session_state.analysis_start_time
                    minutes = int(elapsed // 60)
                    seconds = int(elapsed % 60)
                    time_str = f"{minutes:02d}:{seconds:02d}"
                    st.write(f"Time: Analysis time: {time_str} | Status: {st.session_state.analysis_status}")
            
            with col2:
                if st.button("Stop Stop Analysis", key="stop_btn"):
                    st.session_state.stop_analysis = True
                    st.session_state.is_analyzing = False
                    st.warning("Analysis stopped by user")
                    
            with col3:
                st.write("Processing Processing...")
        
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
    
    if run_all:
        st.session_state.analysis_status = f"Running comprehensive motif detection on {seq_name}..."
        motifs = all_motifs(seq)
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
        
        color = PUBLICATION_COLORS.get(motif_class, "#666666")
        
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
    # Premium header with modern glassmorphism design
    st.markdown("""
    <div style='text-align: center; padding: 40px; background: linear-gradient(135deg, rgba(255, 255, 255, 0.25) 0%, rgba(255, 255, 255, 0.1) 100%); backdrop-filter: blur(20px); border-radius: 24px; margin-bottom: 40px; border: 1px solid rgba(255, 255, 255, 0.2); box-shadow: 0 16px 64px rgba(31, 38, 135, 0.25); position: relative; overflow: hidden;'>
        <div style='position: absolute; top: 0; left: 0; right: 0; bottom: 0; background: linear-gradient(135deg, rgba(102, 126, 234, 0.1) 0%, rgba(118, 75, 162, 0.1) 100%); pointer-events: none;'></div>
        <h1 style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; font-family: Inter, sans-serif; font-weight: 900; margin-bottom: 20px; font-size: 3.5rem; letter-spacing: -0.05em; position: relative; z-index: 1;'>
            NBDFinder: Advanced Non-B DNA Analysis Platform
        </h1>
        <p style='background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; font-size: 1.5rem; font-weight: 600; margin-top: 16px; position: relative; z-index: 1;'>
            Next-Generation Computational Framework for Genome-Wide Motif Discovery
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Create main content layout
    left, right = st.columns([1.2, 1])
    
    with left:
        # Enhanced image with caption
        st.image("nbdcircle.JPG", use_container_width=True, caption="Non-B DNA structural diversity: From canonical B-form to complex alternative conformations")
        

    
    with right:
        st.markdown("""
        <div style='font-family: Inter, sans-serif; font-size: 1rem; color: #2d3748; line-height: 1.7; padding: 32px; background: linear-gradient(135deg, rgba(255, 255, 255, 0.95) 0%, rgba(255, 255, 255, 0.8) 100%); backdrop-filter: blur(20px); border-radius: 24px; box-shadow: 0 16px 64px rgba(31, 38, 135, 0.15); border: 1px solid rgba(255, 255, 255, 0.2); position: relative; overflow: hidden;'>
        
        <div style='position: absolute; top: 0; left: 0; right: 0; bottom: 0; background: linear-gradient(135deg, rgba(102, 126, 234, 0.03) 0%, rgba(118, 75, 162, 0.03) 100%); pointer-events: none;'></div>
                    
        <div style='margin-bottom: 32px; position: relative; z-index: 1;'>
            <h3 style='margin-top: 0; margin-bottom: 16px; font-size: 1.875rem; font-weight: 800;'>Premium DNA Analysis Platform</h3>
            <p style='margin-bottom: 20px; font-size: 1.125rem; font-weight: 400;'>
                <strong style='background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;'>
                Non-canonical DNA structures
                </strong> play crucial roles in genome organization, gene regulation, and disease pathogenesis. Our platform provides cutting-edge computational tools for comprehensive structural analysis.
            </p>
        </div>
        
        <div style='margin-bottom: 32px; position: relative; z-index: 1;'>
            <h4 style='margin-bottom: 20px; font-size: 1.5rem; font-weight: 700;'>Advanced Motif Detection Suite</h4>
            <div style='background: linear-gradient(135deg, rgba(102, 126, 234, 0.08) 0%, rgba(118, 75, 162, 0.08) 100%); padding: 24px; border-radius: 16px; border-left: 4px solid #667eea; box-shadow: 0 8px 32px rgba(102, 126, 234, 0.1);'>
                <div style='margin-bottom: 16px;'>
                    <span style='font-weight: 700; font-size: 1rem; background: linear-gradient(135deg, #1E40AF 0%, #3B82F6 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;'>G-Quadruplex Family:</span> 
                    <span style='font-size: 0.95rem; color: #4a5568;'>Canonical G4, Relaxed G4, Bulged G4, Bipartite G4, Multimeric G4, Imperfect G4</span>
                </div>
                <div style='margin-bottom: 16px;'>
                    <span style='font-weight: 700; font-size: 1rem; background: linear-gradient(135deg, #7C3AED 0%, #8B5CF6 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;'>Triplex Structures:</span> 
                    <span style='font-size: 0.95rem; color: #4a5568;'>G-Triplex, Triplex DNA, i-Motif</span>
                </div>
                <div style='margin-bottom: 16px;'>
                    <span style='font-weight: 700; font-size: 1rem; background: linear-gradient(135deg, #DC2626 0%, #EF4444 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;'>Helix Deviations:</span> 
                    <span style='font-size: 0.95rem; color: #4a5568;'>Z-DNA, eGZ (Extruded-G), Curved DNA, AC-Motif</span>
                </div>
                <div style='margin-bottom: 16px;'>
                    <span style='font-weight: 700; font-size: 1rem; background: linear-gradient(135deg, #059669 0%, #10B981 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;'>Junction/Repeat:</span> 
                    <span style='font-size: 0.95rem; color: #4a5568;'>Slipped DNA, Cruciform, Sticky DNA, R-Loop</span>
                </div>
                <div style='margin-bottom: 0px;'>
                    <span style='font-weight: 700; font-size: 1rem; background: linear-gradient(135deg, #06B6D4 0%, #0891B2 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text;'>Advanced Analysis:</span> 
                    <span style='font-size: 0.95rem; color: #4a5568;'>Hybrid Motifs, Non-B DNA Clusters</span>
                </div>
            </div>
        </div>
        
        <div style='position: relative; z-index: 1;'>
            <h4 style='margin-bottom: 20px; font-size: 1.25rem; font-weight: 700;'>Premium Features</h4>
            <ul style='padding-left: 0; margin-bottom: 0; list-style: none;'>
                <li style='margin-bottom: 12px; font-size: 1rem; padding: 8px 0; display: flex; align-items: center;'>
                    <span style='width: 8px; height: 8px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border-radius: 50%; margin-right: 12px; flex-shrink: 0;'></span>
                    <span><strong style='color: #667eea;'>Multiple Algorithms:</strong> Advanced detection methods optimized for different motif types</span>
                </li>
                <li style='margin-bottom: 12px; font-size: 1rem; padding: 8px 0; display: flex; align-items: center;'>
                    <span style='width: 8px; height: 8px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border-radius: 50%; margin-right: 12px; flex-shrink: 0;'></span>
                    <span><strong style='color: #667eea;'>Premium Interface:</strong> Intuitive, responsive browser-based platform</span>
                </li>
                <li style='margin-bottom: 12px; font-size: 1rem; padding: 8px 0; display: flex; align-items: center;'>
                    <span style='width: 8px; height: 8px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border-radius: 50%; margin-right: 12px; flex-shrink: 0;'></span>
                    <span><strong style='color: #667eea;'>Interactive Visualizations:</strong> Publication-quality charts and statistical analysis</span>
                </li>
                <li style='margin-bottom: 0px; font-size: 1rem; padding: 8px 0; display: flex; align-items: center;'>
                    <span style='width: 8px; height: 8px; background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); border-radius: 50%; margin-right: 12px; flex-shrink: 0;'></span>
                    <span><strong style='color: #667eea;'>Flexible Input:</strong> Support for individual sequences and multi-FASTA files</span>
                </li>
            </ul>
        </div>
        
        </div>
# Simplified working version for demonstration
st.markdown("""
---
<div style='font-size: 1.05rem; color: #1e293b; margin-top: 30px; text-align: left; font-family:Inter,Arial;'>
<b>Developed by Dr. Venkata Rajesh Yella</b><br>
Premium NBDFinder Web Application with Top 1 Percentile Design
</div>
""", unsafe_allow_html=True)
    


# ---------- UPLOAD & ANALYZE ----------
with tab_pages["Upload & Analyze"]:
    # Premium page header for Upload & Analyze
    st.markdown("""
    <div style='text-align: center; padding: 32px; background: linear-gradient(135deg, rgba(255, 255, 255, 0.25) 0%, rgba(255, 255, 255, 0.1) 100%); backdrop-filter: blur(20px); border-radius: 24px; margin-bottom: 32px; border: 1px solid rgba(255, 255, 255, 0.2); box-shadow: 0 16px 64px rgba(31, 38, 135, 0.25); position: relative; overflow: hidden;'>
        <div style='position: absolute; top: 0; left: 0; right: 0; bottom: 0; background: linear-gradient(135deg, rgba(67, 172, 254, 0.1) 0%, rgba(0, 242, 254, 0.1) 100%); pointer-events: none;'></div>
        <h2 style='background: linear-gradient(135deg, #4facfe 0%, #00f2fe 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; font-family: Inter, sans-serif; font-weight: 900; margin-bottom: 16px; font-size: 2.75rem; letter-spacing: -0.025em; position: relative; z-index: 1;'>
            Advanced Sequence Analysis
        </h2>
        <p style='background: linear-gradient(135deg, #43e97b 0%, #38f9d7 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; background-clip: text; font-size: 1.25rem; font-weight: 600; margin-top: 8px; position: relative; z-index: 1;'>
            Premium Non-B DNA Structure Analysis Pipeline
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    # Enhanced Motif class selection with scientific categorization
    st.markdown("""
    <div style='background: linear-gradient(135deg, #f0f9ff 0%, #e3f2fd 50%, #f8fdff 100%); border-radius: 16px; padding: 25px; margin-bottom: 30px; border: 3px solid #1565c0; box-shadow: 0 6px 20px rgba(21, 101, 192, 0.15);'>
        <h3 style='color: #0d47a1; margin-top: 0; margin-bottom: 18px; font-size: 1.5rem; font-weight: 800;'>Select Motif Classes for Analysis</h3>
        <p style='margin-bottom: 15px; color: #555; font-size: 1.15rem; line-height: 1.6;'>Select specific Non-B DNA motif classes for targeted analysis. Our comprehensive detection suite covers all major structural categories validated by experimental studies.</p>
    </div>
    """, unsafe_allow_html=True)
    # Main multiselect with enhanced help
    selected_motifs = st.multiselect(
        "Select specific motif classes for analysis:", 
        MOTIF_ORDER, 
        default=MOTIF_ORDER,
        help="""
        **Motif Categories & Clinical Significance:**
        
        G-quadruplex-related: Therapeutic targets in cancer, telomere biology
        G-Triplex: Immunoglobulin diversification, antibody engineering
        i-motif related: pH-sensing, metabolic regulation
        Helix deviations: Gene regulation, chromatin structure
        Repeat/junction: Genetic instability, disease mechanisms
        Hybrid: Complex regulatory networks
        Non-B DNA Clusters: Mutational hotspots, evolutionary breakpoints
        """
    )
    st.session_state.selected_motifs = selected_motifs if selected_motifs else MOTIF_ORDER

    # Enhanced input method selection
    st.markdown("""
    <div style='background: linear-gradient(135deg, #f0fdf4 0%, #f8fdff 50%, #f0f9ff 100%); border-radius: 16px; padding: 25px; margin: 30px 0; border: 3px solid #059669; box-shadow: 0 6px 20px rgba(5, 150, 105, 0.15);'>
        <h3 style='color: #059669; margin-top: 0; margin-bottom: 18px; font-size: 1.5rem; font-weight: 800;'>Input Method Selection</h3>
        <p style='margin-bottom: 0; color: #555; font-size: 1.15rem; line-height: 1.6;'>Choose your preferred method for sequence input. All formats support both single sequences and batch processing.</p>
    </div>
    """, unsafe_allow_html=True)
    st.markdown('<p class="input-method-title">Input Method:</p>', unsafe_allow_html=True)
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

    # --- File upload ---
    if input_method == "Upload FASTA / Multi-FASTA File":
        fasta_file = st.file_uploader("Drag and drop FASTA/multi-FASTA file here", type=["fa", "fasta", "txt"])
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

    # --- Paste sequence ---
    elif input_method == "Paste Sequence(s)":
        # Show example format in an expandable section
        with st.expander("View Single FASTA or Multiple FASTA Format Example"):
            st.markdown("**Multi-FASTA Format Example:**")
            st.code(EXAMPLE_MULTI_FASTA, language='text')
            st.markdown("**Tips:**")
            st.markdown("- Each sequence starts with `>` followed by sequence name")
            st.markdown("- DNA sequence follows on the next line(s)")
            st.markdown("- Multiple sequences can be pasted at once")
        
        seq_input = st.text_area("Paste FASTA or raw sequence(s)", height=150)
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
        
        # Display example options with descriptions
        example_choices = list(EXAMPLE_SEQUENCES.keys())
        
        selected_example = st.selectbox(
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

    # --- Enhanced NCBI Query with famous examples ---
    elif input_method == "NCBI Fetch":
        st.markdown("**- Fetch sequences from NCBI database:**")
        
        # Show famous examples
        with st.expander("▸ Famous genes/sequences with Non-B DNA motifs"):
            st.markdown("**Select from famous examples:**")
            
            # Create a searchable selectbox for famous examples
            example_options = [""] + [f"{gene} ({accession})" for gene, accession in FAMOUS_NCBI_EXAMPLES.items()]
            selected_example = st.selectbox(
                "Choose a famous gene/sequence:",
                options=example_options,
                help="Select a well-known gene or sequence with Non-B DNA motifs"
            )
            
            if selected_example:
                # Extract accession from the selected example
                accession = selected_example.split("(")[1].rstrip(")")
                st.session_state.ncbi_query = accession
                st.success(f"Loaded Selected: {selected_example}")
            
            # Alternative: Show examples as formatted text for easy reference
            st.markdown("**Available examples:**")
            examples_text = ""
            for gene, accession in FAMOUS_NCBI_EXAMPLES.items():
                examples_text += f"- **{gene}**: `{accession}`\n"
            st.markdown(examples_text)
        
        # NCBI query input
        ncbi_query = st.text_input(
            "Enter NCBI Query:", 
            value=st.session_state.get('ncbi_query', ''),
            placeholder="Gene name, accession number, or search term",
            help="Examples: BRCA1, NM_007294.3, 'human telomerase'"
        )
        
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
    if st.session_state.seqs and st.button("Analyze Sequences"):
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

# ---------- RESULTS ----------
with tab_pages["Results"]:
    st.markdown('<h2>Analysis Results and Visualization</h2>', unsafe_allow_html=True)
    if not st.session_state.results:
        st.info("No analysis results. Please run motif analysis first.")
    else:
        # Enhanced summary table with key metrics only
        st.subheader("- Analysis Summary")
        summary_cols = ["Sequence Name", "Length (bp)", "GC %", "Motif Count", "Motif Coverage (%)", "Top Motifs"]
        st.dataframe(st.session_state.summary_df[summary_cols], use_container_width=True)
        
        seq_idx = 0
        if len(st.session_state.seqs) > 1:
            seq_idx = st.selectbox("Choose Sequence for Detailed Analysis:", range(len(st.session_state.seqs)), format_func=lambda i: st.session_state.names[i])
        
        # Extract motifs from results tuple (name, motifs)
        _, motifs = st.session_state.results[seq_idx]
        if not motifs:
            st.warning("No motifs detected for this sequence.")
        else:
            df = pd.DataFrame(motifs)
            
            # Enhanced motif table with essential columns only
            st.markdown(f"<h3>- Detailed Motifs for <b>{st.session_state.names[seq_idx]}</b></h3>", unsafe_allow_html=True)
            
            # Check available columns and adapt
            available_columns = df.columns.tolist()
            
            # Key columns for display (adapting to actual column names)
            essential_columns = []
            column_mapping = {
                'Class': ['Class'],
                'Subtype': ['Subtype'],
                'Start': ['Start'],
                'End': ['End'],
                'Length': ['Length'],
                'Score': ['Score']
            }
            
            # Add columns that exist
            for display_name, possible_names in column_mapping.items():
                for col_name in possible_names:
                    if col_name in available_columns:
                        essential_columns.append(col_name)
                        break
            
            # Add predicted sequence preview (first 50 chars)
            if 'Sequence' in df.columns:
                df['Sequence Preview'] = df['Sequence'].apply(lambda x: str(x).replace('\n', '')[:50] + '...' if len(str(x).replace('\n', '')) > 50 else str(x).replace('\n', ''))
                essential_columns.append('Sequence Preview')
            
            # Add serial number starting from 1
            df['S.No'] = range(1, len(df) + 1)
            essential_columns = ['S.No'] + essential_columns
            
            # Only use columns that actually exist
            final_columns = [col for col in essential_columns if col in df.columns]
            display_df = df[final_columns].copy()
            
            # Format numeric columns
            if 'Score' in display_df.columns:
                display_df['Score'] = display_df['Score'].apply(lambda x: f"{x:.2f}" if isinstance(x, (int, float)) else str(x))
            
            # Add normalized confidence scoring and significance
            if 'Score' in display_df.columns:
                def get_normalized_confidence_score(score, motif_class):
                    """
                    Normalize scores to a unified 1.0-3.0 scale where:
                    1.0 = Minimum viable threshold (Moderate confidence)
                    2.0 = High confidence threshold 
                    3.0 = Optimal/Validated confidence
                    """
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
                    """Convert normalized score to confidence label with enhanced styling"""
                    if normalized_score >= 2.5:
                        return "Green Optimal (>=2.5)"
                    elif normalized_score >= 2.0:
                        return "High (>=2.0)"
                    elif normalized_score >= 1.0:
                        return "Yellow Moderate (>=1.0)"
                    else:
                        return "Red Low (<1.0)"
                
                # Add normalized confidence columns
                display_df['Normalized Score'] = display_df.apply(
                    lambda row: f"{get_normalized_confidence_score(row.get('Score', 0), row.get('Class', '')):.2f}",
                    axis=1
                )
                
                display_df['Confidence Level'] = display_df.apply(
                    lambda row: get_confidence_label(
                        get_normalized_confidence_score(row.get('Score', 0), row.get('Class', ''))
                    ),
                    axis=1
                )
                
            # Enhanced dataframe styling and display
            st.dataframe(
                display_df, 
                use_container_width=True, 
                height=400,
                hide_index=True  # Hide the default pandas index to avoid confusion
            )
            
            # Enhanced visualizations - Comprehensive Overview
            st.markdown('<h4>- Comprehensive Motif Analysis - All 19 Types</h4>', unsafe_allow_html=True)
            
            # Use new comprehensive visualization
            sequence_length = len(st.session_state.seqs[seq_idx])
            comprehensive_fig = create_comprehensive_motif_overview(
                motifs, sequence_length, st.session_state.names[seq_idx]
            )
            st.plotly_chart(comprehensive_fig, use_container_width=True)
            
            # Additional detailed visualizations in columns
            col1, col2 = st.columns(2)
            
            with col1:
                st.markdown('<h4>- Detection Summary</h4>', unsafe_allow_html=True)
                
                # Summary statistics 
                if len(motifs) > 0:
                    total_motifs = len(motifs)
                    unique_types = len(set(m.get('Class', '') for m in motifs))
                    sequence_length = len(st.session_state.seqs[seq_idx])
                    total_coverage = sum(m.get('Length', 0) for m in motifs)
                    coverage_percent = (total_coverage / sequence_length * 100) if sequence_length > 0 else 0
                    
                    st.metric("Total Motifs Detected", f"{total_motifs}")
                    st.metric("Unique Motif Types", f"{unique_types} / 19")
                    st.metric("Sequence Coverage", f"{coverage_percent:.1f}%")
                    
                    # Detection rate
                    detection_rate = (unique_types / 19) * 100
                    if detection_rate >= 50:
                        st.success(f"High diversity: {detection_rate:.1f}% of motif types detected")
                    elif detection_rate >= 25:
                        st.info(f"Moderate diversity: {detection_rate:.1f}% of motif types detected")
                    else:
                        st.warning(f"Limited diversity: {detection_rate:.1f}% of motif types detected")
                else:
                    st.info("No motifs detected in this sequence.")
            
            with col2:
                st.markdown('<h4>- Motif Categories Overview</h4>', unsafe_allow_html=True)
                
                # Category-based analysis
                if len(motifs) > 0:
                    # Group motifs by categories
                    category_counts = {
                        "G-quadruplex-related": 0,
                        "G-Triplex": 0, 
                        "i-motif related": 0,
                        "Helix deviations": 0,
                        "Repeat/junction": 0,
                        "Hybrid": 0,
                        "Non-B DNA Clusters": 0
                    }
                    
                    for motif in motifs:
                        motif_class = motif.get('Class', '')
                        # Handle special case for eGZ
                        if motif_class == "Z-DNA" and motif.get("Subclass", "") == "eGZ (Extruded-G)":
                            motif_class = "eGZ (Extruded-G)"
                        
                        # Map to categories
                        if motif_class in ["Bipartite G4", "Bulged G4", "Canonical G4", "Imperfect G4", "Multimeric G4", "Relaxed G4"]:
                            category_counts["G-quadruplex-related"] += 1
                        elif motif_class == "G-Triplex":
                            category_counts["G-Triplex"] += 1
                        elif motif_class in ["AC-Motif", "i-Motif"]:
                            category_counts["i-motif related"] += 1
                        elif motif_class in ["Curved DNA", "eGZ (Extruded-G)", "Z-DNA"]:
                            category_counts["Helix deviations"] += 1
                        elif motif_class in ["Cruciform", "R-Loop", "Slipped DNA", "Sticky DNA", "Triplex DNA"]:
                            category_counts["Repeat/junction"] += 1
                        elif motif_class == "Hybrid":
                            category_counts["Hybrid"] += 1
                        elif motif_class == "Non-B DNA Clusters":
                            category_counts["Non-B DNA Clusters"] += 1
                    
                    # Display category breakdown
                    for category, count in category_counts.items():
                        if count > 0:
                            st.write(f"**{category}**: {count} motifs")
                        else:
                            st.write(f"**{category}**: Not found")
                else:
                    st.info("No motifs available for category analysis.")
            
            # Enhanced interactive motif map
            st.markdown('<h4>Interactive Motif Map</h4>', unsafe_allow_html=True)
            
            # Create the enhanced visualization
            interactive_fig = create_enhanced_motif_visualization(motifs, st.session_state.names[seq_idx], len(st.session_state.seqs[seq_idx]))
            if interactive_fig:
                st.plotly_chart(interactive_fig, use_container_width=True)
            
            # Motif density heatmap
            if len(motifs) > 0:
                st.markdown('<h4>▦ Motif Density Heatmap</h4>', unsafe_allow_html=True)
                
                # Create density array
                seq_len = len(st.session_state.seqs[seq_idx])
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
            
            # ========== ENHANCED ADVANCED VISUALIZATIONS ==========
            if ADVANCED_VIZ_AVAILABLE and len(motifs) > 0:
                st.markdown("---")
                st.markdown('<h3>🚀 Advanced Analysis Dashboard</h3>', unsafe_allow_html=True)
                st.markdown("*Publication-quality visualizations and comprehensive analysis*")
                
                # Create enhanced dashboard
                motifs_df = pd.DataFrame(motifs)
                seq_length = len(st.session_state.seqs[seq_idx])
                
                try:
                    enhanced_plots = create_enhanced_dashboard(motifs_df, seq_length)
                    
                    # Create tabs for different analysis types
                    viz_tabs = st.tabs([
                        "Comprehensive Analysis", 
                        "Clinical Significance", 
                        "🤖 ML Predictions",
                        "🗺️ Genomic Map",
                        "🌟 Cluster Analysis",
                        "🔗 3D Structure"
                    ])
                    
                    with viz_tabs[0]:
                        st.markdown("### Comprehensive Motif Distribution Analysis")
                        if 'distribution' in enhanced_plots:
                            st.plotly_chart(enhanced_plots['distribution'], use_container_width=True)
                            
                            # Add export button
                            if st.button("🔽 Export Distribution Plot", key="export_dist"):
                                try:
                                    img_bytes = export_publication_figure(
                                        enhanced_plots['distribution'], 
                                        "motif_distribution", 
                                        format='png'
                                    )
                                    st.download_button(
                                        label="Download PNG",
                                        data=img_bytes,
                                        file_name="motif_distribution.png",
                                        mime="image/png"
                                    )
                                except Exception as e:
                                    st.error(f"Export failed: {e}")
                    
                    with viz_tabs[1]:
                        st.markdown("### Clinical Significance Analysis")
                        if 'clinical' in enhanced_plots:
                            st.plotly_chart(enhanced_plots['clinical'], use_container_width=True)
                            
                            # Display clinical summary
                            disease_motifs = motifs_df[motifs_df['Class'] == 'Disease-Associated Motif']
                            if not disease_motifs.empty:
                                st.markdown("#### Clinical Summary")
                                col1, col2, col3 = st.columns(3)
                                
                                with col1:
                                    pathogenic_count = disease_motifs[
                                        disease_motifs.get('Clinical_Significance', '').str.contains('Pathogenic', na=False)
                                    ].shape[0]
                                    st.metric("Pathogenic Variants", pathogenic_count)
                                
                                with col2:
                                    avg_risk = disease_motifs.get('Risk_Score', [0]).astype(float).mean()
                                    st.metric("Average Risk Score", f"{avg_risk:.1f}")
                                
                                with col3:
                                    diseases = disease_motifs.get('Disease_Name', []).nunique()
                                    st.metric("Disease Categories", diseases)
                    
                    with viz_tabs[2]:
                        st.markdown("### Machine Learning Enhanced Predictions")
                        if 'ml_predictions' in enhanced_plots:
                            st.plotly_chart(enhanced_plots['ml_predictions'], use_container_width=True)
                            
                            # ML performance metrics
                            ml_motifs = motifs_df.dropna(subset=['ML_Probability'])
                            if not ml_motifs.empty:
                                st.markdown("#### ML Performance Metrics")
                                col1, col2, col3, col4 = st.columns(4)
                                
                                with col1:
                                    avg_prob = ml_motifs['ML_Probability'].astype(float).mean()
                                    st.metric("Avg ML Probability", f"{avg_prob:.3f}")
                                
                                with col2:
                                    avg_conf = ml_motifs['ML_Confidence'].astype(float).mean()
                                    st.metric("Avg Confidence", f"{avg_conf:.3f}")
                                
                                with col3:
                                    high_conf = (ml_motifs['ML_Confidence'].astype(float) > 0.8).sum()
                                    st.metric("High Confidence", high_conf)
                                
                                with col4:
                                    enhanced_count = ml_motifs['Enhanced_Score'].notna().sum()
                                    st.metric("ML Enhanced", enhanced_count)
                    
                    with viz_tabs[3]:
                        st.markdown("### Interactive Genomic Map")
                        if 'genomic_map' in enhanced_plots:
                            st.plotly_chart(enhanced_plots['genomic_map'], use_container_width=True)
                            st.markdown("*Click and drag to zoom, hover for details*")
                    
                    with viz_tabs[4]:
                        st.markdown("### Advanced Cluster Analysis")
                        if 'cluster_analysis' in enhanced_plots:
                            st.plotly_chart(enhanced_plots['cluster_analysis'], use_container_width=True)
                            
                            # Cluster statistics
                            cluster_motifs = motifs_df[motifs_df['Class'].str.contains('Cluster', na=False)]
                            if not cluster_motifs.empty:
                                st.markdown("#### Cluster Statistics")
                                col1, col2 = st.columns(2)
                                
                                with col1:
                                    st.metric("Cluster Regions", len(cluster_motifs))
                                
                                with col2:
                                    avg_diversity = cluster_motifs.get('Diversity_Index', [0]).astype(float).mean()
                                    st.metric("Avg Diversity Index", f"{avg_diversity:.3f}")
                    
                    with viz_tabs[5]:
                        st.markdown("### 3D Structure Visualization")
                        if '3d_structure' in enhanced_plots:
                            st.plotly_chart(enhanced_plots['3d_structure'], use_container_width=True)
                            st.markdown("*Drag to rotate, scroll to zoom*")
                            st.info("This is a simplified 3D representation. For detailed structural analysis, consider molecular modeling software.")
                
                except Exception as e:
                    st.error(f"Advanced visualization error: {e}")
                    st.info("Using standard visualizations instead.")
            
            # ========== ENHANCED SUMMARY STATISTICS ==========
            st.markdown("---")
            st.markdown('<h3>Enhanced Summary Statistics</h3>', unsafe_allow_html=True)
            
            # Create comprehensive statistics
            stats_cols = st.columns(4)
            
            with stats_cols[0]:
                st.metric(
                    "Total Motifs",
                    len(motifs),
                    delta=None
                )
            
            with stats_cols[1]:
                unique_classes = df['Class'].nunique() if 'Class' in df.columns else 0
                st.metric(
                    "Motif Classes",
                    unique_classes
                )
            
            with stats_cols[2]:
                if 'Score' in df.columns:
                    avg_score = pd.to_numeric(df['Score'], errors='coerce').mean()
                    st.metric(
                        "Average Score",
                        f"{avg_score:.1f}" if not pd.isna(avg_score) else "N/A"
                    )
                else:
                    st.metric("Average Score", "N/A")
            
            with stats_cols[3]:
                seq_len = len(st.session_state.seqs[seq_idx])
                total_coverage = sum(m.get('Length', 0) for m in motifs)
                coverage_pct = (total_coverage / seq_len) * 100 if seq_len > 0 else 0
                st.metric(
                    "Coverage",
                    f"{coverage_pct:.1f}%"
                )
            
            # Enhanced feature breakdown
            if any('ML_Probability' in str(m) for m in motifs):
                st.markdown("#### 🤖 Machine Learning Features")
                ml_stats_cols = st.columns(3)
                
                ml_motifs = [m for m in motifs if 'ML_Probability' in m]
                
                with ml_stats_cols[0]:
                    avg_ml_prob = np.mean([float(m.get('ML_Probability', 0)) for m in ml_motifs])
                    st.metric("Avg ML Probability", f"{avg_ml_prob:.3f}")
                
                with ml_stats_cols[1]:
                    high_conf_count = sum(1 for m in ml_motifs if float(m.get('ML_Confidence', 0)) > 0.8)
                    st.metric("High Confidence Motifs", high_conf_count)
                
                with ml_stats_cols[2]:
                    enhanced_count = sum(1 for m in ml_motifs if 'Enhanced_Score' in m)
                    st.metric("ML Enhanced", enhanced_count)
            
            # Disease analysis summary
            disease_motifs = [m for m in motifs if m.get('Class') == 'Disease-Associated Motif']
            if disease_motifs:
                st.markdown("#### Clinical Analysis")
                disease_stats_cols = st.columns(4)
                
                with disease_stats_cols[0]:
                    st.metric("Disease Motifs", len(disease_motifs))
                
                with disease_stats_cols[1]:
                    pathogenic_count = sum(1 for m in disease_motifs if 'Pathogenic' in m.get('Clinical_Significance', ''))
                    st.metric("Pathogenic", pathogenic_count)
                
                with disease_stats_cols[2]:
                    avg_risk = np.mean([float(m.get('Risk_Score', 0)) for m in disease_motifs])
                    st.metric("Avg Risk Score", f"{avg_risk:.1f}")
                
                with disease_stats_cols[3]:
                    unique_diseases = len(set(m.get('Disease_Name', 'Unknown') for m in disease_motifs))
                    st.metric("Disease Types", unique_diseases)

# ---------- ADVANCED DISEASE RESULTS ----------
with tab_pages["Advanced Disease Results"]:
    st.markdown("""
    <div style='text-align: center; background: linear-gradient(135deg, #ffebee 0%, #f3e5f5 50%, #e8f5e8 100%); padding: 20px; border-radius: 15px; margin-bottom: 25px; border: 2px solid #e91e63;'>
        <h2 style='color: #ad1457; margin-bottom: 10px;'>Advanced Disease-Associated Non-B DNA Motifs Analysis</h2>
        <p style='color: #6a1b9a; font-size: 1.1rem; margin: 0;'>
            Comprehensive clinical analysis and disease association assessment
        </p>
    </div>
    """, unsafe_allow_html=True)
    
    if not st.session_state.results:
        st.info("Please run motif analysis first to view disease-associated results.")
    else:
        # Aggregate disease motifs from all sequences
        all_disease_motifs = []
        disease_sequences = []
        
        for seq_name, motifs in st.session_state.results:
            # Filter for disease-associated motifs
            disease_motifs = [m for m in motifs if 
                            m.get('Class') == 'Disease-Associated Motif' or
                            any(keyword in str(m.get('Subtype', '')).lower() for keyword in 
                                ['repeat', 'expansion', 'disease', 'pathogenic', 'clinical'])]
            
            if disease_motifs:
                all_disease_motifs.extend(disease_motifs)
                disease_sequences.append((seq_name, disease_motifs))
        
        if not all_disease_motifs:
            st.warning("No disease-associated motifs detected in the analyzed sequences.")
            st.info("""
            **Disease-Associated Motifs Include:**
            - GAA/TTC repeats (Friedreich Ataxia)
            - CGG repeats (Fragile X Syndrome) 
            - CAG repeats (Huntington Disease)
            - G4C2 repeats (ALS/FTD)
            - Other pathogenic repeat expansions
            """)
        else:
            # Disease Analysis Dashboard
            st.markdown("### Disease Motif Analysis Dashboard")
            
            disease_cols = st.columns(4)
            with disease_cols[0]:
                st.metric("Total Disease Motifs", len(all_disease_motifs))
            with disease_cols[1]:
                affected_seqs = len(disease_sequences)
                st.metric("Affected Sequences", affected_seqs)
            with disease_cols[2]:
                risk_scores = [float(m.get('Risk_Score', 0)) for m in all_disease_motifs if m.get('Risk_Score')]
                avg_risk = sum(risk_scores) / len(risk_scores) if risk_scores else 0
                st.metric("Warning: Average Risk Score", f"{avg_risk:.2f}")
            with disease_cols[3]:
                pathogenic_count = len([m for m in all_disease_motifs 
                                     if str(m.get('Clinical_Significance', '')).startswith('Pathogenic')])
                st.metric("🚨 Pathogenic Motifs", pathogenic_count)
            
            # Clinical Significance Distribution
            st.markdown("### 🏥 Clinical Significance Analysis")
            clinical_significance = {}
            for motif in all_disease_motifs:
                sig = motif.get('Clinical_Significance', 'Unknown')
                clinical_significance[sig] = clinical_significance.get(sig, 0) + 1
            
            if clinical_significance:
                import plotly.graph_objects as go
                fig_clinical = go.Figure(data=[
                    go.Bar(x=list(clinical_significance.keys()), 
                          y=list(clinical_significance.values()),
                          marker_color=['#d32f2f' if 'Pathogenic' in k else 
                                       '#ff9800' if 'VUS' in k else 
                                       '#4caf50' if 'Benign' in k else '#9e9e9e' 
                                       for k in clinical_significance.keys()])
                ])
                fig_clinical.update_layout(
                    title="Clinical Significance Distribution",
                    xaxis_title="Clinical Classification",
                    yaxis_title="Number of Motifs",
                    height=400
                )
                st.plotly_chart(fig_clinical, use_container_width=True)
            
            # Disease Type Analysis
            st.markdown("### 🦠 Disease Association Analysis")
            disease_types = {}
            for motif in all_disease_motifs:
                disease = motif.get('Disease_Name', 'Unknown')
                if disease != 'Unknown':
                    disease_types[disease] = disease_types.get(disease, 0) + 1
            
            if disease_types:
                fig_diseases = go.Figure(data=[
                    go.Pie(labels=list(disease_types.keys()), 
                          values=list(disease_types.values()),
                          hole=0.4)
                ])
                fig_diseases.update_layout(
                    title="Disease Association Distribution",
                    height=500
                )
                st.plotly_chart(fig_diseases, use_container_width=True)
            
            # Detailed Disease Motifs Table
            st.markdown("### Detailed Disease Motifs Analysis")
            
            # Create comprehensive disease motifs dataframe
            disease_df = pd.DataFrame(all_disease_motifs)
            
            # Select relevant columns for disease analysis
            disease_columns = ['Sequence Name', 'Subtype', 'Start', 'End', 'Length', 
                             'Repeat_Count', 'Clinical_Significance', 'Disease_Name', 
                             'Risk_Score', 'Pathogenic_Threshold']
            
            # Filter to available columns
            available_disease_columns = [col for col in disease_columns if col in disease_df.columns]
            
            if available_disease_columns:
                st.dataframe(disease_df[available_disease_columns], use_container_width=True)
            else:
                # Fallback display
                basic_columns = ['Subtype', 'Start', 'End', 'Length', 'Score']
                display_columns = [col for col in basic_columns if col in disease_df.columns]
                st.dataframe(disease_df[display_columns], use_container_width=True)
            
            # Therapeutic Implications
            st.markdown("### 💊 Therapeutic Implications")
            
            therapeutic_targets = set()
            for motif in all_disease_motifs:
                target = motif.get('Therapeutic_Target', '')
                if target and target != 'None identified':
                    therapeutic_targets.add(target)
            
            if therapeutic_targets:
                st.success("**Potential Therapeutic Targets Identified:**")
                for target in therapeutic_targets:
                    st.write(f"- {target}")
            else:
                st.info("No specific therapeutic targets identified in current analysis.")
            
            # Genetic Counseling Recommendations
            counseling_notes = set()
            for motif in all_disease_motifs:
                note = motif.get('Genetic_Counseling', '')
                if note:
                    counseling_notes.add(note)
            
            if counseling_notes:
                st.markdown("### 🧑‍⚕️ Genetic Counseling Recommendations")
                for note in counseling_notes:
                    st.info(f"**Recommendation:** {note}")

# ---------- DOWNLOAD ----------
with tab_pages["Download"]:
    st.header("Export Data")
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
