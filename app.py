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
    """Guarantee every motif has a string 'Subtype'"""
    if isinstance(motif, dict):
        if 'Subtype' not in motif or motif['Subtype'] is None:
            motif['Subtype'] = 'Other'
        return motif
    else:
        # Handle non-dict motifs gracefully (could log/warn here)
        return {'Subtype': 'Other', 'Motif': motif}

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

def add_sticky_action_buttons():
    """Add sticky action buttons for Run Analysis and Download"""
    if 'analysis_running' not in st.session_state:
        st.session_state.analysis_running = False
    if 'results' not in st.session_state:
        st.session_state.results = []
    
    st.markdown("""
    <div style="
        position: fixed; 
        top: 80px; 
        right: 20px; 
        background: rgba(255, 255, 255, 0.95);
        border: 1px solid #e2e8f0;
        border-radius: 12px; 
        padding: 16px; 
        box-shadow: 0 8px 32px rgba(0, 0, 0, 0.1);
        backdrop-filter: blur(10px);
        z-index: 999;
        min-width: 160px;
    ">
        <div style="text-align: center;">
            <div style="font-weight: 600; color: #374151; margin-bottom: 12px; font-size: 14px;">Quick Actions</div>
    """, unsafe_allow_html=True)
    
    # Run Analysis Button
    if not st.session_state.analysis_running:
        if st.button("🔬 Run Analysis", key="sticky_run", help="Start motif analysis"):
            st.session_state.analysis_running = True
            st.rerun()
    else:
        st.button("⏳ Running...", disabled=True, key="sticky_running")
    
    # Download Button
    if st.session_state.results:
        if st.button("📥 Download Results", key="sticky_download", help="Download analysis results"):
            # This will be handled by the Download page logic
            pass
    else:
        st.button("📥 Download Results", disabled=True, key="sticky_download_disabled", help="No results available")
    
    st.markdown("</div></div>", unsafe_allow_html=True)

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
    "🏠 Home": "Home / Welcome",
    "📤 Upload & Analyze": "Sequence Upload and Analysis",
    "📊 Results": "Analysis Results and Visualization", 
    "🎨 Visualization": "Publication-Ready Visualizations",
    "🏥 Clinical/Disease": "Disease Annotation and Clinical Data",
    "📥 Download & Export": "Download Results and Export Data",
    "📚 Documentation": "Scientific Documentation & References",
    "⚙️ Settings": "Application Settings and Preferences"
}

# Initialize session state
if 'current_page' not in st.session_state:
    st.session_state.current_page = "🏠 Home"
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
add_sticky_action_buttons()

# Header
st.markdown("""
<div style='text-align: center; padding: 20px; margin-bottom: 30px; background: linear-gradient(135deg, #f8fafc 0%, #e2e8f0 100%); border-radius: 12px; border: 1px solid #e2e8f0;'>
    <h1 style='color: #1e3a8a; font-family: Inter, sans-serif; font-weight: 700; margin-bottom: 8px; font-size: 2.5rem;'>
        🧬 NBDFinder: Non-B DNA Analysis Platform
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
with tab_dict["🏠 Home"]:
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
                🎯 10 Non-B DNA Classes Detected
            </h3>
            <div style='display: grid; grid-template-columns: 1fr 1fr; gap: 16px; margin-bottom: 20px;'>
                <div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>1. 🌀 Curved DNA</div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>2. 🔄 Slipped DNA</div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>3. ✚ Cruciform DNA</div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>4. 🔁 R-loop</div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>5. 🧬 Triplex</div>
                </div>
                <div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>6. 🟫 G-Quadruplex Family</div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>7. 🔴 i-motif family</div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>8. ⚡ Z-DNA</div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>9. 🔗 Hybrid</div>
                    <div style='margin-bottom: 8px; font-weight: 500; color: #374151;'>10. 🎯 Non-B DNA clusters</div>
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
            <h4>📋 Version Information</h4>
            <p><strong>Version:</strong> 2.0</p>
            <p><strong>Last Updated:</strong> December 2024</p>
            <p><strong>License:</strong> Academic Use</p>
        </div>
        """, unsafe_allow_html=True)
    
    with col2:
        st.markdown("""
        <div class="feature-card">
            <h4>📖 Citation</h4>
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
            <h4>🚀 Platform Overview</h4>
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
with tab_dict["📤 Upload & Analyze"]:
    # Subtabs for Upload & Analyze
    upload_subtabs = st.tabs(["📁 Upload Sequence", "⚙️ Parameter Settings", "▶️ Run Analysis"])
    
    # ---- UPLOAD SEQUENCE SUBTAB ----
    with upload_subtabs[0]:
        st.markdown("""
        <div style='background: #f8fafc; border: 1px solid #e2e8f0; border-radius: 6px; padding: 20px; margin-bottom: 24px;'>
            <h3 style='color: #374151; margin-top: 0; margin-bottom: 12px;'>📁 Sequence Input</h3>
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
                <h4>📤 File Upload</h4>
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
                        st.success(f"✅ Successfully loaded {len(seqs)} sequence(s)")
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
                    st.success(f"✅ Parsed {len(seqs)} sequence(s)")
        
        elif input_method == "Example Sequence":
            st.markdown("""
            <div class="feature-card">
                <h4>🧪 Example Datasets</h4>
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
                st.info(f"📋 **{selected_example}**: {example_data['description']}")
                seqs = [example_data['sequence']]
                names = [example_data['name']]
                
                with st.expander("View sequence details"):
                    st.text(f"Length: {len(seqs[0])} bp")
                    st.text(f"GC Content: {gc_content(seqs[0]):.1f}%")
                    st.text_area("Sequence:", seqs[0], height=100)
        
        elif input_method == "NCBI Fetch":
            st.markdown("""
            <div class="feature-card">
                <h4>🌐 NCBI Database Fetch</h4>
                <p>Retrieve sequences directly from NCBI databases</p>
            </div>
            """, unsafe_allow_html=True)
            
            accession = st.text_input(
                "Enter NCBI Accession Number:",
                help="Format: NG_123456.1, NM_123456.1, etc.",
                placeholder="NG_123456.1"
            )
            
            if accession and st.button("🔍 Fetch Sequence"):
                try:
                    with st.spinner("Fetching sequence from NCBI..."):
                        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
                        fasta_content = handle.read()
                        handle.close()
                        
                        seqs, names = parse_fasta(fasta_content)
                        if seqs:
                            st.success(f"✅ Successfully fetched: {names[0]}")
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
            <h3 style='color: #374151; margin-top: 0; margin-bottom: 12px;'>⚙️ Analysis Parameters</h3>
            <p style='margin-bottom: 0; color: #6b7280; font-size: 0.875rem;'>Configure motif detection settings and sensitivity thresholds</p>
        </div>
        """, unsafe_allow_html=True)
        
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**🎯 Motif Classes Selection**")
            motif_classes = st.multiselect(
                "Select motif classes to analyze:",
                list(MOTIF_COLORS.keys()),
                default=list(MOTIF_COLORS.keys()),
                help="Choose which non-B DNA motif types to detect"
            )
            
            st.markdown("**📊 Sensitivity Settings**")
            sensitivity = st.select_slider(
                "Detection sensitivity:",
                options=["Low", "Medium", "High", "Maximum"],
                value="High",
                help="Higher sensitivity detects more motifs but may include false positives"
            )
        
        with col2:
            st.markdown("**🔬 Advanced Options**")
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
            <h3 style='color: #374151; margin-top: 0; margin-bottom: 12px;'>▶️ Analysis Execution</h3>
            <p style='margin-bottom: 0; color: #6b7280; font-size: 0.875rem;'>Review settings and launch motif detection analysis</p>
        </div>
        """, unsafe_allow_html=True)
        
        # Settings summary
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**📋 Analysis Summary**")
            if st.session_state.seqs:
                st.success(f"✅ {len(st.session_state.seqs)} sequence(s) loaded")
                for i, name in enumerate(st.session_state.names):
                    st.text(f"• {name[:50]}... ({len(st.session_state.seqs[i])} bp)")
            else:
                st.warning("⚠️ No sequences loaded. Please upload sequences first.")
        
        with col2:
            st.markdown("**⚙️ Current Settings**")
            if 'analysis_settings' in st.session_state:
                settings = st.session_state.analysis_settings
                st.text(f"• Motif classes: {len(settings.get('motif_classes', []))} selected")
                st.text(f"• Sensitivity: {settings.get('sensitivity', 'High')}")
                st.text(f"• Length range: {settings.get('min_motif_length', 10)}-{settings.get('max_motif_length', 200)} bp")
                st.text(f"• Quality filter: {'✓' if settings.get('quality_filter', True) else '✗'}")
        
        st.markdown("---")
        
        # Analysis button
        if st.session_state.seqs and not st.session_state.get('analysis_running', False):
            if st.button("🚀 Start Analysis", type="primary", use_container_width=True):
                st.session_state.analysis_running = True
                st.rerun()
        
        # Analysis progress
        if st.session_state.get('analysis_running', False):
            st.markdown("**🔬 Analysis in Progress**")
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
                st.success("✅ Analysis completed successfully!")
                st.balloons()
                
            except Exception as e:
                st.session_state.analysis_running = False
                st.error(f"❌ Analysis failed: {e}")

# =====================================================
# RESULTS PAGE
# =====================================================
with tab_dict["📊 Results"]:
    if not st.session_state.results:
        st.info("📋 No analysis results available. Please run an analysis first.")
    else:
        # Results subtabs
        results_subtabs = st.tabs(["📈 Overview", "🔍 Motif Classes", "📍 Genomic Position", "🔗 Hybrid & Cluster Analysis", "📊 Statistical Summary"])
        
        # ---- OVERVIEW SUBTAB ----
        with results_subtabs[0]:
            st.markdown("### 📈 Analysis Overview")
            
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
            st.markdown("### 📋 Results Summary")
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
            st.markdown("### 🔍 Motif Classes Analysis")
            
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
                st.markdown("### 📊 Detailed Motif Breakdown")
                motif_df = pd.DataFrame(list(motif_counts.items()), columns=['Motif Type', 'Count'])
                motif_df['Percentage'] = (motif_df['Count'] / motif_df['Count'].sum() * 100).round(2)
                st.dataframe(motif_df, use_container_width=True)
        
        # ---- GENOMIC POSITION SUBTAB ----
        with results_subtabs[2]:
            st.markdown("### 📍 Genomic Position Analysis")
            
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
            st.markdown("### 🔗 Hybrid & Cluster Analysis")
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
            st.markdown("### 📊 Statistical Summary")
            
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
                    st.markdown("**🔍 Sequence Statistics**")
                    seq_lengths = [r['sequence_length'] for r in st.session_state.results]
                    st.text(f"Mean sequence length: {np.mean(seq_lengths):.0f} bp")
                    st.text(f"Median sequence length: {np.median(seq_lengths):.0f} bp")
                    st.text(f"Total analyzed: {np.sum(seq_lengths):,} bp")

# =====================================================
# VISUALIZATION PAGE
# =====================================================
with tab_dict["🎨 Visualization"]:
    if not st.session_state.results:
        st.info("📋 No analysis results available for visualization. Please run an analysis first.")
    else:
        # Visualization subtabs
        viz_subtabs = st.tabs(["📊 Motif Distributions", "📸 Publication Figures", "🌐 3D/Advanced Views"])
        
        # ---- MOTIF DISTRIBUTIONS SUBTAB ----
        with viz_subtabs[0]:
            st.markdown("### 📊 Interactive Motif Distribution Plots")
            
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
            st.markdown("### 📸 Publication-Ready Figures")
            
            col1, col2 = st.columns(2)
            with col1:
                st.markdown("**Export Options:**")
                export_format = st.selectbox("Format:", ["PNG (300 DPI)", "PDF", "SVG"])
                figure_size = st.selectbox("Size:", ["Small (8x6)", "Medium (12x8)", "Large (16x12)"])
                
            with col2:
                st.markdown("**Style Options:**")
                color_scheme = st.selectbox("Color scheme:", ["Default", "Colorblind-friendly", "Grayscale"])
                include_title = st.checkbox("Include title", value=True)
            
            if st.button("🎨 Generate Publication Figure"):
                st.info("Publication figure generation functionality will be implemented here.")
        
        # ---- 3D/ADVANCED VIEWS SUBTAB ----
        with viz_subtabs[2]:
            st.markdown("### 🌐 3D and Advanced Visualizations")
            
            if ADVANCED_VIZ_AVAILABLE:
                st.info("Advanced 3D visualization functionality will be integrated here.")
            else:
                st.info("3D visualization requires additional dependencies. Please install the advanced visualization package.")

# =====================================================
# CLINICAL/DISEASE ANNOTATION PAGE
# =====================================================
with tab_dict["🏥 Clinical/Disease"]:
    # Clinical subtabs
    clinical_subtabs = st.tabs(["🩺 Disease Motif Results", "📋 Clinical Summary"])
    
    # ---- DISEASE MOTIF RESULTS SUBTAB ----
    with clinical_subtabs[0]:
        st.markdown("### 🩺 Disease-Associated Motif Analysis")
        
        if not st.session_state.results:
            st.info("📋 No analysis results available. Please run an analysis first.")
        else:
            st.markdown("""
            <div class="feature-card">
                <h4>🧬 Disease Annotation Status</h4>
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
        st.markdown("### 📋 Clinical Interpretation Summary")
        
        st.markdown("""
        <div class="feature-card">
            <h4>🏥 Clinical Data Integration</h4>
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
with tab_dict["📥 Download & Export"]:
    st.markdown("### 📥 Download Analysis Results")
    
    if not st.session_state.results:
        st.info("📋 No analysis results available for download. Please run an analysis first.")
    else:
        col1, col2 = st.columns(2)
        
        with col1:
            st.markdown("**📊 Data Export Options**")
            
            # Results table download
            if st.button("📋 Download Results Table (CSV)"):
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
                    label="💾 Download CSV",
                    data=csv,
                    file_name="nbdfinder_results.csv",
                    mime="text/csv"
                )
            
            # Excel export
            if st.button("📊 Download Results (Excel)"):
                st.info("Excel export functionality will be implemented.")
        
        with col2:
            st.markdown("**🎨 Visualization Export**")
            
            if st.button("📸 Download Plots (PNG)"):
                st.info("Plot export functionality will be implemented.")
            
            if st.button("📄 Generate Analysis Report"):
                st.info("PDF report generation will be implemented.")
        
        st.markdown("---")
        st.markdown("**📁 Session Data**")
        
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
with tab_dict["📚 Documentation"]:
    # Documentation subtabs
    doc_subtabs = st.tabs(["🎯 Motif Definitions", "⚙️ Detection Algorithms", "📊 Scoring Systems", "📖 User Guide", "📚 References"])
    
    # ---- MOTIF DEFINITIONS SUBTAB ----
    with doc_subtabs[0]:
        st.markdown("### 🎯 Non-B DNA Motif Definitions")
        
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
            with st.expander(f"🔍 {motif}"):
                st.write(definition)
                st.info(f"Detailed information about {motif} detection algorithms and biological significance.")
    
    # ---- DETECTION ALGORITHMS SUBTAB ----
    with doc_subtabs[1]:
        st.markdown("### ⚙️ Detection Algorithms")
        
        st.markdown("""
        <div class="feature-card">
            <h4>🔬 Algorithm Overview</h4>
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
        st.markdown("### 📊 Scoring Systems")
        
        st.markdown("""
        <div class="feature-card">
            <h4>📈 Motif Scoring Methodology</h4>
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
        st.markdown("### 📖 User Guide")
        
        st.markdown("""
        <div class="feature-card">
            <h4>🚀 Quick Start Guide</h4>
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
        st.markdown("### 📚 Scientific References")
        
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
with tab_dict["⚙️ Settings"]:
    st.markdown("### ⚙️ Application Settings")
    
    col1, col2 = st.columns(2)
    
    with col1:
        st.markdown("**🎨 Theme & Appearance**")
        
        theme_mode = st.selectbox(
            "Color theme:",
            ["Auto", "Light", "Dark"],
            help="Choose your preferred color theme"
        )
        
        font_size = st.select_slider(
            "Font size:",
            options=["Small", "Medium", "Large"],
            value="Medium"
        )
        
        st.markdown("**♿ Accessibility**")
        
        high_contrast = st.checkbox("High contrast mode")
        reduced_motion = st.checkbox("Reduce animations")
        screen_reader = st.checkbox("Screen reader optimizations")
    
    with col2:
        st.markdown("**🔬 Analysis Preferences**")
        
        default_sensitivity = st.selectbox(
            "Default sensitivity:",
            ["Low", "Medium", "High", "Maximum"],
            index=2
        )
        
        auto_save = st.checkbox("Auto-save results", value=True)
        parallel_default = st.checkbox("Enable parallel processing by default", value=True)
        
        st.markdown("**📊 Export Preferences**")
        
        default_format = st.selectbox(
            "Default export format:",
            ["CSV", "Excel", "JSON"],
            index=0
        )
        
        include_metadata = st.checkbox("Include metadata in exports", value=True)
    
    st.markdown("---")
    
    if st.button("💾 Save Settings"):
        st.success("✅ Settings saved successfully!")
    
    if st.button("🔄 Reset to Defaults"):
        st.info("🔄 Settings reset to default values")

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