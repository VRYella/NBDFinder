"""
NBDFinder - Optimized Non-B DNA Analysis Platform
================================================

Streamlined for complete mode analysis with publication-quality results.
Focus on scientific accuracy and performance without decorative elements.

Author: Dr. Venkata Rajesh Yella
License: Academic Use
"""

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

# Essential configuration
st.set_page_config(
    page_title="NBDFinder - Non-B DNA Analysis Platform",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Scientific styling - clean and professional
st.markdown("""
<style>
body, [data-testid="stAppViewContainer"] {
    background: #ffffff;
    font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif;
    color: #1f2937;
    line-height: 1.6;
}

.stTabs [data-baseweb="tab-list"] {
    background: #f9fafb;
    border-bottom: 2px solid #e5e7eb;
    gap: 0;
}

.stTabs [data-baseweb="tab"] {
    background: transparent;
    border: none;
    color: #6b7280;
    font-weight: 500;
    padding: 12px 24px;
}

.stTabs [aria-selected="true"] {
    background: #ffffff;
    color: #1f2937;
    border-bottom: 2px solid #3b82f6;
}

.metric-container {
    background: #f9fafb;
    border: 1px solid #e5e7eb;
    border-radius: 6px;
    padding: 16px;
    margin: 8px 0;
}

.analysis-section {
    background: #ffffff;
    border: 1px solid #e5e7eb;
    border-radius: 8px;
    padding: 20px;
    margin: 16px 0;
}
</style>
""", unsafe_allow_html=True)

# Example sequences for testing
EXAMPLE_FASTA = """>Example Sequence - G4-rich region
ATCGATCGATCGAAAATTTTATTTAAATTTAAATTTGGGTTAGGGTTAGGGTTAGGGCCCCCTCCCCCTCCCCCTCCCC
ATCGATCGCGCGCGCGATCGCACACACACAGCTGCTGCTGCTTGGGAAAGGGGAAGGGTTAGGGAAAGGGGTTT
GGGTTTAGGGGGGAGGGGCTGCTGCTGCATGCGGGAAGGGAGGGTAGAGGGTCCGGTAGGAACCCCTAACCCCTAA
GAAAGAAGAAGAAGAAGAAGAAAGGAAGGAAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGAGGG
CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC"""

def ensure_subtype(motif):
    """Ensure every motif has a valid subtype"""
    if isinstance(motif, dict):
        if 'Subtype' not in motif or motif['Subtype'] is None:
            motif['Subtype'] = motif.get('Class', 'Unknown')
        return motif
    return {'Subtype': 'Unknown', 'Motif': str(motif)}

def analyze_sequence_complete(seq, seq_name, selected_motifs):
    """
    Complete mode analysis with hybrid and cluster detection.
    Ensures all overlapping classes are detected as hybrids (2+)
    and clusters use 3+ motifs in 100nt window.
    """
    if not seq or len(seq) < 10:
        return []
    
    # Clean sequence
    seq = parse_fasta(seq) if '>' in seq else seq.upper().replace(' ', '').replace('U', 'T')
    
    # Complete analysis mode - no fast mode shortcuts
    motif_results = all_motifs(
        seq, 
        nonoverlap=False, 
        report_hotspots=True,  # Always include clusters
        sequence_name=seq_name, 
        fast_mode=False  # Always complete mode
    )
    
    # Ensure all motifs have proper subtype
    motif_results = [ensure_subtype(m) for m in motif_results]
    
    # Map UI names to actual class names for filtering
    ui_to_class_mapping = {
        "Canonical G4": "G-Quadruplex Family",
        "Relaxed G4": "G-Quadruplex Family", 
        "Bulged G4": "G-Quadruplex Family",
        "Imperfect G4": "G-Quadruplex Family",
        "Bipartite G4": "G-Quadruplex Family",
        "Multimeric G4": "G-Quadruplex Family",
        "G-Triplex": "G-Quadruplex Family",
        "i-Motif": "i-motif family",
        "AC-Motif": "i-motif family",
        "Curved DNA": "Curved DNA",
        "Cruciform": "Cruciform DNA",
        "Z-DNA": "Z-DNA",
        "eGZ (Extruded-G)": "Z-DNA",
        "R-Loop": "R-loop",
        "Triplex DNA": "Triplex",
        "Slipped DNA": "Slipped DNA",
        "Sticky DNA": "Triplex",
        "Hybrid": "Hybrid",
        "Non-B DNA Clusters": "Non-B DNA Clusters"
    }
    
    # Filter by selected motifs if specified
    if selected_motifs and len(selected_motifs) > 0:
        # Map UI names to actual class names
        actual_classes = set()
        for ui_name in selected_motifs:
            if ui_name in ui_to_class_mapping:
                actual_classes.add(ui_to_class_mapping[ui_name])
        
        if actual_classes:
            motif_results = [m for m in motif_results if m.get('Class', '') in actual_classes]
    
    return motif_results

def create_results_summary(motifs):
    """Create scientific summary of results"""
    if not motifs:
        return "No motifs detected."
    
    # Count by class
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    total_motifs = len(motifs)
    
    # Hybrid and cluster analysis
    hybrid_count = class_counts.get('Hybrid', 0)
    cluster_count = class_counts.get('Non-B DNA Clusters', 0)
    
    summary = f"""
    **Analysis Summary:**
    - Total motifs detected: {total_motifs}
    - Motif classes: {len(class_counts)}
    - Hybrid structures: {hybrid_count}
    - Cluster regions: {cluster_count}
    """
    
    if hybrid_count > 0:
        hybrid_motifs = [m for m in motifs if m.get('Class') == 'Hybrid']
        class_combos = set()
        for h in hybrid_motifs:
            if 'MotifClasses' in h:
                combo = tuple(sorted(h['MotifClasses']))
                class_combos.add(combo)
        summary += f"\n    - Hybrid class combinations: {len(class_combos)}"
    
    if cluster_count > 0:
        cluster_motifs = [m for m in motifs if m.get('Class') == 'Non-B DNA Clusters']
        avg_density = np.mean([m.get('MotifCount', 0) for m in cluster_motifs])
        summary += f"\n    - Average cluster density: {avg_density:.1f} motifs per region"
    
    return summary

def create_class_distribution_chart(motifs):
    """Create publication-quality class distribution chart"""
    if not motifs:
        return None
    
    class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
    
    # Scientific color palette (colorblind-friendly)
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    fig = go.Figure(data=[
        go.Bar(
            x=list(class_counts.keys()),
            y=list(class_counts.values()),
            marker_color=colors[:len(class_counts)],
            text=list(class_counts.values()),
            textposition='auto',
        )
    ])
    
    fig.update_layout(
        title="Non-B DNA Motif Class Distribution",
        xaxis_title="Motif Class",
        yaxis_title="Count",
        font=dict(family="Arial", size=12),
        plot_bgcolor='white',
        paper_bgcolor='white',
        height=400
    )
    
    fig.update_xaxis(tickangle=45)
    
    return fig

def create_position_plot(motifs, seq_len):
    """Create genomic position plot"""
    if not motifs:
        return None
    
    fig = go.Figure()
    
    # Group by class for color coding
    class_groups = {}
    for motif in motifs:
        cls = motif.get('Class', 'Unknown')
        if cls not in class_groups:
            class_groups[cls] = []
        class_groups[cls].append(motif)
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    y_pos = 0
    for i, (cls, motifs_in_class) in enumerate(class_groups.items()):
        for motif in motifs_in_class:
            start = motif.get('Start', 0)
            end = motif.get('End', 0)
            
            fig.add_trace(go.Scatter(
                x=[start, end],
                y=[y_pos, y_pos],
                mode='lines',
                line=dict(color=colors[i % len(colors)], width=4),
                name=cls,
                showlegend=(motif == motifs_in_class[0]),  # Only show legend once per class
                hovertemplate=f"<b>{cls}</b><br>Position: {start}-{end}<br>Length: {end-start+1}<extra></extra>"
            ))
        y_pos += 1
    
    fig.update_layout(
        title="Genomic Position Map",
        xaxis_title="Position (bp)",
        yaxis_title="Motif Class",
        font=dict(family="Arial", size=12),
        plot_bgcolor='white',
        paper_bgcolor='white',
        height=400,
        yaxis=dict(showticklabels=False)
    )
    
    return fig

def main():
    """Main application"""
    
    # Initialize session state
    if 'motif_results' not in st.session_state:
        st.session_state.motif_results = []
    if 'analyzed_sequences' not in st.session_state:
        st.session_state.analyzed_sequences = {}
    
    # Main navigation
    tab1, tab2, tab3, tab4, tab5 = st.tabs([
        "Home", "Upload & Analyze", "Results", "Download", "Documentation"
    ])
    
    with tab1:
        st.title("NBDFinder: Non-B DNA Analysis Platform")
        
        col1, col2 = st.columns([2, 1])
        
        with col1:
            st.markdown("""
            ### Complete Mode Analysis Platform
            
            Advanced computational tool for detecting and analyzing non-canonical DNA structures.
            Optimized for publication-quality results with comprehensive motif detection.
            
            **Key Features:**
            - **10 Non-B DNA Classes**: Comprehensive detection suite
            - **Hybrid Detection**: All overlapping classes (2+) automatically detected
            - **Cluster Analysis**: Regions with 3+ motifs in 100nt windows
            - **Publication Quality**: Export-ready visualizations and data
            - **Scientific Accuracy**: Validated algorithms with statistical rigor
            """)
        
        with col2:
            st.markdown("""
            **Detected Classes:**
            1. Curved DNA
            2. Slipped DNA  
            3. Cruciform DNA
            4. R-loop
            5. Triplex
            6. G-Quadruplex Family
            7. i-motif family
            8. Z-DNA
            9. Hybrid Structures
            10. Non-B DNA Clusters
            """)
    
    with tab2:
        st.header("Sequence Analysis")
        
        # Motif selection
        st.subheader("Motif Classes for Analysis")
        
        col1, col2, col3 = st.columns([1, 1, 1])
        
        with col1:
            g4_family = st.multiselect(
                "G-Quadruplex Family",
                ["Canonical G4", "Relaxed G4", "Bulged G4", "Imperfect G4", 
                 "Bipartite G4", "Multimeric G4", "G-Triplex"],
                default=["Canonical G4", "Relaxed G4", "G-Triplex"]
            )
            
            alternative = st.multiselect(
                "Alternative Conformations",
                ["Curved DNA", "Cruciform", "Z-DNA", "eGZ (Extruded-G)"],
                default=["Curved DNA", "Z-DNA"]
            )
        
        with col2:
            imotif_family = st.multiselect(
                "i-Motif Family",
                ["i-Motif", "AC-Motif"],
                default=["i-Motif"]
            )
            
            complex_structures = st.multiselect(
                "Complex Structures", 
                ["R-Loop", "Triplex DNA", "Slipped DNA", "Sticky DNA"],
                default=["R-Loop", "Triplex DNA"]
            )
        
        with col3:
            composite = st.multiselect(
                "Composite Analysis",
                ["Hybrid", "Non-B DNA Clusters"],
                default=["Hybrid", "Non-B DNA Clusters"]
            )
        
        # Combine selected motifs
        selected_motifs = g4_family + imotif_family + alternative + complex_structures + composite
        
        if st.button("Select All Classes"):
            selected_motifs = [
                "Canonical G4", "Relaxed G4", "Bulged G4", "Imperfect G4", 
                "Bipartite G4", "Multimeric G4", "G-Triplex", "i-Motif", "AC-Motif",
                "Curved DNA", "Cruciform", "Z-DNA", "eGZ (Extruded-G)",
                "R-Loop", "Triplex DNA", "Slipped DNA", "Sticky DNA",
                "Hybrid", "Non-B DNA Clusters"
            ]
            st.rerun()
        
        st.markdown(f"**Selected classes:** {len(selected_motifs)}")
        
        # Input method
        st.subheader("Input Method")
        
        input_method = st.radio(
            "Choose input method:",
            ["Upload FASTA File", "Paste Sequence", "Example Sequence", "NCBI Fetch"],
            horizontal=True
        )
        
        sequences = []
        sequence_names = []
        
        if input_method == "Upload FASTA File":
            uploaded_file = st.file_uploader(
                "Upload FASTA file",
                type=['fa', 'fasta', 'txt'],
                help="Single or multi-FASTA format"
            )
            
            if uploaded_file:
                content = uploaded_file.read().decode('utf-8')
                if '>' in content:
                    # Multi-FASTA
                    current_seq = ""
                    current_name = ""
                    for line in content.split('\n'):
                        line = line.strip()
                        if line.startswith('>'):
                            if current_seq and current_name:
                                sequences.append(current_seq)
                                sequence_names.append(current_name)
                            current_name = line[1:].split()[0]
                            current_seq = ""
                        else:
                            current_seq += line
                    if current_seq and current_name:
                        sequences.append(current_seq)
                        sequence_names.append(current_name)
                else:
                    sequences = [content]
                    sequence_names = [uploaded_file.name]
        
        elif input_method == "Paste Sequence":
            seq_input = st.text_area(
                "Paste your DNA sequence(s)",
                height=200,
                help="FASTA format or raw sequence"
            )
            
            if seq_input:
                if '>' in seq_input:
                    # Parse as FASTA
                    current_seq = ""
                    current_name = ""
                    for line in seq_input.split('\n'):
                        line = line.strip()
                        if line.startswith('>'):
                            if current_seq and current_name:
                                sequences.append(current_seq)
                                sequence_names.append(current_name)
                            current_name = line[1:].split()[0]
                            current_seq = ""
                        else:
                            current_seq += line
                    if current_seq and current_name:
                        sequences.append(current_seq)
                        sequence_names.append(current_name)
                else:
                    sequences = [seq_input]
                    sequence_names = ["Input_Sequence"]
        
        elif input_method == "Example Sequence":
            if st.button("Load Example"):
                seq = parse_fasta(EXAMPLE_FASTA)
                sequences = [seq]
                sequence_names = ["Example_G4_Rich"]
        
        elif input_method == "NCBI Fetch":
            accession = st.text_input("Enter NCBI accession number or gene name:")
            if accession and st.button("Fetch"):
                try:
                    Entrez.email = "user@example.com"
                    handle = Entrez.esearch(db="nucleotide", term=accession)
                    search_results = Entrez.read(handle)
                    handle.close()
                    
                    if search_results['IdList']:
                        seq_id = search_results['IdList'][0]
                        handle = Entrez.efetch(db="nucleotide", id=seq_id, rettype="fasta")
                        fasta_content = handle.read()
                        handle.close()
                        
                        seq = parse_fasta(fasta_content)
                        sequences = [seq]
                        sequence_names = [accession]
                        st.success(f"Successfully fetched {len(seq):,} bp")
                    else:
                        st.error("No sequences found for this identifier")
                except Exception as e:
                    st.error(f"NCBI fetch failed: {str(e)}")
        
        # Analysis button
        if sequences and st.button("Analyze Sequences", type="primary"):
            with st.spinner("Analyzing sequences..."):
                st.session_state.motif_results = []
                st.session_state.analyzed_sequences = {}
                
                for seq, seq_name in zip(sequences, sequence_names):
                    if len(seq) >= 10:
                        results = analyze_sequence_complete(seq, seq_name, selected_motifs)
                        st.session_state.motif_results.extend(results)
                        st.session_state.analyzed_sequences[seq_name] = {
                            'sequence': seq,
                            'length': len(seq),
                            'gc_content': gc_content(seq),
                            'motif_count': len(results)
                        }
                
                st.success(f"Analysis complete! Found {len(st.session_state.motif_results)} motifs across {len(sequences)} sequence(s).")
    
    with tab3:
        st.header("Analysis Results")
        
        if not st.session_state.motif_results:
            st.info("No analysis results available. Please analyze sequences first.")
            return
        
        # Results organization by analysis type
        result_view = st.selectbox(
            "View Results By:",
            ["Overview", "Motif Classes", "Genomic Position", "Hybrid Analysis", "Cluster Analysis", "Statistical Summary"]
        )
        
        motifs = st.session_state.motif_results
        
        if result_view == "Overview":
            st.subheader("Analysis Overview")
            
            col1, col2 = st.columns([1, 1])
            
            with col1:
                st.markdown(create_results_summary(motifs))
                
                # Sequence statistics
                st.subheader("Sequence Statistics")
                seq_data = []
                for seq_name, info in st.session_state.analyzed_sequences.items():
                    seq_data.append({
                        'Sequence': seq_name,
                        'Length (bp)': f"{info['length']:,}",
                        'GC Content (%)': f"{info['gc_content']:.1f}",
                        'Motifs Found': info['motif_count']
                    })
                
                if seq_data:
                    st.dataframe(pd.DataFrame(seq_data), use_container_width=True)
            
            with col2:
                # Class distribution chart
                fig = create_class_distribution_chart(motifs)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
        
        elif result_view == "Motif Classes":
            st.subheader("Motif Classes Analysis")
            
            class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
            
            for cls in sorted(class_counts.keys()):
                with st.expander(f"{cls} ({class_counts[cls]} motifs)"):
                    class_motifs = [m for m in motifs if m.get('Class') == cls]
                    
                    # Create dataframe for this class
                    class_data = []
                    for motif in class_motifs:
                        class_data.append({
                            'Sequence': motif.get('Sequence Name', ''),
                            'Subtype': motif.get('Subtype', ''),
                            'Start': motif.get('Start', 0),
                            'End': motif.get('End', 0),
                            'Length': motif.get('Length', 0),
                            'Score': f"{motif.get('Score', 0):.3f}"
                        })
                    
                    if class_data:
                        st.dataframe(pd.DataFrame(class_data), use_container_width=True)
        
        elif result_view == "Genomic Position":
            st.subheader("Genomic Position Map")
            
            if st.session_state.analyzed_sequences:
                seq_name = st.selectbox("Select sequence:", list(st.session_state.analyzed_sequences.keys()))
                seq_motifs = [m for m in motifs if m.get('Sequence Name') == seq_name]
                seq_len = st.session_state.analyzed_sequences[seq_name]['length']
                
                fig = create_position_plot(seq_motifs, seq_len)
                if fig:
                    st.plotly_chart(fig, use_container_width=True)
        
        elif result_view == "Hybrid Analysis":
            st.subheader("Hybrid Structure Analysis")
            
            hybrid_motifs = [m for m in motifs if m.get('Class') == 'Hybrid']
            
            if hybrid_motifs:
                st.markdown(f"**Found {len(hybrid_motifs)} hybrid structures**")
                
                hybrid_data = []
                for h in hybrid_motifs:
                    classes = h.get('MotifClasses', [])
                    class_breakdown = h.get('ClassBreakdown', {})
                    
                    hybrid_data.append({
                        'Sequence': h.get('Sequence Name', ''),
                        'Position': f"{h.get('Start', 0)}-{h.get('End', 0)}",
                        'Length': h.get('Length', 0),
                        'Classes Involved': len(classes),
                        'Class Types': ', '.join(classes),
                        'Class Counts': ', '.join([f"{k}({v})" for k, v in class_breakdown.items()]),
                        'Score': f"{h.get('Score', 0):.3f}"
                    })
                
                st.dataframe(pd.DataFrame(hybrid_data), use_container_width=True)
                
                # Hybrid class combination analysis
                st.subheader("Class Combination Frequency")
                combo_counts = Counter()
                for h in hybrid_motifs:
                    classes = tuple(sorted(h.get('MotifClasses', [])))
                    combo_counts[classes] += 1
                
                combo_data = []
                for combo, count in combo_counts.most_common():
                    combo_data.append({
                        'Class Combination': ' + '.join(combo),
                        'Frequency': count,
                        'Percentage': f"{count/len(hybrid_motifs)*100:.1f}%"
                    })
                
                st.dataframe(pd.DataFrame(combo_data), use_container_width=True)
            else:
                st.info("No hybrid structures detected.")
        
        elif result_view == "Cluster Analysis":
            st.subheader("Non-B DNA Cluster Analysis")
            
            cluster_motifs = [m for m in motifs if m.get('Class') == 'Non-B DNA Clusters']
            
            if cluster_motifs:
                st.markdown(f"**Found {len(cluster_motifs)} cluster regions**")
                
                cluster_data = []
                for c in cluster_motifs:
                    cluster_data.append({
                        'Sequence': c.get('Sequence Name', ''),
                        'Position': f"{c.get('Start', 0)}-{c.get('End', 0)}",
                        'Length': c.get('Length', 0),
                        'Motif Count': c.get('MotifCount', 0),
                        'Class Diversity': c.get('ClassDiversity', 0),
                        'Complexity': c.get('ComplexityLevel', 'Unknown'),
                        'Score': f"{c.get('Score', 0):.3f}"
                    })
                
                st.dataframe(pd.DataFrame(cluster_data), use_container_width=True)
                
                # Cluster statistics
                st.subheader("Cluster Statistics")
                col1, col2, col3 = st.columns(3)
                
                with col1:
                    avg_length = np.mean([c.get('Length', 0) for c in cluster_motifs])
                    st.metric("Average Length", f"{avg_length:.1f} bp")
                
                with col2:
                    avg_density = np.mean([c.get('MotifCount', 0) for c in cluster_motifs])
                    st.metric("Average Density", f"{avg_density:.1f} motifs")
                
                with col3:
                    avg_diversity = np.mean([c.get('ClassDiversity', 0) for c in cluster_motifs])
                    st.metric("Average Class Diversity", f"{avg_diversity:.1f}")
            else:
                st.info("No cluster regions detected.")
        
        elif result_view == "Statistical Summary":
            st.subheader("Statistical Summary")
            
            # Overall statistics
            col1, col2, col3, col4 = st.columns(4)
            
            with col1:
                st.metric("Total Motifs", len(motifs))
            
            with col2:
                unique_classes = len(set(m.get('Class', 'Unknown') for m in motifs))
                st.metric("Unique Classes", unique_classes)
            
            with col3:
                hybrid_count = len([m for m in motifs if m.get('Class') == 'Hybrid'])
                st.metric("Hybrid Structures", hybrid_count)
            
            with col4:
                cluster_count = len([m for m in motifs if m.get('Class') == 'Non-B DNA Clusters'])
                st.metric("Cluster Regions", cluster_count)
            
            # Length distribution
            lengths = [m.get('Length', 0) for m in motifs if m.get('Length', 0) > 0]
            if lengths:
                st.subheader("Motif Length Distribution")
                
                fig = go.Figure(data=[go.Histogram(x=lengths, nbinsx=20)])
                fig.update_layout(
                    title="Motif Length Distribution",
                    xaxis_title="Length (bp)",
                    yaxis_title="Frequency",
                    font=dict(family="Arial", size=12),
                    plot_bgcolor='white',
                    paper_bgcolor='white',
                    height=400
                )
                st.plotly_chart(fig, use_container_width=True)
                
                # Length statistics
                col1, col2, col3, col4 = st.columns(4)
                with col1:
                    st.metric("Mean Length", f"{np.mean(lengths):.1f} bp")
                with col2:
                    st.metric("Median Length", f"{np.median(lengths):.1f} bp")
                with col3:
                    st.metric("Min Length", f"{min(lengths)} bp")
                with col4:
                    st.metric("Max Length", f"{max(lengths)} bp")
    
    with tab4:
        st.header("Download Results")
        
        if not st.session_state.motif_results:
            st.info("No results available for download. Please analyze sequences first.")
            return
        
        motifs = st.session_state.motif_results
        
        # Create downloadable dataframe
        download_data = []
        for motif in motifs:
            download_data.append({
                'Sequence_Name': motif.get('Sequence Name', ''),
                'Class': motif.get('Class', ''),
                'Subtype': motif.get('Subtype', ''),
                'Start': motif.get('Start', 0),
                'End': motif.get('End', 0),
                'Length': motif.get('Length', 0),
                'Score': motif.get('Score', 0),
                'Sequence': motif.get('Sequence', ''),
                'Additional_Info': str(motif.get('Arms/Repeat Unit/Copies', ''))
            })
        
        df = pd.DataFrame(download_data)
        
        # Download options
        col1, col2, col3 = st.columns(3)
        
        with col1:
            csv_data = df.to_csv(index=False)
            st.download_button(
                label="Download CSV",
                data=csv_data,
                file_name="nbdfinder_results.csv",
                mime="text/csv"
            )
        
        with col2:
            excel_buffer = io.BytesIO()
            df.to_excel(excel_buffer, index=False, engine='openpyxl')
            excel_data = excel_buffer.getvalue()
            
            st.download_button(
                label="Download Excel",
                data=excel_data,
                file_name="nbdfinder_results.xlsx",
                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
            )
        
        with col3:
            # Export summary report
            summary_report = f"""NBDFinder Analysis Report
Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}

ANALYSIS SUMMARY:
- Total motifs detected: {len(motifs)}
- Unique classes: {len(set(m.get('Class', 'Unknown') for m in motifs))}
- Sequences analyzed: {len(st.session_state.analyzed_sequences)}

CLASS DISTRIBUTION:
"""
            class_counts = Counter(m.get('Class', 'Unknown') for m in motifs)
            for cls, count in class_counts.most_common():
                summary_report += f"- {cls}: {count}\n"
            
            st.download_button(
                label="Download Report",
                data=summary_report,
                file_name="nbdfinder_summary.txt",
                mime="text/plain"
            )
        
        # Preview data
        st.subheader("Data Preview")
        st.dataframe(df, use_container_width=True)
    
    with tab5:
        st.header("Documentation")
        
        st.markdown("""
        ## NBDFinder Complete Mode Analysis
        
        ### Overview
        NBDFinder is optimized for complete mode analysis with focus on scientific accuracy and publication-quality results.
        
        ### Key Features
        
        **Hybrid Structure Detection:**
        - Automatically detects all overlapping motif classes (2 or more)
        - Provides detailed class breakdown and interaction scores
        - Identifies structural complexity levels
        
        **Cluster Analysis:**
        - Detects regions with 3+ motifs in 100bp windows
        - Calculates class diversity and complexity metrics
        - Provides statistical scoring for hotspot significance
        
        **Publication Quality:**
        - Export-ready visualizations and data tables
        - Statistical rigor in all calculations
        - Standardized output formats for scientific reporting
        
        ### Algorithm Details
        
        **Complete Mode Processing:**
        1. Full motif detection across all 10 classes
        2. Hybrid identification for overlapping structures
        3. Cluster detection using sliding window approach
        4. Statistical validation and scoring
        
        **Quality Assurance:**
        - Validated against experimental datasets
        - Peer-reviewed algorithmic approaches
        - Standardized output formatting
        
        ### Citation
        Please cite NBDFinder in your publications:
        
        > Yella VR. NBDFinder: A comprehensive platform for non-B DNA structure prediction and analysis. 
        > Bioinformatics. 2024.
        
        ### Contact
        Dr. Venkata Rajesh Yella  
        Email: yvrajesh_bt@kluniversity.in  
        GitHub: [VRYella](https://github.com/VRYella)
        """)

if __name__ == "__main__":
    main()