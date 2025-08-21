"""
Advanced Visualization Module for NBDFinder
===========================================

This module provides publication-quality visualizations and interactive
features for the NBDFinder web interface, including 3D structure models,
advanced plotting, and real-time analysis capabilities.

Key Features:
- Publication-quality figure generation
- Interactive 3D structure visualization
- Advanced statistical plots
- Real-time collaborative analysis
- Accessibility optimized design
- Professional color schemes

Author: Dr. Venkata Rajesh Yella
Updated: 2024 with latest visualization techniques
License: Academic Use
"""

import streamlit as st
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from typing import Dict, List, Optional, Tuple
import io
import base64

# Professional color schemes for accessibility
PROFESSIONAL_COLORS = {
    'primary': '#1565c0',
    'secondary': '#2e7d32', 
    'accent': '#d32f2f',
    'warning': '#f57c00',
    'success': '#388e3c',
    'info': '#1976d2',
    'disease': '#c62828',
    'ml': '#7b1fa2',
    'cluster': '#f57c00',
    'hybrid': '#5e35b1'
}

# Colorblind-friendly palette
COLORBLIND_PALETTE = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
    '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf'
]

class AdvancedVisualizer:
    """Advanced visualization engine for NBDFinder"""
    
    def __init__(self):
        self.color_palette = COLORBLIND_PALETTE
        self.professional_colors = PROFESSIONAL_COLORS
        
    def create_publication_quality_plot(self, data: pd.DataFrame, plot_type: str, **kwargs) -> go.Figure:
        """
        Create publication-quality plots with professional styling
        
        Args:
            data: DataFrame with motif data
            plot_type: Type of plot to create
            **kwargs: Additional plotting parameters
            
        Returns:
            Plotly figure object
        """
        if plot_type == "motif_distribution":
            return self._create_motif_distribution_plot(data, **kwargs)
        elif plot_type == "genomic_map":
            return self._create_genomic_map(data, **kwargs)
        elif plot_type == "clinical_analysis":
            return self._create_clinical_analysis_plot(data, **kwargs)
        elif plot_type == "ml_predictions":
            return self._create_ml_predictions_plot(data, **kwargs)
        elif plot_type == "cluster_analysis":
            return self._create_cluster_analysis_plot(data, **kwargs)
        elif plot_type == "3d_structure":
            return self._create_3d_structure_plot(data, **kwargs)
        else:
            raise ValueError(f"Unknown plot type: {plot_type}")
    
    def _create_motif_distribution_plot(self, data: pd.DataFrame, **kwargs) -> go.Figure:
        """Create advanced motif distribution visualization"""
        # Count motifs by class and subtype
        class_counts = data['Class'].value_counts()
        subtype_counts = data['Subtype'].value_counts()
        
        # Create subplot with secondary y-axis
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                "Motif Class Distribution", 
                "Top Motif Subtypes",
                "Score Distribution by Class",
                "Length vs Score Analysis"
            ),
            specs=[
                [{"type": "pie"}, {"type": "bar"}],
                [{"type": "box"}, {"type": "scatter"}]
            ]
        )
        
        # Pie chart for class distribution
        fig.add_trace(
            go.Pie(
                labels=class_counts.index,
                values=class_counts.values,
                hole=0.4,
                marker_colors=self.color_palette[:len(class_counts)],
                textinfo='label+percent',
                textfont_size=12
            ),
            row=1, col=1
        )
        
        # Bar chart for top subtypes
        top_subtypes = subtype_counts.head(10)
        fig.add_trace(
            go.Bar(
                x=top_subtypes.values,
                y=top_subtypes.index,
                orientation='h',
                marker_color=self.professional_colors['secondary'],
                text=top_subtypes.values,
                textposition='outside'
            ),
            row=1, col=2
        )
        
        # Box plot for score distribution by class
        for i, motif_class in enumerate(class_counts.index):
            class_data = data[data['Class'] == motif_class]
            if 'Score' in class_data.columns:
                scores = pd.to_numeric(class_data['Score'], errors='coerce').dropna()
                if len(scores) > 0:
                    fig.add_trace(
                        go.Box(
                            y=scores,
                            name=motif_class,
                            marker_color=self.color_palette[i % len(self.color_palette)]
                        ),
                        row=2, col=1
                    )
        
        # Scatter plot for length vs score
        if 'Length' in data.columns and 'Score' in data.columns:
            lengths = pd.to_numeric(data['Length'], errors='coerce')
            scores = pd.to_numeric(data['Score'], errors='coerce')
            
            fig.add_trace(
                go.Scatter(
                    x=lengths,
                    y=scores,
                    mode='markers',
                    marker=dict(
                        size=8,
                        color=scores,
                        colorscale='viridis',
                        colorbar=dict(title="Score"),
                        opacity=0.7
                    ),
                    text=data['Subtype'],
                    hovertemplate="<b>%{text}</b><br>Length: %{x}<br>Score: %{y}<extra></extra>"
                ),
                row=2, col=2
            )
        
        # Update layout with professional styling
        fig.update_layout(
            title=dict(
                text="<b>Comprehensive Motif Analysis Dashboard</b>",
                x=0.5,
                font=dict(size=20, color=self.professional_colors['primary'])
            ),
            showlegend=False,
            height=800,
            plot_bgcolor='white',
            paper_bgcolor='white',
            font=dict(family="Arial, sans-serif", size=12)
        )
        
        return fig
    
    def _create_genomic_map(self, data: pd.DataFrame, sequence_length: int, **kwargs) -> go.Figure:
        """Create interactive genomic map visualization"""
        fig = go.Figure()
        
        # Color mapping for different motif classes
        class_colors = {
            'G-Quadruplex': self.professional_colors['primary'],
            'Disease-Associated Motif': self.professional_colors['disease'],
            'Hybrid Structures': self.professional_colors['hybrid'],
            'Advanced Non-B DNA Clusters': self.professional_colors['cluster'],
            'Z-DNA': self.professional_colors['secondary'],
            'R-Loop': self.professional_colors['accent'],
            'Triplex': self.professional_colors['warning']
        }
        
        # Create tracks for different motif classes
        y_positions = {}
        y_counter = 0
        
        for motif_class in data['Class'].unique():
            class_data = data[data['Class'] == motif_class]
            y_positions[motif_class] = y_counter
            
            # Add motifs as rectangles
            for _, motif in class_data.iterrows():
                start = int(motif['Start'])
                end = int(motif['End'])
                score = float(motif.get('Score', 0))
                
                # Color intensity based on score
                color = class_colors.get(motif_class, self.color_palette[y_counter % len(self.color_palette)])
                
                fig.add_trace(
                    go.Scatter(
                        x=[start, end, end, start, start],
                        y=[y_counter-0.4, y_counter-0.4, y_counter+0.4, y_counter+0.4, y_counter-0.4],
                        fill='toself',
                        fillcolor=color,
                        line=dict(color=color, width=2),
                        opacity=min(1.0, score / 100.0 + 0.3),
                        name=motif_class,
                        hovertemplate=(
                            f"<b>{motif['Subtype']}</b><br>"
                            f"Position: {start}-{end}<br>"
                            f"Length: {end-start}<br>"
                            f"Score: {score:.1f}<br>"
                            f"Class: {motif_class}<extra></extra>"
                        ),
                        showlegend=True if motif_class not in [trace.name for trace in fig.data] else False
                    )
                )
            
            y_counter += 1
        
        # Add sequence scale
        fig.add_trace(
            go.Scatter(
                x=[0, sequence_length],
                y=[-1, -1],
                mode='lines',
                line=dict(color='black', width=3),
                name='Sequence Scale',
                showlegend=False
            )
        )
        
        # Add scale markers
        scale_intervals = max(1, sequence_length // 10)
        for i in range(0, sequence_length + 1, scale_intervals):
            fig.add_annotation(
                x=i,
                y=-1.2,
                text=str(i),
                showarrow=False,
                font=dict(size=10)
            )
        
        fig.update_layout(
            title=dict(
                text="<b>Interactive Genomic Map of Non-B DNA Structures</b>",
                x=0.5,
                font=dict(size=18, color=self.professional_colors['primary'])
            ),
            xaxis=dict(
                title="Genomic Position (bp)",
                range=[0, sequence_length],
                showgrid=True,
                gridcolor='lightgray'
            ),
            yaxis=dict(
                title="Motif Classes",
                tickvals=list(range(len(y_positions))),
                ticktext=list(y_positions.keys()),
                showgrid=False
            ),
            height=600,
            plot_bgcolor='white',
            paper_bgcolor='white',
            hovermode='closest'
        )
        
        return fig
    
    def _create_clinical_analysis_plot(self, data: pd.DataFrame, **kwargs) -> go.Figure:
        """Create clinical significance analysis visualization"""
        # Filter disease-associated motifs
        disease_data = data[data['Class'] == 'Disease-Associated Motif']
        
        if disease_data.empty:
            # Create placeholder plot
            fig = go.Figure()
            fig.add_annotation(
                text="No disease-associated motifs detected",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            return fig
        
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                "Clinical Significance Distribution",
                "Risk Score Analysis",
                "Repeat Count vs Risk",
                "Disease Categories"
            ),
            specs=[
                [{"type": "pie"}, {"type": "histogram"}],
                [{"type": "scatter"}, {"type": "bar"}]
            ]
        )
        
        # Clinical significance pie chart
        if 'Clinical_Significance' in disease_data.columns:
            clin_sig_counts = disease_data['Clinical_Significance'].value_counts()
            colors = [
                self.professional_colors['disease'] if 'Pathogenic' in sig 
                else self.professional_colors['warning'] if 'VUS' in sig
                else self.professional_colors['success']
                for sig in clin_sig_counts.index
            ]
            
            fig.add_trace(
                go.Pie(
                    labels=clin_sig_counts.index,
                    values=clin_sig_counts.values,
                    marker_colors=colors,
                    hole=0.3
                ),
                row=1, col=1
            )
        
        # Risk score histogram
        if 'Risk_Score' in disease_data.columns:
            risk_scores = pd.to_numeric(disease_data['Risk_Score'], errors='coerce').dropna()
            fig.add_trace(
                go.Histogram(
                    x=risk_scores,
                    nbinsx=20,
                    marker_color=self.professional_colors['disease'],
                    opacity=0.7
                ),
                row=1, col=2
            )
        
        # Repeat count vs risk scatter
        if 'Repeat_Count' in disease_data.columns and 'Risk_Score' in disease_data.columns:
            repeat_counts = pd.to_numeric(disease_data['Repeat_Count'], errors='coerce')
            risk_scores = pd.to_numeric(disease_data['Risk_Score'], errors='coerce')
            
            fig.add_trace(
                go.Scatter(
                    x=repeat_counts,
                    y=risk_scores,
                    mode='markers',
                    marker=dict(
                        size=10,
                        color=risk_scores,
                        colorscale='reds',
                        opacity=0.8
                    ),
                    text=disease_data['Disease_Name'],
                    hovertemplate="<b>%{text}</b><br>Repeats: %{x}<br>Risk: %{y}<extra></extra>"
                ),
                row=2, col=1
            )
        
        # Disease categories bar chart
        if 'Disease_Name' in disease_data.columns:
            disease_counts = disease_data['Disease_Name'].value_counts()
            fig.add_trace(
                go.Bar(
                    x=disease_counts.index,
                    y=disease_counts.values,
                    marker_color=self.professional_colors['disease']
                ),
                row=2, col=2
            )
        
        fig.update_layout(
            title=dict(
                text="<b>Clinical Significance Analysis</b>",
                x=0.5,
                font=dict(size=18, color=self.professional_colors['disease'])
            ),
            height=700,
            showlegend=False
        )
        
        return fig
    
    def _create_ml_predictions_plot(self, data: pd.DataFrame, **kwargs) -> go.Figure:
        """Create machine learning predictions visualization"""
        # Filter ML enhanced motifs
        ml_data = data.dropna(subset=['ML_Probability'])
        
        if ml_data.empty:
            fig = go.Figure()
            fig.add_annotation(
                text="No ML predictions available",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            return fig
        
        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=(
                "ML Probability Distribution",
                "Confidence vs Probability",
                "Enhanced Score Comparison",
                "ML Performance by Class"
            )
        )
        
        # ML probability histogram
        ml_probs = pd.to_numeric(ml_data['ML_Probability'], errors='coerce').dropna()
        fig.add_trace(
            go.Histogram(
                x=ml_probs,
                nbinsx=25,
                marker_color=self.professional_colors['ml'],
                opacity=0.7,
                name="ML Probability"
            ),
            row=1, col=1
        )
        
        # Confidence vs Probability scatter
        if 'ML_Confidence' in ml_data.columns:
            ml_conf = pd.to_numeric(ml_data['ML_Confidence'], errors='coerce')
            fig.add_trace(
                go.Scatter(
                    x=ml_probs,
                    y=ml_conf,
                    mode='markers',
                    marker=dict(
                        size=8,
                        color=ml_probs,
                        colorscale='viridis',
                        opacity=0.7
                    ),
                    text=ml_data['Subtype'],
                    hovertemplate="<b>%{text}</b><br>Probability: %{x:.3f}<br>Confidence: %{y:.3f}<extra></extra>"
                ),
                row=1, col=2
            )
        
        # Enhanced vs Original score comparison
        if 'Enhanced_Score' in ml_data.columns and 'Score' in ml_data.columns:
            original_scores = pd.to_numeric(ml_data['Score'], errors='coerce')
            enhanced_scores = pd.to_numeric(ml_data['Enhanced_Score'], errors='coerce')
            
            fig.add_trace(
                go.Scatter(
                    x=original_scores,
                    y=enhanced_scores,
                    mode='markers',
                    marker=dict(
                        size=8,
                        color=self.professional_colors['ml'],
                        opacity=0.7
                    ),
                    name="Score Enhancement"
                ),
                row=2, col=1
            )
            
            # Add diagonal line for reference
            max_score = max(original_scores.max(), enhanced_scores.max())
            fig.add_trace(
                go.Scatter(
                    x=[0, max_score],
                    y=[0, max_score],
                    mode='lines',
                    line=dict(dash='dash', color='gray'),
                    name="No Change Line"
                ),
                row=2, col=1
            )
        
        # ML performance by class
        class_ml_performance = ml_data.groupby('Class')['ML_Probability'].mean().sort_values(ascending=True)
        fig.add_trace(
            go.Bar(
                x=class_ml_performance.values,
                y=class_ml_performance.index,
                orientation='h',
                marker_color=self.professional_colors['ml']
            ),
            row=2, col=2
        )
        
        fig.update_layout(
            title=dict(
                text="<b>Machine Learning Predictions Analysis</b>",
                x=0.5,
                font=dict(size=18, color=self.professional_colors['ml'])
            ),
            height=700,
            showlegend=False
        )
        
        return fig
    
    def _create_cluster_analysis_plot(self, data: pd.DataFrame, **kwargs) -> go.Figure:
        """Create cluster analysis visualization"""
        # Filter cluster data
        cluster_data = data[data['Class'].str.contains('Cluster', na=False)]
        
        if cluster_data.empty:
            fig = go.Figure()
            fig.add_annotation(
                text="No cluster data available",
                xref="paper", yref="paper",
                x=0.5, y=0.5, showarrow=False,
                font=dict(size=16)
            )
            return fig
        
        # Create cluster visualization
        fig = go.Figure()
        
        for _, cluster in cluster_data.iterrows():
            start = int(cluster['Start'])
            end = int(cluster['End'])
            
            fig.add_trace(
                go.Scatter(
                    x=[start, end],
                    y=[0, 0],
                    mode='lines+markers',
                    line=dict(width=10, color=self.professional_colors['cluster']),
                    marker=dict(size=12),
                    name=f"Cluster {cluster.name}",
                    hovertemplate=f"Cluster: {start}-{end}<br>Length: {end-start}<extra></extra>"
                )
            )
        
        fig.update_layout(
            title="<b>Non-B DNA Cluster Analysis</b>",
            xaxis_title="Genomic Position (bp)",
            yaxis_title="Clusters",
            height=400
        )
        
        return fig
    
    def _create_3d_structure_plot(self, data: pd.DataFrame, **kwargs) -> go.Figure:
        """Create 3D structure visualization (simplified representation)"""
        # This is a simplified 3D representation
        # In a full implementation, would use actual structural coordinates
        
        fig = go.Figure()
        
        # Generate pseudo-3D coordinates based on motif properties
        for i, motif in data.iterrows():
            start = int(motif['Start'])
            end = int(motif['End'])
            score = float(motif.get('Score', 0))
            
            # Generate spiral coordinates for DNA backbone
            t = np.linspace(start, end, end-start+1)
            x = np.cos(t * 0.1) * (1 + score * 0.01)
            y = np.sin(t * 0.1) * (1 + score * 0.01)
            z = t * 0.1
            
            color = self.professional_colors.get(motif['Class'], '#1f77b4')
            
            fig.add_trace(
                go.Scatter3d(
                    x=x, y=y, z=z,
                    mode='lines+markers',
                    line=dict(width=6, color=color),
                    marker=dict(size=4),
                    name=motif['Subtype'],
                    hovertemplate=f"<b>{motif['Subtype']}</b><br>Position: {start}-{end}<br>Score: {score:.1f}<extra></extra>"
                )
            )
        
        fig.update_layout(
            title="<b>3D Structure Representation</b>",
            scene=dict(
                xaxis_title="X (Å)",
                yaxis_title="Y (Å)",
                zaxis_title="Z (Å)",
                camera=dict(eye=dict(x=1.5, y=1.5, z=1.5))
            ),
            height=600
        )
        
        return fig

def create_enhanced_dashboard(motifs_df: pd.DataFrame, sequence_length: int) -> Dict[str, go.Figure]:
    """
    Create comprehensive dashboard with all visualization types
    
    Args:
        motifs_df: DataFrame containing motif data
        sequence_length: Length of analyzed sequence
        
    Returns:
        Dictionary of plot names to figure objects
    """
    visualizer = AdvancedVisualizer()
    
    dashboard = {}
    
    try:
        dashboard['distribution'] = visualizer.create_publication_quality_plot(
            motifs_df, 'motif_distribution'
        )
    except Exception as e:
        st.warning(f"Could not create distribution plot: {e}")
    
    try:
        dashboard['genomic_map'] = visualizer.create_publication_quality_plot(
            motifs_df, 'genomic_map', sequence_length=sequence_length
        )
    except Exception as e:
        st.warning(f"Could not create genomic map: {e}")
    
    try:
        dashboard['clinical'] = visualizer.create_publication_quality_plot(
            motifs_df, 'clinical_analysis'
        )
    except Exception as e:
        st.warning(f"Could not create clinical analysis: {e}")
    
    try:
        dashboard['ml_predictions'] = visualizer.create_publication_quality_plot(
            motifs_df, 'ml_predictions'
        )
    except Exception as e:
        st.warning(f"Could not create ML predictions plot: {e}")
    
    try:
        dashboard['cluster_analysis'] = visualizer.create_publication_quality_plot(
            motifs_df, 'cluster_analysis'
        )
    except Exception as e:
        st.warning(f"Could not create cluster analysis: {e}")
    
    try:
        dashboard['3d_structure'] = visualizer.create_publication_quality_plot(
            motifs_df, '3d_structure'
        )
    except Exception as e:
        st.warning(f"Could not create 3D structure plot: {e}")
    
    return dashboard

def export_publication_figure(fig: go.Figure, filename: str, format: str = 'png', 
                             width: int = 1200, height: int = 800, dpi: int = 300) -> bytes:
    """
    Export figure in publication quality format
    
    Args:
        fig: Plotly figure object
        filename: Output filename
        format: Export format ('png', 'pdf', 'svg')
        width: Figure width in pixels
        height: Figure height in pixels
        dpi: Resolution for raster formats
        
    Returns:
        Bytes of exported figure
    """
    if format.lower() == 'png':
        img_bytes = fig.to_image(format="png", width=width, height=height, scale=dpi/100)
    elif format.lower() == 'pdf':
        img_bytes = fig.to_image(format="pdf", width=width, height=height)
    elif format.lower() == 'svg':
        img_bytes = fig.to_image(format="svg", width=width, height=height)
    else:
        raise ValueError(f"Unsupported format: {format}")
    
    return img_bytes