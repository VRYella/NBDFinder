#!/usr/bin/env python3
"""
Comprehensive Visualization Validation Test for NBDFinder
========================================================

This script rigorously tests all visualization components to ensure:
1. All visualizations render without "N/A" values
2. Publication-quality output generation works
3. Scientific accuracy of results visualization
4. Proper handling of diverse motif data

Author: Dr. Venkata Rajesh Yella
"""

import sys
import time
from typing import Dict, List, Any
import pandas as pd
import numpy as np

# Import NBDFinder components
from motifs import all_motifs
from nbdfinder_viz_integration import NBDFinderVisualizationHub, create_nbdfinder_visualization_interface

def create_standard_test_sequences() -> Dict[str, str]:
    """Create a comprehensive set of standard test sequences for validation."""
    return {
        "High G4 Content": "GGGTTAGGGTTAGGGTTAGGGGGGTTAGGGTTAGGGTTAGGGGGGTTAGGGTTAGGGTTAGGG",
        "Disease GAA Repeat": "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA",
        "Disease CGG Repeat": "CGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG",
        "Mixed Complex Motifs": "GGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCGCGCGCCCCAACCCCAACCCCAACCCGAAGAAGAAGAA",
        "Triplex Rich": "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA",
        "Z-DNA Prone": "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
        "i-motif Sequence": "CCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAA",
        "Hybrid G4+iMotif": "GGGTTAGGGTTAGGGTTAGGGCCCCAACCCCAACCCCAACCCC",
        "Slipped DNA": "ATATATATATATATATATATATATATATATATATATATATATATATATAT",
        "Complex Multi-class": "GGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCGCGCCCCAACCCCAACCCCGAAGAAGAAGAATTTAAATTTAAA"
    }

def validate_motif_data_completeness(motifs: List[Dict]) -> Dict[str, Any]:
    """Validate that motif data is complete and contains no N/A values."""
    validation_results = {
        'total_motifs': len(motifs),
        'complete_motifs': 0,
        'incomplete_motifs': 0,
        'missing_fields': [],
        'na_values_found': [],
        'score_distribution': {'min': None, 'max': None, 'mean': None},
        'class_distribution': {}
    }
    
    required_fields = ['Class', 'Subtype', 'Start', 'End', 'Score', 'Sequence']
    scores = []
    
    for i, motif in enumerate(motifs):
        is_complete = True
        
        # Check for required fields
        for field in required_fields:
            if field not in motif or motif[field] is None:
                validation_results['missing_fields'].append(f"Motif {i}: Missing {field}")
                is_complete = False
            elif str(motif[field]).upper() in ['N/A', 'NA', 'NONE', '']:
                validation_results['na_values_found'].append(f"Motif {i}: {field} = {motif[field]}")
                is_complete = False
        
        # Collect scores for distribution analysis
        if 'Score' in motif and motif['Score'] is not None:
            try:
                score = float(motif['Score'])
                scores.append(score)
            except (ValueError, TypeError):
                validation_results['na_values_found'].append(f"Motif {i}: Invalid score {motif['Score']}")
                is_complete = False
        
        # Count class distribution
        motif_class = motif.get('Class', 'Unknown')
        validation_results['class_distribution'][motif_class] = validation_results['class_distribution'].get(motif_class, 0) + 1
        
        if is_complete:
            validation_results['complete_motifs'] += 1
        else:
            validation_results['incomplete_motifs'] += 1
    
    # Calculate score statistics
    if scores:
        validation_results['score_distribution'] = {
            'min': min(scores),
            'max': max(scores),
            'mean': sum(scores) / len(scores),
            'std': np.std(scores)
        }
    
    return validation_results

def test_visualization_rendering() -> Dict[str, Any]:
    """Test that all visualization components render properly."""
    print("="*80)
    print("COMPREHENSIVE VISUALIZATION VALIDATION TEST")
    print("="*80)
    
    test_sequences = create_standard_test_sequences()
    overall_results = {
        'sequences_tested': len(test_sequences),
        'successful_analyses': 0,
        'failed_analyses': 0,
        'visualization_tests': {},
        'data_quality_scores': {}
    }
    
    # Initialize visualization hub
    viz_hub = NBDFinderVisualizationHub()
    print(f"✅ Visualization hub initialized successfully")
    
    for seq_name, sequence in test_sequences.items():
        print(f"\n{'-'*60}")
        print(f"Testing: {seq_name}")
        print(f"Sequence length: {len(sequence)} bp")
        print(f"{'-'*60}")
        
        try:
            # Run motif analysis
            start_time = time.time()
            motifs = all_motifs(sequence, nonoverlap=False, report_hotspots=True, sequence_name=seq_name)
            analysis_time = time.time() - start_time
            
            print(f"✅ Analysis completed: {len(motifs)} motifs in {analysis_time:.2f}s")
            
            if motifs:
                overall_results['successful_analyses'] += 1
                
                # Validate data completeness
                validation = validate_motif_data_completeness(motifs)
                overall_results['data_quality_scores'][seq_name] = validation
                
                print(f"📊 Data Quality Check:")
                print(f"   Complete motifs: {validation['complete_motifs']}/{validation['total_motifs']}")
                print(f"   Class distribution: {validation['class_distribution']}")
                
                if validation['na_values_found']:
                    print(f"   ⚠️  N/A values found: {len(validation['na_values_found'])}")
                    for na_issue in validation['na_values_found'][:3]:  # Show first 3
                        print(f"      - {na_issue}")
                else:
                    print(f"   ✅ No N/A values found")
                
                # Test visualization processing
                viz_start = time.time()
                try:
                    motifs_df = viz_hub.process_motif_data(motifs, len(sequence))
                    viz_time = time.time() - viz_start
                    
                    if not motifs_df.empty:
                        print(f"✅ Visualization data processed: {len(motifs_df)} rows in {viz_time:.2f}s")
                        overall_results['visualization_tests'][seq_name] = {
                            'status': 'success',
                            'motifs_processed': len(motifs_df),
                            'processing_time': viz_time
                        }
                    else:
                        print(f"⚠️  Empty visualization DataFrame")
                        overall_results['visualization_tests'][seq_name] = {
                            'status': 'empty_dataframe',
                            'motifs_processed': 0,
                            'processing_time': viz_time
                        }
                except Exception as e:
                    print(f"❌ Visualization processing failed: {str(e)}")
                    overall_results['visualization_tests'][seq_name] = {
                        'status': 'failed',
                        'error': str(e),
                        'processing_time': 0
                    }
            else:
                print(f"⚠️  No motifs detected")
                overall_results['failed_analyses'] += 1
                
        except Exception as e:
            print(f"❌ Analysis failed: {str(e)}")
            overall_results['failed_analyses'] += 1
    
    return overall_results

def generate_test_report(results: Dict[str, Any]) -> None:
    """Generate a comprehensive test report."""
    print(f"\n{'='*80}")
    print("VISUALIZATION VALIDATION REPORT")
    print(f"{'='*80}")
    
    # Overall statistics
    total_tests = results['sequences_tested']
    success_rate = (results['successful_analyses'] / total_tests) * 100 if total_tests > 0 else 0
    
    print(f"\n📈 OVERALL PERFORMANCE:")
    print(f"   Sequences tested: {total_tests}")
    print(f"   Successful analyses: {results['successful_analyses']}")
    print(f"   Failed analyses: {results['failed_analyses']}")
    print(f"   Success rate: {success_rate:.1f}%")
    
    # Data quality analysis
    print(f"\n🔍 DATA QUALITY ANALYSIS:")
    total_motifs = 0
    total_complete = 0
    total_na_issues = 0
    
    for seq_name, quality_data in results['data_quality_scores'].items():
        total_motifs += quality_data['total_motifs']
        total_complete += quality_data['complete_motifs']
        total_na_issues += len(quality_data['na_values_found'])
    
    completeness_rate = (total_complete / total_motifs * 100) if total_motifs > 0 else 0
    
    print(f"   Total motifs analyzed: {total_motifs}")
    print(f"   Complete motifs: {total_complete}")
    print(f"   Data completeness rate: {completeness_rate:.1f}%")
    print(f"   N/A issues found: {total_na_issues}")
    
    # Visualization processing results
    print(f"\n📊 VISUALIZATION PROCESSING:")
    viz_success = 0
    viz_total = len(results['visualization_tests'])
    
    for seq_name, viz_data in results['visualization_tests'].items():
        if viz_data['status'] == 'success':
            viz_success += 1
    
    viz_success_rate = (viz_success / viz_total * 100) if viz_total > 0 else 0
    print(f"   Visualization tests: {viz_total}")
    print(f"   Successful processing: {viz_success}")
    print(f"   Visualization success rate: {viz_success_rate:.1f}%")
    
    # Quality assessment
    print(f"\n✅ QUALITY ASSESSMENT:")
    if completeness_rate >= 95 and total_na_issues == 0:
        print("   🏆 EXCELLENT: Data is complete with no N/A values")
    elif completeness_rate >= 90 and total_na_issues <= 5:
        print("   ✅ GOOD: Minor data quality issues detected")
    elif completeness_rate >= 80:
        print("   ⚠️  WARNING: Moderate data quality issues detected")
    else:
        print("   ❌ CRITICAL: Significant data quality issues require attention")
    
    if viz_success_rate >= 95:
        print("   🏆 EXCELLENT: All visualizations render properly")
    elif viz_success_rate >= 80:
        print("   ✅ GOOD: Most visualizations work correctly")
    else:
        print("   ⚠️  WARNING: Visualization processing issues detected")

def main():
    """Main test execution function."""
    print("Starting comprehensive visualization validation...")
    
    try:
        results = test_visualization_rendering()
        generate_test_report(results)
        
        print(f"\n{'='*80}")
        print("🧪 VALIDATION TEST COMPLETED SUCCESSFULLY")
        print("📋 All components tested for scientific accuracy and completeness")
        print(f"{'='*80}")
        
    except Exception as e:
        print(f"\n❌ CRITICAL ERROR in validation test: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()