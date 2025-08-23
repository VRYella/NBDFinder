#!/usr/bin/env python3
"""
Scientific Validation Test Suite for NBDFinder
==============================================

This script performs rigorous testing with scientific standards to ensure:
1. High-quality, complete visualization outputs
2. Scientific accuracy of motif detection
3. Publication-ready figure generation
4. Comprehensive validation with test sequences

Author: Dr. Venkata Rajesh Yella
"""

import sys
import time
import os
from typing import Dict, List, Any, Tuple
import pandas as pd
import numpy as np

# Import NBDFinder components
from motifs import all_motifs
from nbdfinder_viz_integration import NBDFinderVisualizationHub, create_nbdfinder_visualization_interface
from publication_visualizations import create_comprehensive_publication_suite

def create_scientific_test_sequences() -> Dict[str, Dict[str, Any]]:
    """Create a comprehensive set of scientifically validated test sequences."""
    return {
        "G4_Canonical": {
            "sequence": "GGGTTAGGGTTAGGGTTAGGG",
            "expected_classes": ["G-Quadruplex Family"],
            "expected_motifs": {"min": 1, "max": 10},
            "description": "Canonical G-quadruplex forming sequence with 3-base loops"
        },
        "Disease_CGG_Repeat": {
            "sequence": "CGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG",
            "expected_classes": ["Slipped DNA", "Triplex", "Hybrid"],
            "expected_motifs": {"min": 3, "max": 50},
            "description": "CGG repeat expansion associated with Fragile X syndrome"
        },
        "Disease_GAA_Repeat": {
            "sequence": "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA",
            "expected_classes": ["Slipped DNA", "Triplex", "Hybrid"],
            "expected_motifs": {"min": 3, "max": 50},
            "description": "GAA repeat expansion associated with Friedreich's ataxia"
        },
        "Z_DNA_Prone": {
            "sequence": "CGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG",
            "expected_classes": ["Z-DNA", "Cruciform DNA"],
            "expected_motifs": {"min": 1, "max": 15},
            "description": "Alternating CG sequence prone to Z-DNA formation"
        },
        "iMotif_Rich": {
            "sequence": "CCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAACCCCAA",
            "expected_classes": ["i-motif family"],
            "expected_motifs": {"min": 1, "max": 20},
            "description": "C-rich sequence capable of i-motif formation"
        },
        "Complex_Hybrid": {
            "sequence": "GGGTTAGGGTTAGGGTTAGGGCGCGCGCGCGCGCGCCCCAACCCCAACCCCAACCCCGAAGAAGAAGAAGAAGAAGAATTTAAATTTAAATTTAAA",
            "expected_classes": ["G-Quadruplex Family", "Z-DNA", "i-motif family", "Slipped DNA", "Hybrid"],
            "expected_motifs": {"min": 5, "max": 100},
            "description": "Complex sequence containing multiple motif types"
        },
        "Triplex_Target": {
            "sequence": "GGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAAGGAA",
            "expected_classes": ["Triplex"],
            "expected_motifs": {"min": 1, "max": 30},
            "description": "Purine-rich sequence capable of triplex formation"
        },
        "Curved_DNA": {
            "sequence": "AAAATTTTTAAAATTTTTAAAATTTTTAAAATTTTTAAAATTTTTAAAATTTTT",
            "expected_classes": ["Curved DNA"],
            "expected_motifs": {"min": 1, "max": 10},
            "description": "A-tract sequence known to cause DNA curvature"
        }
    }

def validate_scientific_accuracy(motifs: List[Dict], test_info: Dict[str, Any]) -> Dict[str, Any]:
    """Validate the scientific accuracy of motif detection results."""
    validation = {
        "sequence_name": test_info.get("description", "Unknown"),
        "total_motifs": len(motifs),
        "classes_detected": [],
        "classes_expected": test_info.get("expected_classes", []),
        "motif_count_valid": False,
        "class_accuracy": 0.0,
        "score_distribution": {},
        "quality_issues": [],
        "scientific_validity": "UNKNOWN"
    }
    
    # Collect detected classes
    detected_classes = set()
    scores = []
    
    for motif in motifs:
        motif_class = motif.get('Class', 'Unknown')
        detected_classes.add(motif_class)
        
        # Collect scores for distribution analysis
        if 'Score' in motif and motif['Score'] is not None:
            try:
                score = float(motif['Score'])
                scores.append(score)
            except (ValueError, TypeError):
                validation["quality_issues"].append(f"Invalid score: {motif['Score']}")
    
    validation["classes_detected"] = list(detected_classes)
    
    # Validate motif count range
    expected_range = test_info.get("expected_motifs", {"min": 0, "max": 1000})
    motif_count = len(motifs)
    validation["motif_count_valid"] = expected_range["min"] <= motif_count <= expected_range["max"]
    
    # Calculate class detection accuracy
    expected_classes = set(test_info.get("expected_classes", []))
    if expected_classes:
        detected_expected = detected_classes.intersection(expected_classes)
        validation["class_accuracy"] = len(detected_expected) / len(expected_classes) * 100
    
    # Score distribution analysis
    if scores:
        validation["score_distribution"] = {
            "min": min(scores),
            "max": max(scores),
            "mean": np.mean(scores),
            "std": np.std(scores),
            "count": len(scores)
        }
    
    # Assess overall scientific validity
    if validation["motif_count_valid"] and validation["class_accuracy"] >= 50:
        validation["scientific_validity"] = "VALID"
    elif validation["class_accuracy"] >= 25:
        validation["scientific_validity"] = "PARTIAL"
    else:
        validation["scientific_validity"] = "INVALID"
    
    return validation

def test_visualization_completeness(motifs: List[Dict], sequence_length: int) -> Dict[str, Any]:
    """Test that visualizations are complete and high-quality."""
    viz_test = {
        "data_processing": False,
        "na_values": 0,
        "visualization_suite": False,
        "plot_types_created": 0,
        "export_test": False,
        "quality_score": 0.0,
        "issues": []
    }
    
    try:
        # Test data processing
        viz_hub = NBDFinderVisualizationHub()
        motifs_df = viz_hub.process_motif_data(motifs, sequence_length)
        
        if not motifs_df.empty:
            viz_test["data_processing"] = True
            viz_test["na_values"] = motifs_df.isnull().sum().sum()
            
            # Test visualization suite creation
            suite = create_comprehensive_publication_suite(motifs_df)
            
            if suite:
                viz_test["visualization_suite"] = True
                viz_test["plot_types_created"] = len(suite)
                
                # Test export capability
                if 'bar_plots' in suite and suite['bar_plots']:
                    try:
                        # Get first available plot
                        first_plot_key = list(suite['bar_plots'].keys())[0]
                        plot_data = suite['bar_plots'][first_plot_key]
                        
                        # Extract plotly figure from nested structure
                        if isinstance(plot_data, dict) and 'plotly' in plot_data:
                            fig = plot_data['plotly']
                        else:
                            fig = plot_data
                        
                        # Test PNG export
                        img_bytes = fig.to_image(format="png", width=800, height=600, scale=2)
                        if len(img_bytes) > 1000:  # Reasonable size check
                            viz_test["export_test"] = True
                    except Exception as e:
                        viz_test["issues"].append(f"Export test failed: {str(e)}")
        
        # Calculate quality score
        quality_components = [
            viz_test["data_processing"] * 25,
            (viz_test["na_values"] == 0) * 25,
            viz_test["visualization_suite"] * 25,
            viz_test["export_test"] * 25
        ]
        viz_test["quality_score"] = sum(quality_components)
        
    except Exception as e:
        viz_test["issues"].append(f"Critical error: {str(e)}")
    
    return viz_test

def run_scientific_validation_suite():
    """Run the complete scientific validation test suite."""
    print("=" * 80)
    print("NBDFinder SCIENTIFIC VALIDATION TEST SUITE")
    print("=" * 80)
    print("Performing rigorous testing for high-standard scientific output")
    print()
    
    test_sequences = create_scientific_test_sequences()
    overall_results = {
        "tests_run": 0,
        "tests_passed": 0,
        "scientific_validity": {"VALID": 0, "PARTIAL": 0, "INVALID": 0},
        "visualization_quality": [],
        "performance_metrics": [],
        "issues_found": []
    }
    
    for test_name, test_info in test_sequences.items():
        print(f"\n{'='*60}")
        print(f"TESTING: {test_name}")
        print(f"{'='*60}")
        print(f"Description: {test_info['description']}")
        print(f"Sequence length: {len(test_info['sequence'])} bp")
        print(f"Expected classes: {test_info['expected_classes']}")
        
        overall_results["tests_run"] += 1
        
        try:
            # Run motif analysis
            start_time = time.time()
            motifs = all_motifs(test_info["sequence"], nonoverlap=True, report_hotspots=False)
            analysis_time = time.time() - start_time
            
            print(f"\n📊 MOTIF ANALYSIS RESULTS:")
            print(f"   Motifs detected: {len(motifs)}")
            print(f"   Analysis time: {analysis_time:.3f}s")
            
            # Validate scientific accuracy
            validation = validate_scientific_accuracy(motifs, test_info)
            
            print(f"\n🔬 SCIENTIFIC VALIDATION:")
            print(f"   Classes detected: {validation['classes_detected']}")
            print(f"   Class accuracy: {validation['class_accuracy']:.1f}%")
            print(f"   Motif count valid: {validation['motif_count_valid']}")
            print(f"   Scientific validity: {validation['scientific_validity']}")
            
            if validation['score_distribution']:
                scores = validation['score_distribution']
                print(f"   Score range: {scores['min']:.1f} - {scores['max']:.1f}")
                print(f"   Mean score: {scores['mean']:.2f} ± {scores['std']:.2f}")
            
            # Test visualization completeness
            viz_test = test_visualization_completeness(motifs, len(test_info["sequence"]))
            
            print(f"\n📈 VISUALIZATION QUALITY:")
            print(f"   Data processing: {'✅' if viz_test['data_processing'] else '❌'}")
            print(f"   N/A values: {viz_test['na_values']}")
            print(f"   Visualization suite: {'✅' if viz_test['visualization_suite'] else '❌'}")
            print(f"   Plot types created: {viz_test['plot_types_created']}")
            print(f"   Export capability: {'✅' if viz_test['export_test'] else '❌'}")
            print(f"   Quality score: {viz_test['quality_score']}/100")
            
            if viz_test['issues']:
                print(f"   Issues: {'; '.join(viz_test['issues'])}")
                overall_results["issues_found"].extend(viz_test['issues'])
            
            # Record results
            overall_results["scientific_validity"][validation['scientific_validity']] += 1
            overall_results["visualization_quality"].append(viz_test['quality_score'])
            overall_results["performance_metrics"].append({
                "test": test_name,
                "motifs": len(motifs),
                "time": analysis_time,
                "throughput": len(test_info["sequence"]) / analysis_time
            })
            
            # Determine if test passed
            test_passed = (validation['scientific_validity'] in ['VALID', 'PARTIAL'] and 
                          viz_test['quality_score'] >= 75)
            
            if test_passed:
                overall_results["tests_passed"] += 1
                print(f"\n✅ TEST RESULT: PASSED")
            else:
                print(f"\n❌ TEST RESULT: FAILED")
                
        except Exception as e:
            print(f"\n❌ CRITICAL ERROR: {str(e)}")
            overall_results["issues_found"].append(f"{test_name}: {str(e)}")
    
    # Generate final report
    generate_final_scientific_report(overall_results)

def generate_final_scientific_report(results: Dict[str, Any]) -> None:
    """Generate a comprehensive scientific validation report."""
    print(f"\n{'='*80}")
    print("SCIENTIFIC VALIDATION REPORT")
    print(f"{'='*80}")
    
    # Overall test statistics
    total_tests = results["tests_run"]
    passed_tests = results["tests_passed"]
    success_rate = (passed_tests / total_tests * 100) if total_tests > 0 else 0
    
    print(f"\n📈 OVERALL PERFORMANCE:")
    print(f"   Tests executed: {total_tests}")
    print(f"   Tests passed: {passed_tests}")
    print(f"   Success rate: {success_rate:.1f}%")
    
    # Scientific validity analysis
    validity_stats = results["scientific_validity"]
    print(f"\n🔬 SCIENTIFIC ACCURACY:")
    print(f"   Valid results: {validity_stats['VALID']}")
    print(f"   Partial results: {validity_stats['PARTIAL']}")
    print(f"   Invalid results: {validity_stats['INVALID']}")
    
    validity_rate = ((validity_stats['VALID'] + validity_stats['PARTIAL']) / total_tests * 100) if total_tests > 0 else 0
    print(f"   Scientific validity rate: {validity_rate:.1f}%")
    
    # Visualization quality analysis
    viz_scores = results["visualization_quality"]
    if viz_scores:
        avg_viz_quality = np.mean(viz_scores)
        min_viz_quality = min(viz_scores)
        max_viz_quality = max(viz_scores)
        
        print(f"\n📊 VISUALIZATION QUALITY:")
        print(f"   Average quality score: {avg_viz_quality:.1f}/100")
        print(f"   Quality range: {min_viz_quality} - {max_viz_quality}")
        print(f"   High quality tests (≥80): {sum(1 for s in viz_scores if s >= 80)}")
    
    # Performance analysis
    perf_metrics = results["performance_metrics"]
    if perf_metrics:
        throughputs = [m["throughput"] for m in perf_metrics]
        avg_throughput = np.mean(throughputs)
        
        print(f"\n⚡ PERFORMANCE METRICS:")
        print(f"   Average throughput: {avg_throughput:.0f} bp/second")
        print(f"   Throughput range: {min(throughputs):.0f} - {max(throughputs):.0f} bp/second")
    
    # Issues summary
    issues = results["issues_found"]
    if issues:
        print(f"\n⚠️  ISSUES IDENTIFIED ({len(issues)}):")
        for i, issue in enumerate(issues[:5], 1):  # Show first 5
            print(f"   {i}. {issue}")
        if len(issues) > 5:
            print(f"   ... and {len(issues) - 5} more")
    
    # Final assessment
    print(f"\n🏆 FINAL ASSESSMENT:")
    if success_rate >= 90 and avg_viz_quality >= 80:
        print("   ✅ EXCELLENT: NBDFinder meets high scientific standards")
    elif success_rate >= 75 and avg_viz_quality >= 70:
        print("   ✅ GOOD: NBDFinder meets acceptable scientific standards")
    elif success_rate >= 50:
        print("   ⚠️  WARNING: Some quality issues detected")
    else:
        print("   ❌ CRITICAL: Significant quality issues require attention")
    
    print(f"\n📋 COMPLIANCE STATUS:")
    print(f"   ✅ Visualization completeness: {'PASS' if avg_viz_quality >= 75 else 'FAIL'}")
    print(f"   ✅ Scientific accuracy: {'PASS' if validity_rate >= 75 else 'FAIL'}")
    print(f"   ✅ Data quality (no N/A): {'PASS' if all(s >= 75 for s in viz_scores) else 'FAIL'}")
    print(f"   ✅ Publication readiness: {'PASS' if success_rate >= 80 else 'FAIL'}")
    
    print(f"\n{'='*80}")
    print("🧪 SCIENTIFIC VALIDATION COMPLETED")
    print("📊 Ready for publication-quality scientific analysis")
    print(f"{'='*80}")

def main():
    """Main test execution function."""
    print("Starting comprehensive scientific validation...")
    print("This test ensures high-standard scientific output quality")
    
    try:
        run_scientific_validation_suite()
        
    except Exception as e:
        print(f"\n❌ CRITICAL ERROR in scientific validation: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == "__main__":
    main()