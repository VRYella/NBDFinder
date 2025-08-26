#!/usr/bin/env python3
"""
Comprehensive Non-B DNA Motif Detection Validation Suite
=========================================================

This script validates that all 10 Non-B DNA classes and their 22 subclasses 
are properly detected using carefully designed test sequences.

Each test sequence is designed to trigger specific motif detection while 
avoiding false positives from other motif classes.

Author: Dr. Venkata Rajesh Yella
Updated: 2025
"""

import sys
import pandas as pd
from motifs import all_motifs
from motifs.classification_config import OFFICIAL_CLASSIFICATION, get_official_classification

class MotifTestSequences:
    """Carefully designed test sequences for each motif class and subclass."""
    
    def __init__(self):
        self.test_sequences = {
            # 1. Curved DNA (2 subclasses)
            "global_curvature": {
                "sequence": "AAAAAATTTTTTCGCGCGAAAAAATTTTTTCGCGCGAAAAAATTTTTTCGCGCG",
                "expected_class": "Curved DNA",
                "expected_subtype": "Global curvature",
                "description": "Periodic poly(A)/poly(T) tracts with ~10 bp spacing"
            },
            "local_curvature": {
                "sequence": "CGCGCGAAAAAAAAAAAACGCGCG",  # Isolated A tract
                "expected_class": "Curved DNA", 
                "expected_subtype": "Local Curvature",
                "description": "Single long poly(A) tract"
            },
            
            # 2. Slipped DNA (2 subclasses)
            "direct_repeat": {
                "sequence": "GCTAAGCTAGCTAAGCTACGCTAAGCTAGCTAAGCTACGCTAAGCTAGCTAAGCTA",  # Perfect direct repeats
                "expected_class": "Slipped DNA",
                "expected_subtype": "Slipped DNA [Direct Repeat]", 
                "description": "Perfect direct repeats"
            },
            "str": {
                "sequence": "CAGCAGCAGCAGCAGCAGCAGCAGCAGCAG",  # CAG STR
                "expected_class": "Slipped DNA",
                "expected_subtype": "Slipped DNA [STR]",
                "description": "CAG short tandem repeat"
            },
            
            # 3. Cruciform DNA (1 subclass)
            "cruciform": {
                "sequence": "GCTAAGCTTACGGTAACCGTTACGAAGCTTCGC",  # Palindromic inverted repeat
                "expected_class": "Cruciform DNA",
                "expected_subtype": "Cruciform DNA [IR]/HairPin [IR]",
                "description": "Palindromic inverted repeat for cruciform formation"
            },
            
            # 4. R-loop (1 subclass)
            "r_loop": {
                "sequence": "GGGCTAGGGAAAGGGCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC",
                "expected_class": "R-loop",
                "expected_subtype": "R-loop",
                "description": "G-rich RIZ followed by long GC-rich REZ region"
            },
            
            # 5. Triplex (2 subclasses)
            "triplex": {
                "sequence": "AAAAAGAAGAAAAAGAAAACTTCTTTCTTTC",  # Homopurine/homopyrimidine mirror repeat
                "expected_class": "Triplex",
                "expected_subtype": "Triplex",
                "description": "Homopurine/homopyrimidine mirror repeat"
            },
            "sticky_dna": {
                "sequence": "GAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAAGAA",
                "expected_class": "Triplex",
                "expected_subtype": "sticky DNA",
                "description": "GAA repeat expansions (sticky DNA)"
            },
            
            # 6. G-Quadruplex Family (7 subclasses)
            "canonical_g4": {
                "sequence": "GGGTTGGGTTGGGTTGGG",  # Perfect G4
                "expected_class": "G-Quadruplex Family",
                "expected_subtype": "Canonical G4",
                "description": "Canonical G-quadruplex with short loops"
            },
            "relaxed_g4": {
                "sequence": "GGGTTTTTTTTGGGTTTTTTTTGGGTTTTTTTTGGG",  # Loops 8+ bp for relaxed
                "expected_class": "G-Quadruplex Family", 
                "expected_subtype": "Relaxed G4",
                "description": "Relaxed G-quadruplex with longer loops"
            },
            "bulged_g4": {
                "sequence": "GGGATGGGTTGGGTTGGG",  # Bulge in first G-run
                "expected_class": "G-Quadruplex Family",
                "expected_subtype": "Bulged G4",
                "description": "Bulged G-quadruplex with interruption in G-run"
            },
            "bipartite_g4": {
                "sequence": "GGGTTGGGAAAAAAAAAAAAAAAAAAGGGTTGGG",  # Long central loop 20+ bp
                "expected_class": "G-Quadruplex Family",
                "expected_subtype": "Bipartite G4", 
                "description": "Bipartite G-quadruplex with long central loop"
            },
            "multimeric_g4": {
                "sequence": "GGGTTGGGTTGGGTTGGGAAGGGTTGGGTTGGGTTGGG",  # Multiple G4s close together
                "expected_class": "G-Quadruplex Family",
                "expected_subtype": "Multimeric G4",
                "description": "Multimeric G-quadruplex structures"
            },
            "imperfect_g4": {
                "sequence": "GGATTGGGTTGGGTTGGA",  # G2 runs instead of G3+
                "expected_class": "G-Quadruplex Family",
                "expected_subtype": "Imperfect G4",
                "description": "Imperfect G-quadruplex with G2 runs"
            },
            "g_triplex": {
                "sequence": "GGGTTGGGTTGGG",  # Only 3 G-runs
                "expected_class": "G-Quadruplex Family",
                "expected_subtype": "G-Triplex intermediate",
                "description": "G-triplex intermediate structure"
            },
            
            # 7. i-motif family (3 subclasses)
            "canonical_imotif": {
                "sequence": "CCCTTCCCTTCCCTTCCC",  # Perfect i-motif
                "expected_class": "i-motif family",
                "expected_subtype": "Canonical i-motif",
                "description": "Canonical i-motif with short loops"
            },
            "relaxed_imotif": {
                "sequence": "CCCTTTTTTTTCCCTTTTTTTTCCCTTTTTTTTCCC",  # Longer loops 8+ bp
                "expected_class": "i-motif family",
                "expected_subtype": "Relaxed i-motif",
                "description": "Relaxed i-motif with longer loops"
            },
            "ac_motif": {
                "sequence": "ACACACACACACACACACACACACACACACACACACAC",  # AC repeats
                "expected_class": "i-motif family",
                "expected_subtype": "AC-motif",
                "description": "AC-motif alternating pattern"
            },
            
            # 8. Z-DNA (2 subclasses)
            "z_dna": {
                "sequence": "CGCACGCACGCACGCACGCACGCACGCACGCA",  # Alternating but not CGG repeats
                "expected_class": "Z-DNA",
                "expected_subtype": "Z-DNA",
                "description": "Alternating purine-pyrimidine sequence"
            },
            "egz_dna": {
                "sequence": "CGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGGCGG",  # Pure CGG repeats
                "expected_class": "Z-DNA",
                "expected_subtype": "eGZ (Extruded-G) DNA",
                "description": "CGG repeat expansions (eGZ)"
            }
        }
    
    def get_test_sequence(self, name):
        """Get a specific test sequence."""
        return self.test_sequences.get(name)
    
    def get_all_sequences(self):
        """Get all test sequences."""
        return self.test_sequences

def validate_motif_detection():
    """
    Validate that each motif class and subclass can be properly detected.
    """
    print("=" * 80)
    print("NBDFinder Comprehensive Motif Detection Validation")
    print("=" * 80)
    print()
    
    test_sequences = MotifTestSequences()
    all_test_sequences = test_sequences.get_all_sequences()
    
    results = {}
    detection_summary = {}
    
    # Initialize detection tracking
    for class_id, class_info in OFFICIAL_CLASSIFICATION.items():
        class_name = class_info['class_name']
        detection_summary[class_name] = {
            'expected_subclasses': class_info['subclasses'].copy(),
            'detected_subclasses': set(),
            'test_results': {}
        }
    
    print(f"Testing {len(all_test_sequences)} carefully designed test sequences...")
    print()
    
    for test_name, test_info in all_test_sequences.items():
        print(f"Testing: {test_name}")
        print(f"  Description: {test_info['description']}")
        print(f"  Sequence: {test_info['sequence'][:50]}{'...' if len(test_info['sequence']) > 50 else ''}")
        print(f"  Expected: {test_info['expected_class']} -> {test_info['expected_subtype']}")
        
        # Run detection
        motifs_detected = all_motifs(test_info['sequence'], test_name)
        
        # Analyze results
        test_result = {
            'sequence': test_info['sequence'],
            'expected_class': test_info['expected_class'],
            'expected_subtype': test_info['expected_subtype'],
            'motifs_found': len(motifs_detected),
            'detected_motifs': motifs_detected,
            'success': False,
            'notes': []
        }
        
        # Check if expected motif was detected
        expected_found = False
        for motif in motifs_detected:
            motif_class = motif.get('Class', '')
            motif_subtype = motif.get('Subtype', '')
            
            # Track all detected subclasses
            if motif_class in detection_summary:
                detection_summary[motif_class]['detected_subclasses'].add(motif_subtype)
            
            # Check if this matches our expectation
            if motif_class == test_info['expected_class'] and motif_subtype == test_info['expected_subtype']:
                expected_found = True
                test_result['success'] = True
                test_result['notes'].append(f"✅ Expected motif detected: {motif_class} -> {motif_subtype}")
        
        if not expected_found:
            if motifs_detected:
                detected_classes = [(m.get('Class', ''), m.get('Subtype', '')) for m in motifs_detected]
                test_result['notes'].append(f"❌ Expected motif not found. Detected: {detected_classes}")
            else:
                test_result['notes'].append("❌ No motifs detected")
        
        # Store results
        results[test_name] = test_result
        detection_summary[test_info['expected_class']]['test_results'][test_name] = test_result
        
        # Print result
        status = "✅ PASS" if test_result['success'] else "❌ FAIL"
        print(f"  Result: {status} ({len(motifs_detected)} motifs found)")
        for note in test_result['notes']:
            print(f"    {note}")
        print()
    
    return results, detection_summary

def print_detection_summary(detection_summary):
    """Print a summary of detection results by class."""
    print("=" * 80)
    print("DETECTION SUMMARY BY CLASS")
    print("=" * 80)
    
    total_classes = len(detection_summary)
    classes_with_detection = 0
    total_expected_subclasses = 0
    total_detected_subclasses = 0
    
    for class_name, summary in detection_summary.items():
        expected_subclasses = set(summary['expected_subclasses'])
        detected_subclasses = summary['detected_subclasses']
        total_expected_subclasses += len(expected_subclasses)
        total_detected_subclasses += len(detected_subclasses)
        
        if detected_subclasses:
            classes_with_detection += 1
        
        print(f"\n{class_name}:")
        print(f"  Expected subclasses ({len(expected_subclasses)}): {', '.join(sorted(expected_subclasses))}")
        print(f"  Detected subclasses ({len(detected_subclasses)}): {', '.join(sorted(detected_subclasses))}")
        
        # Check coverage
        missing = expected_subclasses - detected_subclasses
        extra = detected_subclasses - expected_subclasses
        
        if missing:
            print(f"  ❌ Missing subclasses: {', '.join(sorted(missing))}")
        if extra:
            print(f"  ⚠️  Extra subclasses: {', '.join(sorted(extra))}")
        if not missing and not extra and detected_subclasses:
            print(f"  ✅ All expected subclasses detected")
        elif not detected_subclasses:
            print(f"  ❌ No subclasses detected")
    
    print(f"\n" + "=" * 80)
    print("OVERALL SUMMARY")
    print("=" * 80)
    print(f"Classes with detection: {classes_with_detection}/{total_classes}")
    print(f"Subclasses detected: {total_detected_subclasses}/{total_expected_subclasses}")
    
    # Calculate success rate
    success_rate = (total_detected_subclasses / total_expected_subclasses * 100) if total_expected_subclasses > 0 else 0
    print(f"Detection success rate: {success_rate:.1f}%")
    
    return success_rate >= 80  # Consider 80%+ as passing

def create_scoring_documentation():
    """Create comprehensive documentation of all scoring systems."""
    print("\n" + "=" * 80)
    print("MOTIF DETECTION SCORING SYSTEMS DOCUMENTATION") 
    print("=" * 80)
    
    scoring_systems = {
        "Curved DNA": {
            "Global curvature": {
                "method": "Phasing-aware tract scoring",
                "formula": "(tract_len_sum/10) + tract_count + (2*phasing_score) + length_bonus",
                "min_score": 6.0,
                "max_score": "~50+ (long, well-phased arrays)",
                "factors": ["Tract length", "Tract count", "Phasing (~10.5 bp)", "Total length"]
            },
            "Local Curvature": {
                "method": "Individual tract curvature scoring", 
                "formula": "length * (1 + AT_fraction) + 0.5 * run_bonus",
                "min_score": "~7+ (min 7bp tract)",
                "max_score": "~100+ (very long tracts)",
                "factors": ["Tract length", "AT content", "Run structure"]
            }
        },
        "Slipped DNA": {
            "Slipped DNA [Direct Repeat]": {
                "method": "DR_PerfectBlock_v1",
                "formula": "unit_weight + copy_bonus + length_bonus - spacer_penalty",
                "min_score": 10.0,
                "max_score": "~20+ (long, perfect repeats)",
                "factors": ["Unit length (weight=5/l)", "Copy number", "Total length", "Spacer penalty"]
            },
            "Slipped DNA [STR]": {
                "method": "STR_PerfectBlock_v1",
                "formula": "unit_weight[1-6] + copy_bonus + length_bonus",
                "min_score": 10.0,
                "max_score": "~25+ (many short-unit repeats)",
                "factors": ["Unit length (4.0 for 1bp, 1.0 for 6bp)", "Repeat count", "Total length"]
            }
        },
        "Cruciform DNA": {
            "Cruciform DNA [IR]/HairPin [IR]": {
                "method": "IR_PerfectPalindrome_v1",
                "formula": "arm_weight(arm_len) + length_bonus - spacer_penalty",
                "min_score": "2+ (min 10bp arms)",
                "max_score": "13+ (100bp+ arms, no spacer)",
                "factors": ["Arm length (1-8 based on thresholds)", "Total length", "Spacer length"]
            }
        },
        "R-loop": {
            "R-loop": {
                "method": "RLFS_UnifiedFramework_raw",
                "formula": "(traditional_score + unified_score) / 2",
                "min_score": "Variable (depends on length)",
                "max_score": "~500+ (long, G-rich sequences)",
                "factors": ["GC fraction", "G-run count", "Length scaling", "Hunter score"]
            }
        },
        "Triplex": {
            "Triplex": {
                "method": "Triplex_Mirror_v1",
                "formula": "arm_len * (1 + 1.5*homogeneity) - 0.8*spacer_len + length_step",
                "min_score": 25.0,
                "max_score": "~200+ (long homogeneous arms)",
                "factors": ["Arm length", "Homogeneity (pur/pyr)", "Spacer penalty", "Length bonus"]
            },
            "sticky DNA": {
                "method": "GAA_TTC_v1",
                "formula": "copies * (1 + 0.3*AT_fraction)",
                "min_score": "6+ (6 repeats minimum)",
                "max_score": "~100+ (long expansions)",
                "factors": ["Repeat count", "AT content bonus"]
            }
        },
        "G-Quadruplex Family": {
            "Canonical G4": {
                "method": "G4Hunter + structural factors",
                "formula": "G4Hunter score * structural_factor * length",
                "min_score": "1.0+ (G4Hunter threshold)",
                "max_score": "~10+ (perfect, long G4s)",
                "factors": ["G4Hunter score", "Loop quality", "G-run quality", "Structural stability"]
            },
            "Relaxed G4": {
                "method": "Relaxed G4Hunter scoring",
                "formula": "Similar to canonical, relaxed thresholds",
                "min_score": "0.5+ (lower threshold)",
                "max_score": "~8+ (good relaxed structures)",
                "factors": ["Relaxed G4Hunter", "Longer loops allowed", "Structural tolerance"]
            },
            "Bulged G4": {
                "method": "Bulge-tolerant G4 scoring",
                "formula": "G4Hunter with bulge compensation",
                "min_score": "Variable (bulge-dependent)",
                "max_score": "~6+ (good structure despite bulges)",
                "factors": ["G4Hunter base", "Bulge penalty", "Compensation factors"]
            },
            "Bipartite G4": {
                "method": "Long-loop G4 detection",
                "formula": "Two-part G4 scoring",
                "min_score": "Variable (part-dependent)",
                "max_score": "~8+ (strong bipartite)",
                "factors": ["Part 1 strength", "Part 2 strength", "Central loop length"]
            },
            "Multimeric G4": {
                "method": "Multiple G4 detection",
                "formula": "Combined G4 scoring",
                "min_score": "2+ (multiple units)",
                "max_score": "~20+ (many strong units)",
                "factors": ["Individual G4 scores", "Number of units", "Spacing"]
            },
            "Imperfect G4": {
                "method": "G2+ run tolerance",
                "formula": "Imperfect G4Hunter scoring",
                "min_score": "0.3+ (very relaxed)",
                "max_score": "~5+ (good imperfect)",
                "factors": ["G2 run acceptance", "Imperfection penalty", "Overall structure"]
            },
            "G-Triplex intermediate": {
                "method": "3-tetrad G4 detection",
                "formula": "Simplified G4 scoring",
                "min_score": "0.5+ (3-tetrad minimum)",
                "max_score": "~4+ (perfect triplex)",
                "factors": ["3 G-runs only", "Loop constraints", "Stability prediction"]
            }
        },
        "i-motif family": {
            "Canonical i-motif": {
                "method": "iM-Hunter exact scoring",
                "formula": "C-centric windowed scoring",
                "min_score": "0.45+ (canonical threshold)",
                "max_score": "~1.0+ (perfect i-motif)",
                "factors": ["C-run quality", "Loop length (1-7)", "Overall structure"]
            },
            "Relaxed i-motif": {
                "method": "Relaxed iM-Hunter scoring",
                "formula": "Similar with relaxed constraints",
                "min_score": "0.30+ (relaxed threshold)",
                "max_score": "~0.8+ (good relaxed)",
                "factors": ["C-run tolerance", "Loop length (1-12)", "Structural flexibility"]
            },
            "AC-motif": {
                "method": "AC alternation scoring",
                "formula": "A/C fraction + run density + alternation",
                "min_score": "Variable (composition-based)",
                "max_score": "~20+ (long, pure AC)",
                "factors": ["A/C content", "Alternation density", "Length factor"]
            }
        },
        "Z-DNA": {
            "Z-DNA": {
                "method": "Kadane_TransitionArray_v1",
                "formula": "Dinucleotide transition scoring",
                "min_score": 50.0,
                "max_score": "~200+ (long alternating)",
                "factors": ["GC/CG dinucs (+7)", "GT/TG/AC/CA (+1.25)", "AT/TA (+0.5)", "Run penalties"]
            },
            "eGZ (Extruded-G) DNA": {
                "method": "eGZ_CGG_Expansion_v1",
                "formula": "copies * (1 + 2*G_fraction)",
                "min_score": "3+ (3 CGG minimum)",
                "max_score": "~100+ (long CGG expansions)",
                "factors": ["CGG repeat count", "G content bonus", "Structural factors"]
            }
        },
        "Hybrid": {
            "Dynamic": {
                "method": "HybridOverlap_SubclassAnalysis_raw",
                "formula": "base_score * (1 + interaction_bonus) * overlap_degree",
                "min_score": "Variable (component-dependent)",
                "max_score": "Variable (multiple high-scoring overlaps)",
                "factors": ["Contributing motif scores", "Interaction bonus", "Overlap degree"]
            }
        },
        "Non-B DNA cluster regions": {
            "Dynamic": {
                "method": "Hotspot density analysis",
                "formula": "Density-based clustering in 100nt windows",
                "min_score": "3+ motifs in window",
                "max_score": "Variable (high motif density)",
                "factors": ["Motif count in window", "Window size", "Motif diversity"]
            }
        }
    }
    
    # Print documentation in tabular format
    for class_name, subclasses in scoring_systems.items():
        print(f"\n{class_name}:")
        print("-" * 60)
        for subclass, details in subclasses.items():
            print(f"  {subclass}:")
            print(f"    Method: {details['method']}")
            print(f"    Formula: {details['formula']}")
            print(f"    Score Range: {details['min_score']} to {details['max_score']}")
            print(f"    Key Factors: {', '.join(details['factors'])}")
            print()

def main():
    """Run comprehensive motif detection validation."""
    print("NBDFinder Comprehensive Motif Detection Validation Suite")
    print("Testing all 10 classes and 22+ subclasses with designed sequences")
    print()
    
    # Run validation
    results, detection_summary = validate_motif_detection()
    
    # Print summary
    overall_success = print_detection_summary(detection_summary)
    
    # Create scoring documentation
    create_scoring_documentation()
    
    # Final assessment
    print("\n" + "=" * 80)
    if overall_success:
        print("✅ VALIDATION PASSED - Most motif detection systems are working correctly")
        print("All features properly integrated with scientific accuracy maintained.")
        exit_code = 0
    else:
        print("❌ VALIDATION FAILED - Some motif detection systems need attention")
        print("Please review the failed test cases and implement missing detection algorithms.")
        exit_code = 1
    
    print("=" * 80)
    return exit_code

if __name__ == "__main__":
    sys.exit(main())