#!/usr/bin/env python3
"""
Standard Test Sequences for NBDFinder
=====================================

This module contains curated test sequences from model organisms and pathogenic genomes
for validating NBDFinder's detection capabilities.
"""

# Standard pathogenic genome sequences with known Non-B DNA motifs
PATHOGENIC_SEQUENCES = {
    "Friedreich_Ataxia_GAA": {
        "sequence": "GAA" * 30,  # 90 bp GAA repeat
        "description": "Friedreich Ataxia GAA Repeat Expansion",
        "pathology": "Neurodegeneration",
        "expected_motifs": ["Triplex_DNA", "Slipped_DNA", "Sticky_DNA"],
        "clinical_threshold": "200+ repeats",
        "organism": "Human"
    },
    
    "Fragile_X_CGG": {
        "sequence": "CGG" * 25,  # 75 bp CGG repeat
        "description": "Fragile X Syndrome CGG Repeat",
        "pathology": "Intellectual disability",
        "expected_motifs": ["Z-DNA", "Cruciform"],
        "clinical_threshold": "200+ repeats",
        "organism": "Human"
    },
    
    "Huntington_CAG": {
        "sequence": "CAG" * 20,  # 60 bp CAG repeat
        "description": "Huntington Disease CAG Repeat",
        "pathology": "Neurodegeneration",
        "expected_motifs": ["Slipped_DNA", "Hairpin"],
        "clinical_threshold": "36+ repeats",
        "organism": "Human"
    },
    
    "cMYC_G4_Promoter": {
        "sequence": "TGGGGAGGGTGGGGAGGGTGGGGAAGG",
        "description": "c-MYC Oncogene G4 Forming Sequence",
        "pathology": "Cancer (multiple types)",
        "expected_motifs": ["Canonical_G4", "G-Triplex"],
        "clinical_threshold": "Overexpression in cancer",
        "organism": "Human"
    },
    
    "BRCA1_G4": {
        "sequence": "GGGAGGTGGGGAGGGTGGGGAAGG",
        "description": "BRCA1 G4 Forming Sequence",
        "pathology": "Breast/Ovarian cancer susceptibility",
        "expected_motifs": ["Canonical_G4"],
        "clinical_threshold": "Loss of function mutations",
        "organism": "Human"
    },
    
    "HIV_LTR_G4": {
        "sequence": "GGGAGGGAGGGAGGGGGGCCC",
        "description": "HIV-1 LTR G4 Motif",
        "pathology": "Viral replication regulation",
        "expected_motifs": ["Canonical_G4", "G-Triplex"],
        "clinical_threshold": "Therapeutic target",
        "organism": "HIV-1"
    }
}

# Model organism sequences with well-characterized Non-B DNA structures
MODEL_ORGANISM_SEQUENCES = {
    "Human_Telomeric": {
        "sequence": "TTAGGGTTAGGGTTAGGGTTAGGG",
        "description": "Human Telomeric Repeat",
        "function": "Chromosome protection",
        "expected_motifs": ["Canonical_G4", "G-Triplex"],
        "conservation": "High across mammals",
        "organism": "Homo sapiens"
    },
    
    "Yeast_Telomeric": {
        "sequence": "TGTGGGTGTGGTGTGGGTGTGG",
        "description": "S. cerevisiae Telomeric Sequence",
        "function": "Chromosome stability", 
        "expected_motifs": ["G4", "Bulged_G4"],
        "conservation": "Species-specific",
        "organism": "Saccharomyces cerevisiae"
    },
    
    "Plant_G4": {
        "sequence": "GGGTTGGGTTGGGTTGGGTT",
        "description": "Arabidopsis G4 Motif",
        "function": "Gene regulation",
        "expected_motifs": ["Canonical_G4"],
        "conservation": "Moderate across plants",
        "organism": "Arabidopsis thaliana"
    },
    
    "Drosophila_G4": {
        "sequence": "GGGAGGGAGGGAGGGA",
        "description": "Drosophila Regulatory G4",
        "function": "Developmental gene control",
        "expected_motifs": ["Canonical_G4"],
        "conservation": "High in insects",
        "organism": "Drosophila melanogaster"
    },
    
    "E_coli_rRNA_G4": {
        "sequence": "GGGCGGGGCGGGGCGGG",
        "description": "E. coli Ribosomal G4",
        "function": "rRNA processing",
        "expected_motifs": ["Canonical_G4"],
        "conservation": "High in bacteria",
        "organism": "Escherichia coli"
    },
    
    "Mouse_Immunoglobulin": {
        "sequence": "GGGGTGGGGTGGGGTGGGGAAACCCAAACCCAAACCCAAA",
        "description": "Mouse Immunoglobulin Switch Region",
        "function": "Class switch recombination",
        "expected_motifs": ["G4", "R-Loop"],
        "conservation": "High in mammals",
        "organism": "Mus musculus"
    },
    
    "Zebrafish_Hox": {
        "sequence": "GGGCGGGCGGGCGGGCGGGCGGG",
        "description": "Zebrafish Hox Gene G4",
        "function": "Developmental regulation",
        "expected_motifs": ["Canonical_G4", "Bipartite_G4"],
        "conservation": "High in vertebrates",
        "organism": "Danio rerio"
    }
}

# High-confidence G4 sequences with different formation potentials
G4_REFERENCE_SEQUENCES = {
    "High_Formation_G4": {
        "sequence": "GGGGGGGGGGGGGGGGGGGGG",  # Pure G-tract
        "g4hunter_expected": ">1.5",
        "description": "High G4 Formation Potential",
        "formation_probability": "85-95%",
        "stability": "High"
    },
    
    "Moderate_Formation_G4": {
        "sequence": "GGGGGGGGGGGGGGGCCCCCCCCCCCCCCC",
        "g4hunter_expected": "1.0-1.5", 
        "description": "Moderate G4 Formation Potential",
        "formation_probability": "60-85%",
        "stability": "Moderate"
    },
    
    "Low_Formation_G4": {
        "sequence": "GGGCCCGGGCCCGGGCCCGGG",
        "g4hunter_expected": "<1.0",
        "description": "Low G4 Formation Potential", 
        "formation_probability": "10-60%",
        "stability": "Low"
    },
    
    "Telomeric_G4": {
        "sequence": "GGGGTGGGGTGGGGTGGGG",
        "g4hunter_expected": "0.6-1.2",
        "description": "Telomeric G4 (Literature Validated)",
        "formation_probability": "Experimentally confirmed",
        "stability": "High (with K+ ions)"
    },
    
    "Promoter_G4": {
        "sequence": "GGGCGGGGCGGGGCGGGG",
        "g4hunter_expected": "0.5-1.0",
        "description": "Promoter G4 (CpG Rich)",
        "formation_probability": "Context-dependent",
        "stability": "Moderate to High"
    }
}

def get_test_dataset(dataset_type="all"):
    """
    Get test sequences for validation
    
    Parameters:
    dataset_type (str): Type of dataset to return
        - "pathogenic": Only pathogenic sequences
        - "model_organisms": Only model organism sequences  
        - "g4_reference": Only G4 reference sequences
        - "all": All sequences combined
    
    Returns:
    dict: Selected test sequences
    """
    if dataset_type == "pathogenic":
        return PATHOGENIC_SEQUENCES
    elif dataset_type == "model_organisms":
        return MODEL_ORGANISM_SEQUENCES
    elif dataset_type == "g4_reference":
        return G4_REFERENCE_SEQUENCES
    elif dataset_type == "all":
        all_sequences = {}
        all_sequences.update(PATHOGENIC_SEQUENCES)
        all_sequences.update(MODEL_ORGANISM_SEQUENCES)
        all_sequences.update(G4_REFERENCE_SEQUENCES)
        return all_sequences
    else:
        raise ValueError(f"Unknown dataset_type: {dataset_type}")

def create_validation_report():
    """Create a comprehensive validation report for all test sequences"""
    from motifs import all_motifs, g4hunter_score, get_g4_formation_category, calculate_conservation_score
    
    print("="*80)
    print("NBDFinder Validation Report - Standard Test Sequences")
    print("="*80)
    
    all_sequences = get_test_dataset("all")
    
    for seq_name, seq_data in all_sequences.items():
        print(f"\n{seq_name.upper()}")
        print("-" * len(seq_name))
        
        sequence = seq_data["sequence"]
        print(f"Sequence: {sequence}")
        print(f"Length: {len(sequence)} bp")
        print(f"Description: {seq_data['description']}")
        
        # Run analysis
        g4h_score = g4hunter_score(sequence)
        conservation = calculate_conservation_score(sequence, "G4")
        g4_category = get_g4_formation_category(g4h_score)
        
        print(f"G4Hunter Score: {g4h_score:.3f}")
        print(f"Formation Category: {g4_category['category']}")
        print(f"Conservation Score: {conservation:.3f}")
        
        # Run full motif detection
        motifs = all_motifs(sequence)
        print(f"Total Motifs Detected: {len(motifs)}")
        
        # Show motif breakdown
        motif_counts = {}
        for motif in motifs:
            motif_class = motif.get('Class', 'Unknown')
            motif_counts[motif_class] = motif_counts.get(motif_class, 0) + 1
        
        for motif_class, count in sorted(motif_counts.items()):
            print(f"  {motif_class}: {count}")
        
        # Compare with expected
        if 'expected_motifs' in seq_data:
            expected = set(seq_data['expected_motifs'])
            detected = set(motif_counts.keys())
            overlap = expected.intersection(detected)
            print(f"Expected Motifs: {expected}")
            print(f"Detection Overlap: {overlap}")
            if overlap:
                print(f"✓ Successfully detected {len(overlap)}/{len(expected)} expected motifs")
            else:
                print("⚠ No expected motifs detected")
    
    print("\n" + "="*80)
    print("Validation Report Complete")
    print("="*80)

if __name__ == "__main__":
    create_validation_report()