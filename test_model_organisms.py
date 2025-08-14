#!/usr/bin/env python3
"""
Test Enhanced G4Hunter with High-Scoring Sequences
"""

from motifs import find_gquadruplex, g4hunter_score, get_g4_formation_category, calculate_conservation_score

def test_high_scoring_g4_sequences():
    """Test with sequences designed to have high G4Hunter scores"""
    print("="*70)
    print("Testing Enhanced G4Hunter with High-Scoring Sequences")
    print("="*70)
    
    # Sequences designed for different G4Hunter score ranges
    test_sequences = [
        # High G-content sequences for high scores (≥1.5)
        ("GGGGGGGGGGGGGGGGGGGGG", "Pure G-tract (should be ≥1.5)"),
        ("GGGGGGGGGGGGGGGGGGCCCCCCCCCC", "G-rich with C balance"),
        ("GGGGGGGGGGGGGGGGGG", "G-rich tract"),
        
        # Moderate G4Hunter scores (1.0-1.5)
        ("GGGGGGGGGGGGGGGCCCCCCCCCCCCCCC", "Balanced G/C content"),
        ("GGGGGGGGGGGGCCCCCCCCCCC", "Moderate G/C balance"),
        
        # Low G4Hunter scores (<1.0)
        ("GGGCCCGGGCCCGGGCCCGGG", "Alternating G/C"),
        ("GGGATAGGGCTAGGGCTGGG", "G-rich with A/T dilution"),
        
        # Real G4 forming sequences with known high potential
        ("GGGGTGGGGTGGGGTGGGG", "Telomeric G4 motif"),
        ("GGGGGCGGGGGCGGGGGCGGGG", "CpG island G4"),
        ("GGGGCGGGGCGGGGCGGGG", "Compact G4")
    ]
    
    print(f"{'Sequence':<30} {'Length':<8} {'G4Hunter':<10} {'Category':<25} {'Evidence':<15} {'Conservation':<12}")
    print("-" * 110)
    
    for seq, description in test_sequences:
        # Calculate G4Hunter score
        g4h_score = g4hunter_score(seq)
        
        # Get formation category  
        category = get_g4_formation_category(g4h_score)
        
        # Calculate conservation score
        conservation = calculate_conservation_score(seq, "G4")
        
        # Truncate description for table
        short_desc = description[:28] + ".." if len(description) > 30 else description
        
        print(f"{short_desc:<30} {len(seq):<8} {g4h_score:<10.3f} {category['category'][:24]:<25} {category['experimental_evidence']:<15} {conservation:<12.3f}")
        
        # Show detailed formation info for high scorers
        if g4h_score >= 1.0:
            print(f"    Formation Probability: {category['formation_probability']}")
            print(f"    Stability: {category['stability']}")
            
        # Run motif detection
        motifs = find_gquadruplex(seq)
        if motifs:
            print(f"    Detected {len(motifs)} G4 motif(s)")
            for motif in motifs:
                print(f"      Position {motif['Start']}-{motif['End']}: Score={motif['Score']:.2f}")
        
        print()

def create_model_organism_test_dataset():
    """Create a test dataset with model organism sequences"""
    print("="*70)
    print("Model Organism G4 Motif Analysis")
    print("="*70)
    
    # Model organism sequences with known G4 motifs
    model_sequences = [
        # Human telomeric repeats
        ("TTAGGGTTAGGGTTAGGGTTAGGG", "Human Telomeric Repeat", "Homo sapiens"),
        
        # c-MYC promoter G4
        ("TGGGGAGGGTGGGGAGGGTGGGGAAGG", "c-MYC Promoter G4", "Homo sapiens"),
        
        # VEGF promoter G4
        ("GGGCGGGGGCGGGGGCGGGGGAGG", "VEGF Promoter G4", "Homo sapiens"),
        
        # BCL2 promoter G4
        ("GGGCGCGGGAGGAAGGGGGCGGG", "BCL2 Promoter G4", "Homo sapiens"),
        
        # Yeast telomeric sequence
        ("TGTGGGTGTGGTGTGGGTGTGG", "Yeast Telomeric", "Saccharomyces cerevisiae"),
        
        # Arabidopsis G4 motif
        ("GGGTTGGGTTGGGTTGGGTT", "Plant G4 Motif", "Arabidopsis thaliana"),
        
        # Drosophila G4 sequence
        ("GGGAGGGAGGGAGGGA", "Drosophila G4", "Drosophila melanogaster"),
        
        # E. coli ribosomal G4
        ("GGGCGGGGCGGGGCGGG", "E. coli Ribosomal G4", "Escherichia coli")
    ]
    
    print(f"{'Organism':<25} {'Sequence Type':<20} {'G4Hunter':<10} {'Category':<20} {'Conservation':<12}")
    print("-" * 95)
    
    for seq, seq_type, organism in model_sequences:
        g4h_score = g4hunter_score(seq)
        category = get_g4_formation_category(g4h_score)
        conservation = calculate_conservation_score(seq, "G4")
        
        print(f"{organism[:24]:<25} {seq_type[:19]:<20} {g4h_score:<10.3f} {category['category'][:19]:<20} {conservation:<12.3f}")
        
        # Detailed analysis for each sequence
        motifs = find_gquadruplex(seq)
        if motifs:
            print(f"    Detected {len(motifs)} G4 motif(s)")
            for i, motif in enumerate(motifs):
                print(f"      Motif {i+1}: Pos {motif['Start']}-{motif['End']}, "
                      f"Score={motif['Score']:.2f}, Prob={motif['Formation_Probability']}")
        else:
            print("    No G4 motifs detected above threshold")
        print()

if __name__ == "__main__":
    test_high_scoring_g4_sequences()
    print()
    create_model_organism_test_dataset()
    print("\n" + "="*70)
    print("Enhanced G4Hunter Model Organism Testing Complete!")
    print("="*70)