#!/usr/bin/env python3
"""
Test Enhanced G4Hunter Implementation with Conservation Scores and Experimental Data
"""

from motifs import find_gquadruplex, g4hunter_score, get_g4_formation_category, calculate_conservation_score

def test_enhanced_g4_features():
    """Test the enhanced G4Hunter features with conservation scores and experimental formation data"""
    print("="*70)
    print("Testing Enhanced G4Hunter Implementation")
    print("="*70)
    
    # Test sequences with different G4Hunter scores
    test_sequences = [
        # High formation potential (≥1.5)
        ("GGGTTAGGGTTAGGGTTAGGG", "High Formation Potential (>1.5)"),
        # Moderate formation potential (1.0-1.5)  
        ("GGGTTTTGGGTTTGGGTTTGGG", "Moderate Formation Potential (1.0-1.5)"),
        # Low formation potential (<1.0)
        ("GGGAAAAAGGGAAAAAGGGAAAAAGGG", "Low Formation Potential (<1.0)"),
        # Real pathogenic sequences
        ("GGGTTAGGGTTCGGTTAGGG", "CGG Repeat (Fragile X-related)"),
        ("GGGTTAGGGTTAGGGTTGGGAAAGGG", "Mixed G4 motif")
    ]
    
    for seq, description in test_sequences:
        print(f"\n{description}")
        print(f"Sequence: {seq}")
        print(f"Length: {len(seq)} bp")
        
        # Calculate G4Hunter score
        g4h_score = g4hunter_score(seq)
        print(f"G4Hunter Score: {g4h_score:.3f}")
        
        # Get formation category
        category = get_g4_formation_category(g4h_score)
        print(f"Formation Category: {category['category']}")
        print(f"Threshold: {category['threshold']}")
        print(f"Experimental Evidence: {category['experimental_evidence']}")
        print(f"Formation Probability: {category['formation_probability']}")
        
        # Calculate conservation score
        conservation = calculate_conservation_score(seq, "G4")
        print(f"Conservation Score: {conservation:.3f}")
        
        # Run full G4 detection
        motifs = find_gquadruplex(seq)
        print(f"Detected Motifs: {len(motifs)}")
        
        for motif in motifs:
            print(f"  - Start: {motif['Start']}, End: {motif['End']}")
            print(f"    G4Hunter: {motif['G4Hunter_Mean']:.3f}")
            print(f"    Conservation: {motif['Conservation_Score']:.3f}")
            print(f"    Category: {motif['Formation_Category']}")
            print(f"    Evidence: {motif['Experimental_Evidence']}")
        
        print("-" * 50)

def test_pathogenic_sequences():
    """Test on known pathogenic sequences"""
    print("\n" + "="*70)
    print("Testing on Known Pathogenic Sequences")
    print("="*70)
    
    pathogenic_sequences = [
        # Friedreich ataxia GAA repeats
        ("GAA" * 30, "Friedreich Ataxia GAA Repeats"),
        # Fragile X CGG repeats  
        ("CGG" * 25, "Fragile X CGG Repeats"),
        # Huntington CAG repeats
        ("CAG" * 20, "Huntington CAG Repeats"),
        # c-MYC G4 forming sequence
        ("TGGGGAGGGTGGGGAGGGTGGGGAAGG", "c-MYC Oncogene G4"),
        # BRCA1 G4 forming sequence
        ("GGGAGGTGGGGAGGGTGGGGAAGG", "BRCA1 G4 motif")
    ]
    
    for seq, description in pathogenic_sequences:
        print(f"\n{description}")
        print(f"Sequence: {seq}")
        print(f"Length: {len(seq)} bp")
        
        # Run G4 detection
        motifs = find_gquadruplex(seq)
        print(f"G4 Motifs Detected: {len(motifs)}")
        
        if motifs:
            for i, motif in enumerate(motifs):
                print(f"  Motif {i+1}:")
                print(f"    Position: {motif['Start']}-{motif['End']}")
                print(f"    G4Hunter Score: {motif['G4Hunter_Mean']:.3f}")
                print(f"    Conservation: {motif['Conservation_Score']:.3f}")
                print(f"    Formation Category: {motif['Formation_Category']}")
                print(f"    Formation Probability: {motif['Formation_Probability']}")
                print(f"    Experimental Evidence: {motif['Experimental_Evidence']}")
        else:
            print("  No G4 motifs detected (may be due to threshold settings)")
        
        print("-" * 50)

if __name__ == "__main__":
    test_enhanced_g4_features()
    test_pathogenic_sequences()
    print("\n" + "="*70)
    print("Enhanced G4Hunter Testing Complete!")
    print("="*70)