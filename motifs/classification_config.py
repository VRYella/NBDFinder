"""
NBDFinder Motif Classification Configuration
==========================================

This module defines the exact 10 Non-B DNA classes and 22 subclasses
as specified in the requirements.

Scientific Classification Structure:
1. Curved DNA (2 subclasses)
2. Slipped DNA (2 subclasses) 
3. Cruciform DNA (1 subclass)
4. R-loop (1 subclass)
5. Triplex (2 subclasses)
6. G-Quadruplex Family (7 subclasses)
7. i-motif family (3 subclasses)
8. Z-DNA (2 subclasses)
9. Hybrid (variable subclasses based on overlaps)
10. Non-B DNA cluster regions (variable subclasses)

Author: Dr. Venkata Rajesh Yella
Updated: 2024
"""

# Official 10 Non-B DNA Classes and 22 Subclasses Configuration
OFFICIAL_CLASSIFICATION = {
    1: {
        "class_name": "Curved DNA",
        "subclasses": [
            "Global curvature",
            "Local Curvature"
        ]
    },
    2: {
        "class_name": "Slipped DNA", 
        "subclasses": [
            "Slipped DNA [Direct Repeat]",
            "Slipped DNA [STR]"
        ]
    },
    3: {
        "class_name": "Cruciform DNA",
        "subclasses": [
            "Cruciform DNA [IR]/HairPin [IR]"
        ]
    },
    4: {
        "class_name": "R-loop",
        "subclasses": [
            "R-loop"
        ]
    },
    5: {
        "class_name": "Triplex",
        "subclasses": [
            "Triplex",
            "sticky DNA"
        ]
    },
    6: {
        "class_name": "G-Quadruplex Family",
        "subclasses": [
            "Multimeric G4",
            "Canonical G4", 
            "Relaxed G4",
            "Bulged G4",
            "Bipartite G4",
            "Imperfect G4",
            "G-Triplex intermediate"
        ],
        "priority_order": [
            "Multimeric G4",
            "Canonical G4",
            "Relaxed G4", 
            "Bulged G4",
            "Bipartite G4",
            "Imperfect G4",
            "G-Triplex intermediate"
        ]
    },
    7: {
        "class_name": "i-motif family",
        "subclasses": [
            "Canonical i-motif",
            "Relaxed i-motif", 
            "AC-motif"
        ]
    },
    8: {
        "class_name": "Z-DNA",
        "subclasses": [
            "Z-DNA",
            "eGZ (Extruded-G) DNA"
        ]
    },
    9: {
        "class_name": "Hybrid",
        "subclasses": []  # Dynamic based on overlaps between any two classes
    },
    10: {
        "class_name": "Non-B DNA cluster regions",
        "subclasses": []  # Dynamic based on occurring motifs: any three classes occurring 3+ times in 100 nt
    }
}

# Mapping from current implementation to official classification
LEGACY_TO_OFFICIAL_MAPPING = {
    # Current -> Official
    "Curved_DNA": "Curved DNA",
    "Local_Curved_Strict_PolyA_or_PolyT": "Local Curvature",
    "Global_Curved_PolyA_PolyT": "Global curvature",
    
    "Slipped_DNA": "Slipped DNA",
    "Direct_Repeat": "Slipped DNA [Direct Repeat]",
    "STR": "Slipped DNA [STR]",
    
    "Cruciform": "Cruciform DNA", 
    "Cruciform_IR": "Cruciform DNA [IR]/HairPin [IR]",
    
    "R-Loop": "R-loop",
    "RLFS": "R-loop",
    
    "Mirror_Repeat": "Triplex",
    "Triplex_Motif": "Triplex",
    "Sticky_DNA": "sticky DNA",
    
    "G-Triplex": "G-Quadruplex Family",
    "G_Triplex": "G-Triplex intermediate",
    "Canonical G4": "G-Quadruplex Family",
    "G4_Canonical": "Canonical G4",
    "Relaxed G4": "G-Quadruplex Family", 
    "G4_Relaxed": "Relaxed G4",
    "Bulged G4": "G-Quadruplex Family",
    "G4_Bulged": "Bulged G4",
    "Bipartite G4": "G-Quadruplex Family",
    "G4_Bipartite": "Bipartite G4",
    "Multimeric G4": "G-Quadruplex Family",
    "G4_Multimeric": "Multimeric G4",
    "Imperfect G4": "G-Quadruplex Family", 
    "G4_Imperfect": "Imperfect G4",
    
    "i-Motif": "i-motif family",
    "Canonical_iMotif": "Canonical i-motif",
    "Relaxed_iMotif": "Relaxed i-motif",
    "AC-Motif": "i-motif family",
    "AC_Motif": "AC-motif",
    
    "Z-DNA": "Z-DNA",
    "Z-Seeker_Kadane": "Z-DNA",
    "eGZ": "Z-DNA",
    "CGG_Expansion": "eGZ (Extruded-G) DNA",
    
    "Hybrid": "Hybrid",
    "Non-B_Clusters": "Non-B DNA cluster regions"
}

# G4 family priority enforcement
G4_PRIORITY_ORDER = [
    "Multimeric G4",
    "Canonical G4", 
    "Relaxed G4",
    "Bulged G4",
    "Bipartite G4", 
    "Imperfect G4",
    "G-Triplex intermediate"
]

def get_official_classification(legacy_class, legacy_subtype):
    """
    Convert legacy classification to official 10-class system.
    
    Args:
        legacy_class (str): Current class name
        legacy_subtype (str): Current subtype name
        
    Returns:
        tuple: (official_class, official_subtype)
    """
    # Map legacy class to official class
    official_class = LEGACY_TO_OFFICIAL_MAPPING.get(legacy_class, legacy_class)
    
    # Map legacy subtype to official subtype
    official_subtype = LEGACY_TO_OFFICIAL_MAPPING.get(legacy_subtype, legacy_subtype)
    
    return official_class, official_subtype

def apply_g4_priority_filter(motifs):
    """
    Apply G4 family priority order: Multimeric G4 > Canonical G4 > Relaxed G4 > 
    Bulged G4 > Bipartite G4 > Imperfect G4 > G-Triplex intermediate
    
    Note: Canonical G4 is not screened if G-triplex is shown.
    
    Args:
        motifs (list): List of detected motifs
        
    Returns:
        list: Filtered motifs with G4 priority applied
    """
    # Group G4 motifs by genomic position
    g4_groups = {}
    non_g4_motifs = []
    
    for motif in motifs:
        if motif.get("Class") == "G-Quadruplex Family":
            start = motif.get("Start", 0)
            end = motif.get("End", 0)
            # Group overlapping G4 motifs
            position_key = (start, end)
            if position_key not in g4_groups:
                g4_groups[position_key] = []
            g4_groups[position_key].append(motif)
        else:
            non_g4_motifs.append(motif)
    
    # Apply priority filtering within each group
    filtered_g4_motifs = []
    for position_key, group_motifs in g4_groups.items():
        if len(group_motifs) == 1:
            filtered_g4_motifs.append(group_motifs[0])
        else:
            # Sort by priority order
            group_motifs.sort(key=lambda m: G4_PRIORITY_ORDER.index(
                m.get("Subtype", "G-Triplex intermediate")
                if m.get("Subtype", "") in G4_PRIORITY_ORDER 
                else len(G4_PRIORITY_ORDER)
            ))
            
            # Take highest priority motif
            filtered_g4_motifs.append(group_motifs[0])
    
    return non_g4_motifs + filtered_g4_motifs

def get_class_counts():
    """Return the official class and subclass counts for validation."""
    total_classes = len(OFFICIAL_CLASSIFICATION)
    total_subclasses = sum(len(class_info["subclasses"]) for class_info in OFFICIAL_CLASSIFICATION.values())
    return total_classes, total_subclasses

def validate_classification():
    """Validate that we have exactly 10 classes and 22+ subclasses."""
    classes, subclasses = get_class_counts()
    
    expected_fixed_subclasses = 20  # 2+2+1+1+2+7+3+2 = 20 fixed subclasses
    
    print(f"Official classification validation:")
    print(f"  Classes: {classes} (expected: 10)")
    print(f"  Fixed subclasses: {subclasses} (expected: ≥{expected_fixed_subclasses})")
    print(f"  Variable subclasses: Hybrid + Cluster regions (dynamic)")
    
    return classes == 10 and subclasses >= expected_fixed_subclasses

if __name__ == "__main__":
    validate_classification()