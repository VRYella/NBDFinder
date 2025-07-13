import re
import numpy as np
from collections import defaultdict
from utils import gc_content, reverse_complement, g4hunter_score

# --------------------------
# Core Detection Functions
# --------------------------

def overlapping_finditer(pattern: str, seq: str) -> re.Match:
    """Generator for overlapping regex matches"""
    regex = re.compile(pattern, re.IGNORECASE)
    for match in regex.finditer(seq):
        yield match

def validate_motif(motif: Dict, seq_len: int) -> bool:
    """Validate motif coordinates and sequence"""
    return (0 < motif['Start'] <= motif['End'] <= seq_len and
            motif['Length'] == (motif['End'] - motif['Start'] + 1) and
            re.match("^[ATGC]+$", motif['Sequence'], re.IGNORECASE))

# --------------------------
# 1. CURVED DNA (nBSTAPR logic)
# --------------------------

def find_curved_dna(seq: str) -> List[Dict]:
    """
    Detect:
    - APRs with 3+ A-tracts (3-11bp) with 10.5bp spacing
    - T-tracts on reverse strand
    - Local curved DNA (A7/C7) non-overlapping with APRs
    """
    # APR detection (nBSTAPR logic)
    apr_pattern = r"(?=(A{3,11}[ATGC]{2,10}){3,})"
    tpr_pattern = r"(?=(T{3,11}[ATGC]{2,10}){3,})"
    
    curved_motifs = []
    apr_coords = set()
    
    # Process APRs and TPRs
    for pattern, subtype in [(apr_pattern, "A-Phased_Repeat"), 
                           (tpr_pattern, "T-Phased_Repeat")]:
        for m in overlapping_finditer(pattern, seq):
            seq_frag = m.group()
            a_tracts = re.findall(r"[AT]{3,11}", seq_frag)
            spacing = np.mean([m.start(i+1)-m.start(i) for i in range(len(a_tracts)-1)])
            
            score = min(1.0, len(a_tracts)/3 + (0.5 if 9 <= spacing <= 11 else 0))
            start, end = m.start()+1, m.end()
            
            curved_motifs.append({
                "Class": "Curved_DNA",
                "Subtype": subtype,
                "Start": start,
                "End": end,
                "Length": end-start+1,
                "Sequence": seq_frag,
                "ScoreMethod": "nBSTAPR",
                "Score": f"{score:.2f}"
            })
            apr_coords.update(range(start, end+1))
    
    # Local curved DNA (non-overlapping)
    for m in overlapping_finditer(r"(A{7,}|T{7,})", seq):
        start, end = m.start()+1, m.end()
        if not any(p in apr_coords for p in range(start, end+1)):
            curved_motifs.append({
                "Class": "Curved_DNA",
                "Subtype": "Local_Curve",
                "Start": start,
                "End": end,
                "Length": end-start+1,
                "Sequence": m.group(),
                "ScoreMethod": "Length_Score",
                "Score": f"{min(1.0, len(m.group())/10):.2f}"
            })
    
    return curved_motifs

# --------------------------
# 2. Z-DNA (Z-Seeker scoring)
# --------------------------

def find_zdna(seq: str) -> List[Dict]:
    """Z-DNA detection with Z-Seeker scoring (alternating purine-pyrimidine)"""
    zdna_motifs = []
    for m in overlapping_finditer(r"(([CG][GC]|[TA][AC]){6,})", seq):
        seq_frag = m.group(1)
        if len(seq_frag) < 10: continue
        
        # Z-Seeker scoring (Ho et al. 2010)
        cg_pairs = len(re.findall(r"CG|GC", seq_frag))
        score = min(1.0, cg_pairs/5 + len(seq_frag)/20)
        
        zdna_motifs.append({
            "Class": "Z-DNA",
            "Subtype": "Alternating_PuPy",
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(seq_frag),
            "Sequence": seq_frag,
            "ScoreMethod": "Z-Seeker",
            "Score": f"{score:.2f}"
        })
    return zdna_motifs

# --------------------------
# 3. SLIPPED DNA (nBST logic)
# --------------------------

def find_slipped_dna(seq: str) -> List[Dict]:
    """Direct repeats with 10-300bp units and ≤100bp spacers"""
    slipped_motifs = []
    for m in overlapping_finditer(r"(.{10,300})(.{0,100})\1", seq):
        unit, spacer = m.group(1), m.group(2)
        if len(spacer) > 100: continue
        
        score = min(1.0, len(unit)/150 + (1 - len(spacer)/100))
        slipped_motifs.append({
            "Class": "Slipped_DNA",
            "Subtype": f"DR_{len(unit)}bp",
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(unit)*2 + len(spacer),
            "Sequence": unit + spacer + unit,
            "ScoreMethod": "nBST_DR",
            "Score": f"{score:.2f}"
        })
    return slipped_motifs

# --------------------------
# 4. CRUCIFORM (nBSTAPR logic)
# --------------------------

def find_cruciform(seq: str) -> List[Dict]:
    """Inverted repeats with ≥6bp arms and ≤100bp loops"""
    cruciforms = []
    for m in overlapping_finditer(r"([ATGC]{6,}).{0,100}\1", seq):
        arm = m.group(1)
        loop = m.group(2) if m.group(2) else ""
        
        # nBSTAPR scoring
        at_content = (arm.count('A') + arm.count('T'))/len(arm)
        score = min(1.0, len(arm)/12 + at_content*0.5)
        
        cruciforms.append({
            "Class": "Cruciform",
            "Subtype": f"Hairpin_{len(arm)}bp",
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(arm)*2 + len(loop),
            "Sequence": arm + loop + arm,
            "ScoreMethod": "nBSTAPR",
            "Score": f"{score:.2f}"
        })
    return cruciforms

# --------------------------
# 5. TRIPLEX & STICKY DNA
# --------------------------

def find_triplex(seq: str) -> List[Dict]:
    """Mirror repeats ≥10bp with ≤100bp loops + sticky DNA (GAA/TTC)"""
    triplex_motifs = []
    
    # Mirror repeats (nBST logic)
    for m in overlapping_finditer(r"([CT]{10,}).{0,100}\1", seq):
        arm = m.group(1)
        score = min(1.0, len(arm)/15 + arm.count('C')/len(arm))
        triplex_motifs.append({
            "Class": "Triplex_DNA",
            "Subtype": "Mirror_Repeat",
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(arm)*2 + len(m.group(2)),
            "Sequence": arm + (m.group(2) if m.group(2) else "") + arm,
            "ScoreMethod": "nBST_Mirror",
            "Score": f"{score:.2f}"
        })
    
    # Sticky DNA (GAA/TTC repeats)
    for m in overlapping_finditer(r"(GAA){3,}|(TTC){3,}", seq):
        triplex_motifs.append({
            "Class": "Triplex_DNA",
            "Subtype": "Sticky_DNA",
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(m.group()),
            "Sequence": m.group(),
            "ScoreMethod": "Length_Score",
            "Score": f"{min(1.0, len(m.group())/9):.2f}"
        })
    
    return triplex_motifs

# --------------------------
# 6. G-TRIPLEX (non-G4)
# --------------------------

def find_gtriplex(seq: str) -> List[Dict]:
    """G-triplex not part of G-quadruplex"""
    g4_coords = [(m.start(), m.end()) for m in 
                overlapping_finditer(r"(G{3,}\w{1,7}){4,}", seq)]
    
    gtriplex_motifs = []
    for m in overlapping_finditer(r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,})", seq):
        start, end = m.start()+1, m.end()
        if any(s <= end and e >= start for s,e in g4_coords): continue
        
        g_runs = re.findall(r"G{3,}", m.group())
        score = min(1.0, sum(len(r) for r in g_runs)/9)
        gtriplex_motifs.append({
            "Class": "G-Triplex",
            "Subtype": "Isolated_G3",
            "Start": start,
            "End": end,
            "Length": len(m.group()),
            "Sequence": m.group(),
            "ScoreMethod": "G3_Score",
            "Score": f"{score:.2f}"
        })
    return gtriplex_motifs

# --------------------------
# 7. G-QUADRUPLEX (G4Hunter)
# --------------------------

def find_gquadruplex(seq: str) -> List[Dict]:
    """Hierarchical G4 detection with G4Hunter scoring"""
    g4_motifs = []
    
    # Detection order: Multimeric > Bipartite > Canonical > Relaxed > Bulged
    patterns = [
        (r"(G{3,}\w{1,12}){6,}", "Multimeric_G4", 1.2),  # Bonus for multimers
        (r"(G{3,}\w{1,30}){6,}", "Bipartite_G4", 1.1),    # Bipartite bonus
        (r"(G{3,}\w{1,7}){4,}", "Canonical_G4", 1.0),
        (r"(G{3,}\w{1,12}){4,}", "Relaxed_G4", 0.8),
        (r"(G{3,}\w{0,3}){4,}", "Bulged_G4", 0.7)
    ]
    
    detected_coords = set()
    for pattern, subtype, weight in patterns:
        for m in overlapping_finditer(pattern, seq):
            start, end = m.start()+1, m.end()
            if any(s <= end and e >= start for s,e in detected_coords): continue
            
            seq_frag = m.group()
            score = g4hunter_score(seq_frag) * weight
            if score >= (0.7 if "Bulged" in subtype else 1.0):
                g4_motifs.append({
                    "Class": "G4",
                    "Subtype": subtype,
                    "Start": start,
                    "End": end,
                    "Length": len(seq_frag),
                    "Sequence": seq_frag,
                    "ScoreMethod": "G4Hunter",
                    "Score": f"{score:.2f}"
                })
                detected_coords.add((start, end))
    
    return sorted(g4_motifs, key=lambda x: -float(x['Score']))

# --------------------------
# 8. I-MOTIF (G4Hunter)
# --------------------------

def find_imotif(seq: str) -> List[Dict]:
    """C-rich i-Motifs with 1-7nt loops"""
    imotifs = []
    for m in overlapping_finditer(r"(C{3,}\w{0,7}C{3,}\w{0,7}C{3,}\w{0,7}C{3,})", seq):
        seq_frag = m.group()
        c_runs = re.findall(r"C{3,}", seq_frag)
        score = min(1.0, sum(len(r) for r in c_runs)/12)
        imotifs.append({
            "Class": "i-Motif",
            "Subtype": "C_Quadruplex",
            "Start": m.start()+1,
            "End": m.end(),
            "Length": len(seq_frag),
            "Sequence": seq_frag,
            "ScoreMethod": "G4Hunter_Adapted",
            "Score": f"{score:.2f}"
        })
    return imotifs

# --------------------------
# 9. HYBRID MOTIFS
# --------------------------

def find_hybrids(seq: str, motifs: List[Dict]) -> List[Dict]:
    """Overlapping motifs of different classes"""
    class_map = defaultdict(list)
    for m in motifs:
        class_map[(m['Start'], m['End'])].append(m['Class'])
    
    hybrids = []
    for (start, end), classes in class_map.items():
        if len(set(classes)) >= 2:
            seq_frag = seq[start-1:end]
            rev_comp = reverse_complement(seq_frag)
            
            # Check both strands
            for s in [seq_frag, rev_comp]:
                hybrids.append({
                    "Class": "Hybrid",
                    "Subtype": "+".join(sorted(set(classes))),
                    "Start": start,
                    "End": end,
                    "Length": end-start+1,
                    "Sequence": s,
                    "ScoreMethod": "Multi_Class",
                    "Score": f"{min(1.0, len(set(classes))/3):.2f}"
                })
    return hybrids

# --------------------------
# 10. NON-B DNA CLUSTERS
# --------------------------

def find_clusters(motifs: List[Dict], seq_len: int) -> List[Dict]:
    """Regions with ≥3 non-B motifs in 100bp"""
    density = np.zeros(seq_len)
    for m in motifs:
        density[m['Start']-1:m['End']] += 1
    
    clusters = []
    for i in range(0, seq_len - 100 + 1):
        window = density[i:i+100]
        if np.sum(window >= 1) >= 3:  # At least 3 motifs
            motifs_in_window = [m for m in motifs 
                              if m['Start'] <= i+100 and m['End'] >= i+1]
            
            # Calculate metrics
            total_motif_bp = sum(m['Length'] for m in motifs_in_window)
            unique_types = len(set(m['Class'] for m in motifs_in_window))
            
            clusters.append({
                "Start": i+1,
                "End": i+100,
                "MotifCount": len(motifs_in_window),
                "TypeDiversity": unique_types,
                "Coverage": f"{total_motif_bp/100:.2%}",
                "Score": f"{min(1.0, len(motifs_in_window)/5 + unique_types/3):.2f}"
            })
    
    return clusters

# --------------------------
# Main Detection Pipeline
# --------------------------

def all_motifs(seq: str) -> List[Dict]:
    """Execute all detectors in specified order"""
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    
    seq = seq.upper()
    detectors = [
        find_curved_dna, find_zdna, find_slipped_dna,
        find_cruciform, find_triplex, find_gtriplex,
        find_gquadruplex, find_imotif
    ]
    
    # Run primary detectors
    motifs = []
    for detector in detectors:
        try:
            motifs.extend(detector(seq))
        except Exception as e:
            print(f"Error in {detector.__name__}: {e}")
    
    # Add hybrid detection
    motifs.extend(find_hybrids(seq, motifs))
    
    # Add clusters
    motifs.extend(find_clusters(motifs, len(seq)))
    
    return [m for m in motifs if validate_motif(m, len(seq))]
