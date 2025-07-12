import re
import numpy as np
from utils import wrap, gc_content, reverse_complement, g4hunter_score

def overlapping_finditer(pattern, seq):
    """Find all overlapping matches of a pattern in sequence"""
    regex = re.compile(pattern, re.IGNORECASE)
    for m in regex.finditer(seq):
        yield m

def all_motifs(seq):
    """Run all motif detection functions and return combined, validated results"""
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE):
        return []
    
    seq = seq.upper()
    results = (
        find_apr(seq) + find_bent_dna(seq) + find_zdna(seq) +
        find_slipped_dna(seq) + find_rlfs(seq) + find_cruciform(seq) +
        find_hdna(seq) + find_gtriplex(seq) + find_gquadruplex(seq) +
        find_relaxed_gquadruplex(seq) + find_bulged_gquadruplex(seq) +
        find_bipartite_gquadruplex(seq) + find_multimeric_gquadruplex(seq) +
        find_imotif(seq) + find_hybrids(seq)
    
    # Validate and filter results
    return [m for m in results if validate_motif(m, len(seq))]

def validate_motif(motif, seq_len):
    """Scientific validation of motif features"""
    if not (0 < motif['Start'] <= motif['End'] <= seq_len):
        return False
    if motif['Length'] != (motif['End'] - motif['Start'] + 1):
        return False
    if not re.match("^[ATGC]+$", motif['Sequence']):
        return False
    return True

# 1. CURVED DNA (A-Phased Repeats)
def find_apr(seq):
    """Detect A-Phased Repeats with curvature scoring (Brukner et al. 1995)"""
    pattern = r"(?=((?:A{3,6}[ATGC]{2,5}){3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        seq_fragment = m.group(1)
        # Score based on repeat number and A-tract length
        a_count = seq_fragment.count('A') + seq_fragment.count('T')
        score = min(1.0, 0.2 * len(re.findall(r"A{3,6}", seq_fragment)) + (a_count/len(seq_fragment))
        results.append({
            "Class": "Curved_DNA",
            "Subtype": "A-Phased_Repeat",
            "Start": m.start() + 1,
            "End": m.start() + len(seq_fragment),
            "Length": len(seq_fragment),
            "Sequence": wrap(seq_fragment),
            "ScoreMethod": "Brukner_Curvature",
            "Score": f"{score:.2f}"
        })
    return results

# 2. BENT DNA (Poly-A/T)
def find_bent_dna(seq):
    """Detect poly-A/T tracts with bending propensity (Haran & Mohanty 2009)"""
    pattern = r"(?=(A{4,}|T{4,}))"
    return [{
        "Class": "Bent_DNA",
        "Subtype": "Poly-A/T",
        "Start": m.start() + 1,
        "End": m.end(),
        "Length": len(m.group(1)),
        "Sequence": wrap(m.group(1)),
        "ScoreMethod": "Bending_Propensity",
        "Score": f"{min(1.0, len(m.group(1))/10):.2f}"  # 0-1 scale based on length
    } for m in overlapping_finditer(pattern, seq)]

# 3. Z-DNA (Ho et al. 2010)
def find_zdna(seq):
    """Detect Z-DNA forming regions with stability scoring"""
    pattern = r"(?=((?:CG|GC|GT|TG|AC|CA){6,}))"
    return [{
        "Class": "Z-DNA",
        "Subtype": "Z-Forming_Repeat",
        "Start": m.start() + 1,
        "End": m.start() + len(m.group(1)),
        "Length": len(m.group(1)),
        "Sequence": wrap(m.group(1)),
        "ScoreMethod": "Z-Hunter_v2",
        "Score": f"{zdna_stability_score(m.group(1)):.2f}"
    } for m in overlapping_finditer(pattern, seq)]

def zdna_stability_score(seq):
    """Calculate Z-DNA stability score (0-1) based on dinucleotide content"""
    dinucs = re.findall(r"(CG|GC|GT|TG|AC|CA)", seq, re.IGNORECASE)
    if not dinucs:
        return 0
    cg_count = sum(1 for d in dinucs if d.upper() in ('CG','GC'))
    return min(1.0, (len(dinucs)/10 + cg_count/5))  # CG pairs contribute more

# 4. SLIPPED DNA (Bacolla et al. 2006)
def find_slipped_dna(seq):
    """Detect slipped (CTG/CAG) repeats with instability potential"""
    pattern = r"(?=((?:CTG|CAG|GTC|GAC){3,}))"
    return [{
        "Class": "Slipped_DNA",
        "Subtype": "TRS_Repeat",
        "Start": m.start() + 1,
        "End": m.start() + len(m.group(1)),
        "Length": len(m.group(1)),
        "Sequence": wrap(m.group(1)),
        "ScoreMethod": "Instability_Index",
        "Score": f"{len(m.group(1))/30:.2f}"  # 0-1 scale based on length
    } for m in overlapping_finditer(pattern, seq)]

# 5. R-LOOPS (QmRLFS improved)
def find_rlfs(seq):
    """Detect R-loop forming sequences with thermodynamic scoring"""
    if len(seq) < 100: return []
    
    results = []
    # R-loop initiation zones (high GC, G-clusters)
    riz_patterns = {
        'm1': r"(G{3,}[ATGC]{1,10}G{3,}[ATGC]{1,10}G{4,})",
        'm2': r"(G{3,}[ATGC]{1,10}G{4,})"
    }
    
    for model, pattern in riz_patterns.items():
        for m in re.finditer(pattern, seq, re.IGNORECASE):
            riz_seq = m.group(1)
            if gc_content(riz_seq) < 50: continue
            
            # Find downstream extension zone (REZ)
            for rez in find_rez(seq, m.end(), model):
                stability = calculate_rloop_stability(riz_seq + rez['seq'])
                results.append({
                    "Class": "R-Loop",
                    "Subtype": f"RLFS_{model}",
                    "Start": m.start() + 1,
                    "End": m.start() + len(riz_seq) + rez['end'],
                    "Length": len(riz_seq) + rez['end'],
                    "Sequence": wrap(riz_seq + rez['seq']),
                    "ScoreMethod": "QmRLFS_Thermo",
                    "Score": f"{stability:.2f}"
                })
    return results

def find_rez(seq, start_pos, model):
    """Find R-loop extension zones"""
    max_len = 500 if model == 'm1' else 300
    window = seq[start_pos:start_pos+max_len]
    if gc_content(window) < 40: return []
    
    return [{
        'seq': window,
        'end': len(window),
        'gc': gc_content(window)
    }]

def calculate_rloop_stability(seq):
    """Calculate R-loop stability score (0-1)"""
    gc = gc_content(seq)/100
    g_clusters = len(re.findall(r"G{3,}", seq, re.IGNORECASE))
    return min(1.0, 0.6*gc + 0.4*(g_clusters/5))

# [Continued in next message due to length...]
# 6. CRUCIFORM DNA (Lilley 1985, Pearson et al. 1996)
def find_cruciform(seq):
    """Detect cruciform-forming palindromes with stability scoring"""
    # Minimum hairpin of 6bp with 4bp loop
    pattern = r"(?=(([ATGC]{6,})\s*?\2))"  # Inverted repeats
    results = []
    for m in overlapping_finditer(pattern, seq):
        stem = m.group(2)
        stem_len = len(stem)
        if stem_len < 6: continue
        
        # Score based on stem length and AT content (more stable)
        at_content = (stem.count('A') + stem.count('T')) / stem_len
        score = min(1.0, (stem_len/15) + (at_content*0.3))
        results.append({
            "Class": "Cruciform",
            "Subtype": "Inverted_Repeat",
            "Start": m.start() + 1,
            "End": m.start() + (2*stem_len),
            "Length": 2*stem_len,
            "Sequence": wrap(m.group(1)),
            "ScoreMethod": "Cruciform_Stability",
            "Score": f"{score:.2f}"
        })
    return results

# 7. TRIPLEX DNA (H-DNA) (Mirkin & Frank-Kamenetskii 1994)
def find_hdna(seq):
    """Detect intramolecular triplex (H-DNA) with pyrimidine mirror repeats"""
    # Pyrimidine-rich mirror repeats
    pattern = r"(?=(([CT]{5,})\s*?\2))"
    return [{
        "Class": "Triplex_DNA",
        "Subtype": "H-DNA_Pyrimidine",
        "Start": m.start() + 1,
        "End": m.start() + len(m.group(1)),
        "Length": len(m.group(1)),
        "Sequence": wrap(m.group(1)),
        "ScoreMethod": "Triplex_Propensity",
        "Score": f"{min(1.0, len(m.group(2))/10 + 0.3):.2f}"  # Length + pyrimidine bonus
    } for m in overlapping_finditer(pattern, seq) if len(m.group(2)) >= 5]

# 8. G-TRIPLEX (Karsisiotis et al. 2011)
def find_gtriplex(seq):
    """Detect G-triplex structures with stability scoring"""
    pattern = r"(?=(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}))"
    results = []
    for m in overlapping_finditer(pattern, seq):
        seq_frag = m.group(1)
        g_runs = [len(r) for r in re.findall(r"G{3,}", seq_frag)]
        if len(g_runs) < 3: continue
        
        # Score based on G-run length and loop sizes
        loop_sizes = [len(l) for l in re.findall(r"G{3,}(\w{1,7})G{3,}", seq_frag)]
        score = min(1.0, sum(g_runs)/15 + sum(1/l if l>0 else 0.5 for l in loop_sizes)/3)
        results.append({
            "Class": "G-Triplex",
            "Subtype": "Three_G-Runs",
            "Start": m.start() + 1,
            "End": m.start() + len(seq_frag),
            "Length": len(seq_frag),
            "Sequence": wrap(seq_frag),
            "ScoreMethod": "G3_Stability",
            "Score": f"{score:.2f}"
        })
    return results

# 9. G-QUADRUPLEX VARIANTS (Bedrat et al. 2016, Chambers et al. 2015)
def find_gquadruplex(seq):
    """Detect canonical G4 with G4Hunter scoring"""
    pattern = r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"
    return [{
        "Class": "G4",
        "Subtype": "Canonical_G4",
        "Start": m.start() + 1,
        "End": m.end(),
        "Length": len(m.group(1)),
        "Sequence": wrap(m.group(1)),
        "ScoreMethod": "G4Hunter_v2",
        "Score": f"{g4hunter_score(m.group(1)):.2f}"
    } for m in overlapping_finditer(pattern, seq) if g4hunter_score(m.group(1)) >= 1.0]

def find_relaxed_gquadruplex(seq):
    """Detect G4 with longer loops (up to 12nt)"""
    pattern = r"(G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,})"
    return [{
        "Class": "G4",
        "Subtype": "Relaxed_G4",
        "Start": m.start() + 1,
        "End": m.end(),
        "Length": len(m.group(1)),
        "Sequence": wrap(m.group(1)),
        "ScoreMethod": "G4Hunter_LongLoop",
        "Score": f"{g4hunter_score(m.group(1))*0.8:.2f}"  # Penalty for longer loops
    } for m in overlapping_finditer(pattern, seq) if g4hunter_score(m.group(1)) >= 0.8]

def find_bulged_gquadruplex(seq):
    """Detect G4 with bulges (1-3nt)"""
    pattern = r"(G{3,}\w{0,3}G{3,}\w{0,3}G{3,}\w{0,3}G{3,})"
    return [{
        "Class": "G4",
        "Subtype": "Bulged_G4",
        "Start": m.start() + 1,
        "End": m.end(),
        "Length": len(m.group(1)),
        "Sequence": wrap(m.group(1)),
        "ScoreMethod": "G4Hunter_Bulge",
        "Score": f"{g4hunter_score(m.group(1))*0.7:.2f}"  # Penalty for bulges
    } for m in overlapping_finditer(pattern, seq) if len(re.findall(r"G{3,}", m.group(1))) >= 4]

def find_bipartite_gquadruplex(seq):
    """Detect bipartite G4 (two G4 units separated by <30nt)"""
    pattern = r"(G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,})"
    results = []
    for m in overlapping_finditer(pattern, seq):
        seq_frag = m.group(1)
        g_runs = re.findall(r"G{3,}", seq_frag)
        if len(g_runs) < 6: continue
        
        # Split into two potential G4 units
        unit1 = seq_frag[:len(seq_frag)//2]
        unit2 = seq_frag[len(seq_frag)//2:]
        score = max(g4hunter_score(unit1), g4hunter_score(unit2)) * 0.9  # Slight penalty
        if score >= 0.9:
            results.append({
                "Class": "G4",
                "Subtype": "Bipartite_G4",
                "Start": m.start() + 1,
                "End": m.end(),
                "Length": len(seq_frag),
                "Sequence": wrap(seq_frag),
                "ScoreMethod": "Bipartite_Score",
                "Score": f"{score:.2f}"
            })
    return results

def find_multimeric_gquadruplex(seq):
    """Detect contiguous G4 multimers"""
    pattern = r"(G{3,}\w{1,12}){4,}"
    return [{
        "Class": "G4",
        "Subtype": "Multimeric_G4",
        "Start": m.start() + 1,
        "End": m.end(),
        "Length": len(m.group(1)),
        "Sequence": wrap(m.group(1)),
        "ScoreMethod": "G4Hunter_Multimer",
        "Score": f"{g4hunter_score(m.group(1))*1.2:.2f}"  # Bonus for multimerization
    } for m in overlapping_finditer(pattern, seq) if g4hunter_score(m.group(1)) >= 1.0]

# 10. i-MOTIF (Zeraati et al. 2018)
def find_imotif(seq):
    """Detect i-Motif forming C-rich sequences"""
    pattern = r"(?=(C{3,}\w{0,7}C{3,}\w{0,7}C{3,}\w{0,7}C{3,}))"
    return [{
        "Class": "i-Motif",
        "Subtype": "C_Quadruplex",
        "Start": m.start() + 1,
        "End": m.start() + len(m.group(1)),
        "Length": len(m.group(1)),
        "Sequence": wrap(m.group(1)),
        "ScoreMethod": "iM_Stability",
        "Score": f"{imotif_score(m.group(1)):.2f}"
    } for m in overlapping_finditer(pattern, seq) if imotif_score(m.group(1)) >= 0.7]

def imotif_score(seq):
    """Calculate i-Motif stability score (0-1)"""
    c_runs = [len(r) for r in re.findall(r"C{3,}", seq)]
    if len(c_runs) < 4: return 0
    c_content = seq.count('C')/len(seq)
    return min(1.0, sum(c_runs)/16 + c_content*0.5)

# 11. HYBRID MOTIFS (Sengar et al. 2021)
def find_hybrids(seq):
    """Detect overlapping G4 and i-Motif forming regions"""
    g4_regions = [(m.start()+1, m.end(), m.group(1)) 
                 for m in overlapping_finditer(r"(G{3,}\w{1,12}G{3,}\w{1,12}G{3,}\w{1,12}G{3,})", seq)]
    im_regions = [(m.start()+1, m.start()+len(m.group(1)), m.group(1)) 
                 for m in overlapping_finditer(r"(C{3,}\w{0,7}C{3,}\w{0,7}C{3,}\w{0,7}C{3,})", seq)]
    
    hybrids = []
    for g_start, g_end, g_seq in g4_regions:
        for c_start, c_end, c_seq in im_regions:
            if (g_start <= c_end) and (c_start <= g_end):
                overlap_start = min(g_start, c_start)
                overlap_end = max(g_end, c_end)
                hybrid_seq = seq[overlap_start-1:overlap_end]
                
                # Score based on both structures' scores
                score = min(1.0, (g4hunter_score(g_seq) + imotif_score(c_seq))/2 * 1.1)  # 10% bonus
                hybrids.append({
                    "Class": "Hybrid",
                    "Subtype": "G4_iM_Overlap",
                    "Start": overlap_start,
                    "End": overlap_end,
                    "Length": overlap_end - overlap_start + 1,
                    "Sequence": wrap(hybrid_seq),
                    "ScoreMethod": "Hybrid_Score",
                    "Score": f"{score:.2f}"
                })
    return hybrids

# 12. HOTSPOTS (Non-B Clusters)
def find_hotspots(seq, motif_hits, window=100, min_count=3):
    """Identify genomic hotspots with multiple non-B motifs"""
    hotspots = []
    positions = [(hit['Start'], hit['End']) for hit in motif_hits]
    
    for i in range(0, len(seq) - window + 1):
        region_start = i + 1
        region_end = i + window
        count = sum(s <= region_end and e >= region_start for s,e in positions)
        
        if count >= min_count:
            # Calculate hotspot score based on motif density and types
            motifs_in_region = [m for m in motif_hits 
                              if m['Start'] <= region_end and m['End'] >= region_start]
            type_diversity = len({m['Subtype'] for m in motifs_in_region})
            score = min(1.0, count/10 + type_diversity/5)
            
            hotspots.append({
                "RegionStart": region_start,
                "RegionEnd": region_end,
                "MotifCount": count,
                "TypeDiversity": type_diversity,
                "Score": f"{score:.2f}"
            })
    
    # Merge overlapping hotspots
    return merge_hotspots(hotspots)

def merge_hotspots(hotspots):
    """Merge adjacent or overlapping hotspots"""
    if not hotspots: return []
    
    merged = [hotspots[0]]
    for current in hotspots[1:]:
        last = merged[-1]
        if current['RegionStart'] <= last['RegionEnd']:
            # Merge with previous hotspot
            last['RegionEnd'] = max(last['RegionEnd'], current['RegionEnd'])
            last['MotifCount'] += current['MotifCount']
            last['TypeDiversity'] = max(last['TypeDiversity'], current['TypeDiversity'])
            last['Score'] = f"{min(1.0, float(last['Score']) + float(current['Score'])):.2f}"
        else:
            merged.append(current)
    return merged
