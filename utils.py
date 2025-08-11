# ===========================
# Motif Aggregation & Reporting
# ===========================

import re
from motifs import (
    find_sticky_dna, find_curved_DNA, find_zdna, find_egz_motif,
    find_slipped_dna, find_rlfs, find_cruciform, find_hdna,
    find_gtriplex, find_gquadruplex, find_relaxed_gquadruplex,
    find_bulged_gquadruplex, find_bipartite_gquadruplex,
    find_multimeric_gquadruplex, find_imotif, find_ac_motifs
)

# Validate motif structure
def validate_motif(motif, seq_length):
    req=["Class","Subtype","Start","End","Length","Sequence"]
    if not all(k in motif for k in req): return False
    if not (1<=motif["Start"]<=motif["End"]<=seq_length): return False
    if len(motif["Sequence"].replace('\n',''))==0: return False
    return True

# Select best non-overlapping motifs of each class
def select_best_nonoverlapping_motifs(motifs, motif_priority=None):
    if motif_priority is None:
        motif_priority=['Multimeric_G4','Bipartite_G4','Dimeric_G4','Canonical_G4',
                       'Relaxed_G4','Non_canonical_G4','Bulged_G4','Three_G-Runs']
    subtype_rank={s:i for i,s in enumerate(motif_priority)}
    def motif_key(m):
        rank=subtype_rank.get(m.get('Subtype'),len(subtype_rank))
        try: score=float(m.get('Score',0))
        except: score=0.0
        length=m.get('Length',0)
        return (m.get('Class',''),rank,-score,-length)
    sorted_motifs=sorted(motifs,key=motif_key)
    selected=[]; occupied_per_class={}
    for m in sorted_motifs:
        cls=m.get('Class','Other'); region=set(range(m['Start'],m['End']+1))
        if cls not in occupied_per_class: occupied_per_class[cls]=set()
        if occupied_per_class[cls].isdisjoint(region):
            selected.append(m); occupied_per_class[cls].update(region)
    return selected

# Hybrid motif discovery (overlap of two motif classes)
def find_hybrids(motifs, seq):
    events=[(m['Start'],'start',i) for i,m in enumerate(motifs)]+[(m['End']+1,'end',i) for i,m in enumerate(motifs)]
    events.sort(); active=set(); region_start=None; results=[]
    for pos,typ,idx in events:
        if typ=='start': active.add(idx); 
        elif typ=='end': active.discard(idx)
        if typ=='start' and len(active)==2: region_start=pos
        if typ=='end' and len(active)==1 and region_start is not None:
            region_end=pos-1; involved_idxs=list(active)
            involved_classes={motifs[i]['Class'] for i in involved_idxs}
            if len(involved_classes)>=2:
                region_motifs=[motifs[i] for i in involved_idxs]
                results.append({"Class":"Hybrid","Subtype":"_".join(sorted(involved_classes))+"_Overlap",
                                "Start":region_start,"End":region_end,"Length":region_end-region_start+1,
                                "MotifClasses":sorted(involved_classes),"ContributingMotifs":region_motifs,
                                "ScoreMethod":"HybridOverlap","Score":f"{min(1.0,len(involved_classes)/5+len(region_motifs)/10):.2f}",
                                "Sequence":seq[region_start-1:region_end]})
            region_start=None
    return results

# Hotspot regions reporting
def find_hotspots(motif_hits, seq_len, window=100, min_count=3):
    hotspots=[]; positions=[(m['Start'],m['End']) for m in motif_hits]
    for i in range(0,seq_len-window+1):
        r1,r2=i+1,i+window
        count=sum(s<=r2 and e>=r1 for s,e in positions)
        if count>=min_count:
            motifs_in_region=[m for m in motif_hits if m['Start']<=r2 and m['End']>=r1]
            type_div=len({m['Subtype'] for m in motifs_in_region})
            seq_region=motif_hits[0]['Sequence'] if motif_hits else ""
            hotspots.append({"Class":"Non-B DNA Clusters","Subtype":"Hotspot","Start":r1,"End":r2,
                            "Length":r2-r1+1,"Sequence":"",
                            "ScoreMethod":"Hotspot","Score":f"{min(1.0,count/10+type_div/5):.2f}",
                            "MotifCount":count,"TypeDiversity":type_div})
    return merge_hotspots(hotspots)

def merge_hotspots(hotspots):
    if not hotspots: return []
    merged=[hotspots[0]]
    for curr in hotspots[1:]:
        last=merged[-1]
        if curr['Start']<=last['End']:
            last['End']=max(last['End'],curr['End'])
            last['Length']=last['End']-last['Start']+1
            last['MotifCount']+=curr['MotifCount']
            last['TypeDiversity']=max(last['TypeDiversity'],curr['TypeDiversity'])
            last['Score']=f"{min(1.0,float(last['Score'])+float(curr['Score'])):.2f}"
        else: merged.append(curr)
    return merged

# Basic sequence and motif statistics
def get_basic_stats(seq, motifs=None):
    seq=seq.upper(); length=len(seq); gc=sum(1 for c in seq if c in 'GC')/max(1,length)*100 if length else 0
    at=(seq.count('A')+seq.count('T'))/length*100 if length else 0
    stats={"Length":length,"GC%":round(gc,2),"AT%":round(at,2),
           "A":seq.count('A'),"T":seq.count('T'),"G":seq.count('G'),"C":seq.count('C')}
    if motifs is not None:
        covered=set()
        for m in motifs: covered.update(range(m['Start'],m['End']))
        coverage_pct=(len(covered)/length*100) if length else 0
        stats["Motif Coverage %"]=round(coverage_pct,2)
    return stats

# Aggregation: all motifs in sequence
def all_motifs(seq, nonoverlap=False, report_hotspots=False):
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE): return []
    seq=seq.upper()
    motif_list=(find_sticky_dna(seq)+find_curved_DNA(seq)+find_zdna(seq)+find_egz_motif(seq)+
                find_slipped_dna(seq)+find_rlfs(seq)+find_cruciform(seq)+find_hdna(seq)+find_gtriplex(seq)+
                find_gquadruplex(seq)+find_relaxed_gquadruplex(seq)+find_bulged_gquadruplex(seq)+
                find_bipartite_gquadruplex(seq)+find_multimeric_gquadruplex(seq)+find_imotif(seq)+find_ac_motifs(seq))
    motif_list=[m for m in motif_list if validate_motif(m,len(seq))]
    motif_list+=find_hybrids(motif_list,seq)
    if nonoverlap: motif_list=select_best_nonoverlapping_motifs(motif_list)
    if report_hotspots: motif_list+=find_hotspots(motif_list,len(seq))
    return motif_list
