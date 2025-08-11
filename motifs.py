# ===========================
# Motif Definitions & Scoring
# ===========================

import re; import numpy as np

# PolyA/PolyT tract finder
def find_polyA_polyT_tracts(seq, min_len=7):
    results = []; i=0; n=len(seq)
    while i<n:
        if seq[i] in 'AT':
            s=i; while i<n and seq[i]==seq[s]: i+=1
            if i-s>=min_len: results.append((s,i-1,seq[s:i]))
        else: i+=1
    return results

# Curvature scoring as tract length
def curvature_score(seq): return len(seq)

# Reverse complement
def reverse_complement(seq): return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

# GC content (%)
def gc_content(seq): gc=seq.count('G')+seq.count('C'); return (gc/max(1,len(seq)))*100 if seq else 0

# Overlapping regex search
def overlapping_finditer(pattern, seq):
    regex=re.compile(pattern,re.IGNORECASE); pos=0
    while pos<len(seq):
        m=regex.search(seq,pos)
        if not m: break
        yield m; pos=m.start()+1

# Z-DNA seeker scoring (weighted walk)
def zdna_seeker_scoring_array(seq, GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25,
        consecutive_AT_scoring=(0.5,0.5,0.5,0.5,0.0,0.0,-5.0,-100.0), mismatch_penalty_type="linear",
        mismatch_penalty_starting_value=3, mismatch_penalty_linear_delta=3, cadence_reward=0.0):
    scoring=np.empty(len(seq)-1,dtype=float); mismatches=0; at_count=0
    for i in range(len(seq)-1):
        t=seq[i:i+2].upper()
        if t in ("GC","CG"): scoring[i]=GC_weight; mismatches=0; at_count=0
        elif t in ("GT","TG"): scoring[i]=GT_weight; mismatches=0; at_count=0
        elif t in ("AC","CA"): scoring[i]=AC_weight; mismatches=0; at_count=0
        elif t in ("AT","TA"):
            w=AT_weight+(consecutive_AT_scoring[at_count] if at_count<len(consecutive_AT_scoring) else consecutive_AT_scoring[-1])
            scoring[i]=w; at_count+=1; mismatches=0
        else:
            mismatches+=1; at_count=0
            if mismatch_penalty_type=="exponential": scoring[i]=-mismatch_penalty_starting_value**mismatches if mismatches<15 else -32000
            elif mismatch_penalty_type=="linear": scoring[i]=-mismatch_penalty_starting_value-mismatch_penalty_linear_delta*(mismatches-1)
            else: scoring[i]=-10
        if t in ("GC","CG","GT","TG","AC","CA","AT","TA"): scoring[i]+=cadence_reward
    return scoring

# G4Hunter score (G/C balance)
def g4hunter_score(seq):
    scores=[1 if c=='G' else -1 if c=='C' else 0 for c in seq.upper()]
    return np.mean(scores) if scores else 0

# Purine/pyrimidine fraction
def purine_fraction(seq): return (seq.count('A')+seq.count('G'))/max(1,len(seq))
def pyrimidine_fraction(seq): return (seq.count('C')+seq.count('T'))/max(1,len(seq))

# i-Motif scoring
def imotif_score(seq):
    c_runs=[len(r) for r in re.findall(r"C{3,}",seq)]
    c_fraction=seq.count('C')/len(seq) if seq else 0
    if len(c_runs)<4: return 0
    c_run_spans=[m.span() for m in re.finditer(r"C{3,}",seq)]
    loops=[c_run_spans[i+1][0]-c_run_spans[i][1] for i in range(len(c_run_spans)-1)]
    loop_score=sum(1/(l+1) for l in loops)/max(1,len(loops)) if loops else 0.5
    return min(1.0,sum(c_runs)/16+c_fraction*0.5+loop_score*0.3)

# Motif definitions below (strictly motif-finding and scoring only, no aggregation/reporting):

def find_global_curved_polyA_polyT(seq,min_tract_len=3,min_repeats=3,min_spacing=8,max_spacing=12,min_score=6):
    tracts=find_polyA_polyT_tracts(seq,min_tract_len); results=[]; apr_regions=[]
    for i in range(len(tracts)-min_repeats+1):
        group=[tracts[i]]
        for j in range(1,min_repeats):
            p=((tracts[i+j-1][0]+tracts[i+j-1][1])//2); c=((tracts[i+j][0]+tracts[i+j][1])//2)
            spacing=c-p
            if min_spacing<=spacing<=max_spacing: group.append(tracts[i+j])
            else: break
        if len(group)>=min_repeats:
            mseq=seq[group[0][0]:group[-1][1]+1]; score=curvature_score(mseq)
            if score>=min_score:
                motif={"Class":"Curved_DNA","Subtype":"Global_Curved_Strict_PolyA_or_PolyT",
                       "Start":group[0][0]+1,"End":group[-1][1]+1,"Length":group[-1][1]-group[0][0]+1,
                       "Sequence":mseq,"ScoreMethod":"Strict: PolyA/PolyT curvature tract score","Score":score}
                results.append(motif); apr_regions.append((motif["Start"],motif["End"]))
    return results,apr_regions

def find_local_curved_polyA_polyT(seq,apr_regions,min_len=7):
    results=[]; tracts=find_polyA_polyT_tracts(seq,min_len)
    for s,e,tract_seq in tracts:
        S,E=s+1,e+1
        if not any(r1<=S<=r2 or r1<=E<=r2 for r1,r2 in apr_regions):
            results.append({"Class":"Curved_DNA","Subtype":"Local_Curved_Strict_PolyA_or_PolyT",
                            "Start":S,"End":E,"Length":len(tract_seq),"Sequence":tract_seq,
                            "ScoreMethod":"Strict: PolyA/PolyT tract length","Score":len(tract_seq)})
    return results

def find_curved_DNA(seq):
    global_results,apr_regions=find_global_curved_polyA_polyT(seq)
    local_results=find_local_curved_polyA_polyT(seq,apr_regions)
    return global_results+local_results

def find_zdna(seq,threshold=50,drop_threshold=50,GC_weight=7.0,AT_weight=0.5,GT_weight=1.25,AC_weight=1.25,
              consecutive_AT_scoring=(0.5,0.5,0.5,0.5,0.0,0.0,-5.0,-100.0),mismatch_penalty_type="linear",
              mismatch_penalty_starting_value=3,mismatch_penalty_linear_delta=3,cadence_reward=0.0):
    seq=seq.upper(); motifs=[]
    if len(seq)<12: return []
    scoring=zdna_seeker_scoring_array(seq,GC_weight,AT_weight,GT_weight,AC_weight,consecutive_AT_scoring,
                                      mismatch_penalty_type,mismatch_penalty_starting_value,
                                      mismatch_penalty_linear_delta,cadence_reward)
    start=0; max_here=scoring[0]; curr_max=0; cand=None; end=1
    for i in range(1,len(scoring)):
        num=scoring[i]
        if num>=max_here+num: start=i; end=i+1; max_here=num
        else: max_here+=num; end=i+1
        if max_here>=threshold and (cand is None or curr_max<max_here): cand=(start,end,max_here); curr_max=max_here
        if cand and (max_here<0 or curr_max-max_here>=drop_threshold):
            s,e,score=cand
            motifs.append({"Class":"Z-DNA","Subtype":"Z-Seeker","Start":s+1,"End":e+1,"Length":e-s+1,
                           "Sequence":seq[s:e+1],"ScoreMethod":"Z-Seeker Weighted","Score":f"{score:.2f}"})
            cand=None; max_here=curr_max=0
    if cand:
        s,e,score=cand
        motifs.append({"Class":"Z-DNA","Subtype":"Z-Seeker","Start":s+1,"End":e+1,"Length":e-s+1,
                       "Sequence":seq[s:e+1],"ScoreMethod":"Z-Seeker Weighted","Score":f"{score:.2f}"})
    return motifs

def find_egz_motif(seq):
    pattern=re.compile(r'(CGG){4,}',re.IGNORECASE); results=[]
    for m in pattern.finditer(seq):
        mseq=m.group(0); n_repeats=len(mseq)//3; score=min(1.0,n_repeats/12)
        results.append({"Family":"Double-stranded","Class":"Z-DNA","Subclass":"eGZ (extruded-G)",
                        "Start":m.start()+1,"End":m.end(),"Length":len(mseq),"Sequence":mseq,
                        "ScoreMethod":"Repeat_normalized","Score":f"{score:.2f}","CGG_Repeats":n_repeats})
    return results

def find_slipped_dna(seq):
    results=[]; min_len_dr=10; max_len_dr=300
    for i in range(len(seq)-min_len_dr*2+1):
        for l in range(min_len_dr,min(max_len_dr+1,(len(seq)-i)//2+1)):
            repeat=seq[i:i+l]
            if seq[i+l:i+2*l]==repeat:
                results.append({"Class":"Slipped_DNA","Subtype":"Direct_Repeat","Start":i+1,"End":i+2*l,
                                "Length":2*l,"Sequence":repeat+repeat,"ScoreMethod":"nBST_DR",
                                "Score":f"{min(1.0,l/300):.2f}"})
    min_unit=1; max_unit=6; min_reps=5; min_len_str=15; i=0
    while i<len(seq)-min_unit*min_reps+1:
        found=False
        for unit in range(min_unit,max_unit+1):
            if i+unit*min_reps>len(seq): continue
            repeat_unit=seq[i:i+unit]
            if 'n' in repeat_unit.lower(): continue
            reps=1
            while (i+reps*unit+unit<=len(seq) and seq[i+reps*unit:i+(reps+1)*unit]==repeat_unit): reps+=1
            if reps>=min_reps and reps*unit>=min_len_str:
                remainder=0; rs=i+reps*unit; re_idx=rs
                while (re_idx<len(seq) and seq[re_idx]==repeat_unit[re_idx%unit]):
                    remainder+=1; re_idx+=1
                results.append({"Class":"Slipped_DNA","Subtype":"STR","Start":i+1,"End":i+reps*unit+remainder,
                                "Length":reps*unit+remainder,"Unit":repeat_unit,"Copies":reps,
                                "Sequence":seq[i:i+reps*unit+remainder],"ScoreMethod":"nBST_STR",
                                "Score":f"{min(1.0,reps/20):.2f}"})
                i=i+reps*unit+remainder-1; found=True; break
        if not found: i+=1
    return results

RLFS_MODELS={"m1":r"G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}","m2":r"G{4,}(?:[ATGC]{1,10}?G{4,}){1,}"}

def find_rlfs(seq,models=("m1","m2")):
    if len(seq)<100: return []
    results=[]
    for model_name in models:
        pattern=RLFS_MODELS[model_name]
        for m in re.finditer(pattern,seq,re.IGNORECASE):
            riz_seq=m.group(0)
            if gc_content(riz_seq)<50: continue
            rez=find_rez_max(seq,m.end())
            if rez:
                rez_seq=rez['seq']; stability=min(1.0,0.6*(gc_content(riz_seq+rez_seq)/100)+0.4*(len(re.findall(r"G{3,}",riz_seq+rez_seq))/5))
                results.append({"Class":"R-Loop","Subtype":f"RLFS_{model_name}","Start":m.start()+1,"End":m.start()+len(riz_seq)+rez['end'],
                                "Length":len(riz_seq)+rez['end'],"Sequence":riz_seq+rez_seq,"ScoreMethod":"QmRLFS_Thermo",
                                "Score":f"{stability:.2f}"})
    return results

def find_rez_max(seq,start_pos,max_len=2000,step=100,min_gc=40):
    max_window=""
    for win_start in range(start_pos,min(len(seq),start_pos+max_len),step):
        win_end=min(win_start+step,len(seq))
        window=seq[win_start:win_end]
        if gc_content(window)>=min_gc and len(window)>len(max_window): max_window=window
    return {'seq':max_window,'end':len(max_window)} if max_window else None

def find_cruciform(seq):
    results=[]
    for i in range(len(seq)-2*10):
        for arm_len in range(10,min(101,(len(seq)-i)//2)):
            for spacer_len in range(0,4):
                arm=seq[i:i+arm_len]; rev_arm=reverse_complement(arm)
                mid=i+arm_len+spacer_len
                if mid+arm_len>len(seq): continue
                candidate=seq[mid:mid+arm_len]
                if candidate==rev_arm:
                    full=seq[i:mid+arm_len]
                    score=min(1.0,(arm_len/100)+((arm.count('A')+arm.count('T'))/arm_len*0.3))
                    results.append({"Class":"Cruciform","Subtype":f"Inverted_Repeat_spacer{spacer_len}",
                                    "Start":i+1,"End":mid+arm_len,"Length":len(full),"Sequence":full,
                                    "ScoreMethod":"nBST_IR","Score":f"{score:.2f}"})
    return results

def find_hdna(seq):
    results=[]; n=len(seq)
    for rep_len in range(10,min(101,n//2)):
        for spacer in range(0,9):
            pattern=re.compile(rf"(?=(([ATGC]{{{rep_len}}})[ATGC]{{{spacer}}}\2))",re.IGNORECASE)
            for m in pattern.finditer(seq):
                repeat=m.group(2); mirror_start=m.start(); mirror_end=mirror_start+2*rep_len+spacer
                if mirror_end>n: continue
                full_seq=seq[mirror_start:mirror_end]
                pur_frac=purine_fraction(full_seq); pyr_frac=pyrimidine_fraction(full_seq)
                is_triplex=(pur_frac>=0.9 or pyr_frac>=0.9)
                results.append({"Class":"Triplex_DNA" if is_triplex else "Mirror_Repeat","Subtype":"Triplex_Motif" if is_triplex else "Mirror_Repeat",
                                "Start":mirror_start+1,"End":mirror_end,"Length":len(full_seq),"Spacer":spacer,
                                "Sequence":full_seq,"PurineFrac":round(pur_frac,2),"PyrimidineFrac":round(pyr_frac,2)})
    return results

def find_sticky_dna(seq):
    motifs=[]; seq=seq.replace('\n','').replace(' ','').upper()
    pattern=r"(?:GAA){59,}|(?:TTC){59,}"
    for m in re.finditer(pattern,seq):
        repeat_count=len(m.group())//3
        motifs.append({"Class":"Sticky_DNA","Subtype":"GAA_TTC_Repeat","Start":m.start()+1,"End":m.end(),
                       "Length":len(m.group()),"RepeatCount":repeat_count,"Sequence":m.group(),
                       "ScoreMethod":"Sakamoto1999","Score":f"{min(1.0,repeat_count/270):.2f}"})
    return motifs

def find_multimeric_gquadruplex(seq):
    results=[]; pattern=r"(G{3,}\w{1,12}){4,}"
    for m in overlapping_finditer(pattern,seq):
        mseq=m.group(0)
        if g4hunter_score(mseq)>=1.0:
            results.append({"Class":"G4","Subtype":"Multimeric_G4","Start":m.start()+1,"End":m.end(),
                            "Length":len(mseq),"Sequence":mseq,"ScoreMethod":"G4Hunter_Multimer",
                            "Score":f"{g4hunter_score(mseq)*1.2:.2f}"})
    return results

def find_bipartite_gquadruplex(seq):
    results=[]; pattern=r"(G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{10,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,})"
    for m in overlapping_finditer(pattern,seq):
        mseq=m.group(1)
        if len(re.findall(r"G{3,}",mseq))<8: continue
        half=len(mseq)//2; unit1,unit2=mseq[:half],mseq[half:]
        score=max(g4hunter_score(unit1),g4hunter_score(unit2))*0.9
        if score>=0.9:
            results.append({"Class":"G4","Subtype":"Bipartite_G4","Start":m.start()+1,"End":m.end(),
                            "Length":len(mseq),"Sequence":mseq,"ScoreMethod":"Bipartite_Score",
                            "Score":f"{score:.2f}"})
    return results

def find_gquadruplex(seq):
    pattern=r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        mseq=m.group(1); score=g4hunter_score(mseq)
        if score>=1.2:
            results.append({"Class":"G4","Subtype":"Canonical_G4","Start":m.start()+1,"End":m.end(),
                            "Length":len(mseq),"Sequence":mseq,"ScoreMethod":"G4Hunter_v2",
                            "Score":f"{score:.2f}"})
    return results

def find_relaxed_gquadruplex(seq):
    pattern=r"(G{3,}\w{8,12}G{3,}\w{8,12}G{3,}\w{8,12}G{3,})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        mseq=m.group(1); score=g4hunter_score(mseq)
        if score>=0.8:
            results.append({"Class":"G4","Subtype":"Relaxed_G4","Start":m.start()+1,"End":m.end(),
                            "Length":len(mseq),"Sequence":mseq,"ScoreMethod":"G4Hunter_LongLoop",
                            "Score":f"{score*0.8:.2f}"})
    return results

def find_bulged_gquadruplex(seq):
    pattern=r"(G{3,}\w{0,3}G{3,}\w{0,3}G{3,}\w{0,3}G{3,})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        mseq=m.group(1)
        if len(re.findall(r"G{3,}",mseq))>=4:
            score=g4hunter_score(mseq)
            results.append({"Class":"G4","Subtype":"Bulged_G4","Start":m.start()+1,"End":m.end(),
                            "Length":len(mseq),"Sequence":mseq,"ScoreMethod":"G4Hunter_Bulge",
                            "Score":f"{score*0.7:.2f}"})
    return results

def find_imperfect_gquadruplex(seq):
    pattern=r"(G{2,3}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"\
            r"|(G{3,}\w{1,7}G{2,3}\w{1,7}G{3,}\w{1,7}G{3,})"\
            r"|(G{3,}\w{1,7}G{3,}\w{1,7}G{2,3}\w{1,7}G{3,})"\
            r"|(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{2,3})"
    results=[]
    for m in overlapping_finditer(pattern,seq):
        mseq=m.group(0); score=g4hunter_score(mseq)
        if score>=1.2:
            results.append({"Class":"G4","Subtype":"Imperfect_G4","Start":m.start()+1,"End":m.end(),
                            "Length":len(mseq),"Sequence":mseq,"ScoreMethod":"G4Hunter_Imperfect",
                            "Score":f"{score:.2f}"})
    return results

def find_gtriplex(seq):
    pattern=r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        mseq=m.group(1)
        g_runs=[len(r) for r in re.findall(r"G{3,}",mseq)]
        if len(g_runs)<3: continue
        loops=[len(l) for l in re.findall(r"G{3,}(\w{1,7})G{3,}",mseq)]
        score=min(1.0,sum(g_runs)/15+sum(1/l if l>0 else 0.5 for l in loops)/3)
        results.append({"Class":"G-Triplex","Subtype":"Three_G-Runs","Start":m.start()+1,"End":m.end(),
                        "Length":len(mseq),"Sequence":mseq,"ScoreMethod":"G3_Stability",
                        "Score":f"{score:.2f}"})
    return results

def find_imotif(seq):
    results=[]; pattern=r"(?=(C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}))"
    for m in overlapping_finditer(pattern,seq):
        mseq=m.group(1); score=imotif_score(mseq)
        if score>=0.7:
            c_run_spans=[mt.span() for mt in re.finditer(r"C{3,}",mseq)]
            loops=[c_run_spans[i+1][0]-c_run_spans[i][1] for i in range(len(c_run_spans)-1)]
            if loops and all(1<=l<=7 for l in loops): subtype="Canonical_iMotif"
            elif loops and any(8<=l<=12 for l in loops): subtype="LongLoop_iMotif"
            else: subtype="Other_iMotif"
            results.append({"Class":"i-Motif","Subtype":subtype,"Start":m.start()+1,"End":m.start()+len(mseq),
                            "Length":len(mseq),"Sequence":mseq,"ScoreMethod":"iM_G4HunterStyle",
                            "Score":f"{score:.2f}"})
    return results

def find_ac_motifs(seq):
    pattern=re.compile(
        r"(?=(?:A{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}|"
        r"C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}A{3}))",re.IGNORECASE)
    results=[]
    for m in pattern.finditer(seq):
        mseq=m.group(0).upper()
        results.append({"Class":"AC-Motif","Subtype":"Consensus","Start":m.start()+1,"End":m.start()+len(mseq),
                        "Length":len(mseq),"Sequence":mseq,"ScoreMethod":"PatternMatch","Score":"1.0"})
    return results
