import re; import numpy as np

# [1] FASTA Parser
def parse_fasta(fasta_str): return "".join([l.strip() for l in fasta_str.split('\n') if not l.startswith(">")]).upper().replace(" ", "").replace("U", "T")

# [2] Sequence Wrapper
def wrap(seq, width=60): return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))

# [3] GC Content
def gc_content(seq): return (seq.count('G')+seq.count('C'))/max(1,len(seq))*100 if seq else 0

# [4] Reverse Complement
def reverse_complement(seq): return seq.translate(str.maketrans("ATGC","TACG"))[::-1]

# [5] Overlapping regex finditer
def overlapping_finditer(pattern, seq): regex=re.compile(pattern,re.IGNORECASE); pos=0; 
    while pos<len(seq): m=regex.search(seq,pos); 
        if not m: break; yield m; pos=m.start()+1

# [6] PolyA/PolyT tract finder
def find_polyA_polyT_tracts(seq, min_len=7): results=[]; i=0; n=len(seq)
    while i<n:
        if seq[i]=='A': s=i; while i<n and seq[i]=='A': i+=1; 
            if i-s>=min_len: results.append((s,i-1,seq[s:i]))
        elif seq[i]=='T': s=i; while i<n and seq[i]=='T': i+=1; 
            if i-s>=min_len: results.append((s,i-1,seq[s:i]))
        else: i+=1
    return results

# [7] Curvature score (basic: length)
def curvature_score(seq): return len(seq)/10

# [8] Global curved DNA (strict PolyA/PolyT arrays)
def find_global_curved_polyA_polyT(seq, min_tract_len=3, min_repeats=3, min_spacing=8, max_spacing=12, min_score=1):
    tracts=find_polyA_polyT_tracts(seq,min_tract_len); results=[]; apr_regions=[]
    for i in range(len(tracts)-min_repeats+1):
        group=[tracts[i]]
        for j in range(1,min_repeats):
            pc=(tracts[i+j-1][0]+tracts[i+j-1][1])//2; cc=(tracts[i+j][0]+tracts[i+j][1])//2; sp=cc-pc
            if min_spacing<=sp<=max_spacing: group.append(tracts[i+j]); else: break
        if len(group)>=min_repeats:
            motif_seq=seq[group[0][0]:group[-1][1]+1]
            score=curvature_score(motif_seq)
            if score>=min_score:
                motif={"Class":"Curved_DNA","Subtype":"Global_Curved_Strict_PolyA_or_PolyT","Start":group[0][0]+1,"End":group[-1][1]+1,
                    "Length":group[-1][1]-group[0][0]+1,"Sequence":wrap(motif_seq),
                    "Score":round(score,2),"Arm":motif_seq}
                results.append(motif); apr_regions.append((motif["Start"],motif["End"]))
    return results,apr_regions

# [9] Local curved DNA (strict PolyA/PolyT tracts)
def find_local_curved_polyA_polyT(seq, apr_regions, min_len=7): results=[]; tracts=find_polyA_polyT_tracts(seq,min_len)
    for s,e,tract_seq in tracts:
        st=s+1; en=e+1
        if not any(rs<=st<=re or rs<=en<=re for rs,re in apr_regions):
            results.append({"Class":"Curved_DNA","Subtype":"Local_Curved_Strict_PolyA_or_PolyT",
                "Start":st,"End":en,"Length":len(tract_seq),"Sequence":wrap(tract_seq),"Score":round(curvature_score(tract_seq),2),"Arm":tract_seq})
    return results

# [10] Curved DNA main
def find_curved_DNA(seq): g,apr=find_global_curved_polyA_polyT(seq); l=find_local_curved_polyA_polyT(seq,apr); return g+l

# [11] Slipped DNA (direct repeats & STR)
def find_slipped_dna(seq): results=[]; n=len(seq)
    # Direct repeats
    for i in range(n-20):
        for l in range(10,min(301,(n-i)//2+1)):
            rep=seq[i:i+l]
            if seq[i+l:i+2*l]==rep:
                score=min(1.0,l/300)
                results.append({"Class":"Slipped_DNA","Subtype":"Direct_Repeat","Start":i+1,"End":i+2*l,"Length":2*l,
                    "Sequence":wrap(rep+rep),"Score":score,"Unit":rep,"Copies":2})
    # STR
    min_unit,max_unit,min_reps,min_len=1,6,5,15;i=0
    while i<n-min_unit*min_reps+1:
        found=False
        for u in range(min_unit,max_unit+1):
            if i+u*min_reps>n: continue
            unit=seq[i:i+u]
            if 'n' in unit.lower(): continue
            reps=1
            while i+reps*u+u<=n and seq[i+reps*u:i+(reps+1)*u]==unit: reps+=1
            if reps>=min_reps and reps*u>=min_len:
                rem=0; rs=i+reps*u; re_idx=rs
                while re_idx<n and seq[re_idx]==unit[re_idx%u]: rem+=1; re_idx+=1
                score=min(1.0,reps/20)
                results.append({"Class":"Slipped_DNA","Subtype":"STR","Start":i+1,"End":i+reps*u+rem,"Length":reps*u+rem,
                    "Sequence":wrap(seq[i:i+reps*u+rem]),"Score":score,"Unit":unit,"Copies":reps})
                i=i+reps*u+rem-1; found=True; break
        if not found: i+=1
    return results

# [12] Cruciform DNA (inverted repeats)
def find_cruciform(seq): results=[]
    for i in range(len(seq)-20):
        for arm_len in range(10,min(101,(len(seq)-i)//2)):
            for spacer_len in range(0,4):
                arm=seq[i:i+arm_len]; rev_arm=reverse_complement(arm); mid=i+arm_len+spacer_len
                if mid+arm_len>len(seq): continue
                cand=seq[mid:mid+arm_len]
                if cand==rev_arm:
                    full=seq[i:mid+arm_len]; score=min(1.0,(arm_len/100)+((arm.count('A')+arm.count('T'))/arm_len*0.3))
                    results.append({"Class":"Cruciform","Subtype":f"Inverted_Repeat_spacer{spacer_len}","Start":i+1,"End":mid+arm_len,
                        "Length":len(full),"Sequence":wrap(full),"Score":round(score,2),"Arm":arm,"Spacer":spacer_len})
    return results

# [13] Triplex DNA & Mirror repeats (purine/pyrimidine fraction, consensus scoring)
def purine_fraction(seq): return (seq.count('A')+seq.count('G'))/max(1,len(seq))
def pyrimidine_fraction(seq): return (seq.count('C')+seq.count('T'))/max(1,len(seq))
def find_hdna(seq): results=[]; n=len(seq)
    for rep_len in range(10,min(101,n//2)):
        for spacer in range(0,9):
            pattern=re.compile(rf"(?=(([ATGC]{{{rep_len}}})[ATGC]{{{spacer}}}\2))",re.IGNORECASE)
            for m in pattern.finditer(seq):
                repeat=m.group(2); ms=m.start(); me=ms+2*rep_len+spacer; 
                if me>n: continue
                full_seq=seq[ms:me]; pur=purine_fraction(full_seq); pyr=pyrimidine_fraction(full_seq)
                is_triplex=pur>=0.9 or pyr>=0.9
                score=1.0 if is_triplex else 0.5
                results.append({"Class":"Triplex_DNA" if is_triplex else "Mirror_Repeat",
                    "Subtype":"Triplex_Motif" if is_triplex else "Mirror_Repeat",
                    "Start":ms+1,"End":me,"Length":len(full_seq),"Sequence":wrap(full_seq),
                    "Score":score,"Arm":repeat,"Spacer":spacer})
    return results

# [14] Sticky DNA (GAA/TTC repeats, consensus scoring)
def find_sticky_dna(seq): motifs=[]; seq=seq.replace('\n','').replace(' ','').upper(); pattern=r"(?:GAA){59,}|(?:TTC){59,}"
    for m in re.finditer(pattern,seq):
        rc=len(m.group())//3; motifs.append({"Class":"Sticky_DNA","Subtype":"GAA_TTC_Repeat",
            "Start":m.start()+1,"End":m.end(),"Length":len(m.group()),"Sequence":wrap(m.group()),
            "Score":1.0,"Unit":m.group()[:3],"Copies":rc})
    return motifs

# [15] i-Motif (C-rich, consensus scoring)
def find_imotif(seq): results=[]; pattern=r"(?=(C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}))"
    for m in overlapping_finditer(pattern,seq):
        motif_seq=m.group(1)
        c_runs=[len(r) for r in re.findall(r"C{3,}",motif_seq)]; c_frac=motif_seq.count('C')/len(motif_seq) if motif_seq else 0
        score=1.0 if len(c_runs)>=4 and c_frac>0.6 else 0.5
        c_run_spans=[match.span() for match in re.finditer(r"C{3,}",motif_seq)]
        loops=[c_run_spans[i+1][0]-c_run_spans[i][1] for i in range(len(c_run_spans)-1)]
        subtype="Canonical_iMotif" if loops and all(1<=l<=7 for l in loops) else "Other_iMotif"
        results.append({"Class":"i-Motif","Subtype":subtype,"Start":m.start()+1,"End":m.start()+len(motif_seq),
            "Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Arm":motif_seq})
    return results

# [16] AC motifs (consensus scoring)
def find_ac_motifs(seq): pattern=re.compile(r"(?=(?:A{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}|C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}A{3}))",re.IGNORECASE)
    results=[]; 
    for m in pattern.finditer(seq):
        motif_seq=m.group(0).upper(); results.append({"Class":"AC-Motif","Subtype":"Consensus",
            "Start":m.start()+1,"End":m.start()+len(motif_seq),"Length":len(motif_seq),
            "Sequence":wrap(motif_seq),"Score":1.0,"Arm":motif_seq})
    return results

# [17] Hybrids (overlap of motif classes)
def find_hybrids(motifs,seq): events=[]
    for idx,m in enumerate(motifs):
        events.append((m['Start'],'start',idx)); events.append((m['End']+1,'end',idx))
    events.sort(); active=set(); region_start=None; results=[]
    for pos,typ,idx in events:
        if typ=='start': active.add(idx); 
            if len(active)==2: region_start=pos
        elif typ=='end': 
            if len(active)==2: region_end=pos-1; involved_idxs=list(active); involved_classes={motifs[i]['Class'] for i in involved_idxs}
                if len(involved_classes)>=2:
                    region_motifs=[motifs[i] for i in involved_idxs]
                    results.append({"Class":"Hybrid","Subtype":"_".join(sorted(involved_classes))+"_Overlap",
                        "Start":region_start,"End":region_end,"Length":region_end-region_start+1,
                        "MotifClasses":sorted(involved_classes),
                        "Sequence":seq[region_start-1:region_end],
                        "Score":1.0,"ContributingRegion":region_motifs})
            active.discard(idx)
    return results

# [18] Hotspot finder (clusters of motifs)
def find_hotspots(motif_hits,seq_len,window=100,min_count=3): hotspots=[]; pos=[(hit['Start'],hit['End']) for hit in motif_hits]
    for i in range(0,seq_len-window+1):
        rs, re = i+1, i+window
        count=sum(s<=re and e>=rs for s,e in pos)
        if count>=min_count:
            motifs_in_region=[m for m in motif_hits if m['Start']<=re and m['End']>=rs]
            type_div=len({m['Subtype'] for m in motifs_in_region})
            hotspots.append({"Class":"Non-B DNA Clusters","Subtype":"Hotspot",
                "Start":rs,"End":re,"Length":re-rs+1,"Sequence":"", "Score":1.0,
                "MotifCount":count,"TypeDiversity":type_div})
    return merge_hotspots(hotspots)
def merge_hotspots(hotspots):
    if not hotspots: return []
    merged=[hotspots[0]]
    for cur in hotspots[1:]:
        last=merged[-1]
        if cur['Start']<=last['End']:
            last['End']=max(last['End'],cur['End']); last['Length']=last['End']-last['Start']+1
            last['MotifCount']+=cur['MotifCount']; last['TypeDiversity']=max(last['TypeDiversity'],cur['TypeDiversity'])
            last['Score']=1.0
        else: merged.append(cur)
    return merged

# [19] Validate motif structure
def validate_motif(motif,seq_length): 
    req=["Class","Subtype","Start","End","Length","Sequence"]
    if not all(k in motif for k in req): return False
    if not (1<=motif["Start"]<=motif["End"]<=seq_length): return False
    if len(motif["Sequence"].replace('\n',''))==0: return False
    return True

# [20] Main motif finder and results formatter
def all_motifs(seq, nonoverlap=False, report_hotspots=False, seq_name="Sequence1"):
    if not seq or not re.match("^[ATGC]+$", seq, re.IGNORECASE): return []
    seq=seq.upper()
    motif_list=(find_sticky_dna(seq) + find_curved_DNA(seq) + find_slipped_dna(seq) + find_cruciform(seq) +
        find_hdna(seq) + find_gtriplex(seq) + find_imotif(seq) + find_ac_motifs(seq))
    motif_list=[m for m in motif_list if validate_motif(m,len(seq))]
    motif_list+=find_hybrids(motif_list,seq)
    if nonoverlap: motif_list=select_best_nonoverlapping_motifs(motif_list)
    if report_hotspots: motif_list+=find_hotspots(motif_list,len(seq))
    # [21] Format for download/output
    out=[]
    for m in motif_list:
        out.append({
            "Sequence Name": seq_name,
            "Class": m.get("Class",""),
            "Subtype": m.get("Subtype",""),
            "Start:End": f"{m.get('Start','')}:{m.get('End','')}",
            "Motif Length": m.get("Length",""),
            "Sequence": m.get("Sequence",""),
            "Score": m.get("Score",1.0),
            "Motif Classes": m.get("MotifClasses",m.get("Class","")),
            "Arm/Stem/Contributing region/Unit{Copies}": m.get("Arm",m.get("Unit",""))+(f"{{{m.get('Copies','')}}" if 'Copies' in m else ""),
            "Spacer": m.get("Spacer","")
        })
    return out

# [22] Nonoverlapping motif selection (priority)
def select_best_nonoverlapping_motifs(motifs, motif_priority=None):
    if motif_priority is None: motif_priority=['Sticky_DNA','Cruciform','Curved_DNA','Triplex_DNA','Mirror_Repeat','Slipped_DNA','i-Motif','AC-Motif']
    subtype_rank={s:i for i,s in enumerate(motif_priority)}
    def motif_key(m): rank=subtype_rank.get(m.get('Class'),len(subtype_rank)); sc=float(m.get('Score',1.0)); l=m.get('Length',0); return (rank,-sc,-l)
    sorted_motifs=sorted(motifs, key=motif_key)
    selected=[]; occ=set()
    for m in sorted_motifs:
        reg=set(range(m['Start'],m['End']+1))
        if occ.isdisjoint(reg): selected.append(m); occ.update(reg)
    return selected
