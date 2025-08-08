# =================== Non-B DNA Motif Detection & Annotation ===================
import re; import numpy as np; import random

# --- Sequence Utilities ---
def parse_fasta(fasta_str: str) -> str: return "".join([line.strip() for line in fasta_str.split('\n') if not line.startswith(">")]).upper().replace(" ", "").replace("U", "T")
def wrap(seq: str, width: int = 60) -> str: return "\n".join(seq[i:i+width] for i in range(0, len(seq), width))
def gc_content(seq: str) -> float: gc = seq.count('G') + seq.count('C'); return (gc / max(1, len(seq))) * 100 if seq else 0
def reverse_complement(seq: str) -> str: return seq.translate(str.maketrans("ATGC", "TACG"))[::-1]
def is_palindrome(seq: str) -> bool: return seq == reverse_complement(seq)
def overlapping_finditer(pattern, seq): regex = re.compile(pattern, re.IGNORECASE); pos = 0; 
    while pos < len(seq): m = regex.search(seq, pos); 
        if not m: break; yield m; pos = m.start() + 1

# --- k-mer Conservation Score (Altschul & Erickson 1985; Clamp 2003) ---
def kmer_conservation_score(seq, motif_seq, motif_start, motif_end, k=8, n_shuffles=1000):
    if len(motif_seq) < k: k = max(4, len(motif_seq)//2)
    motif_kmer = motif_seq[:k]; observed_count = seq.count(motif_kmer)
    if observed_count == 0: return {'conservation_score':0.0,'p_value':1.0,'significance':'not significant'}
    dinucs = [seq[i:i+2] for i in range(len(seq)-1)]; null_counts = []
    for _ in range(n_shuffles):
        shuffled_dinucs = dinucs.copy(); random.shuffle(shuffled_dinucs)
        shuffled_seq = shuffled_dinucs[0][0] + ''.join(d[1] for d in shuffled_dinucs) if shuffled_dinucs else ""
        null_counts.append(shuffled_seq.count(motif_kmer))
    expected = np.mean(null_counts); conservation_score = np.log2(observed_count/(expected+1e-10))
    p_value = sum(1 for x in null_counts if x >= observed_count)/n_shuffles
    if p_value<0.001 and conservation_score>2: significance='very high'
    elif p_value<0.01 and conservation_score>1.5: significance='high'
    elif p_value<0.05 and conservation_score>1: significance='medium'
    elif p_value<0.1 and conservation_score>0.5: significance='low'
    else: significance='not significant'
    return {'conservation_score':round(conservation_score,3),'p_value':round(p_value,4),'significance':significance}

# --- Curved DNA Motif Detection ---
def find_polyA_polyT_tracts(seq: str, min_len: int = 7) -> list:
    results=[]; i=0; n=len(seq)
    while i<n:
        if seq[i]=='A': start=i; while i<n and seq[i]=='A': i+=1; 
            if i-start>=min_len: results.append((start,i-1,seq[start:i]))
        elif seq[i]=='T': start=i; while i<n and seq[i]=='T': i+=1; 
            if i-start>=min_len: results.append((start,i-1,seq[start:i]))
        else: i+=1
    return results

def curved_dna_scoring(motif_seq, a_tract_positions):
    base_score=1.0; length_score=min(3.0,np.mean([len(pos[2]) for pos in a_tract_positions if len(pos)>2])/5.0) if a_tract_positions else len(motif_seq)/20.0
    phasing_score=1.0
    if len(a_tract_positions)>1:
        spacings=[(a_tract_positions[i+1][0]+a_tract_positions[i+1][1])/2-(a_tract_positions[i][0]+a_tract_positions[i][1])/2 for i in range(len(a_tract_positions)-1)]
        if spacings: phasing_deviation=np.mean([abs(s-10.5) for s in spacings]); phasing_score=max(0.5,2.0-phasing_deviation/5.0)
    at_content=(motif_seq.count('A')+motif_seq.count('T'))/len(motif_seq); at_bonus=1.0+(at_content-0.5)*0.5
    final_score=base_score+length_score*phasing_score*at_bonus
    return max(1.0,round(final_score,2))

def find_global_curved_polyA_polyT(seq: str, sequence_name='Unknown', min_tract_len=3, min_repeats=3, min_spacing=8, max_spacing=12, min_score=2) -> tuple:
    tracts=find_polyA_polyT_tracts(seq,min_tract_len); results=[]; apr_regions=[]
    for i in range(len(tracts)-min_repeats+1):
        group=[tracts[i]]
        for j in range(1,min_repeats):
            prev_center=(tracts[i+j-1][0]+tracts[i+j-1][1])//2; curr_center=(tracts[i+j][0]+tracts[i+j][1])//2
            spacing=curr_center-prev_center
            if min_spacing<=spacing<=max_spacing: group.append(tracts[i+j])
            else: break
        if len(group)>=min_repeats:
            motif_seq=seq[group[0][0]:group[-1][1]+1]; score=curved_dna_scoring(motif_seq,group)
            conservation=kmer_conservation_score(seq,motif_seq,group[0][0],group[-1][1]+1)
            significance='very high' if score>=4 else 'high' if score>=3 else 'medium' if score>=2 else 'low'
            motif={"SequenceName":sequence_name,"Class":"Curved_DNA","Subtype":"Global_Curved_Strict_PolyA_or_PolyT","Start":group[0][0]+1,"End":group[-1][1]+1,"Length":group[-1][1]-group[0][0]+1,"Sequence":wrap(motif_seq),"Score":score,"Significance":significance,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":"; ".join([tract[2] for tract in group]),"Spacer_Linker":"8-12bp"}
            results.append(motif); apr_regions.append((motif["Start"],motif["End"]))
    return results, apr_regions

def find_local_curved_polyA_polyT(seq: str, apr_regions: list, sequence_name='Unknown', min_len: int = 7) -> list:
    results=[]; tracts=find_polyA_polyT_tracts(seq,min_len)
    for start,end,tract_seq in tracts:
        s,e=start+1,end+1
        if not any(r_start<=s<=r_end or r_start<=e<=r_end for r_start,r_end in apr_regions):
            score=curved_dna_scoring(tract_seq,[(start,end,tract_seq)])
            conservation=kmer_conservation_score(seq,tract_seq,start,end+1)
            significance='high' if score>=3 else 'medium' if score>=2 else 'low'
            results.append({"SequenceName":sequence_name,"Class":"Curved_DNA","Subtype":"Local_Curved_Strict_PolyA_or_PolyT","Start":s,"End":e,"Length":len(tract_seq),"Sequence":wrap(tract_seq),"Score":score,"Significance":significance,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":tract_seq,"Spacer_Linker":"None"})
    return results

def find_curved_DNA(seq: str, sequence_name='Unknown') -> list:
    global_results,apr_regions=find_global_curved_polyA_polyT(seq,sequence_name)
    local_results=find_local_curved_polyA_polyT(seq,apr_regions,sequence_name)
    return global_results+local_results

# --- Z-DNA Motif Detection ---
def zdna_seeker_scoring_array(seq, GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25, consecutive_AT_scoring=(0.5,0.5,0.5,0.5,0.0,0.0,-5.0,-100.0), mismatch_penalty_type="linear", mismatch_penalty_starting_value=3, mismatch_penalty_linear_delta=3, cadence_reward=0.0):
    scoring_array=np.empty(len(seq)-1,dtype=float); mismatches_counter=0; consecutive_AT_counter=0
    for i in range(len(seq)-1):
        t=seq[i:i+2].upper()
        if t in ("GC","CG"): scoring_array[i]=GC_weight; mismatches_counter=0; consecutive_AT_counter=0
        elif t in ("GT","TG"): scoring_array[i]=GT_weight; mismatches_counter=0; consecutive_AT_counter=0
        elif t in ("AC","CA"): scoring_array[i]=AC_weight; mismatches_counter=0; consecutive_AT_counter=0
        elif t in ("AT","TA"):
            adjusted_weight=AT_weight+(consecutive_AT_scoring[consecutive_AT_counter] if consecutive_AT_counter<len(consecutive_AT_scoring) else consecutive_AT_scoring[-1])
            scoring_array[i]=adjusted_weight; consecutive_AT_counter+=1; mismatches_counter=0
        else:
            mismatches_counter+=1; consecutive_AT_counter=0
            if mismatch_penalty_type=="exponential": scoring_array[i]=-mismatch_penalty_starting_value**mismatches_counter if mismatches_counter<15 else -32000
            elif mismatch_penalty_type=="linear": scoring_array[i]=-mismatch_penalty_starting_value-mismatch_penalty_linear_delta*(mismatches_counter-1)
            else: scoring_array[i]=-10
        if t in ("GC","CG","GT","TG","AC","CA","AT","TA"): scoring_array[i]+=cadence_reward
    return scoring_array

def find_zdna(seq, sequence_name='Unknown', threshold=50, drop_threshold=50, GC_weight=7.0, AT_weight=0.5, GT_weight=1.25, AC_weight=1.25, consecutive_AT_scoring=(0.5,0.5,0.5,0.5,0.0,0.0,-5.0,-100.0), mismatch_penalty_type="linear", mismatch_penalty_starting_value=3, mismatch_penalty_linear_delta=3, cadence_reward=0.0):
    seq=seq.upper(); motifs=[]; scoring=zdna_seeker_scoring_array(seq,GC_weight,AT_weight,GT_weight,AC_weight,consecutive_AT_scoring,mismatch_penalty_type,mismatch_penalty_starting_value,mismatch_penalty_linear_delta,cadence_reward)
    start_idx=0; max_ending_here=scoring[0]; current_max=0; candidate=None; end_idx=1
    for i in range(1,len(scoring)):
        num=scoring[i]
        if num>=max_ending_here+num: start_idx=i; end_idx=i+1; max_ending_here=num
        else: max_ending_here+=num; end_idx=i+1
        if max_ending_here>=threshold and (candidate is None or current_max<max_ending_here):
            candidate=(start_idx,end_idx,max_ending_here); current_max=max_ending_here
        if candidate and (max_ending_here<0 or current_max-max_ending_here>=drop_threshold):
            s,e,score=candidate; motif_seq=seq[s:e+1]; conservation=kmer_conservation_score(seq,motif_seq,s,e+1)
            significance='very high' if score>=150 else 'high' if score>=100 else 'medium' if score>=75 else 'low'
            motifs.append({"SequenceName":sequence_name,"Class":"Z-DNA","Subtype":"Z-Seeker","Start":s+1,"End":e+1,"Length":e-s+1,"Sequence":wrap(motif_seq),"Score":round(score,2),"Significance":significance,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":"None"}); candidate=None; max_ending_here=current_max=0
    if candidate:
        s,e,score=candidate; motif_seq=seq[s:e+1]; conservation=kmer_conservation_score(seq,motif_seq,s,e+1)
        significance='very high' if score>=150 else 'high' if score>=100 else 'medium' if score>=75 else 'low'
        motifs.append({"SequenceName":sequence_name,"Class":"Z-DNA","Subtype":"Z-Seeker","Start":s+1,"End":e+1,"Length":e-s+1,"Sequence":wrap(motif_seq),"Score":round(score,2),"Significance":significance,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":"None"})
    return motifs

# --- eGZ DNA Motif Detection ---
def find_egz_motif(seq, sequence_name='Unknown'):
    pattern=re.compile(r'(CGG){4,}',re.IGNORECASE); results=[]
    for m in pattern.finditer(seq):
        motif_seq=m.group(0); n_repeats=len(motif_seq)//3
        score=1.0 if n_repeats>=100 else 0.8 if n_repeats>=55 else 0.5 if n_repeats>=30 else 0.3
        conservation=kmer_conservation_score(seq,motif_seq,m.start(),m.end())
        sig='very high' if score==1.0 else 'high' if score==0.8 else 'medium' if score==0.5 else 'low'
        results.append({"SequenceName":sequence_name,"Class":"Z-DNA","Subtype":"eGZ_extruded_G","Start":m.start()+1,"End":m.end(),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":sig,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":f"CGG x {n_repeats}","Spacer_Linker":"None"})
    return results

# --- Slipped DNA Motif Detection ---
def slipped_dna_scoring(repeat_length, copy_number, motif_type='direct'):
    base_score=1.0
    if motif_type=='STR':
        if copy_number>=200: instability_score=5.0
        elif copy_number>=100: instability_score=4.0
        elif copy_number>=50: instability_score=3.0
        elif copy_number>=20: instability_score=2.0
        else: instability_score=1.0+(copy_number/20.0)
    else:
        length_factor=min(3.0,repeat_length/50.0); instability_score=1.0+length_factor
    final_score=base_score*instability_score
    return max(1.0,round(final_score,2))
def find_slipped_dna(seq, sequence_name='Unknown'):
    results=[]; min_len_dr=10; max_len_dr=300
    for i in range(len(seq)-min_len_dr*2+1):
        for l in range(min_len_dr,min(max_len_dr+1,(len(seq)-i)//2+1)):
            repeat=seq[i:i+l]
            if seq[i+l:i+2*l]==repeat:
                motif_seq=repeat+repeat
                score=slipped_dna_scoring(l,2,'direct')
                conservation=kmer_conservation_score(seq,motif_seq,i,i+2*l)
                sig='very high' if l>=100 else 'high' if l>=50 else 'medium' if l>=25 else 'low'
                results.append({"SequenceName":sequence_name,"Class":"Slipped_DNA","Subtype":"Direct_Repeat","Start":i+1,"End":i+2*l,"Length":2*l,"Sequence":wrap(motif_seq),"Score":score,"Significance":sig,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":repeat,"Spacer_Linker":"None"})
    min_unit_str=1; max_unit_str=6; min_reps_str=5; min_len_str=15; i=0
    while i<len(seq)-min_unit_str*min_reps_str+1:
        found=False
        for unit in range(min_unit_str,max_unit_str+1):
            if i+unit*min_reps_str>len(seq): continue
            repeat_unit=seq[i:i+unit]
            if 'n' in repeat_unit.lower(): continue
            reps=1
            while(i+reps*unit+unit<=len(seq) and seq[i+reps*unit:i+(reps+1)*unit]==repeat_unit): reps+=1
            if reps>=min_reps_str and reps*unit>=min_len_str:
                remainder=0; rs=i+reps*unit; re_idx=rs
                while(re_idx<len(seq) and seq[re_idx]==repeat_unit[re_idx%unit]): remainder+=1; re_idx+=1
                motif_seq=seq[i:i+reps*unit+remainder]
                score=slipped_dna_scoring(unit,reps,'STR')
                conservation=kmer_conservation_score(seq,motif_seq,i,i+reps*unit+remainder)
                sig='very high' if reps>=100 else 'high' if reps>=50 else 'medium' if reps>=20 else 'low'
                results.append({"SequenceName":sequence_name,"Class":"Slipped_DNA","Subtype":"STR","Start":i+1,"End":i+reps*unit+remainder,"Length":reps*unit+remainder,"Sequence":wrap(motif_seq),"Score":score,"Significance":sig,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":f"{repeat_unit} x {reps}","Spacer_Linker":f"{remainder}bp remainder" if remainder>0 else "None"}); i=i+reps*unit+remainder-1; found=True; break
        if not found: i+=1
    return results

# --- Cruciform DNA Motif Detection ---
def cruciform_scoring(arm_seq,spacer_length,total_length):
    base_score=1.0; arm_length=len(arm_seq); length_score=min(3.0,arm_length/15.0)
    at_score=1.0+(arm_seq.count('A')+arm_seq.count('T'))/len(arm_seq)*0.8
    spacer_penalty=max(0.7,1.0-spacer_length*0.1)
    gc_pairs=arm_seq.count('G')+arm_seq.count('C'); thermo_bonus=1.0+(gc_pairs/len(arm_seq))*0.3
    final_score=base_score+length_score*at_score*spacer_penalty*thermo_bonus
    return max(1.0,round(final_score,2))
def find_cruciform(seq, sequence_name='Unknown'):
    results=[]
    for i in range(len(seq)-2*10):
        for arm_len in range(10,min(101,(len(seq)-i)//2)):
            for spacer_len in range(0,4):
                arm=seq[i:i+arm_len]; rev_arm=reverse_complement(arm)
                mid=i+arm_len+spacer_len
                if mid+arm_len>len(seq): continue
                candidate=seq[mid:mid+arm_len]
                if candidate==rev_arm:
                    full_seq=seq[i:mid+arm_len]
                    score=cruciform_scoring(arm,spacer_len,len(full_seq))
                    conservation=kmer_conservation_score(seq,full_seq,i,mid+arm_len)
                    sig='very high' if score>=4 else 'high' if score>=3 else 'medium' if score>=2 else 'low'
                    results.append({"SequenceName":sequence_name,"Class":"Cruciform_DNA","Subtype":"Palindrome_Arm","Start":i+1,"End":mid+arm_len,"Length":len(full_seq),"Sequence":wrap(full_seq),"Score":score,"Significance":sig,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":f"{arm}|{rev_arm}","Spacer_Linker":f"{spacer_len}bp"})
    return results

# --- Triplex DNA/Mirror Repeat Motif Detection ---
def triplex_scoring(motif_seq,purine_frac,pyrimidine_frac,repeat_length):
    base_score=1.0; homog=max(purine_frac,pyrimidine_frac)
    homog_score=3.0 if homog>=0.95 else 2.5 if homog>=0.90 else 2.0 if homog>=0.80 else 1.0+(homog-0.5)*2
    length_score=min(2.5,repeat_length/20.0); ph_factor=1.2 if pyrimidine_frac>purine_frac else 1.0
    final_score=base_score+homog_score*length_score*ph_factor
    return max(1.0,round(final_score,2))
def find_hdna(seq, sequence_name='Unknown'):
    results=[]; n=len(seq)
    for rep_len in range(10,min(101,n//2)):
        for spacer in range(0,9):
            pattern=re.compile(rf"(?=(([ATGC]{{{rep_len}}})[ATGC]{{{spacer}}}\2))",re.IGNORECASE)
            for m in pattern.finditer(seq):
                repeat=m.group(2); mirror_start=m.start(); mirror_end=mirror_start+2*rep_len+spacer
                if mirror_end>n: continue
                full_seq=seq[mirror_start:mirror_end]
                pur_frac=(full_seq.count('A')+full_seq.count('G'))/max(1,len(full_seq))
                pyr_frac=(full_seq.count('C')+full_seq.count('T'))/max(1,len(full_seq))
                is_triplex=(pur_frac>=0.9 or pyr_frac>=0.9)
                score=triplex_scoring(full_seq,pur_frac,pyr_frac,len(full_seq)) if is_triplex else 1.0
                conservation=kmer_conservation_score(seq,full_seq,mirror_start,mirror_end)
                sig='high' if score>=3 else 'medium' if score>=2 else 'low'
                results.append({"SequenceName":sequence_name,"Class":"Triplex_DNA" if is_triplex else "Mirror_Repeat","Subtype":"Triplex_Motif" if is_triplex else "Mirror_Repeat","Start":mirror_start+1,"End":mirror_end,"Length":len(full_seq),"Sequence":wrap(full_seq),"Score":score,"Significance":sig,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":repeat,"Spacer_Linker":spacer})
    return results

# --- Sticky DNA Motif Detection ---
def sticky_dna_scoring(repeat_count,motif_type='GAA'):
    base_score=1.0
    if repeat_count>=1000: pathogenicity_score=6.0
    elif repeat_count>=700: pathogenicity_score=5.0
    elif repeat_count>=400: pathogenicity_score=4.0
    elif repeat_count>=200: pathogenicity_score=3.5
    elif repeat_count>=100: pathogenicity_score=3.0
    elif repeat_count>=59: pathogenicity_score=2.0
    else: pathogenicity_score=1.0+(repeat_count/59.0)
    final_score=base_score+pathogenicity_score
    return max(1.0,round(final_score,2))
def find_sticky_dna(seq, sequence_name='Unknown'):
    motifs=[]; pattern=r"(?:GAA){59,}|(?:TTC){59,}"
    for m in re.finditer(pattern,seq):
        motif_seq=m.group(); repeat_count=len(motif_seq)//3
        score=sticky_dna_scoring(repeat_count)
        conservation=kmer_conservation_score(seq,motif_seq,m.start(),m.end())
        sig='very high' if score>=6 else 'high' if score>=5 else 'medium' if score>=3 else 'low'
        motifs.append({"SequenceName":sequence_name,"Class":"Sticky_DNA","Subtype":"GAA_TTC_Repeat","Start":m.start()+1,"End":m.end(),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":sig,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":f"{motif_seq[:3]} x {repeat_count}","Spacer_Linker":"None"})
    return motifs

# --- R-loop Motif Detection (Ginno et al 2012, QmRLFS scoring) ---
RLFS_MODELS={"m1":r"G{3,}[ATGC]{1,10}?G{3,}(?:[ATGC]{1,10}?G{3,}){1,}","m2":r"G{4,}(?:[ATGC]{1,10}?G{4,}){1,}"}
def find_rlfs(seq, sequence_name='Unknown', models=("m1","m2")):
    if len(seq)<100: return []
    results=[]
    for model_name in models:
        pattern=RLFS_MODELS[model_name]
        for m in re.finditer(pattern,seq,re.IGNORECASE):
            riz_seq=m.group(0)
            if gc_content(riz_seq)<50: continue
            rez=find_rez_max(seq,m.end())
            if rez:
                rez_seq=rez['seq']
                stability=min(1.0,0.6*(gc_content(riz_seq+rez_seq)/100)+0.4*(len(re.findall(r"G{3,}",riz_seq+rez_seq))/5))
                conservation=kmer_conservation_score(seq,riz_seq+rez_seq,m.start(),m.start()+len(riz_seq)+rez['end'])
                sig='high' if stability>=0.8 else 'medium' if stability>=0.6 else 'low'
                results.append({"SequenceName":sequence_name,"Class":"R-Loop","Subtype":f"RLFS_{model_name}","Start":m.start()+1,"End":m.start()+len(riz_seq)+rez['end'],"Length":len(riz_seq)+rez['end'],"Sequence":wrap(riz_seq+rez_seq),"Score":round(stability,2),"Significance":sig,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":riz_seq,"Spacer_Linker":f"{len(rez_seq)}bp downstream"})
    return results
def find_rez_max(seq,start_pos,max_len=2000,step=100,min_gc=40):
    max_window=""
    for win_start in range(start_pos,min(len(seq),start_pos+max_len),step):
        win_end=min(win_start+step,len(seq))
        window=seq[win_start:win_end]
        if gc_content(window)>=min_gc and len(window)>len(max_window): max_window=window
    if max_window: return {'seq':max_window,'end':len(max_window)}
    return None

# --- G-triplex Motif Detection ---
def gtriplex_scoring(g_runs,loop_lengths,total_length):
    base_score=1.0; g_score=sum(min(4,run_len) for run_len in g_runs)*0.3
    loop_penalty=sum(max(0,l-3)*0.2 for l in loop_lengths); length_bonus=min(1.5,total_length/40.0)
    final_score=base_score+g_score+length_bonus-loop_penalty
    return max(1.0,round(final_score,2))
def find_gtriplex(seq,sequence_name='Unknown'):
    pattern=r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        motif_seq=m.group(1)
        g_runs=[len(r) for r in re.findall(r"G{3,}",motif_seq)]
        loops=[len(l) for l in re.findall(r"G{3,}(\w{1,7})G{3,}",motif_seq)]
        score=gtriplex_scoring(g_runs,loops,len(motif_seq))
        conservation=kmer_conservation_score(seq,motif_seq,m.start(),m.end())
        sig='high' if score>=3 else 'medium' if score>=2 else 'low'
        results.append({"SequenceName":sequence_name,"Class":"G-Triplex","Subtype":"Three_G-Runs","Start":m.start()+1,"End":m.end(),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":sig,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":f"G runs: {g_runs}","Spacer_Linker":f"Loops: {loops}"})
    return results

# --- G-Quadruplex Motif Detection (Canonical, Relaxed, Bulged, Imperfect, Bipartite, Multimeric) ---
def g4hunter_score(seq):
    scores=[]; seq=seq.upper()
    for c in seq:
        if c=='G': scores.append(1)
        elif c=='C': scores.append(-1)
        else: scores.append(0)
    return np.mean(scores) if scores else 0

def find_gquadruplex(seq,sequence_name='Unknown'):
    pattern=r"(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        motif_seq=m.group(1); score=g4hunter_score(motif_seq)
        if score>=1.2:
            results.append({"SequenceName":sequence_name,"Class":"G4","Subtype":"Canonical_G4","Start":m.start()+1,"End":m.end(),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":"high" if score>=1.5 else "medium" if score>=1.2 else "low","Conservation_Score":None,"Conservation_P_Value":None,"Conservation_Significance":None,"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":"None"})
    return results

def find_relaxed_gquadruplex(seq,sequence_name='Unknown'):
    pattern=r"(G{3,}\w{8,12}G{3,}\w{8,12}G{3,}\w{8,12}G{3,})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        motif_seq=m.group(1); score=g4hunter_score(motif_seq)
        if score>=0.8:
            results.append({"SequenceName":sequence_name,"Class":"G4","Subtype":"Relaxed_G4","Start":m.start()+1,"End":m.end(),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":"high" if score>=1.2 else "medium" if score>=0.8 else "low","Conservation_Score":None,"Conservation_P_Value":None,"Conservation_Significance":None,"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":"LongLoop"})
    return results

def find_bulged_gquadruplex(seq,sequence_name='Unknown'):
    pattern=r"(G{3,}\w{0,3}G{3,}\w{0,3}G{3,}\w{0,3}G{3,})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        motif_seq=m.group(1)
        if len(re.findall(r"G{3,}",motif_seq))>=4:
            score=g4hunter_score(motif_seq)
            results.append({"SequenceName":sequence_name,"Class":"G4","Subtype":"Bulged_G4","Start":m.start()+1,"End":m.end(),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":"high" if score>=1.2 else "medium" if score>=1.0 else "low","Conservation_Score":None,"Conservation_P_Value":None,"Conservation_Significance":None,"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":"Bulge"})
    return results

def find_imperfect_gquadruplex(seq,sequence_name='Unknown'):
    pattern=r"(G{2,3}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{3,})"\
              r"|(G{3,}\w{1,7}G{2,3}\w{1,7}G{3,}\w{1,7}G{3,})"\
              r"|(G{3,}\w{1,7}G{3,}\w{1,7}G{2,3}\w{1,7}G{3,})"\
              r"|(G{3,}\w{1,7}G{3,}\w{1,7}G{3,}\w{1,7}G{2,3})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        motif_seq=m.group(0); score=g4hunter_score(motif_seq)
        if score>=1.2:
            results.append({"SequenceName":sequence_name,"Class":"G4","Subtype":"Imperfect_G4","Start":m.start()+1,"End":m.end(),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":"high" if score>=1.5 else "medium" if score>=1.2 else "low","Conservation_Score":None,"Conservation_P_Value":None,"Conservation_Significance":None,"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":"Imperfect"})
    return results

def find_bipartite_gquadruplex(seq,sequence_name='Unknown'):
    pattern=r"(G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{10,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,}\w{1,30}G{3,})"; results=[]
    for m in overlapping_finditer(pattern,seq):
        motif_seq=m.group(1)
        if len(re.findall(r"G{3,}",motif_seq))<8: continue
        half=len(motif_seq)//2; unit1,unit2=motif_seq[:half],motif_seq[half:]
        score=max(g4hunter_score(unit1),g4hunter_score(unit2))*0.9
        if score>=0.9:
            results.append({"SequenceName":sequence_name,"Class":"G4","Subtype":"Bipartite_G4","Start":m.start()+1,"End":m.end(),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":"high" if score>=1.2 else "medium" if score>=0.9 else "low","Conservation_Score":None,"Conservation_P_Value":None,"Conservation_Significance":None,"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":"Bipartite"})
    return results

def find_multimeric_gquadruplex(seq,sequence_name='Unknown'):
    pattern=r"(G{3,}\w{1,12}){4,}"; results=[]
    for m in overlapping_finditer(pattern,seq):
        motif_seq=m.group(0); score=g4hunter_score(motif_seq)
        if score>=1.0:
            results.append({"SequenceName":sequence_name,"Class":"G4","Subtype":"Multimeric_G4","Start":m.start()+1,"End":m.end(),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score*1.2,"Significance":"high" if score>=1.2 else "medium" if score>=1.0 else "low","Conservation_Score":None,"Conservation_P_Value":None,"Conservation_Significance":None,"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":"Multimeric"})
    return results

# --- i-Motif Motif Detection ---
def imotif_score(seq):
    c_runs=[len(r) for r in re.findall(r"C{3,}",seq)]
    c_fraction=seq.count('C')/len(seq) if seq else 0
    if len(c_runs)<4: return 0
    c_run_spans=[match.span() for match in re.finditer(r"C{3,}",seq)]
    loops=[c_run_spans[i+1][0]-c_run_spans[i][1] for i in range(len(c_run_spans)-1)]
    loop_score=sum(1/(l+1) for l in loops)/max(1,len(loops)) if loops else 0.5
    return min(1.0,sum(c_runs)/16+c_fraction*0.5+loop_score*0.3)
def find_imotif(seq,sequence_name='Unknown'):
    pattern=r"(?=(C{3,}\w{1,12}C{3,}\w{1,12}C{3,}\w{1,12}C{3,}))"; results=[]
    for m in overlapping_finditer(pattern,seq):
        motif_seq=m.group(1); score=imotif_score(motif_seq)
        if score>=0.7:
            c_run_spans=[match.span() for match in re.finditer(r"C{3,}",motif_seq)]
            loops=[c_run_spans[i+1][0]-c_run_spans[i][1] for i in range(len(c_run_spans)-1)]
            subtype="Canonical_iMotif" if all(1<=l<=7 for l in loops) else "LongLoop_iMotif" if any(8<=l<=12 for l in loops) else "Other_iMotif"
            results.append({"SequenceName":sequence_name,"Class":"i-Motif","Subtype":subtype,"Start":m.start()+1,"End":m.start()+len(motif_seq),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":"high" if score>=0.9 else "medium" if score>=0.7 else "low","Conservation_Score":None,"Conservation_P_Value":None,"Conservation_Significance":None,"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":f"Loops:{loops}"})
    return results

# --- AC-motif Detection ---
def ac_motif_scoring(motif_seq):
    base_score=1.0; length_score=min(2.0,len(motif_seq)/30.0)
    a_clusters=len(re.findall(r'A{3,}',motif_seq)); c_clusters=len(re.findall(r'C{3,}',motif_seq))
    pattern_score=min(2.0,(a_clusters+c_clusters)/4.0)
    ac_content=(motif_seq.count('A')+motif_seq.count('C'))/len(motif_seq)
    balance_score=1.0+(1.0-abs(0.5-ac_content))*0.5
    final_score=base_score+length_score+pattern_score*balance_score
    return max(1.0,round(final_score,2))
def find_ac_motifs(seq,sequence_name='Unknown'):
    pattern=re.compile(r"(?=(?:A{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}|C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}C{3}[ACGT]{4,6}A{3}))",re.IGNORECASE)
    results=[]
    for m in pattern.finditer(seq):
        motif_seq=m.group(0).upper()
        score=ac_motif_scoring(motif_seq)
        conservation=kmer_conservation_score(seq,motif_seq,m.start(),m.start()+len(motif_seq))
        sig='high' if score>=3 else 'medium' if score>=2 else 'low'
        results.append({"SequenceName":sequence_name,"Class":"AC-Motif","Subtype":"Consensus","Start":m.start()+1,"End":m.start()+len(motif_seq),"Length":len(motif_seq),"Sequence":wrap(motif_seq),"Score":score,"Significance":sig,"Conservation_Score":conservation['conservation_score'],"Conservation_P_Value":conservation['p_value'],"Conservation_Significance":conservation['significance'],"Repeat_Arm_Contributing_Region":motif_seq,"Spacer_Linker":"Pattern"})
    return results

# --- Hybrid Motif Detection (Overlap regions) ---
def hybrid_scoring(involved_classes,region_length,motif_count):
    # Score hybrid regions based on class diversity, density and length
    base_score=1.0; diversity_score=min(3.0,len(involved_classes)*0.8)
    density_score=min(2.0,motif_count/3.0); length_factor=min(1.5,region_length/50.0)
    final_score=base_score+diversity_score+density_score*length_factor
    return max(1.0,round(final_score,2))

def find_hybrids(motifs,seq,sequence_name='Unknown'):
    # Find hybrid/overlap regions among motifs and score them
    events=[]; results=[]
    for idx,m in enumerate(motifs):
        events.append((m['Start'],'start',idx)); events.append((m['End']+1,'end',idx))
    events.sort(); active=set(); region_start=None
    for pos,typ,idx in events:
        if typ=='start':
            active.add(idx)
            if len(active)==2: region_start=pos
        elif typ=='end':
            if len(active)==2 and region_start is not None:
                region_end=pos-1
                region_motifs=[motifs[i] for i in active]
                involved_classes=list(set([m['Class'] for m in region_motifs]))
                motif_count=len(region_motifs)
                region_length=region_end-region_start+1
                score=hybrid_scoring(involved_classes,region_length,motif_count)
                results.append({'Sequence':sequence_name,'Start':region_start,'End':region_end,
                                'Classes':involved_classes,'Score':score,'Motifs':region_motifs})
                region_start=None
            active.discard(idx)
    return results

# --- [Other motif detection blocks as in your code above...] ---
# [All motif finders: Z-DNA, eGZ, Slipped DNA, Cruciform, Triplex/Mirror, Sticky, RLFS, G-triplex, G-quadruplex (all), i-Motif, AC-motif]
# [Already included above in your message; not repeated here for brevity.]

# --- Master Finder ---
def find_all_nb_motifs(seq, sequence_name='Unknown'):
    # Aggregate all motif types into a single list
    motifs=[]
    motifs+=find_curved_DNA(seq,sequence_name)
    motifs+=find_zdna(seq,sequence_name)
    motifs+=find_egz_motif(seq,sequence_name)
    motifs+=find_slipped_dna(seq,sequence_name)
    motifs+=find_cruciform(seq,sequence_name)
    motifs+=find_hdna(seq,sequence_name)
    motifs+=find_sticky_dna(seq,sequence_name)
    motifs+=find_rlfs(seq,sequence_name)
    motifs+=find_gtriplex(seq,sequence_name)
    motifs+=find_gquadruplex(seq,sequence_name)
    motifs+=find_relaxed_gquadruplex(seq,sequence_name)
    motifs+=find_bulged_gquadruplex(seq,sequence_name)
    motifs+=find_imperfect_gquadruplex(seq,sequence_name)
    motifs+=find_bipartite_gquadruplex(seq,sequence_name)
    motifs+=find_multimeric_gquadruplex(seq,sequence_name)
    motifs+=find_imotif(seq,sequence_name)
    motifs+=find_ac_motifs(seq,sequence_name)
    return motifs

# --- Example Usage ---
# fasta_str = ">seq1\nATGC...."
# seq = parse_fasta(fasta_str)
# motifs = find_all_nb_motifs(seq,"seq1")
# hybrids = find_hybrids(motifs,seq,"seq1")

# --- End of Pipeline ---
