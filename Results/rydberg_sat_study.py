"""Exact small-system stabilizer SAT/SMT + probability-weighted decoder study.

No need for a pip-installed SAT package. A native libz3 (where available) runs
an actual Boolean SMT-LIB satisfiability formulation; enumeration independently
certifies *all* [[4,1]] stabilizer groups and their exact Pauli ML success.

Physical data: 2026-11-13 revised Evered-inspired Rydberg channel.
Toy extension: two disjoint and independent CZ-pair operations on 4 atoms.
The computational Pauli map is trace-decreasing; loss and residual Rydberg
are reported outside the Pauli channel, not silently renormalized away.
"""
import ctypes, itertools, json, re, time
from pathlib import Path
import numpy as np

ROOT=Path(__file__).resolve().parent
N=4; R=3
PAULIS='IXYZ'

def pauli_to_bits(s):
    n=len(s); x=z=0
    for i,p in enumerate(s):
        if p in 'XY': x |= 1<<i
        if p in 'YZ': z |= 1<<i
    return x | (z<<n)

def bits_to_pauli(q,n=N):
    return ''.join({0:'I',1:'X',2:'Z',3:'Y'}[((q>>i)&1)+2*((q>>(n+i))&1)] for i in range(n))

def symp(a,b,n=N):
    mask=(1<<n)-1
    return (((a&mask)&(b>>n)) ^ ((a>>n)&(b&mask))).bit_count()&1

def span(basis):
    return tuple(sorted({0,*[__import__('functools').reduce(lambda a,b:a^b, (basis[j] for j in range(len(basis)) if mask>>j&1),0) for mask in range(1,1<<len(basis))]}))

def enumerate_stabilizers(n=N,r=R):
    """All rank-r totally isotropic binary symplectic subspaces; n=4,r=3."""
    assert (n,r)==(4,3), 'This exact finite-search experiment is for [[4,1]].'
    seen={}
    for a in range(1,1<<(2*n)):
        for b in range(a+1,1<<(2*n)):
            if symp(a,b,n): continue
            ab=a^b
            for c in range(b+1,1<<(2*n)):
                if c==ab or symp(a,c,n) or symp(b,c,n):continue
                group=tuple(sorted((0,a,b,c,ab,a^c,b^c,ab^c)))
                if group not in seen:seen[group]=(a,b,c)
    return seen

def pair_probabilities(npz_file):
    if str(npz_file).endswith('.json'):
        q={str(k):float(v) for k,v in json.loads(Path(npz_file).read_text()).items()}
    else:
        d=np.load(npz_file)
        q={str(k):float(v) for k,v in zip(d['Pauli_labels'],d['Pauli_weights'])}
    assert len(q)==16
    assert abs(sum(q.values())-.9984648622673)<5e-8
    weights=np.zeros(256)
    for p,pp in q.items():
        for s,ps in q.items():
            weights[pauli_to_bits(p+s)] = pp*ps
    assert np.isclose(weights.sum(),sum(q.values())**2,atol=1e-12)
    return q,weights

def stabilizer_syndrome(e,basis,n=N):
    return sum(symp(e,g,n)<<i for i,g in enumerate(basis))

def correctable_error_set(group,basis, errors):
    # Includes I. E_a^dag E_b is in the stabilizer or has nonzero syndrome.
    G=set(group)
    for a,b in itertools.combinations(errors,2):
        d=a^b
        if d not in G and stabilizer_syndrome(d,basis)==0:return False
    return True

def ml_success(group,basis,weights):
    """Exact maximum-likelihood correction of a known Pauli error channel.
    Group same-syndrome errors into stabilizer cosets. For each syndrome,
    take the most likely logical coset (includes degeneracy)."""
    seen=set(); max_coset=np.zeros(1<<len(basis)); winners={}
    for e in range(len(weights)):
        if e in seen:continue
        coset=tuple(e^s for s in group)
        seen.update(coset)
        mass=float(sum(weights[x] for x in coset))
        syn=stabilizer_syndrome(e,basis)
        if mass>max_coset[syn]:
            max_coset[syn]=mass
            winners[syn]=min(coset)
    return float(max_coset.sum()),winners

def detection_weight(group,basis,weights):
    G=set(group)
    return float(sum(p for e,p in enumerate(weights) if e in G or stabilizer_syndrome(e,basis)!=0))

def smt(smt_hard_errors,seed=(),minimum_distance=2):
    """Genuine SAT/SMT formula: arbitrary general Pauli generators.
    Def 3+1 stabilizers with n=4 and 3 generators; GF(2) commutation,
    independence, KL pair constraints including degeneracy, optional seeds.
    """
    n=N;r=R
    var=lambda k,t,i: f'g{k}_{t}{i}'
    opxor=lambda xs: ('false' if not xs else xs[0] if len(xs)==1 else f"(xor {' '.join(xs)})")
    opand=lambda xs: ('true' if not xs else xs[0] if len(xs)==1 else f"(and {' '.join(xs)})")
    opor=lambda xs: ('false' if not xs else xs[0] if len(xs)==1 else f"(or {' '.join(xs)})")
    L=['(set-logic QF_UF)','(set-option :produce-models true)']
    # QF_UF includes booleans and uninterpreted sorts; pure propositional equations
    for k in range(r):
        for i in range(n):
            for t in 'xz':L.append(f'(declare-const {var(k,t,i)} Bool)')
    for k in range(r):
        for i in range(n):
            for t in 'xz':
                if k<len(seed):
                    p=pauli_to_bits(seed[k]); expected=(p>>(i+(n if t=='z' else 0)))&1
                    L.append(f'(assert (= {var(k,t,i)} {str(bool(expected)).lower()}))')
    for k in range(r):
        for l in range(k+1,r):
            terms=[]
            for i in range(n):
                terms += [f'(and {var(k,"x",i)} {var(l,"z",i)})', f'(and {var(k,"z",i)} {var(l,"x",i)})']
            L.append(f'(assert (not {opxor(terms)}))')
    # generators are linearly independent, exact number of logical qubits
    for mask in range(1,1<<r):
        terms=[]
        for i in range(n):
            for t in 'xz':terms.append(opxor([var(k,t,i) for k in range(r) if mask>>k&1]))
        L.append(f'(assert {opor(terms)})')
    def synd(b):
        xx=b&((1<<n)-1);zz=b>>n
        bits=[]
        for k in range(r):
            bits.append(opxor([var(k,'z',i) for i in range(n) if (xx>>i)&1]
                             +[var(k,'x',i) for i in range(n) if (zz>>i)&1]))
        return opor(bits)
    def isin(b):
        terms=[]
        for mask in range(1<<r):
            checks=[]
            for i in range(n):
                for t in 'xz':
                    v=opxor([var(k,t,i) for k in range(r) if (mask>>k)&1])
                    bit=(b>>(i+(n if t=='z' else 0)))&1
                    checks.append(v if bit else f'(not {v})')
            terms.append(opand(checks))
        return opor(terms)
    for a,b in itertools.combinations(smt_hard_errors,2):
        d=a^b
        L.append(f'(assert (or {synd(d)} {isin(d)}))')
    # No nontrivial logical operator of weight < minimum_distance.
    for v in range(1,1<<(2*n)):
        if sum(t != 'I' for t in bits_to_pauli(v,n)) < minimum_distance:
            L.append(f'(assert (or {synd(v)} {isin(v)}))')
    L+=['(check-sat)','(get-model)']
    return '\n'.join(L)+'\n'

def run_z3(smt_text):
    try:
        lib=ctypes.CDLL('libz3.so.4')
    except OSError:
        return 'not-installed',None
    lib.Z3_mk_config.restype=ctypes.c_void_p
    lib.Z3_mk_context.argtypes=[ctypes.c_void_p];lib.Z3_mk_context.restype=ctypes.c_void_p
    lib.Z3_eval_smtlib2_string.argtypes=[ctypes.c_void_p,ctypes.c_char_p];lib.Z3_eval_smtlib2_string.restype=ctypes.c_char_p
    cfg=lib.Z3_mk_config();ctx=lib.Z3_mk_context(cfg)
    out=lib.Z3_eval_smtlib2_string(ctx,smt_text.replace('(get-model)','').encode()).decode()
    status=out.splitlines()[0]
    if status!='sat':return status,None
    ctx=lib.Z3_mk_context(cfg)
    out=lib.Z3_eval_smtlib2_string(ctx,smt_text.encode()).decode()
    bits={m.group(1):m.group(2)=='true' for m in re.finditer(r'\(define-fun\s+(g\d+_[xz]\d+)\s*\(\)\s*Bool\s*(true|false)\s*\)',out)}
    if len(bits)!=2*N*R:raise ValueError(f'Z3 model parse failed: {len(bits)} assigned; {out[:400]}')
    gens=[]
    for k in range(R):
        x=z=0
        for i in range(N):
            if bits[f'g{k}_x{i}']:x |= 1<<i
            if bits[f'g{k}_z{i}']:z |= 1<<i
        gens.append(x|(z<<N))
    return status,gens

def distance(group,basis,n=N):
    group=set(group)
    return min(sum(c!='I' for c in bits_to_pauli(e,n))
               for e in range(1,1<<(2*n)) if e not in group and stabilizer_syndrome(e,basis,n)==0)

def study(source_file,folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    t0=time.perf_counter()
    q,w=pair_probabilities(source_file)
    allcodes=enumerate_stabilizers()
    assert len(allcodes)==11475,len(allcodes)
    hard_strings=('IIII','ZIII','IZII','IIZI','IIIZ')
    hard_errors=[pauli_to_bits(x) for x in hard_strings]
    seeds=[(),('ZZII',),('ZZII','IIZZ')]
    saved={}
    for seed in seeds:
        k='unseeded' if not seed else ('seed_ZZII' if len(seed)==1 else 'seed_ZZII_IIZZ')
        seedints=set(map(pauli_to_bits,seed))
        subset=[(G,basis) for G,basis in allcodes.items() if seedints.issubset(G)]
        feasible=[(G,basis) for G,basis in subset if correctable_error_set(G,basis,hard_errors) and distance(G,basis)>=2]
        scored=[]
        for group,basis in feasible:
            succ,winners=ml_success(group,basis,w)
            # tuple deterministically breaks ties
            scored.append((succ,group,basis,winners))
        scored.sort(key=lambda x:(-round(x[0],16),x[1]))
        if scored:
            success,G,basis,winners=scored[0]
            det=detection_weight(G,basis,w)
            representative=tuple(bits_to_pauli(x) for x in basis)
            code={"generators":representative,"distance":distance(G,basis),"group":[bits_to_pauli(x) for x in G],
                  "decoder_coset_representatives":{str(s):bits_to_pauli(e) for s,e in sorted(winners.items())},
                  "success_unconditional":success,"failure_unconditional":w.sum()-success,
                  "logical_failure_conditional_on_computational_survival":(w.sum()-success)/w.sum(),
                  "detected_or_stabilizer_prob":det}
        else:code=None
        saved[k]={"seed":list(seed),"total_stabilizers_considered":len(subset),"hard_kl_feasible":len(feasible),"best":code}
        # independently solve general stabilizer constraints with Z3, using seed and KL conditions
        form=smt(hard_errors,seed)
        (folder/f'{k}.smt2').write_text(form)
        status,model=run_z3(form)
        saved[k]['smt_status']=status
        saved[k]['smt_generators']=list(map(bits_to_pauli,model)) if model else None
        if model:
            assert len(span(model))==8
            assert correctable_error_set(span(model),model,hard_errors)
            assert distance(span(model),model)>=2
            assert seedints.issubset(set(span(model)))
        if status!='not-installed':
            assert (status=='sat')==bool(feasible)
    # Second mode: no mandatory correction set; retain stabilizer distance >=2
    # and optimize the exact likelihood of successful recovery.
    saved['weighted_only_d2']={}
    for seed in seeds:
        label='unseeded' if not seed else ('seed_ZZII' if len(seed)==1 else 'seed_ZZII_IIZZ')
        seedints=set(map(pauli_to_bits,seed))
        subset=[(G,basis) for G,basis in allcodes.items() if seedints.issubset(G) and distance(G,basis)>=2]
        scored=[(ml_success(G,basis,w)[0],G,basis) for G,basis in subset]
        scored.sort(key=lambda x:(-round(x[0],16),x[1]))
        succ,G,basis=scored[0]
        form=smt([pauli_to_bits('IIII')],seed,minimum_distance=2)
        (folder/f'weighted_{label}.smt2').write_text(form)
        status,model=run_z3(form)
        assert status in ('sat','not-installed')
        assert (model is None) or distance(span(model),model)>=2
        saved['weighted_only_d2'][label]={
            'seed':list(seed), 'distance2_candidates':len(subset),
            'smt_status':status, 'optimal_generators':[bits_to_pauli(g) for g in basis],
            'optimal_distance':distance(G,basis),
            'success_unconditional':succ,
            'logical_failure_conditional_on_computational_survival':float((w.sum()-succ)/w.sum()),
            'corrects_all_single_Z':correctable_error_set(G,basis,hard_errors),
        }
    saved['setup']={"n":4,"k":1,"r":3,"hard_error_set":list(hard_strings),
                    "noise_source":str(source_file),"pair_probability_total":sum(q.values()),
                    "two_pair_computational_survival":float(w.sum()),
                    "two_pair_outside_computational":float(1-w.sum()),
                    "assumption":"Two independent pair channels (atoms 1-2 and 3-4), identical pulse; not a measured four-atom channel.",
                    "source_pauli_weights":q,"n_stabilizer_groups":len(allcodes),"minimum_distance":2,
                    "runtime_seconds":time.perf_counter()-t0,
                    "seed_interpretation":"The 3-level Rydberg-only drive has exact Q_i=diag(1,-1,-1) symmetries (when hyperfine drive is off); on the computational |0>,|1> subspace, Q_1 Q_2 maps to ZZII and Q_3 Q_4 maps to IIZZ. Jumps can break these symmetries; the additional discovered generators need not commute with H."}
    (folder/'results.json').write_text(json.dumps(saved,indent=2))
    print(json.dumps({k:{'N':v['total_stabilizers_considered'],'feasible':v['hard_kl_feasible'],'status':v['smt_status'],'best':v['best']} for k,v in saved.items() if k not in ('setup','weighted_only_d2')},indent=2))
    print('Weighted-only:',json.dumps(saved['weighted_only_d2'],indent=2))
    print('Time',time.perf_counter()-t0)
    return saved

if __name__=='__main__':
    source=ROOT/'rydberg_pair_pauli_weights.json'
    if not source.exists(): source=ROOT/'calibrated_channel_kraus.npz'
    study(source,ROOT/'SAT_results')