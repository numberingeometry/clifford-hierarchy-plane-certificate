"""Exact hierarchy recursion for the gate U_{E6} of (C18) over Gaussian dyadic rationals.

Labels are x + (z << 5), with qubits c, t, d0, d1, d2 from least to most significant. An operator is
(control_xor_shift, four (matrix_id, coherent_i_phase) blocks). Memo keys remove only a scalar common to
all four blocks. No floating-point arithmetic enters a decision.
"""
from __future__ import annotations

import ctypes
import hashlib
import json
import os
from pathlib import Path
import time


def need(value, message):
    if not value:
        raise ValueError(message)


def _normalize(power, values):
    values = tuple(values)
    while power and all(not (a & 1 or b & 1) for a, b in values):
        values = tuple((a//2,b//2) for a,b in values)
        power -= 1
    return power, values


def _phase(matrix, exponent):
    power, data = matrix
    exponent %= 4
    if exponent == 0:
        return matrix
    if exponent == 1:
        return power, tuple((-b,a) for a,b in data)
    if exponent == 2:
        return power, tuple((-a,-b) for a,b in data)
    return power, tuple((b,-a) for a,b in data)


def _multiply_dense(left,right):
    le,left=left;re,right=right;out=[]
    # Sparse rows substantially reduce Gaussian operations on these Cliffords.
    rows=[tuple((k,*left[8*r+k]) for k in range(8) if left[8*r+k]!=(0,0)) for r in range(8)]
    for row in rows:
        for c in range(8):
            a=b=0
            for k,lr,li in row:
                rr,ri=right[8*k+c]
                a += lr*rr-li*ri;b += lr*ri+li*rr
            out.append((a,b))
    return _normalize(le+re,out)


def _adjoint_dense(matrix):
    e,data=matrix
    return e,tuple((data[8*c+r][0],-data[8*c+r][1]) for r in range(8) for c in range(8))


def _omega(a,b,n=3):
    mask=(1<<n)-1
    return (((a&mask)&(b>>n)).bit_count()+((b&mask)&(a>>n)).bit_count())&1


def _pauli_product(a,b,n=3):
    mask=(1<<n)-1;ax=a&mask;az=a>>n;bx=b&mask;bz=b>>n
    out=a^b
    exponent=((ax&az).bit_count()+(bx&bz).bit_count()-((ax^bx)&(az^bz)).bit_count()+2*(az&bx).bit_count())%4
    return out,exponent


def _pauli_dense(label):
    need(type(label) is int and 0<=label<64,'payload Pauli label outside domain')
    x=label&7;z=label>>3;unit=((1,0),(0,1),(-1,0),(0,-1))[(x&z).bit_count()%4]
    values=[(0,0)]*64
    for c in range(8):
        sign=-1 if (z&c).bit_count()&1 else 1
        values[8*(c^x)+c]=(sign*unit[0],sign*unit[1])
    return 0,tuple(values)


_MATRICES=[]
_ACTIONS=[]  # all 64 exact Hermitian Pauli images: (label, i exponent 0 or 2)
_MATRIX_IDS={}
_PRODUCTS={}
_INVERSES={}
_LEFT_PAULIS={}
_PAULI_IDS={}


def _action_from_generators(generators):
    out=[]
    for label in range(64):
        x=label&7;z=label>>3;image=0;exponent=(x&z).bit_count()
        for j in range(6):
            if label>>j&1:
                nxt,sign=generators[j]
                image,extra=_pauli_product(image,nxt)
                exponent += sign+extra
        need(exponent%4 in (0,2),'Clifford action lost Hermitian signs')
        out.append((image,exponent%4))
    need(len({v for v,_ in out})==64,'Clifford action is not bijective')
    return tuple(out)


def _intern(matrix,action):
    power,data=_normalize(*matrix)
    need(type(power) is int and power>=0 and len(data)==64 and all(type(a) is int and type(b) is int for a,b in data),'invalid Gaussian-dyadic payload matrix')
    first=next((v for v in data if v!=(0,0)),None)
    need(first is not None,'zero Clifford matrix')
    a,b=first;rotations=((a,b),(-b,a),(-a,-b),(b,-a))
    k=max(range(4),key=lambda j:rotations[j])
    normalized=_phase((power,data),k)
    ident=_MATRIX_IDS.get(normalized)
    if ident is None:
        ident=len(_MATRICES);_MATRIX_IDS[normalized]=ident
        _MATRICES.append(normalized);_ACTIONS.append(action)
    else:
        need(_ACTIONS[ident]==action,'exact matrix collided with a different signed Clifford action')
    return ident,(-k)%4


def _actual(ref):
    ident,phase=ref
    need(type(ident) is int and 0<=ident<len(_MATRICES) and type(phase) is int and 0<=phase<4,'invalid coherent Clifford reference')
    return _phase(_MATRICES[ident],phase)


def _product(a,b):
    ai,ap=a;bi,bp=b;key=(ai,bi)
    out=_PRODUCTS.get(key)
    if out is None:
        action=tuple((_ACTIONS[ai][label][0],(sign+_ACTIONS[ai][label][1])%4) for label,sign in _ACTIONS[bi])
        out=_intern(_multiply_dense(_MATRICES[ai],_MATRICES[bi]),action)
        _PRODUCTS[key]=out
    ident,phase=out
    return ident,(phase+ap+bp)%4


def _inverse(ref):
    ident,phase=ref;out=_INVERSES.get(ident)
    if out is None:
        action=[None]*64
        for label,(image,sign) in enumerate(_ACTIONS[ident]):
            action[image]=(label,sign)
        out=_intern(_adjoint_dense(_MATRICES[ident]),tuple(action))
        _INVERSES[ident]=out
    return out[0],(out[1]-phase)%4


def _left_pauli(label,ref):
    ident,phase=ref;key=(label,ident);out=_LEFT_PAULIS.get(key)
    if out is None:
        power,data=_MATRICES[ident];x=label&7;z=label>>3
        unit=((1,0),(0,1),(-1,0),(0,-1))[(x&z).bit_count()%4]
        values=[(0,0)]*64
        for r in range(8):
            sign=-1 if (z&r).bit_count()&1 else 1
            ur,ui=sign*unit[0],sign*unit[1]
            for c in range(8):
                a,b=data[8*r+c]
                values[8*(r^x)+c]=(ur*a-ui*b,ur*b+ui*a)
        action=tuple((image,(sign+2*_omega(label,image))%4) for image,sign in _ACTIONS[ident])
        out=_intern((power,tuple(values)),action);_LEFT_PAULIS[key]=out
    return out[0],(out[1]+phase)%4


def _initialize():
    if _MATRICES:
        return
    for label in range(64):
        action=tuple((v,2*_omega(label,v)) for v in range(64))
        ref=_intern(_pauli_dense(label),action)
        _PAULI_IDS[ref[0]]=label


_initialize()


def _root(label):
    _,data=_pauli_dense(label);out=[]
    for index,(a,b) in enumerate(data):
        diagonal=int(index//8==index%8)
        out.append((a+b+diagonal,b-a+diagonal))
    gens=[]
    for q in (1,2,4,8,16,32):
        if _omega(label,q):
            image,exponent=_pauli_product(label,q)
            exponent=(exponent-1)%4
            need(exponent in (0,2),'root Clifford sign malformed')
            gens.append((image,exponent))
        else:
            gens.append((q,0))
    ref=_intern(_normalize(1,out),_action_from_generators(gens))
    need(_multiply_dense(_actual(ref),_adjoint_dense(_actual(ref)))==_pauli_dense(0),'root exact unitarity')
    need(_multiply_dense(_actual(ref),_actual(ref))==_pauli_dense(label),'root exact square')
    return ref


def _word(labels):
    result=_intern(_pauli_dense(0),_ACTIONS[_MATRIX_IDS[_pauli_dense(0)]])
    for label in labels:
        result=_product(result,_root(label))
    return result


def gates():
    identity=_intern(_pauli_dense(0),tuple((v,0) for v in range(64)))
    words={'E6':((1,24,2,48,6),(8,24))}
    out={}
    for name,(aw,bw) in words.items():
        a=_word(aw);b=_word(bw)
        out[name]=(0,(identity,a,b,_product(a,b)))
        for ref in out[name][1]:
            need(_multiply_dense(_actual(ref),_adjoint_dense(_actual(ref)))==_pauli_dense(0),'gate block exact unitarity')
    return out


def validate_operator(op):
    need(type(op) is tuple and len(op)==2,'operator must be a shift/block tuple')
    shift,blocks=op
    need(type(shift) is int and 0<=shift<4 and type(blocks) is tuple and len(blocks)==4,'operator control domain')
    for ref in blocks:
        need(type(ref) is tuple and len(ref)==2,'malformed block reference')
        _actual(ref)


def canonical(op):
    validate_operator(op)
    shift,blocks=op;scalar=blocks[0][1]
    return shift,tuple((ident,(phase-scalar)%4) for ident,phase in blocks)


def conjugate_pauli(op,label):
    need(type(label) is int and 0<=label<1024,'five-qubit Pauli label outside domain')
    shift,blocks=op;x=label&31;z=label>>5
    alpha=x&3;beta=z&3;q=(x>>2)+((z>>2)<<3);out=[]
    for c in range(4):
        source=c^shift;a=blocks[source^alpha];b=blocks[source]
        payload,sign=_ACTIONS[a[0]][q]
        pair=_product((a[0],0),_inverse((b[0],0)))
        ref=_left_pauli(payload,pair)
        exponent=(alpha&beta).bit_count()+2*((source&beta).bit_count()&1)+a[1]-b[1]+sign
        out.append((ref[0],(ref[1]+exponent)%4))
    return alpha,tuple(out)


def to_dense(op):
    validate_operator(op);shift,blocks=op
    dense=tuple(_actual(ref) for ref in blocks);power=max(e for e,_ in dense);values=[(0,0)]*1024
    for c,(e,data) in enumerate(dense):
        for r in range(8):
            for s in range(8):
                a,b=data[8*r+s]
                values[32*(4*r+(c^shift))+4*s+c]=(a<<(power-e),b<<(power-e))
    return _normalize(power,values)


def pauli_operator(label):
    need(type(label) is int and 0<=label<1024,'five-qubit Pauli label outside domain')
    x=label&31;z=label>>5;alpha=x&3;beta=z&3;q=(x>>2)+((z>>2)<<3)
    action=tuple((v,2*_omega(q,v)) for v in range(64));ref=_intern(_pauli_dense(q),action)
    return alpha,tuple((ref[0],(ref[1]+(alpha&beta).bit_count()+2*((c&beta).bit_count()&1))%4) for c in range(4))


def is_pauli(op):
    shift,blocks=op;label=_PAULI_IDS.get(blocks[0][0])
    if label is None or any(ident!=blocks[0][0] for ident,_ in blocks):
        return False
    phases=tuple((p-blocks[0][1])%4 for _,p in blocks)
    return phases[1] in (0,2) and phases[2] in (0,2) and phases[3]==(phases[1]+phases[2])%4



def dense_is_pauli(matrix):
    """Independent exact full-matrix Pauli recognizer; no interned IDs/actions."""
    power,data=matrix
    need(type(power) is int and power>=0 and len(data)==1024,'malformed dense five-qubit matrix')
    first=[r for r in range(32) if data[32*r]!=(0,0)]
    if len(first)!=1:
        return False
    shift=first[0];a,b=data[32*shift]
    if a*a+b*b != 1<<(2*power):
        return False
    beta=0
    for j in range(5):
        value=data[32*((1<<j)^shift)+(1<<j)]
        if value==(-a,-b):
            beta |= 1<<j
        elif value!=(a,b):
            return False
    for r in range(32):
        for c in range(32):
            sign=-1 if (beta&c).bit_count()&1 else 1
            expected=(sign*a,sign*b) if r==c^shift else (0,0)
            if data[32*r+c]!=expected:
                return False
    return True



def audit_exact_blocks(output_dir=None):
    """Bind every optimized action table to direct exact matrix arithmetic.

    Checks every registered block's unitarity and all 64 stored Pauli images.
    This creates no new interned blocks and accepts no cached audit verdict.
    Invoke after both hierarchy and fresh support-table construction.
    """
    started=time.perf_counter();last=started;count=len(_MATRICES);checks=0
    identity=_pauli_dense(0)
    for ident in range(count):
        matrix=_MATRICES[ident];adjoint=_adjoint_dense(matrix)
        need(_multiply_dense(matrix,adjoint)==identity,'interned Clifford block is not exactly unitary')
        power,data=matrix
        for label in range(64):
            x=label&7;z=label>>3
            ur,ui=((1,0),(0,1),(-1,0),(0,-1))[(x&z).bit_count()%4]
            # Direct C*P: signed column permutation, with no stored action table.
            values=[]
            for r in range(8):
                for c in range(8):
                    a,b=data[8*r+(c^x)];sign=-1 if (z&c).bit_count()&1 else 1
                    values.append((sign*(a*ur-b*ui),sign*(a*ui+b*ur)))
            actual=_multiply_dense(_normalize(power,values),adjoint)
            image,phase=_ACTIONS[ident][label]
            need(type(image) is int and 0<=image<64 and type(phase) is int and phase in (0,2),
                 'stored Clifford action has malformed Pauli/sign fields')
            expected=_phase(_pauli_dense(image),phase)
            need(actual==expected,'stored Clifford action does not match its exact matrix')
            checks += 1
        now=time.perf_counter()
        if now-last>=10:
            last=now
            progress={'status':'BLOCK_ALGEBRA_AUDIT_RUNNING','accepted':False,
                      'blocks_completed':ident+1,'total_blocks':count,
                      'signed_Pauli_images_checked':checks,'seconds':now-started}
            if output_dir is not None:
                Path(output_dir).mkdir(parents=True,exist_ok=True)
            _write_json(None if output_dir is None else Path(output_dir)/'block-algebra-audit-progress.json',progress)
            print(json.dumps(progress),flush=True)
    need(len(_MATRICES)==count and checks==64*count,'block algebra audit domain changed or incomplete')
    report={'status':'EXACT_BLOCK_AUDIT_PASS','exact':True,
            'blocks_checked':count,'all_registered_blocks_covered':True,
            'unitarity_blocks_checked':count,'signed_Pauli_images_checked':checks,
            'Pauli_action_domain_per_block':64,'seconds':time.perf_counter()-started}
    if output_dir is not None:
        Path(output_dir).mkdir(parents=True,exist_ok=True)
    _write_json(None if output_dir is None else Path(output_dir)/'block-algebra-audit.json',report)
    return report


def _cache_counts():
    return {'exact_clifford_blocks':len(_MATRICES),'block_product_cache_entries':len(_PRODUCTS),
            'block_inverse_cache_entries':len(_INVERSES),'left_Pauli_cache_entries':len(_LEFT_PAULIS)}


def _upper_counts(engine):
    visits=dict(engine.visits);checks=dict(engine.image_checks)
    need(visits.get(5)==1 and checks.get(5)==1024,'upper C5 domain incomplete')
    need(visits.get(4)==1024 and checks.get(4)==1024*1024,'upper C4 domain incomplete')
    for level in (2,3):
        need(type(visits.get(level)) is int and visits[level]>0 and checks.get(level)==10*visits[level],
             'accepted upper generator-domain coverage incomplete')
    return {'distinct_memoized_operators_by_level':visits,'computed_conjugations_by_level':checks,
            'every_distinct_upper_C4_node_checked_all_1024_Paulis':True,
            'every_upper_C3_and_C2_node_checked_all_ten_generators':True}


def _write_json(path,value):
    if path is not None:
        Path(path).write_text(json.dumps(value,indent=2)+'\n',encoding='utf-8')


def _memory():
    if os.name=='nt':
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_=[('cb',wintypes.DWORD),('PageFaultCount',wintypes.DWORD),('PeakWorkingSetSize',ctypes.c_size_t),('WorkingSetSize',ctypes.c_size_t),('QuotaPeakPagedPoolUsage',ctypes.c_size_t),('QuotaPagedPoolUsage',ctypes.c_size_t),('QuotaPeakNonPagedPoolUsage',ctypes.c_size_t),('QuotaNonPagedPoolUsage',ctypes.c_size_t),('PagefileUsage',ctypes.c_size_t),('PeakPagefileUsage',ctypes.c_size_t)]
        counter=Counters();counter.cb=ctypes.sizeof(counter)
        kernel=ctypes.WinDLL('kernel32');kernel.GetCurrentProcess.restype=wintypes.HANDLE
        psapi=ctypes.WinDLL('psapi');psapi.GetProcessMemoryInfo.argtypes=[wintypes.HANDLE,ctypes.POINTER(Counters),wintypes.DWORD]
        if psapi.GetProcessMemoryInfo(kernel.GetCurrentProcess(),ctypes.byref(counter),counter.cb):
            return {'peak_working_set_bytes':counter.PeakWorkingSetSize,'working_set_bytes':counter.WorkingSetSize}
    return {}


class Engine:
    def __init__(self,progress=None):
        self.memo={};self.visits={};self.image_checks={};self.progress=progress;self.last_progress=0
        self.generators=tuple(1<<j for j in range(10))

    def member(self,op,level):
        need(type(level) is int and 1<=level<=5,'hierarchy level outside E2 domain')
        key=(canonical(op),level)
        if key not in self.memo:
            self.memo[key]=self._member(key[0],level)
        return self.memo[key]

    def _member(self,op,level):
        self.visits[level]=self.visits.get(level,0)+1
        now=time.monotonic()
        if self.progress is not None and now-self.last_progress>=10:
            self.last_progress=now;self.progress(self)
        if level==1:
            return is_pauli(op)
        labels=self.generators if level in (2,3) else range(1024)
        for label in labels:
            self.image_checks[level]=self.image_checks.get(level,0)+1
            if not self.member(conjugate_pauli(op,label),level-1):
                return False
        return True

    def failure(self,op,level):
        need(not self.member(op,level),'requested lower-level failure is false')
        if level==1:
            power,values=to_dense(op)
            return {'level':1,'non_pauli_dense':{'denominator_power':power,'gaussian_integer_pairs':[list(v) for v in values]}}
        labels=self.generators if level in (2,3) else range(1024)
        for label in labels:
            image=conjugate_pauli(op,label)
            if not self.member(image,level-1):
                return {'level':level,'pauli_label':label,'next':self.failure(image,level-1)}
        raise ValueError('missing lower-level failure chain')


def verify_failure(op,chain,expected_level=4):
    """Independent chain traversal: no membership memo or stored booleans."""
    need(type(chain) is dict and type(expected_level) is int and 1<=expected_level<=4,'malformed failure-chain record or expected level')
    level=chain.get('level')
    need(type(level) is int and level==expected_level,'failure-chain level changed or skipped')
    if level==1:
        power,values=to_dense(op);leaf=chain.get('non_pauli_dense',{})
        need(type(leaf) is dict,'malformed failure-leaf record')
        need(type(leaf.get('denominator_power')) is int and leaf['denominator_power']>=0,'malformed failure-leaf denominator')
        pairs=leaf.get('gaussian_integer_pairs')
        need(type(pairs) is list and len(pairs)==1024 and all(type(v) is list and len(v)==2 and all(type(a) is int for a in v) for v in pairs),'malformed failure-leaf Gaussian coefficients')
        need(leaf['denominator_power']==power and pairs==[list(v) for v in values],'failure-chain leaf operator changed')
        need(not dense_is_pauli((power,values)),'failure-chain leaf is Pauli')
        return True
    label=chain.get('pauli_label')
    need(type(label) is int and 0<=label<1024,'malformed failure-chain Pauli label')
    return verify_failure(conjugate_pauli(op,label),chain.get('next',{}),expected_level-1)


def run_hierarchy(output_dir=None):
    started=time.perf_counter();results={};directory=None if output_dir is None else Path(output_dir)
    if directory is not None:
        directory.mkdir(parents=True,exist_ok=True)
    for name,gate in gates().items():
        gate_start=time.perf_counter()
        def progress(engine):
            record={'status':'RUNNING','accepted':False,'gate':name,'seconds':time.perf_counter()-gate_start,'distinct_operators_by_level':engine.visits,'computed_conjugations_by_level':engine.image_checks,'memo_entries':len(engine.memo),**_cache_counts(),**_memory()}
            _write_json(None if directory is None else directory/f'{name}-hierarchy-progress.json',record)
            print(json.dumps(record),flush=True)
        engine=Engine(progress)
        need(engine.member(gate,5),f'{name} upper C5 membership failed')
        upper_coverage=_upper_counts(engine)
        need(not engine.member(gate,4),f'{name} lower C4 failure missing')
        chain=engine.failure(gate,4);verify_failure(gate,chain,4)
        need(engine.image_checks.get(5)==1024,f'{name} root upper domain incomplete')
        power,dense=to_dense(gate)
        result={'in_C5':True,'in_C4':False,'exact':True,'upper_coverage':upper_coverage,'unitarity_checked_exactly':True,'root_lifts_square_to_declared_Paulis':True,'root_words':{'E6':[[1,24,2,48,6],[8,24]]}[name],'control_block_order':['I','A','B','AB'],'complete_root_Pauli_count':1024,'full_Pauli_domain_at_levels_at_least_four':1024,'generator_count_at_levels_two_and_three':10,'lower_failure_chain':chain,'lower_failure_chain_independently_verified':True,'exact_gate_dense_sha256':hashlib.sha256(repr((power,dense)).encode('ascii')).hexdigest(),'distinct_memoized_operators_by_level':dict(engine.visits),'computed_conjugations_by_level':dict(engine.image_checks),'memo_entries':len(engine.memo),'seconds':time.perf_counter()-gate_start,'memory':_memory()}
        results[name]=result
        _write_json(None if directory is None else directory/f'{name}-hierarchy.json',result)
        print(json.dumps({'status':'GATE_HIERARCHY_PASS','gate':name,'seconds':result['seconds'],'visits':engine.visits,'image_checks':engine.image_checks,'memo_entries':len(engine.memo),**_cache_counts(),'memory':result['memory']}),flush=True)
    report={'status':'HIERARCHY_PASS','exact':True,'field':'Q(i)','conductor':4,'entry_ring':'Z[i,1/2]','coordinate_order':'c,t,d0,d1,d2 least-to-most significant','memoization':'full exact control-shift and coherent Clifford-block tuples; only one common Gaussian-unit scalar removed','results':results,'seconds':time.perf_counter()-started,'memory':_memory(),**_cache_counts(),'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    _write_json(None if directory is None else directory/'hierarchy.json',report)
    return report


if __name__=='__main__':
    import sys
    report=run_hierarchy(sys.argv[1] if len(sys.argv)>1 else None)
    print(json.dumps({k:v for k,v in report.items() if k!='results'},indent=2))
