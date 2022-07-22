from functools import reduce
from sys import stdout

import numpy
import scipy
from scipy import linalg

import pyscf
from pyscf import tools
from pyscf import gto, scf, lo
from pyscf import mcscf

def analyze(ci, cas=None, tol=1e-3, nstate=10):
    assert ci is not None

    print('** Largest CI components **')
    print('  [alpha occ-orbitals] [beta occ-orbitals]            CI coefficient')
    res = cas.fcisolver.large_ci(
        ci, cas.ncas, cas.nelecas, tol, return_strs=False
        )
    
    aa_list = []
    bb_list = []
    cc_list = []
    
    for c, ia, ib in res:
        aa = [a+cas.ncore for a in ia]
        bb = [b+cas.ncore for b in ib]

        cc = c**2
        aa_list.append(aa)
        bb_list.append(bb)
        cc_list.append(cc)
    
    idx_list = numpy.argsort(cc_list)[::-1][:nstate]
    
    cc_sum = 0.0
    for idx in idx_list:
        print('  %-20s %-30s %12.4e'%(aa_list[idx], bb_list[idx], cc_list[idx]))
        cc_sum += cc_list[idx]

    print("cc_sum = %6.4f"%cc_sum)

def check_locality(ovlp_ao, c_ao_lo, idx_a, idx_b, tol=1e-1):
    nao, nlo = c_ao_lo.shape
    ct_s_c = reduce(numpy.dot, [c_ao_lo.T, ovlp_ao, c_ao_lo])
    
    assert numpy.linalg.norm(ct_s_c - numpy.eye(nlo)) < tol ** 4

    from scipy.linalg import sqrtm
    w  = numpy.einsum("mn,np->mp", sqrtm(ovlp_ao), c_ao_lo)
    w2 = w**2

    lo_idx_a = []
    lo_idx_b = []

    for p in range(nlo):
        w2_a = sum(w2[idx_a, p])
        w2_b = sum(w2[idx_b, p])

        if w2_a > 0.5 + tol:
            assert w2_a > w2_b
            lo_idx_a.append(p)
        elif w2_b > 0.5 + tol:
            assert w2_b > w2_a
            lo_idx_b.append(p)
        else:
            print("w2_a = %6.4f, w2_b = %6.4" % (w2_a, w2_b))
    
    return lo_idx_a, lo_idx_b


def build_pio(mf, c_ao_lo=None, frag_a_label=None, frag_b_label=None, tol=1e-2):
    assert c_ao_lo is not None
    assert frag_a_label is not None
    assert frag_b_label is not None
    
    mol   = mf.mol
    ao_idx_a = mol.search_ao_label(frag_a_label)
    ao_idx_b = mol.search_ao_label(frag_b_label)

    ovlp_ao = mf.get_ovlp()

    lo_idx_a, lo_idx_b = check_locality(ovlp_ao, c_ao_lo, ao_idx_a, ao_idx_b, tol=0.2)

    assert mf.converged
    dm_ao   = mf.make_rdm1()    
    dm_lo   = reduce(numpy.dot, (c_ao_lo.T, ovlp_ao, dm_ao, ovlp_ao, c_ao_lo))

    pp_mns_p = numpy.dot(dm_lo/2, dm_lo/2) - dm_lo/2
    assert numpy.linalg.norm(pp_mns_p) < 1e-8

    dm_aa = dm_lo[numpy.ix_(lo_idx_a, lo_idx_a)]
    dm_bb = dm_lo[numpy.ix_(lo_idx_b, lo_idx_b)]
    dm_ab = dm_lo[numpy.ix_(lo_idx_a, lo_idx_b)]

    ua, sa, vha = linalg.svd(dm_aa)

    cor_idx_a = []
    act_idx_a = []
    ext_idx_a = []

    for p, sp in enumerate(sa):
        if abs(sp - 2.0) < tol:
            cor_idx_a.append(p)
        elif abs(sp) < tol:
            ext_idx_a.append(p)
        else:
            act_idx_a.append(p)

    ub, sb, vhb = linalg.svd(dm_bb)

    cor_idx_b = []
    act_idx_b = []
    ext_idx_b = []

    for p, sp in enumerate(sb):
        if abs(sp - 2.0) < tol:
            cor_idx_b.append(p)
        elif abs(sp) < tol:
            ext_idx_b.append(p)
        else:
            act_idx_b.append(p)

    cas_nmo  = 0
    cas_nelec = 0.0

    cas_nmo  += len(sa[act_idx_a])
    cas_nelec += sum(sa[act_idx_a])

    cas_nmo  += len(sb[act_idx_b])
    cas_nelec += sum(sb[act_idx_b])

    cas_nelec  = int(numpy.round(cas_nelec))

    nao, nmo = mf.mo_coeff.shape
    coeff    = numpy.empty([nao, nmo])

    cor_list = []
    cas_list = []
    ext_list = []

    c_ao_lo_a = c_ao_lo[:, lo_idx_a]
    c_pio_a   = reduce(numpy.dot, [c_ao_lo_a, ua])

    c_ao_lo_b = c_ao_lo[:, lo_idx_b]
    c_pio_b   = reduce(numpy.dot, [c_ao_lo_b, ub])

    count = 0

    for p in cor_idx_a:
        pp = count
        cor_list.append(pp)
        coeff[:, pp] = c_pio_a[:, p]
        count += 1

    for p in cor_idx_b:
        pp = count
        cor_list.append(pp)
        coeff[:, pp] = c_pio_b[:, p]
        count += 1

    for p in act_idx_a:
        pp = count
        cas_list.append(pp)
        coeff[:, pp] = c_pio_a[:, p]
        count += 1

    for p in act_idx_b:
        pp = count
        cas_list.append(pp)
        coeff[:, pp] = c_pio_b[:, p]
        count += 1
        
    for p in ext_idx_a:
        pp = count
        ext_list.append(pp)
        coeff[:, pp] = c_pio_a[:, p]
        count += 1

    for p in ext_idx_b:
        pp = count
        ext_list.append(pp)
        coeff[:, pp] = c_pio_b[:, p]
        count += 1

    print("cor_list = %s" % cor_list)
    print("cas_list = %s" % cas_list)
    print("ext_list = %s" % ext_list)

    assert count == nmo

    assert numpy.linalg.norm(
        reduce(numpy.dot, [coeff.T, ovlp_ao, coeff]) - numpy.eye(nmo)
    )

    with open("./cube_list.log", "w") as f:
        for pp in cas_list:
            print(f"Processing {pp}")
            file_path = f'cas-{pp}'
            tools.cubegen.orbital(mol, "./cube/"+file_path+".cube", coeff[:, pp])
            f.write(file_path+"\n")

    return coeff, cas_list, (cas_nmo, cas_nelec)

if __name__ == "__main__":
    mol = gto.M(atom="""
    C1    0.0000000    0.0000000   -0.6654946
    H1    0.0000000    0.9374306   -1.2121600
    H1   -0.0000000   -0.9374306   -1.2121600
    C2   -0.0000000    0.0000000    0.6654946
    H2    0.0000000    0.9374306    1.2121600
    H2   -0.0000000   -0.9374306    1.2121600
    """, basis='cc-pvdz')

    mol.build()

    mf = scf.RHF(mol).run()
    pm = pyscf.lo.PM(mol, mf.mo_coeff, mf)
    pm.pop_method = 'iao'
    c_ao_lo = pm.kernel()

    coeff, cas_list, cas_info = build_pio(
        mf, c_ao_lo=c_ao_lo,
        frag_a_label=["C1", "H1"],
        frag_b_label=["C2", "H2"],
        tol=1e-2
        )
    cas_nmo, cas_nelec = cas_info

    cas   = mcscf.CASCI(mf, cas_nmo, cas_nelec)
    coeff = cas.sort_mo(cas_list, mo_coeff=coeff, base=0)
    cas.kernel(coeff)
    cas.verbose = 0

    analyze(cas.ci, cas=cas, tol=1e-4, nstate=10)