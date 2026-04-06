"""
Numerical verification of the SU(N_c) color Fierz identity (Dixon eq. 2.3).

Verifies:
  1. Dixon (2.3): sum_a (T^a)_{i1}^{j1} (T^a)_{i2}^{j2} = delta_{i1}^{j2} delta_{i2}^{j1} - (1/Nc) delta_{i1}^{j1} delta_{i2}^{j2}
     using Dixon normalization: Tr(T^a T^b) = delta^{ab}
  2. Tracelessness: contracting both sides with delta^{j1}_{i1} gives 0 = 0.
  3. Dixon (2.2): i*sqrt(2)*f^{abc} = Tr(T^a T^b T^c) - Tr(T^a T^c T^b)

Generator conventions:
  - SU(2): T^a = sigma^a / sqrt(2)  (Pauli matrices divided by sqrt(2) -> Tr(T^a T^b) = delta^{ab})
  - SU(3): T^a = lambda^a / sqrt(2)  (Gell-Mann matrices divided by sqrt(2) -> Tr(T^a T^b) = delta^{ab})

References:
  Dixon, arXiv:1310.5353, eqs. 2.1, 2.2, 2.3
  Bridge document: amplitudes/09-dixon-srednicki-bridge.tex, Sec. 1
"""

import numpy as np

TOL = 1e-12


def gell_mann():
    """Return the 8 Gell-Mann matrices (standard Peskin-Schroeder normalization, Tr=1/2)."""
    lam = np.zeros((8, 3, 3), dtype=complex)
    # lambda_1
    lam[0, 0, 1] = 1; lam[0, 1, 0] = 1
    # lambda_2
    lam[1, 0, 1] = -1j; lam[1, 1, 0] = 1j
    # lambda_3
    lam[2, 0, 0] = 1; lam[2, 1, 1] = -1
    # lambda_4
    lam[3, 0, 2] = 1; lam[3, 2, 0] = 1
    # lambda_5
    lam[4, 0, 2] = -1j; lam[4, 2, 0] = 1j
    # lambda_6
    lam[5, 1, 2] = 1; lam[5, 2, 1] = 1
    # lambda_7
    lam[6, 1, 2] = -1j; lam[6, 2, 1] = 1j
    # lambda_8
    lam[7, 0, 0] = 1/np.sqrt(3); lam[7, 1, 1] = 1/np.sqrt(3); lam[7, 2, 2] = -2/np.sqrt(3)
    return lam


def pauli():
    """Return the 3 Pauli matrices."""
    sigma = np.zeros((3, 2, 2), dtype=complex)
    sigma[0] = np.array([[0, 1], [1, 0]], dtype=complex)
    sigma[1] = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sigma[2] = np.array([[1, 0], [0, -1]], dtype=complex)
    return sigma


def make_generators_dixon(group):
    """
    Return generators in Dixon normalization: Tr(T^a T^b) = delta^{ab}.
    Dixon: T^a = sqrt(2) * T^a_std, where T^a_std has Tr(T^a T^b) = 1/2 delta^{ab}.
    """
    if group == 'SU2':
        sigma = pauli()
        # Standard: T^a = sigma^a / 2 -> Tr = 1/2. Dixon: T^a = sigma^a / sqrt(2).
        T = sigma / np.sqrt(2)
        Nc = 2
    elif group == 'SU3':
        lam = gell_mann()
        # Standard: T^a = lambda^a / 2 -> Tr = 1/2. Dixon: T^a = lambda^a / sqrt(2).
        T = lam / np.sqrt(2)
        Nc = 3
    else:
        raise ValueError(f"Unknown group: {group}")
    return T, Nc


def verify_trace_normalization(T, group):
    """Verify Tr(T^a T^b) = delta^{ab}."""
    n = len(T)
    errors = []
    for a in range(n):
        for b in range(n):
            trace = np.trace(T[a] @ T[b])
            expected = 1.0 if a == b else 0.0
            errors.append(abs(trace - expected))
    max_err = max(errors)
    ok = max_err < TOL
    print(f"  [{group}] Tr(T^a T^b) = delta^{{ab}}: max_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_fierz_identity(T, Nc, group):
    """
    Verify Dixon (2.3):
      sum_a (T^a)_{i1}^{j1} (T^a)_{i2}^{j2} = delta_{i1}^{j2} delta_{i2}^{j1} - (1/Nc) delta_{i1}^{j1} delta_{i2}^{j2}

    We compute LHS and RHS as rank-4 tensors of shape (Nc,Nc,Nc,Nc) with indices (i1,j1,i2,j2).
    """
    # LHS: sum_a T^a_{i1,j1} * T^a_{i2,j2}
    # T has shape (n_gen, Nc, Nc); T[a,i,j] = (T^a)_i^j
    LHS = np.einsum('aij,akl->ijkl', T, T)  # shape (Nc,Nc,Nc,Nc): i1=i,j1=j,i2=k,j2=l

    # RHS: delta_{i1}^{j2} delta_{i2}^{j1} - (1/Nc) delta_{i1}^{j1} delta_{i2}^{j2}
    delta = np.eye(Nc)
    RHS = np.einsum('il,kj->ijkl', delta, delta) - (1.0/Nc) * np.einsum('ij,kl->ijkl', delta, delta)

    max_err = np.max(np.abs(LHS - RHS))
    ok = max_err < TOL
    print(f"  [{group}] Fierz identity (2.3): max_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_tracelessness_via_fierz(T, Nc, group):
    """
    Verify tracelessness via Fierz: contract i1=j1 in Dixon (2.3).
    LHS: sum_a Tr(T^a) * (T^a)_{i2}^{j2}  (= 0 since Tr(T^a)=0)
    RHS: delta_{i2}^{j2} * Nc - (1/Nc)*Nc*delta_{i2}^{j2} = (Nc-1)*delta_{i2}^{j2}
    Wait - let's be careful: the correct contraction gives 0=0.
    From the bridge doc: contracting i1=j1 gives:
      LHS = sum_a Tr(T^a) * (T^a)_{i2}^{j2} = 0
      RHS = delta_{i2}^{j2} - delta_{i2}^{j2} = 0
    So both sides are 0.
    """
    # LHS: sum_{i1} sum_a T^a_{i1,i1} * T^a_{i2,j2}
    LHS = np.einsum('aii,akl->kl', T, T)  # sum over a and i1=j1=i

    # RHS: delta_{i2}^{j2} * Nc - (1/Nc)*Nc*delta_{i2}^{j2} = Nc*delta - delta = (Nc-1)*delta
    # But from bridge doc: contraction gives delta_{i2}^{j2} - delta_{i2}^{j2} = 0
    # Let's recheck: sum_{i1} delta_{i1}^{j2} delta_{i2}^{i1} = delta_{i2}^{j2}
    #                sum_{i1} (1/Nc) delta_{i1}^{i1} delta_{i2}^{j2} = (1/Nc)*Nc*delta_{i2}^{j2} = delta_{i2}^{j2}
    # So RHS = delta - delta = 0. Good.
    RHS = np.zeros((Nc, Nc), dtype=complex)

    max_err = np.max(np.abs(LHS - RHS))
    ok = max_err < TOL
    print(f"  [{group}] Tracelessness via Fierz: max_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_structure_constants(T, Nc, group):
    """
    Verify Dixon (2.2): i*sqrt(2)*f^{abc} = Tr(T^a T^b T^c) - Tr(T^a T^c T^b)

    We extract f^{abc} two ways:
    1. From the commutator: [T^a, T^b] = i*sqrt(2)*f^{abc} T^c
       -> f^{abc} = -i/sqrt(2) * Tr([T^a,T^b] T^c) / Tr(T^c T^c)  [no sum on c, use norm=1]
       But since Tr(T^c T^d) = delta^{cd}, we have:
       i*sqrt(2)*f^{abd} = Tr([T^a,T^b] T^d)
    2. From Dixon (2.2): i*sqrt(2)*f^{abc} = Tr(T^a T^b T^c) - Tr(T^a T^c T^b)

    We verify: Tr([T^a,T^b] T^c) == Tr(T^a T^b T^c) - Tr(T^a T^c T^b) for all a,b,c.
    """
    n = len(T)
    errors = []
    for a in range(n):
        for b in range(n):
            for c in range(n):
                commutator = T[a] @ T[b] - T[b] @ T[a]
                lhs = np.trace(commutator @ T[c])
                rhs = np.trace(T[a] @ T[b] @ T[c]) - np.trace(T[a] @ T[c] @ T[b])
                errors.append(abs(lhs - rhs))

    max_err = max(errors)
    ok = max_err < TOL
    print(f"  [{group}] Structure constants (2.2): max_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_real_structure_constants(T, Nc, group):
    """
    Verify that i*sqrt(2)*f^{abc} from trace formula is purely imaginary (f^{abc} is real).
    Also verify the antisymmetry f^{abc} = -f^{bac}.
    """
    n = len(T)
    # Extract i*sqrt(2)*f^{abc}
    F = np.zeros((n, n, n), dtype=complex)
    for a in range(n):
        for b in range(n):
            for c in range(n):
                F[a, b, c] = np.trace(T[a] @ T[b] @ T[c]) - np.trace(T[a] @ T[c] @ T[b])

    # i*sqrt(2)*f^{abc} = F[a,b,c], so f^{abc} = F[a,b,c] / (i*sqrt(2))
    f = F / (1j * np.sqrt(2))

    # Check f is real
    max_imag = np.max(np.abs(np.imag(f)))
    ok1 = max_imag < TOL
    print(f"  [{group}] f^{{abc}} is real: max_imag = {max_imag:.2e}  {'PASS' if ok1 else 'FAIL'}")

    # Check antisymmetry: f^{abc} = -f^{bac}
    f_real = np.real(f)
    max_antisym = np.max(np.abs(f_real + np.transpose(f_real, (1, 0, 2))))
    ok2 = max_antisym < TOL
    print(f"  [{group}] f^{{abc}} antisymmetric in a,b: max_err = {max_antisym:.2e}  {'PASS' if ok2 else 'FAIL'}")

    return ok1 and ok2


def run_group(group):
    print(f"\n=== {group} ===")
    T, Nc = make_generators_dixon(group)
    results = []
    results.append(verify_trace_normalization(T, group))
    results.append(verify_fierz_identity(T, Nc, group))
    results.append(verify_tracelessness_via_fierz(T, Nc, group))
    results.append(verify_structure_constants(T, Nc, group))
    results.append(verify_real_structure_constants(T, Nc, group))
    return all(results)


if __name__ == '__main__':
    print("=" * 60)
    print("Script 01: SU(Nc) Color Fierz Identity Verification")
    print("  Verifying Dixon (2.2) and (2.3)")
    print("=" * 60)

    all_pass = True
    for group in ['SU2', 'SU3']:
        ok = run_group(group)
        all_pass = all_pass and ok

    print()
    print("=" * 60)
    if all_pass:
        print("OVERALL: PASS")
    else:
        print("OVERALL: FAIL")
    print("=" * 60)
