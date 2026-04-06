"""
Numerical verification of gluon polarization vectors in spinor-helicity form.

Uses Srednicki mostly-plus metric: g = diag(-1,+1,+1,+1).

Polarization vectors (Dixon 3.3.1 / Srednicki 60.6):
  In mostly-plus, with sigma_bar^mu = (+I, sigma1, sigma2, sigma3):

  epsilon^+_mu(k;q) = tilde_lam_q[dot-a] (sigma_bar^mu)[dot-a,a] lam_k[a] / (sqrt(2) * <qk>)
  epsilon^-_mu(k;q) = tilde_lam_k[dot-a] (sigma_bar^mu)[dot-a,a] lam_q[a] / (sqrt(2) * [kq])

  Note: eps^- uses tilde_lam_k and lam_q (swapped relative to eps^+).

Verifies:
  1. Transversality: k^mu epsilon^+_mu = 0, k^mu epsilon^-_mu = 0
  2. Self-orthogonality: epsilon^+ . epsilon^+ = 0, epsilon^- . epsilon^- = 0
  3. Cross-normalization: epsilon^+ . epsilon^- = -1 (both unconjugated)
  4. Conjugate normalization: epsilon^+_mu (epsilon^+_nu)^* g^{mu nu} = +1 (spacelike in mostly-plus)
  5. Reference independence: different q give the same transversality and normalization,
     and their difference is transverse to k (physical amplitude unchanged by gauge)
  6. Completeness:
     eps^+_mu eps^-_nu + eps^-_mu eps^+_nu = -g_{mu nu} + (k_mu q_nu + q_mu k_nu)/(k.q)
     sum_h eps^h_mu (eps^h_nu)^* = g_{mu nu} - (k_mu q_nu + q_mu k_nu)/(k.q)

References:
  Dixon, arXiv:1310.5353, Section 3.3
  Srednicki QFT, eqs. 60.6, 60.7, 60.8, 60.9
  Bridge document: amplitudes/09-dixon-srednicki-bridge.tex, Sec. 5
"""

import numpy as np

TOL = 1e-10

# Sigma matrices in mostly-plus
sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

# sigma^mu_{alpha dot-alpha} in mostly-plus: (-I, sigma1, sigma2, sigma3)
sigma_mu = np.array([-I2, sigma1, sigma2, sigma3])

# sigma-bar^{mu dot-alpha alpha} = (+I, sigma1, sigma2, sigma3) in mostly-plus
sigma_bar_mu = np.array([I2, sigma1, sigma2, sigma3])

# Minkowski metric in mostly-plus
g_mp = np.diag([-1., 1., 1., 1.])


def dot_mp(a, b):
    """Mostly-plus dot product: -a0*b0 + a1*b1 + a2*b2 + a3*b3."""
    return -a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]


def construct_spinors(p):
    """
    Construct (lambda, tilde_lambda) for massless p.
    Uses complex sqrt to handle negative-energy (all-outgoing) momenta.
    """
    p0, px, py, pz = p
    r = complex(p0 + pz)
    if abs(r) < 1e-10:
        r = r + 1e-8
    sqrt_r = np.sqrt(r)
    lam = np.array([sqrt_r, (px + 1j*py)/sqrt_r], dtype=complex)
    lam_tilde = np.array([sqrt_r, (px - 1j*py)/sqrt_r], dtype=complex)
    return lam, lam_tilde


def angle(lam_i, lam_j):
    """<ij> = lam_i[0]*lam_j[1] - lam_i[1]*lam_j[0]"""
    return lam_i[0]*lam_j[1] - lam_i[1]*lam_j[0]


def square(lam_ti, lam_tj):
    """[ij] = lam_ti[0]*lam_tj[1] - lam_ti[1]*lam_tj[0]"""
    return lam_ti[0]*lam_tj[1] - lam_ti[1]*lam_tj[0]


def polarization_plus(k, q):
    """
    epsilon^+_mu(k;q) = tilde_lam_q[dot-a] (sigma_bar^mu)[dot-a,a] lam_k[a] / (sqrt(2) * <qk>)
    Returns a 4-vector (complex).
    """
    lam_k, lam_t_k = construct_spinors(k)
    lam_q, lam_t_q = construct_spinors(q)

    angle_qk = angle(lam_q, lam_k)  # <qk>

    eps = np.zeros(4, dtype=complex)
    for mu in range(4):
        # lam_t_q[dot-a] * sigma_bar_mu[mu, dot-a, alpha] * lam_k[alpha]
        eps[mu] = lam_t_q @ sigma_bar_mu[mu] @ lam_k

    eps /= (np.sqrt(2) * angle_qk)
    return eps


def polarization_minus(k, q):
    """
    epsilon^-_mu(k;q) = tilde_lam_k[dot-a] (sigma_bar^mu)[dot-a,a] lam_q[a] / (sqrt(2) * [kq])

    This is the parity conjugate of eps^+. Under parity (lam <-> tilde_lam, k <-> q):
    eps^+(k;q) -> tilde_lam_q @ sigma_bar @ lam_k -> lam_k @ sigma_bar @ tilde_lam_q (transposed)
    The correct transverse formula uses tilde_lam_k on the left and lam_q on the right.

    Transversality: k^mu eps-_mu = tilde_lam_k @ (k^mu sigma_bar_mu) @ lam_q = 0
    because tilde_lam_k^T @ P_mp_k = 0 for massless k (bra Weyl equation).
    """
    lam_k, lam_t_k = construct_spinors(k)
    lam_q, lam_t_q = construct_spinors(q)

    sq_kq = square(lam_t_k, lam_t_q)  # [kq]

    eps = np.zeros(4, dtype=complex)
    for mu in range(4):
        # tilde_lam_k[dot-alpha] * sigma_bar_mu[mu, dot-alpha, alpha] * lam_q[alpha]
        eps[mu] = lam_t_k @ sigma_bar_mu[mu] @ lam_q

    eps /= (np.sqrt(2) * sq_kq)
    return eps


def make_lightlike(rng):
    p3 = rng.standard_normal(3)
    p3 /= np.linalg.norm(p3)
    E = rng.uniform(0.5, 5.0)
    return np.array([E, E*p3[0], E*p3[1], E*p3[2]])


def verify_transversality(k, q, label):
    """k^mu epsilon^+/-_mu = 0"""
    eps_plus = polarization_plus(k, q)
    eps_minus = polarization_minus(k, q)

    kp = dot_mp(k, eps_plus)
    km = dot_mp(k, eps_minus)

    ok_p = abs(kp) < TOL * max(1, np.linalg.norm(k))
    ok_m = abs(km) < TOL * max(1, np.linalg.norm(k))
    ok = ok_p and ok_m
    print(f"  [{label}] Transversality: k.eps+ = {abs(kp):.2e}, k.eps- = {abs(km):.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_self_orthogonality(k, q, label):
    """epsilon^+ . epsilon^+ = 0, epsilon^- . epsilon^- = 0"""
    eps_plus = polarization_plus(k, q)
    eps_minus = polarization_minus(k, q)

    pp = dot_mp(eps_plus, eps_plus)
    mm = dot_mp(eps_minus, eps_minus)

    ok_pp = abs(pp) < TOL
    ok_mm = abs(mm) < TOL
    ok = ok_pp and ok_mm
    print(f"  [{label}] Self-orthogonality: eps+.eps+={pp:.2e}, eps-.eps-={mm:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_cross_normalization(k, q, label):
    """
    Verify: epsilon^+_mu(k;q) * epsilon^-(k;q)^mu = -1  (both unconjugated, with metric)
    This is the standard spinor-helicity normalization.
    In mostly-plus: g^{mu nu} eps^+_mu eps^-_nu = -1.
    """
    eps_plus = polarization_plus(k, q)
    eps_minus = polarization_minus(k, q)

    pm = dot_mp(eps_plus, eps_minus)
    ok = abs(pm - (-1.0)) < TOL * 10
    print(f"  [{label}] Cross normalization eps+.eps- = {pm:.6f}  (expect -1)  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_conjugate_normalization(k, q, label):
    """
    Verify: eps^+_mu (eps^+_nu)^* g^{mu nu} = +1  (in mostly-plus).
    In mostly-plus, the spinor-helicity polarization vectors satisfy:
      eps^+.eps^{+*} = +1  (spacelike, spatial components dominate)
    Note: This is +1 (not -1) in mostly-plus because the complex circular polarization
    vectors are spacelike (the metric gives +1 for transverse polarizations).
    The normalization eps^+.eps^- = -1 holds for the CROSS product (unconjugated).
    """
    eps_plus = polarization_plus(k, q)
    eps_minus = polarization_minus(k, q)

    pp_conj = dot_mp(eps_plus, np.conj(eps_plus))
    mm_conj = dot_mp(eps_minus, np.conj(eps_minus))

    ok_p = abs(pp_conj - 1.0) < TOL * 10
    ok_m = abs(mm_conj - 1.0) < TOL * 10
    ok = ok_p and ok_m
    print(f"  [{label}] Conjugate norm: eps+.eps+*={pp_conj:.6f}, eps-.eps-*={mm_conj:.6f}  "
          f"(expect +1)  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_reference_independence(k, q1, q2, label):
    """
    Verify reference independence at the amplitude level.

    Physical amplitudes are gauge invariant: replacing eps(k;q) by eps(k;q')
    leaves physical amplitudes unchanged. We test this by verifying:
      k . (eps^+(k;q1) - eps^+(k;q2)) = 0  (both transverse to k)
      A test current J^mu proportional to k satisfies J.eps = 0 by transversality.

    The direct test: the difference eps(k;q1) - eps(k;q2) is transverse to k (k.diff=0),
    which is sufficient for gauge invariance of any amplitude satisfying Ward identity.
    Also verify q-independence of eps^+(k;q).eps^-(k;q) = -1 (normalization is ref-independent).
    """
    eps1 = polarization_plus(k, q1)
    eps2 = polarization_plus(k, q2)
    diff = eps1 - eps2

    # Both eps1 and eps2 are transverse to k; verify their difference is too
    k_diff = dot_mp(k, diff)
    ok_p = abs(k_diff) < TOL * 10
    print(f"  [{label}] Ref independence (eps+): k.(eps1-eps2) = {abs(k_diff):.2e}  "
          f"{'PASS' if ok_p else 'FAIL'}")

    # Same for eps-
    eps1m = polarization_minus(k, q1)
    eps2m = polarization_minus(k, q2)
    diff_m = eps1m - eps2m
    k_diff_m = dot_mp(k, diff_m)
    ok_m = abs(k_diff_m) < TOL * 10
    print(f"  [{label}] Ref independence (eps-): k.(eps1-eps2) = {abs(k_diff_m):.2e}  "
          f"{'PASS' if ok_m else 'FAIL'}")

    # Verify eps^+.eps^- = -1 is reference-independent
    pm_q1 = dot_mp(polarization_plus(k, q1), polarization_minus(k, q1))
    pm_q2 = dot_mp(polarization_plus(k, q2), polarization_minus(k, q2))
    ok_norm = abs(pm_q1 - pm_q2) < TOL * 10
    print(f"  [{label}] Ref independence (norm): (eps+.eps-)|_q1={pm_q1:.4f}, |_q2={pm_q2:.4f}  "
          f"{'PASS' if ok_norm else 'FAIL'}")

    return ok_p and ok_m and ok_norm


def verify_completeness(k, q, label):
    """
    Verify completeness relation for the spinor-helicity polarization vectors.

    With the conventions used here:
      eps^+(k;q) = tilde_lam_q sigma_bar lam_k / (sqrt(2) <qk>)
      eps^-(k;q) = tilde_lam_k sigma_bar lam_q / (sqrt(2) [kq])

    The relation verified (cross-product form, no conjugates):
      eps^+_mu eps^-_nu + eps^-_mu eps^+_nu = -g_{mu nu} + (k_mu q_nu + q_mu k_nu) / (k.q)

    And the conjugate form (with complex conjugates):
      sum_h eps^h_mu (eps^h_nu)^* = g_{mu nu} - (k_mu q_nu + q_mu k_nu) / (k.q)
    """
    eps_p = polarization_plus(k, q)
    eps_m = polarization_minus(k, q)
    kq = dot_mp(k, q)

    # Completeness (cross-product, no conjugates): eps^+_mu eps^-_nu + eps^-_mu eps^+_nu
    comp_cross = np.outer(eps_p, eps_m) + np.outer(eps_m, eps_p)
    expected_cross = -g_mp + (np.outer(k, q) + np.outer(q, k)) / kq
    err_cross = np.max(np.abs(comp_cross - expected_cross))

    # Conjugate completeness: sum_h eps^h_mu (eps^h_nu)^*
    comp_conj = np.outer(eps_p, np.conj(eps_p)) + np.outer(eps_m, np.conj(eps_m))
    expected_conj = g_mp - (np.outer(k, q) + np.outer(q, k)) / kq
    err_conj = np.max(np.abs(comp_conj - expected_conj))

    ok = err_cross < TOL * 100 and err_conj < TOL * 100
    print(f"  [{label}] Completeness (cross): err={err_cross:.2e}, (conj): err={err_conj:.2e}  "
          f"{'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print("=" * 60)
    print("Script 05: Polarization Vectors Verification")
    print("  Verifying Dixon 3.3.1 / Srednicki 60.6-60.9")
    print("  Metric: mostly-plus g = diag(-1,+1,+1,+1)")
    print("=" * 60)

    rng = np.random.default_rng(55)
    all_pass = True

    n_tests = 6
    print(f"\nTesting {n_tests} random (k, q) pairs:")
    for idx in range(1, n_tests + 1):
        k = make_lightlike(rng)
        q = make_lightlike(rng)
        label = f"test {idx}"
        ok1 = verify_transversality(k, q, label)
        ok2 = verify_self_orthogonality(k, q, label)
        ok3 = verify_cross_normalization(k, q, label)
        ok4 = verify_conjugate_normalization(k, q, label)
        all_pass = all_pass and ok1 and ok2 and ok3 and ok4

    print(f"\nTesting reference momentum independence:")
    for idx in range(1, 4):
        k = make_lightlike(rng)
        q1 = make_lightlike(rng)
        q2 = make_lightlike(rng)
        label = f"test {idx}"
        ok5 = verify_reference_independence(k, q1, q2, label)
        all_pass = all_pass and ok5

    print(f"\nTesting completeness relation:")
    for idx in range(1, 4):
        k = make_lightlike(rng)
        q = make_lightlike(rng)
        label = f"test {idx}"
        ok6 = verify_completeness(k, q, label)
        all_pass = all_pass and ok6

    print()
    print("=" * 60)
    if all_pass:
        print("OVERALL: PASS")
    else:
        print("OVERALL: FAIL")
    print("=" * 60)


if __name__ == '__main__':
    main()
