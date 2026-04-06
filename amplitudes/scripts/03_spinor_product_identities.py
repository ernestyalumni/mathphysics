"""
Numerical verification of spinor product identities.

Uses the all-outgoing convention: sum_i k_i = 0, where some momenta have negative energy
(these correspond to incoming particles in the physical picture).

Verifies for n=4,5 massless momenta satisfying momentum conservation:
  1. Anti-symmetry: <ij> = -<ji>, [ij] = -[ji]     (Dixon 3.13)
  2. Squaring: <ij>[ji] = s_ij = 2 k_i . k_j         (Dixon 3.14 / Srednicki 60.4)
     where s_ij = -2*k_i[0]*k_j[0] + 2*(k_i.k_j spatially) in mostly-plus metric
  3. Momentum conservation: sum_j <ij>[jk] = 0        (Dixon 3.15 / Srednicki 60.18)
  4. Schouten identity: <ij><kl> + <ik><lj> + <il><jk> = 0  (Dixon 3.16 / Srednicki 60.17)

Spinor conventions:
  angle bracket: <ij> = lambda_i[0]*lambda_j[1] - lambda_i[1]*lambda_j[0]
  square bracket: [ij] = tilde_lambda_i[0]*tilde_lambda_j[1] - tilde_lambda_i[1]*tilde_lambda_j[0]

Spinor construction (works for complex momenta via analytic continuation):
  For massless p^mu = (E, px, py, pz) with E > 0 and p0 + p3 > 0:
    lambda_alpha = (sqrt(p0+p3), (p1+i*p2)/sqrt(p0+p3))
    tilde_lambda = (sqrt(p0+p3), (p1-i*p2)/sqrt(p0+p3))
  For negative-energy momenta (incoming particles) p -> -p, then use the formula.

Note: the squaring identity uses 2 k_i . k_j (with mostly-plus metric = -E1E2 + p1.p2)
so for back-to-back pairs k and -k: s_{k,-k} = -2*(k.(-k)) = 2*(E^2 - |p|^2) = 0 (massless).

References:
  Dixon, arXiv:1310.5353, eqs. 3.13, 3.14, 3.15, 3.16
  Srednicki QFT, eqs. 50.16, 50.17, 60.4, 60.17, 60.18
  Bridge document: amplitudes/09-dixon-srednicki-bridge.tex, Sec. 3
"""

import numpy as np

TOL = 1e-8


def make_lightlike_pos_energy(rng):
    """Generate a random massless 4-momentum with positive energy and p0+p3 > 0."""
    while True:
        p3 = rng.standard_normal(3)
        p3 /= np.linalg.norm(p3)
        E = rng.uniform(0.5, 3.0)
        p = np.array([E, E*p3[0], E*p3[1], E*p3[2]])
        if p[0] + p[3] > 0.1 * E:  # ensure p0+p3 > 0
            return p


def make_conserved_momenta_n(n, rng):
    """
    Generate n massless momenta satisfying sum_i k_i = 0.

    For even n: use n/2 back-to-back pairs (k, -k). Note that half the momenta
    will have negative energy. This is the all-outgoing convention where:
    - positive energy = outgoing
    - negative energy = incoming (but treated as outgoing with reversed momentum)

    The spinor formula lambda*tilde_lambda = p_{matrix} uses the absolute value
    and analytic continuation: for k with p0+pz > 0 use formula directly;
    for -k set lambda(-k) = i*lambda(k), tilde_lambda(-k) = i*tilde_lambda(k)
    (overall phase doesn't affect amplitudes, but brackets do get factors of -1).

    Actually: we'll use positive energy momenta only, constructed as ALL-OUTGOING
    by a different method: generate n/2 pairs where pair (k_in, k_out) has
    k_in = k_out for one pair and opposite spatial for others.

    Simplest valid approach: generate n-1 positive-energy massless momenta,
    set the n-th = -(sum of first n-1). Check if it's massless -- it won't be.
    Instead: use RAMBO-like: generate n massless momenta, shift to CM frame.

    Cleanest approach: use n/2 random boosts so each pair sums to zero.
    Actually the cleanest approach for testing the identities is:
    Use 2->2 kinematics for n=4: k1+k2 = k3+k4 where k1,k2 incoming,
    k3,k4 outgoing. Rewrite as k3+k4+(-k1)+(-k2) = 0 (all-outgoing).
    All of k3,k4,(-k1),(-k2) are massless. The negative-energy momenta
    -k1, -k2 have p0 < 0.

    For the spinor formula, it works for any light-like momentum including
    negative energy ones, by using complex spinors. The key is that the
    spinor brackets are well-defined as long as p0+pz != 0.
    """
    if n == 4:
        # Use 4 massless momenta constructed to sum to zero.
        # Method: two pairs of back-to-back momenta in different directions.
        # All momenta are generic (no p0+p3=0 degeneracy).
        # Pair 1: k1 and k2 = -k1
        # Pair 2: k3 and k4 = -k3
        # But -k_i has negative energy. Use complex spinors via sqrt(p0+p3).
        # To avoid the p0+pz=0 issue, tilt all momenta.
        E = rng.uniform(1.0, 5.0)
        # Pick two random non-degenerate light-like momenta
        theta1 = rng.uniform(0.2, np.pi - 0.2)
        phi1 = rng.uniform(0.3, 2*np.pi - 0.3)
        theta2 = rng.uniform(0.2, np.pi - 0.2)
        phi2 = rng.uniform(0.3, 2*np.pi - 0.3)
        E2 = rng.uniform(0.5*E, 1.5*E)
        st1, ct1 = np.sin(theta1), np.cos(theta1)
        sp1, cp1 = np.sin(phi1), np.cos(phi1)
        st2, ct2 = np.sin(theta2), np.cos(theta2)
        sp2, cp2 = np.sin(phi2), np.cos(phi2)
        k1 = np.array([E,  E*st1*cp1,  E*st1*sp1,  E*ct1])
        k2 = np.array([E2, E2*st2*cp2, E2*st2*sp2, E2*ct2])
        # k3 and k4 must be massless and k3+k4 = -(k1+k2)
        P = k1 + k2  # timelike, positive energy
        P2 = -P[0]**2 + P[1]**2 + P[2]**2 + P[3]**2
        M = np.sqrt(-P2)
        # Decay in rest frame of P
        theta_d = rng.uniform(0.2, np.pi - 0.2)
        phi_d = rng.uniform(0, 2*np.pi)
        Ehalf = M / 2
        sp_d, cp_d = np.sin(theta_d), np.cos(theta_d)
        sc_d, cc_d = np.sin(phi_d), np.cos(phi_d)
        k3_rf = Ehalf * np.array([1, sp_d*cc_d, sp_d*sc_d, cp_d])
        k4_rf = Ehalf * np.array([1, -sp_d*cc_d, -sp_d*sc_d, -cp_d])
        # Boost to lab frame
        beta_vec = np.array([P[1], P[2], P[3]]) / P[0]
        beta_mag = np.linalg.norm(beta_vec)
        if beta_mag < 1e-10:
            k3, k4 = k3_rf, k4_rf
        else:
            gamma = P[0] / M
            beta_hat = beta_vec / beta_mag
            def boost_v(k_rf, gamma, beta_hat, beta_mag):
                E_r = k_rf[0]
                p_r = k_rf[1:]
                p_par = np.dot(p_r, beta_hat)
                p_perp = p_r - p_par * beta_hat
                E_lab = gamma * (E_r + beta_mag * p_par)
                p_par_lab = gamma * (p_par + beta_mag * E_r)
                p_lab = p_perp + p_par_lab * beta_hat
                return np.array([E_lab, p_lab[0], p_lab[1], p_lab[2]])
            k3 = boost_v(k3_rf, gamma, beta_hat, beta_mag)
            k4 = boost_v(k4_rf, gamma, beta_hat, beta_mag)
        # All-outgoing: -k1, -k2, k3, k4 sum to zero (since k3+k4=k1+k2)
        # But -k1,-k2 have negative energy. That's fine for complex spinors.
        momenta = [-k1, -k2, k3, k4]
        return momenta

    elif n == 5:
        # Use 3->2 kinematics: 3 incoming, 2 outgoing
        # k1, k2, k3 incoming; k4, k5 outgoing with k4+k5 = k1+k2+k3
        E = rng.uniform(1.0, 3.0)
        # k1, k2, k3: three different light-like momenta
        phi1 = rng.uniform(0, 2*np.pi)
        phi2 = rng.uniform(0, 2*np.pi)
        phi3 = rng.uniform(0, 2*np.pi)
        theta1 = rng.uniform(0.2, np.pi-0.2)
        theta2 = rng.uniform(0.2, np.pi-0.2)
        theta3 = rng.uniform(0.2, np.pi-0.2)
        k1 = E * np.array([1, np.sin(theta1)*np.cos(phi1), np.sin(theta1)*np.sin(phi1), np.cos(theta1)])
        E2 = rng.uniform(0.5*E, 1.5*E)
        k2 = E2 * np.array([1, np.sin(theta2)*np.cos(phi2), np.sin(theta2)*np.sin(phi2), np.cos(theta2)])
        E3 = rng.uniform(0.5*E, 1.5*E)
        k3 = E3 * np.array([1, np.sin(theta3)*np.cos(phi3), np.sin(theta3)*np.sin(phi3), np.cos(theta3)])
        # P = k1+k2+k3 is timelike (positive energy)
        P = k1 + k2 + k3
        # Split P into two light-like momenta k4, k5
        # In P rest frame: M = sqrt(-P^2 in mostly-plus) = sqrt(P0^2 - |P|^2 - using mostly-minus!)
        # Actually in mostly-plus: P^2 = -P0^2 + P1^2 + P2^2 + P3^2 < 0 (timelike incoming)
        P2 = -P[0]**2 + P[1]**2 + P[2]**2 + P[3]**2
        if P2 >= 0:
            # P is spacelike (unlikely for small E ratios). Use default.
            P2 = -4 * E**2  # force timelike

        M = np.sqrt(-P2)  # invariant mass
        # Random decay angle in P rest frame
        psi = rng.uniform(0.1, np.pi-0.1)
        chi = rng.uniform(0, 2*np.pi)
        Ehalf = M / 2
        sp_h, cp_h = np.sin(psi), np.cos(psi)
        sc_h, cc_h = np.sin(chi), np.cos(chi)
        k4_rf = Ehalf * np.array([1, sp_h*cc_h, sp_h*sc_h, cp_h])
        k5_rf = Ehalf * np.array([1, -sp_h*cc_h, -sp_h*sc_h, -cp_h])
        # Boost k4_rf, k5_rf to lab frame (frame where P has its lab-frame value)
        beta_vec = np.array([P[1], P[2], P[3]]) / P[0]
        beta_mag = np.linalg.norm(beta_vec)
        if beta_mag < 1e-10:
            k4 = k4_rf
            k5 = k5_rf
        else:
            gamma = P[0] / M
            beta_hat = beta_vec / beta_mag
            def boost(k_rf):
                E_r = k_rf[0]
                p_r = k_rf[1:]
                p_par = np.dot(p_r, beta_hat)
                p_perp = p_r - p_par * beta_hat
                E_lab = gamma * (E_r + beta_mag * p_par)
                p_par_lab = gamma * (p_par + beta_mag * E_r)
                p_lab = p_perp + p_par_lab * beta_hat
                return np.array([E_lab, p_lab[0], p_lab[1], p_lab[2]])
            k4 = boost(k4_rf)
            k5 = boost(k5_rf)
        # All-outgoing: -k1, -k2, -k3, k4, k5 (3 negative, 2 positive)
        return [-k1, -k2, -k3, k4, k5]
    else:
        raise ValueError(f"n={n} not supported")


def construct_spinors(p):
    """
    Construct (lambda, tilde_lambda) for massless p.
    For p with p0+p3 > 0:
      lambda = (sqrt(p0+p3), (p1+ip2)/sqrt(p0+p3))
      tilde_lambda = (sqrt(p0+p3), (p1-ip2)/sqrt(p0+p3))
    For p with p0+p3 < 0 (negative-energy or negative-z momenta):
      Use |p0+p3| and complex square root. The spinors will be imaginary but
      brackets are still well-defined.
    We use analytic continuation: allow p0+p3 to be any sign.
    """
    p0, px, py, pz = p
    r = p0 + pz
    if abs(r) < 1e-10:
        # Degenerate case: use p0+p3 -> add tiny perturbation
        r = p0 + pz + 1e-8
    sqrt_r = np.sqrt(complex(r))  # complex sqrt handles negative r
    lam = np.array([sqrt_r, (px + 1j*py)/sqrt_r], dtype=complex)
    lam_tilde = np.array([sqrt_r, (px - 1j*py)/sqrt_r], dtype=complex)
    return lam, lam_tilde


def angle(lam_i, lam_j):
    return lam_i[0]*lam_j[1] - lam_i[1]*lam_j[0]


def square(lam_ti, lam_tj):
    return lam_ti[0]*lam_tj[1] - lam_ti[1]*lam_tj[0]


def dot_product_mostly_plus(p, q):
    """Minkowski dot product in mostly-plus: -p0q0 + p1q1 + p2q2 + p3q3"""
    return -p[0]*q[0] + p[1]*q[1] + p[2]*q[2] + p[3]*q[3]


def verify_antisymmetry(momenta, lams, lam_tildes, label):
    n = len(momenta)
    errors_angle, errors_square = [], []
    for i in range(n):
        for j in range(i+1, n):
            aij = angle(lams[i], lams[j])
            aji = angle(lams[j], lams[i])
            errors_angle.append(abs(aij + aji))
            sij = square(lam_tildes[i], lam_tildes[j])
            sji = square(lam_tildes[j], lam_tildes[i])
            errors_square.append(abs(sij + sji))

    max_err = max(max(errors_angle), max(errors_square))
    ok = max_err < TOL
    print(f"  [{label}] Anti-symmetry: max_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_squaring(momenta, lams, lam_tildes, label):
    """Verify <ij>[ji] = s_ij = 2 k_i.k_j."""
    n = len(momenta)
    errors = []
    scale = max(abs(dot_product_mostly_plus(momenta[i], momenta[j]))
                for i in range(n) for j in range(i+1, n) if i != j) + 1
    for i in range(n):
        for j in range(i+1, n):
            aij = angle(lams[i], lams[j])
            sji = square(lam_tildes[j], lam_tildes[i])   # [ji]
            product = aij * sji
            s_ij = 2 * dot_product_mostly_plus(momenta[i], momenta[j])
            errors.append(abs(product - s_ij) / (abs(s_ij) + 1e-10))

    max_err = max(errors)
    ok = max_err < TOL * 100  # relative tolerance
    print(f"  [{label}] Squaring <ij>[ji]=s_ij: max_rel_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_momentum_conservation_sum(momenta, lams, lam_tildes, label):
    """Verify sum_j <ij>[jk] = 0 for all i, k."""
    n = len(momenta)
    errors = []
    for i in range(n):
        for k in range(n):
            total = sum(
                angle(lams[i], lams[j]) * square(lam_tildes[j], lam_tildes[k])
                for j in range(n)
            )
            errors.append(abs(total))

    max_err = max(errors)
    ok = max_err < TOL
    print(f"  [{label}] Mom. conservation sum_j<ij>[jk]=0: max_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_schouten(momenta, lams, lam_tildes, label):
    """Verify <ij><kl> + <ik><lj> + <il><jk> = 0 for distinct i,j,k,l."""
    n = len(momenta)
    errors = []
    for i in range(n):
        for j in range(n):
            if j == i: continue
            for k in range(n):
                if k == i or k == j: continue
                for l in range(n):
                    if l in (i, j, k): continue
                    val = (angle(lams[i], lams[j]) * angle(lams[k], lams[l])
                         + angle(lams[i], lams[k]) * angle(lams[l], lams[j])
                         + angle(lams[i], lams[l]) * angle(lams[j], lams[k]))
                    errors.append(abs(val))

    if not errors:
        print(f"  [{label}] Schouten: not enough particles  SKIP")
        return True

    max_err = max(errors)
    ok = max_err < TOL
    print(f"  [{label}] Schouten identity: max_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def run_n_point(n, rng, trial):
    label = f"n={n}, trial={trial}"
    momenta = make_conserved_momenta_n(n, rng)
    lams, lam_tildes = [], []
    for p in momenta:
        lam, lam_tilde = construct_spinors(p)
        lams.append(lam)
        lam_tildes.append(lam_tilde)

    sum_p = sum(momenta)
    max_p_err = np.max(np.abs(sum_p))
    if max_p_err > 1e-10:
        print(f"  [{label}] WARNING: momentum conservation violated: {max_p_err:.2e}")

    ok1 = verify_antisymmetry(momenta, lams, lam_tildes, label)
    ok2 = verify_squaring(momenta, lams, lam_tildes, label)
    ok3 = verify_momentum_conservation_sum(momenta, lams, lam_tildes, label)
    ok4 = verify_schouten(momenta, lams, lam_tildes, label)
    return ok1 and ok2 and ok3 and ok4


def main():
    print("=" * 60)
    print("Script 03: Spinor Product Identities Verification")
    print("  Verifying Dixon (3.13), (3.14), (3.15), (3.16)")
    print("  Mostly-plus metric g = diag(-1,+1,+1,+1)")
    print("  All-outgoing convention: sum k_i = 0")
    print("=" * 60)

    rng = np.random.default_rng(123)
    all_pass = True

    print("\n--- n=4 ---")
    for trial in range(1, 5):
        ok = run_n_point(4, rng, trial)
        all_pass = all_pass and ok

    print("\n--- n=5 ---")
    for trial in range(1, 5):
        ok = run_n_point(5, rng, trial)
        all_pass = all_pass and ok

    print()
    print("=" * 60)
    if all_pass:
        print("OVERALL: PASS")
    else:
        print("OVERALL: FAIL")
    print("=" * 60)


if __name__ == '__main__':
    main()
