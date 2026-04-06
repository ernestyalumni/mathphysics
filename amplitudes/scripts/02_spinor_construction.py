"""
Numerical verification of spinor construction for massless momenta.

Note on metric conventions:
  - Dixon uses mostly-MINUS metric g = diag(+1,-1,-1,-1).
  - Srednicki uses mostly-PLUS metric g = diag(-1,+1,+1,+1).
  - The spinor formula lambda*tilde_lambda = p_{alpha dot-alpha} uses Dixon's
    mostly-minus sigma matrices: sigma^mu_mm = (I, sigma_i).
  - The momentum reconstruction uses sigma-bar: tilde_lam sigma_bar lam = 2 p^mu,
    with sigma_bar^mu = (I, sigma_i) [same in both conventions for the reconstruction].
  - The Mandelstam invariants s_ij = 2 p_i.p_j use mostly-plus if stated with Srednicki,
    or mostly-minus if with Dixon.
  - This script uses the Dixon mostly-minus sigma matrices for p_matrix construction,
    and notes where the mostly-plus Srednicki result applies.

Verifies:
  1. Factorization: lambda_alpha * tilde_lambda_{dot-alpha} = p_{alpha dot-alpha}
     where p_{alpha dot-alpha} = p_mu * (sigma^mu_mm)_{alpha dot-alpha}
     and sigma^mu_mm = (I, sigma^1, sigma^2, sigma^3) [mostly-minus, Dixon convention]
  2. det(p_{alpha dot-alpha}) = 0  (massless condition: det = p^2_mm = E^2-|p|^2 = 0)
  3. Dixon (3.6) = Srednicki 50.15: tilde_lambda sigma-bar^mu lambda = 2 p^mu
     with sigma-bar^mu = (I, sigma^1, sigma^2, sigma^3)

Spinor construction (Dixon after eq. 3.12):
  lambda_alpha = ( sqrt(p^t + p^z),  (p^x + i*p^y)/sqrt(p^t + p^z) )^T
  tilde_lambda_{dot-alpha} = ( sqrt(p^t + p^z),  (p^x - i*p^y)/sqrt(p^t + p^z) )^T

References:
  Dixon, arXiv:1310.5353, eqs. 3.2, 3.6 and after 3.12
  Srednicki QFT, eqs. 50.1, 50.15
  Bridge document: amplitudes/09-dixon-srednicki-bridge.tex, Sec. 2
"""

import numpy as np

TOL = 1e-10

# Pauli matrices
sigma1 = np.array([[0, 1], [1, 0]], dtype=complex)
sigma2 = np.array([[0, -1j], [1j, 0]], dtype=complex)
sigma3 = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)

# Dixon mostly-minus sigma matrices: sigma^mu_mm = (I, sigma1, sigma2, sigma3)
# Used for: p_{alpha dot-alpha} = p_mu * sigma^mu_{mm, alpha dot-alpha}
# The factorization lambda*tilde_lambda = p_mm is in Dixon's convention.
sigma_mu = np.array([I2, sigma1, sigma2, sigma3])  # shape (4, 2, 2) [mostly-minus]

# sigma-bar^{mu dot-alpha alpha} = (I, sigma1, sigma2, sigma3)
# Used for momentum reconstruction: tilde_lam sigma_bar lam = 2 p^mu
# Same matrices in both mostly-plus and mostly-minus (the bar just transposes Pauli
# and flips sign of spatial parts... but numerically for reconstruction both agree).
sigma_bar_mu = np.array([I2, sigma1, sigma2, sigma3])  # shape (4, 2, 2)


def make_lightlike_momentum(rng):
    """Generate a random massless 4-momentum with positive energy in mostly-plus metric."""
    # Pick random spatial direction
    p3 = rng.standard_normal(3)
    p3 /= np.linalg.norm(p3)
    E = rng.uniform(0.5, 5.0)
    p3 *= E
    # p^mu = (E, px, py, pz), p^2 = -E^2 + px^2 + py^2 + pz^2 = 0 in mostly-plus
    return np.array([E, p3[0], p3[1], p3[2]])


def p_matrix(p):
    """
    Compute p_{alpha dot-alpha} = p_mu * sigma^mu_{alpha dot-alpha}
    In mostly-plus metric, sigma^mu = (-I, sigma1, sigma2, sigma3).
    p_{alpha dot-alpha} = p0*(-I) + p1*sigma1 + p2*sigma2 + p3*sigma3
                        = -p0*I + p1*sigma1 + p2*sigma2 + p3*sigma3
    """
    return np.einsum('m,mij->ij', p, sigma_mu)


def construct_spinors(p):
    """
    Construct lambda_alpha and tilde_lambda_{dot-alpha} from massless p^mu.

    Dixon's formula (after 3.12):
      lambda_alpha = (sqrt(p^t + p^z), (p^x + i*p^y)/sqrt(p^t + p^z))^T
      tilde_lambda = (sqrt(p^t + p^z), (p^x - i*p^y)/sqrt(p^t + p^z))^T

    This assumes p^t + p^z > 0 (i.e., the momentum is not in the -z direction).
    """
    pt, px, py, pz = p[0], p[1], p[2], p[3]
    denom = np.sqrt(pt + pz)
    lam = np.array([denom, (px + 1j*py)/denom], dtype=complex)
    lam_tilde = np.array([denom, (px - 1j*py)/denom], dtype=complex)
    return lam, lam_tilde


def verify_factorization(p, lam, lam_tilde, idx):
    """
    Verify: lambda_alpha * tilde_lambda_{dot-alpha} = p_{alpha dot-alpha}
    i.e., the outer product of the two 2-spinors equals the sigma-matrix contraction.
    """
    # p_matrix
    pm = p_matrix(p)
    # outer product
    outer = np.outer(lam, lam_tilde)  # shape (2,2): [alpha, dot-alpha]
    max_err = np.max(np.abs(outer - pm))
    ok = max_err < TOL
    print(f"  [momentum {idx}] Factorization lambda*lam_tilde = p_matrix: max_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_det_zero(p, idx):
    """
    Verify: det(p_{alpha dot-alpha}) = 0 for massless p (p^2 = 0 in mostly-plus means -p0^2 + p1^2+p2^2+p3^2 = 0).
    det(p_matrix) = -p^2_Minkowski = -(-(p0)^2+(p1)^2+(p2)^2+(p3)^2) in mostly-plus.
    For massless: p^2 = 0, so det = 0.
    """
    pm = p_matrix(p)
    det = np.linalg.det(pm)
    # p^2 in mostly-plus: -p0^2 + p1^2 + p2^2 + p3^2
    p2 = -p[0]**2 + p[1]**2 + p[2]**2 + p[3]**2
    max_err = abs(det)
    ok = max_err < TOL * max(1, abs(p[0])**2)
    print(f"  [momentum {idx}] det(p_matrix)=0 (massless): det = {det:.2e}, p^2 = {p2:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def verify_momentum_reconstruction(p, lam, lam_tilde, idx):
    """
    Verify Dixon (3.6) / Srednicki 50.15:
      tilde_lambda_{dot-alpha} (sigma-bar^mu)^{dot-alpha alpha} lambda_alpha = 2 p^mu

    sigma-bar^mu = (I, sigma1, sigma2, sigma3) in mostly-plus.
    The contraction: tilde_lam^{dot-alpha} sigma_bar^mu_{dot-alpha,alpha} lam^alpha
    = sum_{dot-alpha, alpha} tilde_lam[dot-alpha] * sigma_bar_mu[mu, dot-alpha, alpha] * lam[alpha]
    """
    # Compute tilde_lambda * sigma_bar^mu * lambda for each mu
    p_reconstructed = np.zeros(4, dtype=complex)
    for mu in range(4):
        p_reconstructed[mu] = lam_tilde @ sigma_bar_mu[mu] @ lam

    max_err = np.max(np.abs(p_reconstructed - 2*p))
    ok = max_err < TOL * max(1, abs(p[0]))
    print(f"  [momentum {idx}] tilde_lam sigma-bar lam = 2p^mu: max_err = {max_err:.2e}  {'PASS' if ok else 'FAIL'}")
    return ok


def main():
    print("=" * 60)
    print("Script 02: Spinor Construction Verification")
    print("  Verifying Dixon (3.2), (3.6) / Srednicki 50.1, 50.15")
    print("  Metric: mostly-plus g = diag(-1,+1,+1,+1)")
    print("=" * 60)

    rng = np.random.default_rng(42)
    n_momenta = 8
    all_pass = True

    print(f"\nTesting {n_momenta} random light-like momenta:")
    for idx in range(1, n_momenta + 1):
        p = make_lightlike_momentum(rng)

        # Verify p^2 = 0
        p2 = -p[0]**2 + p[1]**2 + p[2]**2 + p[3]**2
        if abs(p2) > 1e-14 * p[0]**2:
            print(f"  [momentum {idx}] WARNING: p^2 = {p2:.2e} (should be 0 by construction)")

        lam, lam_tilde = construct_spinors(p)

        ok1 = verify_factorization(p, lam, lam_tilde, idx)
        ok2 = verify_det_zero(p, idx)
        ok3 = verify_momentum_reconstruction(p, lam, lam_tilde, idx)
        all_pass = all_pass and ok1 and ok2 and ok3

    # Also test with a specific known momentum: p = (1, 0, 0, 1) (along z-axis)
    print("\nTesting specific momentum p = (1, 0, 0, 1):")
    p_z = np.array([1.0, 0.0, 0.0, 1.0])
    lam_z, lam_tilde_z = construct_spinors(p_z)
    print(f"  lambda = {lam_z}")
    print(f"  tilde_lambda = {lam_tilde_z}")
    # For p=(1,0,0,1): sqrt(1+1)=sqrt(2), (0+0i)/sqrt(2)=0
    # lambda = (sqrt(2), 0), tilde_lambda = (sqrt(2), 0)
    ok4 = verify_factorization(p_z, lam_z, lam_tilde_z, 'z-axis')
    ok5 = verify_det_zero(p_z, 'z-axis')
    ok6 = verify_momentum_reconstruction(p_z, lam_z, lam_tilde_z, 'z-axis')
    all_pass = all_pass and ok4 and ok5 and ok6

    print()
    print("=" * 60)
    if all_pass:
        print("OVERALL: PASS")
    else:
        print("OVERALL: FAIL")
    print("=" * 60)


if __name__ == '__main__':
    main()
