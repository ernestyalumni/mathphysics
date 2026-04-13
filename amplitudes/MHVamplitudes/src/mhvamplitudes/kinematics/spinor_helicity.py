"""
kinematics.spinor_helicity — spinor-helicity bracket operations.

Implements the basic spinor contractions in the SMGA frame (arXiv:2602.12176).

Frame conventions (SMGA eq. 2):
    |i> = lambda_i = (1, z_i)
    |i] = tilde_lambda_i = omega_i * (1, tilde_z_i)

Bracket definitions (SMGA eq. 4):
    <ij> = z_{ij} = z_i - z_j
    [ij] = omega_i * omega_j * tilde_z_{ij}

where z_{ij} = z_i - z_j and tilde_z_{ij} = tilde_z_i - tilde_z_j.

The sign function sg(x) = +1 if x > 0, -1 if x < 0, 0 if x = 0.

Convention comparison: SMGA (this code) vs Srednicki QFT textbook
=================================================================

SPACETIME SIGNATURE — the fundamental difference:
    - Srednicki:  Lorentzian metric eta = diag(-1,+1,+1,+1)  [mostly-plus]
                  (Srednicki eq. 34.1)
    - SMGA:       Klein (2,2) signature eta = diag(-1,+1,+1,-1)
                  (SMGA sec. 2, around eq. 1)
    This is not merely a sign convention — it changes the reality conditions
    on spinors. In Lorentzian (1,3) signature, lambda and tilde-lambda are
    complex conjugates of each other. In Klein (2,2) signature, both lambda
    and tilde-lambda are REAL and independent. The entire SMGA result — that
    single-minus amplitudes are nonzero — depends on this: in Lorentzian
    signature, single-minus tree amplitudes vanish identically by the Parke-
    Taylor/BGK selection rules.

LEVI-CIVITA TENSOR:
    - Srednicki:  epsilon^{12} = +1, epsilon_{12} = -1  (Srednicki eq. 34.5)
    - SMGA:       Same convention (implicit in eq. 4 contraction).
    - This code:  Same. The bracket [a,b] = a[1]*b[0] - a[0]*b[1]
                  corresponds to contracting with epsilon_{ab} where
                  epsilon_{12} = -1 (i.e. epsilon_{21} = +1).
    VERDICT: Epsilon tensor conventions AGREE.

SPINOR BRACKET DEFINITIONS:
    - Srednicki:  <p q> = epsilon^{ab} p_a q_b = p_1 q_2 - p_2 q_1
                  [p q] = epsilon_{dot{a}dot{b}} p^{dot{a}} q^{dot{b}}
                          = p^{dot{1}} q^{dot{2}} - p^{dot{2}} q^{dot{1}}
                  (Srednicki eqs. 34.20, 34.22)
    - SMGA:       <ij> = z_i - z_j, [ij] = omega_i * omega_j * (tz_i - tz_j)
                  Both are real-valued in (2,2) signature.
    - This code:  Both brackets use the contraction a[1]*b[0] - a[0]*b[1].
                  In the SMGA frame with |i> = (1, z_i), this gives
                  <ij> = z_i - z_j, matching SMGA eq. 4.
    VERDICT: Structurally identical contraction. The key difference is that
    Srednicki's brackets are complex (with [ij] = <ij>* for real momenta),
    while SMGA's are both real and independent.

MANDELSTAM VARIABLES:
    - Srednicki:  s_{ij} = -(p_i + p_j)^2 = <ij>[ji]  (Srednicki eq. 34.25)
                  Note the minus sign from mostly-plus metric: p^2 = -m^2.
    - SMGA:       s_{ij} = <ij>[ji] still holds formally, but the momentum
                  dot product p_i . p_j picks up different signs because of
                  the (2,2) metric. In the half-collinear limit the z_i are
                  all equal, so <ij> -> 0 and s_{ij} -> 0, but the ratio
                  s_{ij}/s_{kl} remains well-defined.
    This code does NOT compute Mandelstam variables; it works entirely with
    the sign function sg([ij]) which only depends on the real-valued bracket.

PARKE-TAYLOR FORMULA:
    - Srednicki:  A(1^-,2^-,...,n^+) = <12>^4 / (<12><23>...<n1>)
                  (Srednicki eq. 60.31, complex-valued)
    - SMGA:       PT_{1...n} = 1 / (<12><23>...<n1>)  (SMGA eq. 3)
                  The <12>^4 numerator is stripped off in the "stripped"
                  amplitude convention.
    The BG recursion in this code computes the STRIPPED amplitude, which is
    the reduced quantity after removing the MHV prefactor.

OVERALL SIGN OF AMPLITUDES:
    - Srednicki:  Uses the convention where color-ordered partial amplitudes
                  have a specific overall phase tied to the trace normalization
                  Tr(T^a T^b) = delta^{ab}/2 (Srednicki eq. 25.9).
    - SMGA:       The stripped amplitude A_{1...n} is defined up to the
                  overall sign set by the recursion seed A_3 = sg_{12} and
                  the recursive step. The sign is physical in (2,2) signature
                  because all quantities are real.
    This code follows the SMGA sign convention throughout.
"""

import numpy as np


def sg(x: float) -> int:
    """Sign function: sg(x) = +1 if x > 0, -1 if x < 0, 0 if x = 0."""
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0


def square_bracket(a: np.ndarray, b: np.ndarray) -> float:
    """
    Compute the square bracket [a, b] for two 2-component tilde-lambda vectors.

    SMGA convention: [a, b] = a[1]*b[0] - a[0]*b[1]

    This gives [ij] = omega_i * omega_j * (tilde_z_i - tilde_z_j) when
    a = omega_i * (1, tilde_z_i) and b = omega_j * (1, tilde_z_j).

    Srednicki comparison (eq. 34.22):
        [p q] = epsilon_{dot{a}dot{b}} p^dot{a} q^dot{b}
    With epsilon_{12} = -1 (same convention as here), Srednicki's formula
    gives [p q] = p^{dot{2}} q^{dot{1}} - p^{dot{1}} q^{dot{2}}, which
    matches our a[1]*b[0] - a[0]*b[1] under the identification
    component 0 <-> dotted index 1, component 1 <-> dotted index 2.

    KEY DIFFERENCE: In Srednicki (Lorentzian), [ij] = <ij>* for real momenta,
    so [ij] is complex. In SMGA (Klein), [ij] is real and independent of <ij>.
    This code uses real-valued tilde-lambdas throughout.
    """
    return a[1] * b[0] - a[0] * b[1]


def angle_bracket(lam_i: np.ndarray, lam_j: np.ndarray) -> float:
    """
    Compute the angle bracket <ij> for two 2-component lambda vectors.

    SMGA convention: <ij> = lam_i[1]*lam_j[0] - lam_i[0]*lam_j[1] = z_i - z_j
    when lam_i = (1, z_i) and lam_j = (1, z_j).

    Same contraction convention as square_bracket (epsilon_{12} = -1).

    Srednicki comparison (eq. 34.20):
        <p q> = epsilon^{ab} p_a q_b = p_1 q_2 - p_2 q_1
    With epsilon^{12} = +1 (same convention), this gives the same bilinear
    form. However, Srednicki's angle bracket is complex in Lorentzian
    signature, whereas in SMGA's Klein (2,2) signature this is real-valued.

    NOTE: This function is not used in the Berends-Giele recursion, which
    works entirely with square brackets [ij] of tilde-lambda vectors. It is
    provided for completeness and potential future use with full (non-half-
    collinear) kinematics.
    """
    return lam_i[1] * lam_j[0] - lam_i[0] * lam_j[1]


def construct_spinors(p: np.ndarray):
    """
    Construct spinor-helicity variables from a 4-momentum p = (p0, px, py, pz).

    Uses Srednicki's mostly-plus metric g = diag(-1,+1,+1,+1).
    Returns (lambda, tilde_lambda) as 2-component complex arrays.

    For generic (non-half-collinear) kinematics. Uses the standard
    decomposition via p0 + pz (Srednicki eq. 34.14 / Mangano-Parke).

    IMPORTANT: This function follows SREDNICKI / standard Lorentzian
    conventions, NOT the SMGA Klein (2,2) conventions used by the rest of
    this module. It is provided as a utility for working with generic 4-momenta
    (e.g., for verifying Parke-Taylor in Lorentzian signature, or for
    interfacing with Monte Carlo phase-space generators).

    In Srednicki's Lorentzian signature:
        - p is massless: p^mu p_mu = -p0^2 + px^2 + py^2 + pz^2 = 0
        - lambda and tilde_lambda are complex conjugates for real momenta
        - p_{a dot{a}} = lambda_a * tilde_lambda_{dot{a}}  (Srednicki eq. 34.14)
        - The decomposition uses p0 + pz as the reference component

    In SMGA's Klein (2,2) signature:
        - Both lambda and tilde_lambda are REAL
        - They parametrize as (1, z_i) and omega_i * (1, tilde_z_i)
        - This decomposition does not apply; use make_tilde_lambdas() instead

    This function is NOT used by the Berends-Giele recursion or SMGA formula.
    Those work entirely in the half-collinear frame with real tilde-lambda
    vectors constructed by phase_space.make_tilde_lambdas().
    """
    p0, px, py, pz = p
    r = complex(p0 + pz)
    if abs(r) < 1e-10:
        r = r + 1e-8
    sqrt_r = np.sqrt(r)
    lam = np.array([sqrt_r, (px + 1j * py) / sqrt_r], dtype=complex)
    lam_tilde = np.array([sqrt_r, (px - 1j * py) / sqrt_r], dtype=complex)
    return lam, lam_tilde
