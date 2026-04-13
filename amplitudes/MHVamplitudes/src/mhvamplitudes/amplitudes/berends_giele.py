"""
amplitudes.berends_giele — Berends-Giele recursion for single-minus amplitudes.

Implements the full recursion relation from arXiv:2602.12176 (SMGA paper)
for computing stripped single-minus gluon tree amplitudes A_{1...n}.

The recursion has three layers:
1. Vertex functions V and V-bar (SMGA eq. 7)
2. Preamplitudes A-bar_S (SMGA eq. 5-6)
3. Stripped amplitude A_{1...n} (SMGA eq. 8 / eq. A17)

This implementation works numerically with explicit tilde-lambda vectors
in the half-collinear frame. No symbolic algebra is needed.

Key equations (SMGA paper notation):
    V_{lam1...lamn} = prod_{k=1}^{n-1} sg_{k,k+1} * Theta(-[lam_{1..k}, lam_{k+1..n}] / [lam_k, lam_{k+1}])
    V-bar: same but Theta(+ratio)
    PT-hat = V - V-bar
    A-bar_q = 1,  A-bar_{qp} = 0
    A-bar_{q..p} = -sum_{o.p., A>=3} V_{lam_{S1}..lam_{SA}} * prod A-bar_{Sa}
    A_{1..n} = -sum_{o.p. of (1..n-1), A>=1} PT-hat_{lam_{S1}..lam_{SA}} * prod A-bar_{Sa}

Convention: the LAST particle in the label list is minus-helicity.
To compute A_{1,...,n} with particle 1 minus, use cyclicity:
    A_{1,...,n} = A_{2,...,n,1}
and call stripped_amplitude((2,...,n,1)).

Convention comparison: SMGA Berends-Giele vs Srednicki / standard BG
=====================================================================

WHAT IS "BERENDS-GIELE RECURSION"?

The name "Berends-Giele" refers to a general technique (Berends & Giele,
Nucl. Phys. B306 (1988) 759) for computing tree-level gluon amplitudes
by recursively building off-shell currents. The SMGA paper adapts this
idea to their specific context but with significant differences from the
standard Lorentzian version.

STANDARD BG (Srednicki / Lorentzian):
    - Constructs off-shell currents J^{mu}(1,2,...,k) recursively
    - Uses the Yang-Mills 3-gluon and 4-gluon vertices from the Lagrangian
    - Propagators are i/(P^2 + i*epsilon) with P = sum of momenta
    - The amplitude is obtained by amputating the final propagator
    - Srednicki discusses Yang-Mills vertices in Ch. 72 (eq. 72.7) with
      structure constants f^{abc} normalized via [T^a, T^b] = i*f^{abc}*T^c
    - Srednicki's color-ordering strips f^{abc} down to partial amplitudes
      (Ch. 60), and BG recursion can be formulated for these partial amplitudes

SMGA BG (this code / Klein (2,2)):
    - Does NOT use Feynman-rule vertices from the Lagrangian directly
    - Instead uses "vertex functions" V and V-bar defined combinatorially
      via products of sign functions and Heaviside step functions
    - The key objects are:
        V: product of sg_{k,k+1} * Theta(-ratio)  [SMGA eq. 7]
        V-bar: product of sg_{k,k+1} * Theta(+ratio)
        PT-hat = V - V-bar  (on-shell Parke-Taylor, stripped)
    - These vertex functions encode the piecewise-constant structure of
      amplitudes in (2,2) signature — they are NOT smooth functions of momenta
    - Propagators are absent because the amplitudes are piecewise constant:
      1/P^2 factors cancel against numerator factors in the half-collinear limit
    - The recursion produces the STRIPPED amplitude directly, not an off-shell current

KEY STRUCTURAL DIFFERENCES:

1. VERTEX FUNCTIONS:
    - Srednicki: 3-gluon vertex ~ f^{abc}[g^{mu nu}(k1-k2)^rho + cyclic]
      (eq. 72.7). 4-gluon vertex ~ f^{abe}f^{cde}[g^{mu rho}g^{nu sigma} - ...]
      (eq. 72.8). These are smooth (polynomial) functions of momenta.
    - SMGA: V and V-bar involve step functions Theta and sign functions sg.
      They are piecewise constant with discontinuities at chamber walls.
      There is no 4-gluon vertex analogue; the recursion structure is different.

2. PROPAGATORS:
    - Srednicki: Gluon propagator in Feynman gauge:
      i*delta^{ab}*(-g^{mu nu}) / (k^2 + i*epsilon)  (eq. 72.6)
    - SMGA: No propagators appear. The half-collinear limit is a degenerate
      kinematic regime where s_{ij} -> 0 but ratios remain finite. The
      piecewise-constant result emerges after all propagator poles cancel.

3. COLOR STRUCTURE:
    - Srednicki: Full amplitudes involve color factors f^{abc}; partial
      (color-ordered) amplitudes are obtained after color decomposition.
      Structure constants use [T^a, T^b] = i*f^{abc}*T^c.
    - SMGA: Works directly with color-ordered (partial) amplitudes.
      The recursion is formulated for the stripped amplitude, which is
      the partial amplitude divided by the Parke-Taylor denominator.

4. CYCLICITY AND MINUS-PARTICLE PLACEMENT:
    - Srednicki: Color-ordered amplitudes A(1,2,...,n) are invariant under
      cyclic permutations. The minus-helicity particles can be anywhere.
    - SMGA: Same cyclic symmetry. This code places the minus-helicity
      particle LAST in the label list (the recursion partitions the
      preceding plus-helicity particles). To compute A_{1^-,2^+,...,n^+},
      we use cyclicity to get A_{2^+,...,n^+,1^-}.

5. GAUGE CHOICE / REFERENCE SPINOR:
    - Srednicki: MHV computations often use a reference spinor |q] for
      plus-helicity polarization vectors (Ch. 60). Results are q-independent.
    - SMGA: No explicit gauge choice. The half-collinear frame is a
      specific kinematic limit, not a gauge choice. The vertex functions V
      and V-bar implicitly encode the gauge structure through the Parke-Taylor.

6. REALITY CONDITIONS:
    - Srednicki: In Lorentzian signature, spinor brackets are complex.
      The recursion produces complex amplitudes. For single-minus
      configurations, the result is zero.
    - SMGA: In Klein (2,2) signature, all brackets are real. The recursion
      produces real, piecewise-constant amplitudes taking values in Z/2^{n-2}.
      The sign function sg() and step function Theta() are well-defined
      precisely because the brackets are real.
"""

import numpy as np
from itertools import combinations

from ..kinematics.spinor_helicity import sg, square_bracket


def ordered_partitions(seq, min_parts=1, max_parts=None):
    """
    Generate all ordered partitions of a sequence into consecutive blocks.

    An ordered partition of (a_1, ..., a_n) into k parts is a choice of
    k-1 cut points among the n-1 gaps, producing k consecutive sub-tuples.

    Parameters
    ----------
    seq : tuple
        The sequence to partition.
    min_parts : int
        Minimum number of parts (default 1).
    max_parts : int or None
        Maximum number of parts (default len(seq)).

    Yields
    ------
    tuple of tuples
        Each yielded value is a partition like ((a,b), (c,), (d,e,f)).
    """
    n = len(seq)
    if max_parts is None:
        max_parts = n
    max_parts = min(max_parts, n)

    for k in range(min_parts, max_parts + 1):
        if k == 1:
            yield (tuple(seq),)
            continue
        for cuts in combinations(range(1, n), k - 1):
            parts = []
            prev = 0
            for c in cuts:
                parts.append(tuple(seq[prev:c]))
                prev = c
            parts.append(tuple(seq[prev:]))
            yield tuple(parts)


class BerendsGieleRecursion:
    """
    Numerical computation of single-minus stripped amplitudes via
    the Berends-Giele recursion from arXiv:2602.12176.

    Parameters
    ----------
    tilde_lambdas : dict
        Maps particle label -> 2-component numpy array (tilde-lambda vector).
        In the SMGA frame: tilde_lambda_i = omega_i * (1, tilde_z_i).
    """

    def __init__(self, tilde_lambdas):
        self._tl = dict(tilde_lambdas)
        self._preamplitude_cache = {}

    def partial_sum(self, labels):
        """Compute tilde_lambda_S = sum_{i in S} tilde_lambda_i."""
        result = np.zeros(2)
        for i in labels:
            result = result + self._tl[i]
        return result

    def vertex_V(self, vectors):
        """
        Compute the vertex function V (SMGA eq. 7).

        V_{lam_1...lam_n} = prod_{k=1}^{n-1} sg_{k,k+1}
                            * Theta(-[lam_{1..k}, lam_{k+1..n}] / [lam_k, lam_{k+1}])

        Parameters
        ----------
        vectors : list of np.ndarray
            List of 2-component tilde-lambda vectors (can be block momenta).

        Returns
        -------
        float
            The vertex value (0 or +/-1).

        Srednicki comparison:
            This has no direct analogue in Srednicki. The closest object is
            the color-ordered 3-gluon vertex (Srednicki eq. 72.7, stripped of
            color), but V is defined for an ARBITRARY number of vectors, not
            just 3. V is piecewise constant (0 or +/-1) — it is not a smooth
            function of momenta. The Theta function and sg function are only
            meaningful because all brackets are real in (2,2) signature.
        """
        n = len(vectors)
        if n <= 1:
            return 1.0

        result = 1.0
        for k in range(n - 1):
            brk_adj = square_bracket(vectors[k], vectors[k + 1])
            sg_adj = sg(brk_adj)
            if sg_adj == 0:
                return 0.0

            left_sum = sum(vectors[:k + 1])
            right_sum = sum(vectors[k + 1:])
            brk_lr = square_bracket(left_sum, right_sum)

            if abs(brk_adj) < 1e-15:
                return 0.0
            ratio = brk_lr / brk_adj

            # Theta(-ratio): 1 if ratio < 0, else 0
            theta_val = 1.0 if ratio < -1e-15 else 0.0

            result *= sg_adj * theta_val
            if result == 0.0:
                return 0.0

        return result

    def vertex_Vbar(self, vectors):
        """
        Compute V-bar: same as V but with Theta(+ratio).

        V-bar_{lam_1...lam_n} = prod_{k=1}^{n-1} sg_{k,k+1}
                                * Theta(+[lam_{1..k}, lam_{k+1..n}] / [lam_k, lam_{k+1}])
        """
        n = len(vectors)
        if n <= 1:
            return 1.0

        result = 1.0
        for k in range(n - 1):
            brk_adj = square_bracket(vectors[k], vectors[k + 1])
            sg_adj = sg(brk_adj)
            if sg_adj == 0:
                return 0.0

            left_sum = sum(vectors[:k + 1])
            right_sum = sum(vectors[k + 1:])
            brk_lr = square_bracket(left_sum, right_sum)

            if abs(brk_adj) < 1e-15:
                return 0.0
            ratio = brk_lr / brk_adj

            # Theta(+ratio): 1 if ratio > 0, else 0
            theta_val = 1.0 if ratio > 1e-15 else 0.0

            result *= sg_adj * theta_val
            if result == 0.0:
                return 0.0

        return result

    def pt_hat(self, vectors):
        """
        Compute the on-shell Parke-Taylor factor PT-hat = V - V-bar.

        This is the LSZ-reduced Parke-Taylor factor that appears in the
        recursion for the stripped amplitude (SMGA below eq. 8).

        Srednicki comparison:
            The Parke-Taylor formula in Srednicki is (eq. 60.31):
                A(i^-,j^-,...) = <ij>^4 / (<12><23>...<n1>)
            which is a smooth, complex-valued function of the kinematics.
            PT-hat here is the STRIPPED, ON-SHELL, PIECEWISE-CONSTANT version
            in (2,2) signature: it takes values in {-1, 0, +1}. The
            relationship is: the full Parke-Taylor = PT-hat * product of
            angle brackets (which diverge in the half-collinear limit), and
            the stripped amplitude = full amplitude / full Parke-Taylor.
        """
        return self.vertex_V(vectors) - self.vertex_Vbar(vectors)

    def preamplitude(self, S):
        """
        Compute the preamplitude A-bar_S (SMGA eq. 5-6).

        A-bar_q = 1            (single particle)
        A-bar_{qp} = 0         (two particles)
        A-bar_{q..p} = -sum_{o.p., A>=3} V_{lam_{S1}..lam_{SA}} * prod A-bar_{Sa}

        Parameters
        ----------
        S : tuple
            Ordered tuple of particle labels.

        Returns
        -------
        float
            The preamplitude value.

        Srednicki comparison:
            The preamplitude A-bar_S is the SMGA analogue of the Berends-Giele
            off-shell current J^mu(q,...,p) in Srednicki/standard BG. However:
            - Standard BG current J carries a Lorentz index mu and is a
              complex-valued vector function of momenta
            - A-bar_S is a scalar, real-valued, piecewise-constant function
            - The base cases differ: Srednicki's J^mu(q) = epsilon^mu(q) (the
              polarization vector), while A-bar_q = 1 (dimensionless)
            - A-bar_{qp} = 0 because in the half-collinear limit, the 2-particle
              off-shell current has its propagator pole cancelled. In standard BG,
              J^mu(q,p) is generically nonzero.
            - The recursion sums over partitions into >= 3 blocks (not >= 2 as
              in standard BG), because the "2-block" contribution is absorbed
              into the definition of V.
        """
        S = tuple(S)
        if S in self._preamplitude_cache:
            return self._preamplitude_cache[S]

        n = len(S)
        if n == 1:
            return 1.0
        if n == 2:
            return 0.0

        result = 0.0
        for partition in ordered_partitions(S, min_parts=3):
            # Compute block momenta
            block_vectors = [self.partial_sum(block) for block in partition]

            v = self.vertex_V(block_vectors)
            if v == 0.0:
                continue

            prod_abar = 1.0
            for block in partition:
                ab = self.preamplitude(block)
                if ab == 0.0:
                    prod_abar = 0.0
                    break
                prod_abar *= ab

            if prod_abar == 0.0:
                continue

            result -= v * prod_abar

        self._preamplitude_cache[S] = result
        return result

    def stripped_amplitude(self, labels):
        """
        Compute the stripped amplitude A_{labels} via the full recursion.

        The LAST label in the list is the minus-helicity particle.
        The recursion (SMGA eq. A17) partitions all labels EXCEPT the last:

            A_{1..n} = -sum_{(1..n-1)=S1|..|Sk, k>=1}
                       PT-hat_{lam_{S1}..lam_{Sk}} * prod A-bar_{Sa}

        To compute A_{1,2,...,n} with particle 1 as minus-helicity,
        use cyclicity A_{1,...,n} = A_{2,...,n,1} and call:
            stripped_amplitude((2, 3, ..., n, 1))

        Parameters
        ----------
        labels : tuple
            Ordered particle labels. Last element is minus-helicity.

        Returns
        -------
        float
            The stripped amplitude value.
        """
        labels = tuple(labels)
        n = len(labels)
        assert n >= 3

        # Partition the first n-1 labels (the plus-helicity particles)
        to_partition = labels[:-1]

        result = 0.0
        for partition in ordered_partitions(to_partition, min_parts=1):
            block_vectors = [self.partial_sum(block) for block in partition]

            pt = self.pt_hat(block_vectors)
            if pt == 0.0:
                continue

            prod_abar = 1.0
            for block in partition:
                ab = self.preamplitude(block)
                if ab == 0.0:
                    prod_abar = 0.0
                    break
                prod_abar *= ab

            if prod_abar == 0.0:
                continue

            result -= pt * prod_abar

        return result


def compute_stripped_amplitude(n, omegas, tilde_zs, minus_particle=0):
    """
    Convenience function to compute A_{1,...,n} via Berends-Giele recursion.

    Parameters
    ----------
    n : int
        Number of particles.
    omegas : list of float
        Energy parameters (0-indexed).
    tilde_zs : list of float
        Holomorphic coordinates (0-indexed).
    minus_particle : int
        Index of the minus-helicity particle (default 0, i.e. particle 1).

    Returns
    -------
    float
        The stripped amplitude.
    """
    from ..kinematics.phase_space import make_tilde_lambdas

    tl_vectors = make_tilde_lambdas(omegas, tilde_zs)
    tl_dict = {i: tl_vectors[i] for i in range(n)}

    bg = BerendsGieleRecursion(tl_dict)

    # Build labels with minus particle last (for the recursion convention)
    # Use cyclicity: A_{0,1,...,n-1} = A_{1,...,n-1,0} when minus_particle=0
    all_labels = list(range(n))
    # Rotate so minus_particle is last
    idx = all_labels.index(minus_particle)
    rotated = all_labels[idx + 1:] + all_labels[:idx + 1]

    return bg.stripped_amplitude(tuple(rotated))
