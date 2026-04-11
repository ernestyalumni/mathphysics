"""
spinors.cadabra.declarations — Reusable Cadabra2 setup for spinors.

This makes Cadabra actually useful: proper index declarations,
tensor properties, and common symbolic objects.
"""

import cadabra2
from cadabra2 import Ex, AntiSymmetric, Symmetric


def declare_spinor_indices():
    """Declare all spinor and vector indices with correct properties."""

    # Undotted indices: left-handed (2,1)
    cadabra2.Indices(
        Ex(r"{\alpha, \beta, \gamma, \delta}"), Ex(r"position=fixed"))
    # Dotted indices: right-handed (1,2)
    cadabra2.Indices(Ex(r"{\dal, \dbe, \dga, \dde}"), Ex(r"position=fixed"))
    # Vector indices
    cadabra2.Indices(Ex(r"{\mu, \nu, \rho, \sigma}"), Ex(r"position=free"))

    # Epsilon tensors are antisymmetric
    AntiSymmetric(Ex(r"\epsilon_{\alpha\beta}"))
    AntiSymmetric(Ex(r"\epsilon^{\alpha\beta}"))
    AntiSymmetric(Ex(r"\epsilon_{\dal\dbe}"))
    AntiSymmetric(Ex(r"\epsilon^{\dal\dbe}"))

    # Kronecker delta is symmetric
    Symmetric(Ex(r"\delta^{\alpha}_{\beta}"))
    Symmetric(Ex(r"\delta^{\dal}_{\dbe}"))

    print("Cadabra2 spinor indices and tensor symmetries declared.")
    return True


def get_psi_L():
    """Return symbolic left-handed Weyl spinor ψ_α."""
    return Ex(r"\psi_{\alpha}")


def get_psi_R():
    """Return symbolic right-handed Weyl spinor ψ†^ȧ."""
    return Ex(r"\psidag^{\dal}")


def get_epsilon_lower():
    return Ex(r"\epsilon_{\alpha\beta}")


def get_sigma():
    """Symbolic σ^μ_{aȧ}."""
    return Ex(r"\sigma^{\mu}_{\alpha\dal}")
