from typing import Tuple
import numpy as np

class LorentzMetric:
    """Srednicki metric: diag(-1, +1, +1, +1) — "mostly plus"."""
    def __init__(self):
        self.g = np.diag([-1.0, 1.0, 1.0, 1.0])
        self.g_up = self.g.copy()  # same for this metric

    def __getitem__(self, idx: Tuple[int, int]) -> float:
        return self.g[idx]