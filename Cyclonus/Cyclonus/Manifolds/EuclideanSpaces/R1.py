from sage.manifolds.manifold import Manifold

class R1(object):
    def __init__(self):
        self.M = Manifold(1, 'R1', r'\mathbb{R}', start_index=1)
        self.cartesian_chart = self.M.chart('x')