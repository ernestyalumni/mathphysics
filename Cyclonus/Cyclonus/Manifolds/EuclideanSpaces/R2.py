from sage.manifolds.manifold import Manifold

class R2(object):
    def __init__(self):
        self.M = Manifold(2, 'R2', r'\mathbb{R}^2', start_index=1)
        self.cartesian_chart = self.M.chart('x y')
        self.U_minus = self.M.open_subset(
            'U_minus',
            # f. http://sagemanifolds.obspm.fr/examples/pdf/SM_tutorial.pdf
            # "Introducing a second chart on the manifold" the condition AND
            # written with [] instead of ()
            coord_def={self.cartesian_chart:
                (self.cartesian_chart[1] < 0, self.cartesian_chart[2] != 0)})
        self.spherical_minus_chart = self.U_minus.chart(
            r'r:(0,+oo) ph:[0, 2*pi):\phi')
        self.cartesian_chart_U_minus = self.cartesian_chart.restrict(
            self.U_minus)
        self.transition_spherical_minus_to_cartesian = \
            self.spherical_minus_chart.transition_map(
                self.cartesian_chart_U_minus,
                [
                    self.spherical_minus_chart[1] * cos(self.spherical_minus_chart[2]),
                    self.spherical_minus_chart[1] * sin(self.spherical_minus_chart[2])])
        euclidean_norm = sqrt(
            sum([self.cartesian_chart_U_minus[i[0]]**2 for i in self.M.index_generator(1)]))
        self.transition_spherical_minus_to_cartesian.set_inverse(
            euclidean_norm,
            arccot(
                self.cartesian_chart_U_minus[1],
                self.cartesian_chart_U_minus[2]))

    def equip_with_metric(self):

        self.g = self.M.riemannian_metric('g')
        for i in self.M.index_generator(1):
            self.g[i[0], i[0]] = 1