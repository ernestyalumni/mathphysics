from .SONMatrixRepresentation import SONMatrixRepresentation

import numpy as np

class SO3Matrix(SONMatrixRepresentation):
    # SO(N) is the group of distance-preserving transformations of a Euclidean
    # space of dimension N, here N=N_.
    N_ = 3

    def __init__(self, mat):
        super().__init__(mat)

    @classmethod
    def rotation_about_x(cls, angle_in_radians):
        """
        @brief Form a rotation matrix given an angle in radians about the
        x-axis. Does an active rotation.
        """

        c = np.cos(angle_in_radians)
        s = np.sin(angle_in_radians)

        return cls(np.array([[1., 0., 0.],
            [0., c, -s],
            [0., s, c]]))

    @classmethod
    def inactive_rotation_about_x(cls, angle_in_radians):
        return cls.rotation_about_x(cls, -angle_in_radians)

    @classmethod
    def rotation_about_y(cls, angle_in_radians):
        """
        @brief Form a rotation matrix given an angle in radians about the
        y-axis. Does an active rotation.
        """

        c = np.cos(angle_in_radians)
        s = np.sin(angle_in_radians)

        return cls(np.array([[c, 0., s],
            [0., 1., 0.],
            [-s, 0., c]]))

    @classmethod
    def rotation_about_z(cls, angle_in_radians):
        """
        @brief Form a rotation matrix given an angle in radians about the
        z-axis. Does an active rotation.
        """

        c = np.cos(angle_in_radians)
        s = np.sin(angle_in_radians)

        return cls(np.array([[c, -s, 0.],
            [s, c, 0.],
            [0., 0., 1.]]))
