import numpy as np

class so3:
    # Basis vectors for so(3) as a 3-dim. vector space, in a 3x3 matrix
    # representation.
    # Use convention on wikipedia,
    # https://en.wikipedia.org/wiki/3D_rotation_group
    # for the hat map.

    dim_ = 3

    e_hat_1_ = np.array([[0., 0., 0.],
        [0., 0., -1.],
        [0., 1., 0.]])

    e_hat_2_ = np.array([[0., 0., 1.],
        [0., 0., 0.],
        [-1., 0., 0.]])

    e_hat_3_ = np.array([[0., -1., 0.],
        [1., 0., 0.],
        [0., 0., 0.]])


    @classmethod
    def hat_map(cls, omega):
        """
        @details

        hat_map : \\mathbb{R}^3 \to gl(3, 3) or Mat(3, 3) 3x3 matrices.

        @ref http://asrl.utias.utoronto.ca/~tdb/bib/barfoot_ser17.pdf
        @ref Timothy D. Barfoot, State Estimation For Robotics.
        What Barfoot calls "wedge" operation.

        @ref https://github.com/utiasSTARS/liegroups/blob/master/liegroups/numpy/so3.py
        """
        omega = np.atleast_2d(omega)

        if omega.shape[1] != cls.dim_:
            raise ValueError(
                "input omega must have shape ({},) or (N,{})".format(
                    cls.dim_,
                    cls.dim_))

        # omega.shape[0] yields the number of 3-arrays in the input, e.g.
        # if omega = [1, 2, 3], omega.shape[0] would be 1, and
        # if omega = [[1, 2, 3], [4, 5, 6]], there are 2 arrays of size 3
        # as inputs, and so omega.shape[0] = 2
        result = np.zeros([omega.shape[0], cls.dim_, cls.dim_])

        result[:, 0, 1] = -omega[:, 2]
        result[:, 1, 0] = omega[:, 2]
        result[:, 0, 2] = omega[:, 1]
        result[:, 2, 0] = -omega[:, 1]
        result[:, 1, 2] = -omega[:, 0]
        result[:, 2, 1] = omega[:, 0]

        return np.squeeze(result)


    @classmethod
    def map_to_matrices(cls, omega):
        """
        @brief Map vector in \\mathbb{R}^3 to 3x3 matrix representation.
        """

        omega = np.atleast_2d(omega)

        if omega.shape[1] != cls.dim_:
            raise ValueError(
                "input omega must have shape ({},) or (N,{})".format(
                    cls.dim_,
                    cls.dim_))

        return (
            cls.e_hat_1_ * omega[:, 0] + 
            cls.e_hat_2_ * omega[:, 1] + 
            cls.e_hat_3_ * omega[:, 2])

    @classmethod
    def inverse_hat_map(cls, omega_hat):
        if omega_hat.ndim < 3:
            omega_hat = np.expand_dims(omega_hat, axis=0)

        if omega_hat.shape[1:3] != (cls.dim_, cls.dim_):
            raise ValueError(
                "Omega hat must have shape ({},{}) or (N,{},{}".format(
                    cls.dim_,
                    cls.dim_,
                    cls.dim_,
                    cls.dim_))

        result = np.empty([omega_hat.shape[0], cls.dim_])
        result[:, 0] = omega_hat[:, 2, 1]
        result[:, 1] = omega_hat[:, 0, 2]
        result[:, 2] = omega_hat[:, 1, 0]

        return np.squeeze(result)
