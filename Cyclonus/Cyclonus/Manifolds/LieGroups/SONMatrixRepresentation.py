import numpy as np

class SONMatrixRepresentation:
    """
    @brief Implementation of methods common to SO(N) matrix Lie Groups using
    numpy.
    """
    def __init__(self, mat):
        """
        @brief Create a transformation from a rotation matrix (unsafe, but
        faster).
        """
        self.mat_ = mat
        """Storage for the rotation matrix."""


    def dot(self, multiplicand):
        """
        @brief Multiply another rotation or one or more vectors on the left.

        @ref https://github.com/utiasSTARS/liegroups/blob/master/liegroups/numpy/_base.py
        """
        if isinstance(multiplicand, self.__class__):
            # Compound with another rotation.
            return self.__class__(np.dot(self.mat_, multiplicand.mat_))

        else:
            # View inputs as arrays with at least 2 dims.
            multiplicand = np.atleast_2d(multiplicand)            

            # Transform 1 or more 2-vectors or fail.
            if multiplicand.shape[1] == self.N_:

                return np.squeeze(np.dot(self.mat_, multiplicand.T).T)

            else:

                raise ValueError(
                    "Vector must have shape ({},) or (N,{})".format(
                        self.N_, self.N_))


    def group_operation(self, multiplicand):
        return self.dot(multiplicand)

    def multiply(self, multiplicand):
        return self.dot(multiplicand)