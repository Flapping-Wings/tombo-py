import numpy as np
from scipy.io import loadmat
from numba import njit

import globals as g
from VORTEXm import VORTEXm

@njit(cache=True)
def b_vel_B_by_T_matrix(nXb, nXt, Xb, Xt, RCUT):
    """
    Calculate velocity coefficients at border element nodes (no offset) due to bound vertices

    Parameters
    ----------
    nXb: int
        Number of border elements
    nXt: int
        Total number of elements on the wing
    Xb: ndarray[j, n, i, nwing]
        Border element coordinates j for node n of element i
    Xt: ndarray[j, n, i, nwing]
        Coordinates of the nodes for total elements on the wing.

    Returns
    -------
    VBT: ndarray[j (1 - 3 (x, y, z)), n (4 nodes for each shed element), ibelm, itelm]
        TODO
    """

    cVBT = np.zeros((3, 4, nXb, nXt, g.nwing))
    r = (1, 4, nXb)

    for w in range(g.nwing):
        for i in range(nXt):
            U = np.zeros(r)
            V = np.zeros(r)
            W = np.zeros(r)

            dU, dV, dW = VORTEXm(Xb[0, :, :, w], Xb[1, :, :, w], Xb[2, :, :, w],
                                 Xt[0, 0, i, w], Xt[1, 0, i, w], Xt[2, 0, i, w],
                                 Xt[0, 1, i, w], Xt[1, 1, i, w], Xt[2, 1, i, w],
                                 1.0, RCUT)
            U += dU
            V += dV
            W += dW

            dU, dV, dW = VORTEXm(Xb[0, :, :, w], Xb[1, :, :, w], Xb[2, :, :, w],
                                 Xt[0, 1, i, w], Xt[1, 1, i, w], Xt[2, 1, i, w],
                                 Xt[0, 2, i, w], Xt[1, 2, i, w], Xt[2, 2, i, w],
                                 1.0, RCUT)
            U += dU
            V += dV
            W += dW

            dU, dV, dW = VORTEXm(Xb[0, :, :, w], Xb[1, :, :, w], Xb[2, :, :, w],
                                 Xt[0, 2, i, w], Xt[1, 2, i, w], Xt[2, 2, i, w],
                                 Xt[0, 3, i, w], Xt[1, 3, i, w], Xt[2, 3, i, w],
                                 1.0, RCUT)
            U += dU
            V += dV
            W += dW

            dU, dV, dW = VORTEXm(Xb[0, :, :, w], Xb[1, :, :, w], Xb[2, :, :, w],
                                 Xt[0, 3, i, w], Xt[1, 3, i, w], Xt[2, 3, i, w],
                                 Xt[0, 0, i, w], Xt[1, 0, i, w], Xt[2, 0, i, w],
                                 1.0, RCUT)
            U += dU
            V += dV
            W += dW

            cVBT[0, :, :, i, w] = U
            cVBT[1, :, :, i, w] = V
            cVBT[2, :, :, i, w] = W

    return cVBT
