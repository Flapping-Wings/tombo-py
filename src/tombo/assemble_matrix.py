import numpy as np

def assemble_matrix(MVNs_f, MVNs_r,
                    MVNs_12, MVNs_13, MVNs_14,
                    MVNs_21, MVNs_23, MVNs_24,
                    MVNs_31, MVNs_32, MVNs_34,
                    MVNs_41, MVNs_42, MVNs_43
):
    """
    Assemble 16 sub-matrices.

    A sub-matrix `MVNs_ij` corresponds to 
    `target_wing = i` and `source_wing = j`.

    The sub-matrices along the diagonal (coordinates [i, i])
    are constructed as below:
        `MVNs_11: ndarray[nxt_f, nxt_f] = MVNs_f[nxt_f, nxt_f, 0]`
        `MVNs_22: ndarray[nxt_f, nxt_f] = MVNs_f[nxt_f, nxt_f, 1]`
        `MVNs_33: ndarray[nxt_r, nxt_r] = MVNs_r[nxt_r, nxt_r, 0]`
        `MVNs_44: ndarray[nxt_r, nxt_r] = MVNs_r[nxt_r, nxt_r, 1]`
    
    Parameters
    ----------
    `nxt_f, nxt_r`: ints
        Size of the sub-matrices for front (f) and rear (r) wings
    `MVNs_f`: ndarray[nxt_f, nxt_f, 2]
        Self-influence sub-matrices for front wings
    `MVNs_r`: ndarray[nxt_r, nxt_r, 2]
        Self-influence sub-matrices for rear wings
    TODO
    `MVNs_12`: ndarray[nxt_f, nxt_f]
    `MVNs_13`: ndarray[nxt_f, nxt_r]
    `MVNs_14`: ndarray[nxt_f, nxt_r]
    `MVNs_21`: ndarray[nxt_f, nxt_f]
    `MVNs_23`: ndarray[nxt_f, nxt_r]
    `MVNs_24`: ndarray[nxt_f, nxt_r]
    `MVNs_31`: ndarray[nxt_r, nxt_f]
    `MVNs_32`: ndarray[nxt_r, nxt_f]
    `MVNs_34`: ndarray[nxt_r, nxt_r]
    `MVNs_41`: ndarray[nxt_r, nxt_f]
    `MVNs_42`: ndarray[nxt_r, nxt_f]
    `MVNs_43`: ndarray[nxt_r, nxt_r]
    
    Returns
    -------
    `MVN`: ndarray[2 * (nxt_f + nxt_r), 2 * (nxt_f + nxt_r)]
        Assembled matrix
    """
    MVN = np.block([
        [MVNs_f[..., 0], MVNs_12,        MVNs_13,        MVNs_14       ],
        [MVNs_21,        MVNs_f[..., 1], MVNs_23,        MVNs_24       ],
        [MVNs_31,        MVNs_32,        MVNs_r[..., 0], MVNs_34       ],
        [MVNs_41,        MVNs_42,        MVNs_43,        MVNs_r[..., 1]]
    ])

    return MVN
