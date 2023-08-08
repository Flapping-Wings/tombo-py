import numpy as np
import globals as g

def vel_B_by_T(cVBT, GAM, nXt):
    s = np.shape(cVBT)
    r = (s[0], s[1], s[2], g.nwing)
    VBT = np.zeros(r)
    dVBT = np.zeros(r)

    for w in range(g.nwing):
        for itelem in range(nXt):
            dVBT[..., w] = cVBT[..., itelem, w] * GAM[w, itelem]
            VBT[..., w] = VBT[..., w] + dVBT[..., w]

    return VBT
