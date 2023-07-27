import numpy as np
import globals as g

def nd_data(
        l_f, 
        c_f,
        h_f,
        l_r,
        c_r,
        h_r,
        phiT_,
        phiB_,
        a_,
        beta_,
        delta_,
        gMax_,
        U_,
        Xb_f,
        Xc_f,
        Xb_r,
        Xc_r,
        b_f,
        b_r
):
    """
    Nondimensionalize the input data
    
    Parameters
    ----------
    l_f: float
        Wing span (front)
    c_f: float
        Wing chord length (front)
    h_f: float
        Border element height (front)
    l_r: float
        Wing span (rear)
    c_r: float
        Wing chord length (rear)
    h_r: float
        Border element height (rear)
    phiT_: ndarray
        Top stroke angle in degrees
    phiB_: ndarray
        Bottom stroke angle in degrees
    a_: ndarray
        Rotation axis offset in cm
    beta_: ndarray
        Stroke plane angle in degrees wrt the body axis 
    delta_: float
        Body angle in degrees
    gMax_: ndarray
        Max rotation amplitude in degrees; actual rotation is 2*gMax
    U_: ndarray
        Ambient velocity in (x, y, z); can be interpreted as the flight
        velocity when the wind is calm
    Xb: ndarray[j, n, i]
        Border element coordinates (front)
    Xc: ndarray[j, n, i]
        Center element coordinates (front)
    Xb_r: ndarray[j, n, i]
        Border element coordinates (rear)
    Xc_r: ndarray[j, n, i]
        Center element coordinates (rear)
    b_f: float 
        Front wing location in cm
    b_r: float
        Rear wing location in cm

    Returns
    -------
    l: ndarray
        Wing spans (nondimensional)
    c: ndarray
        Wing chord lengths (nondimensional)
    h: ndarray
        Border element height (nondimensional)
    phiT: ndarray
        Top stroke angle in radians
    phiB: ndarray
        Bottom stroke angle in radians
    a: ndarray
        Rotation axis offset (nondimensional)
    beta: ndarray
        Stroke plane angle in radians wrt the body axis 
    delta: float
        Body angle in radians
    gMax: ndarray
        Max rotation amplitude in radians; actual rotation is 2*gMax
    U: ndarray
        Ambient velocity in (x, y, z); can be interpreted as the flight
        velocity when the wind is calm (nondimensional)
    Xb: ndarray[j, n, i]
        Border element coordinates (front) (nondimensional)
    Xc: ndarray[j, n, i]
        Center element coordinates (front) (nondimensional)
    Xb_r: ndarray[j, n, i]
        Border element coordinates (rear) (nondimensional)
    Xc_r: ndarray[j, n, i]
        Center element coordinates (rear) (nondimensional)
    b_f: float 
        Front wing location (nondimensional)
    b_r: float
        Rear wing location (nondimensional)
    e: ndarray
        Difference between top and bottom stroke angles
    d_: ndarray
        Total stroke length
    v_: ndarray
        Stroke velocity
    """
    # Wing span, chord, and border height
    l_ = np.array([l_f, l_f, l_r, l_r])
    c_ = np.array([c_f, c_f, c_r, c_r])
    h_ = np.array([h_f, h_f, h_r, h_r])

    # Period
    T_ = 1.0 / g.f_
    g.rt = T_[0] / T_

    # Convert angles to radians
    fac = np.pi / 180
    phiT = fac * phiT_
    phiB = fac * phiB_
    beta = fac * beta_
    delta = fac * delta_
    gMax = fac * gMax_

    # Get non-dimensional quantities based on given
    # flight data of actual insect
    dT_ = l_ * phiT
    dB_ = l_ * -phiB
    d_ = dT_ + dB_
    e_ = dT_ - dB_

    a = a_ / d_[0]
    c = c_ / d_[0]
    l = l_ / d_[0]
    h = h_ / d_[0]

    b_f = b_f / d_[0]
    b_r = b_r / d_[0]
    Xb_f = Xb_f / d_[0]
    Xb_r = Xb_r / d_[0]
    Xc_f = Xc_f / d_[0]
    Xc_r = Xc_r / d_[0]

    e = e_ / d_

    # Reference time - use the time for the right front wing
    g.t_ = T_[0] / 2.0

    # Reference velocity - use the right front wing flapping velocity
    v_ = d_ / g.t_

    # Ambient velocity (nondimensional)
    U = U_ / v_[0]

    return l, c, h, phiT, phiB, a, beta, delta, gMax, U, Xb_f, Xc_f, Xb_r, Xc_r, \
        b_f, b_r, e, d_, v_
