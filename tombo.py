import numpy as np
from scipy.io import loadmat

import globals as g
from symmetric_5_sided_mesh import symmetric_5_sided_mesh
from nd_data import nd_data
from wing_total import wing_total
from lr_set_matrix import lr_set_matrix
from wing_m import wing_m
from lr_mass_L2GT import lr_mass_L2GT
from lrs_wing_NVs import lrs_wing_NVs
from n_vel_T_by_W import n_vel_T_by_W
from cross_matrix import cross_matrix
from assemble_matrix import assemble_matrix
from solution import solution
from s_impulse_WT import s_impulse_WT
from b_vel_B_by_T_matrix import b_vel_B_by_T_matrix
from vel_B_by_T import vel_B_by_T
from cross_vel_B_by_T import cross_vel_B_by_T
from assemble_vel_B_by_T import assemble_vel_B_by_T
from add_wake import add_wake
from force_moment import force_moment
from vel_by import vel_by
from plot_function import plot_graphs


def tombo():
    # SETUP
    # -----
    create_directories()

    xb_f, nxb_f, nb_f, xc_f, nxc_f, nc_f, l_f, c_f, h_f = \
        symmetric_5_sided_mesh(1, g.lt_f, g.lr_f, g.bang_f, g.hfactor_f, g.wfactor_f)
    xb_r, nxb_r, nb_r, xc_r, nxc_r, nc_r, l_r, c_r, h_r = \
        symmetric_5_sided_mesh(2, g.lt_r, g.lr_r, g.bang_r, g.hfactor_r, g.wfactor_r)
    
    if g.b_r - g.b_f >= 0.5 * (c_r + c_f):
        print("wing clearance checked")
    else:
        raise ValueError("rear and forward wings interfere")

    l, c, h, phiT, phiB, a, beta, delta, gMax, U, \
        xb_f, xc_f, xb_r, xc_r, b_f, b_r, e, d_, v_, rt = \
        nd_data(l_f, c_f, h_f, l_r, c_r, h_r,
                g.phiT_, g.phiB_, g.a_, g.beta_, g.delta_, g.gMax_, g.U_,
                xb_f, xc_f, xb_r, xc_r, g.b_f, g.b_r, g.f_)

    # Cutoff distance of the extension of a vortex line; velocity evaluation points
    # within this distance from the vortex line and/or its extension are set to zero
    LCUT = 0.1 * h[0]

    # Front right wing
    xc_f, xb_f, xt_f, nxt_f, xC_f, nC_f = \
        wing_total(xb_f, nxb_f, nb_f, xc_f, nxc_f, nc_f)
    # Rear right wing
    xc_r, xb_r, xt_r, nxt_r, xC_r, nC_r = \
        wing_total(xb_r, nxb_r, nb_r, xc_r, nxc_r, nc_r)

    # Wake vortex magnitude array
    GAMw_f = np.zeros((g.nwing, nxb_f))
    GAMw_r = np.zeros((g.nwing, nxb_r))
    # Total wake vortex number
    nxw_f = 0
    nxw_r = 0
    # Wake vortex location array (after convection)
    Xw_f = np.zeros((3, 4, nxb_f * g.nstep, g.nwing))
    Xw_r = np.zeros((3, 4, nxb_r * g.nstep, g.nwing))
    # Shed vortex location array 
    Xs_f = np.zeros((3, 4, nxb_f, g.nwing))
    Xs_r = np.zeros((3, 4, nxb_r, g.nwing))

    if g.nstep > 3:
        # Initialize the linear and angular impulse arrays
        limpa_f = np.zeros((3, g.nstep, g.nwing))
        limpa_r = np.zeros((3, g.nstep, g.nwing))
        aimpa_f = np.zeros((3, g.nstep, g.nwing))
        aimpa_r = np.zeros((3, g.nstep, g.nwing))
        limpw_f = np.zeros((3, g.nstep, g.nwing))
        limpw_r = np.zeros((3, g.nstep, g.nwing))
        aimpw_f = np.zeros((3, g.nstep, g.nwing))
        aimpw_r = np.zeros((3, g.nstep, g.nwing))

    # Normal velocity on the wing due to the wing motion & wake vortices
    Vnc_f = np.zeros((g.nwing, nxt_f))
    Vnc_r = np.zeros((g.nwing, nxt_r))
    Vncw_f = np.zeros((g.nwing, nxt_r))
    Vncw_r = np.zeros((g.nwing, nxt_f))

    # Sub-matrix for the non-penetration condition (self-terms)
    MVNs_f = np.zeros((nxt_f, nxt_f, g.nwing))
    MVNs_r = np.zeros((nxt_r, nxt_r, g.nwing))

    # Velocity value matrices
    VBW_f = np.zeros((3, 4, nxb_f, g.nwing))
    VBW_r = np.zeros((3, 4, nxb_r, g.nwing))
    VWT_f = np.zeros((3, 4, nxb_f * g.nstep, g.nwing))
    VWT_r = np.zeros((3, 4, nxb_r * g.nstep, g.nwing))
    VWW_f = np.zeros((3, 4, nxb_f * g.nstep, g.nwing))
    VWW_r = np.zeros((3, 4, nxb_r * g.nstep, g.nwing))

    # TODO: Document Xc_f/r
    Xc_f = np.zeros((3, 4, nxc_f, 2))
    Xc_r = np.zeros((3, 4, nxc_r, 2))
    # Space-fixed system coords of border elements on the wing
    Xb_f = np.zeros((3, 4, nxb_f, 2))
    Xb_r = np.zeros((3, 4, nxb_r, 2))
    # Global coords of total elements on the wing
    Xt_f = np.zeros((3, 4, nxt_f, 2))
    Xt_r = np.zeros((3, 4, nxt_r, 2))
    # Global coords of the collocation points on the wing
    XC_f = np.zeros((3, nxt_f, 2))
    XC_r = np.zeros((3, nxt_r, 2))
    # Unit normals at the collocation points on the wing
    NC_f = np.zeros((3, nxt_f, 2))
    NC_r = np.zeros((3, nxt_r, 2))

    # TIME MARCH
    # ----------
    MVNs_f = lr_set_matrix(xt_f, nxt_f, xC_f, nC_f, g.RCUT)
    MVNs_r = lr_set_matrix(xt_r, nxt_r, xC_r, nC_r, g.RCUT)

    for istep in range(g.nstep):

        iteration = {}
        
        if g.idebg:
            data = loadmat(f"matlab_data/data{istep + 1}.mat")

        t = istep * g.dt

        # Get wing motion parameters
        phi = np.zeros(g.twing)
        theta = np.zeros(g.twing)
        dph = np.zeros(g.twing)
        dth = np.zeros(g.twing)

        for i in range(g.twing):
            phi[i], theta[i], dph[i], dth[i] = \
                wing_m(g.mpath[i], t, rt[i], g.tau[i], e[i],
                       gMax[i], g.p[i], g.rtOff[i], phiT[i], phiB[i])

        # Get global coordinates of the points on the wing
        for i in range(g.nwing):
            # Front wing
            Xc_f[:, :, :, i], Xb_f[:, :, :, i], Xt_f[:, :, :, i], XC_f[:, :, i], NC_f[:, :, i] = \
                lr_mass_L2GT(i, beta[i], delta, phi[i], theta[i], a[i], U, t, b_f,
                             xc_f, xb_f, xt_f, xC_f, nC_f)
            # Rear wing                                                                       
            Xc_r[:, :, :, i], Xb_r[:, :, :, i], Xt_r[:, :, :, i], XC_r[:, :, i], NC_r[:, :, i] = \
                lr_mass_L2GT(i, beta[i + 2], delta, phi[i + 2], theta[i + 2], a[i + 2], U, t, b_r,
                             xc_r, xb_r, xt_r, xC_r, nC_r)

        if g.idebg:
            print(f"Xc {istep + 1}")
            print(np.allclose(Xc_f, data['Xc_f'], atol=1e-16))
            print(np.allclose(Xc_r, data['Xc_r'], atol=1e-16))

        # Find velocity of the wing
        for i in range(g.nwing):
            # Front wing
            Vnc_f[i, :] = lrs_wing_NVs(0, i, xC_f, XC_f[:, :, i], NC_f[:, :, i], t, theta[i],
                                       phi[i], dph[i], dth[i], a[i], beta[i], U, iteration)
            # Rear wing
            Vnc_r[i, :] = lrs_wing_NVs(1, i, xC_r, XC_r[:, :, i], NC_r[:, :, i], t, theta[i + 2],
                                       phi[i + 2], dph[i + 2], dth[i + 2], a[i + 2], beta[i + 2], U, iteration)

        if g.idebg:
            print(f"Vnc {istep + 1}")
            print(np.allclose(Vnc_f, data['Vnc_f'], atol=1e-16))
            print(np.allclose(Vnc_r, data['Vnc_r'], atol=1e-16))

        # Normal vel on each airfoil by front & rear, right & left wake vortices
        # For each wing, there are 4 wake vortex contributions
        for i in range(g.nwing):
            # Front wing
            Vncw_f[i, :] = n_vel_T_by_W(istep, nxt_f, XC_f[:, :, i], NC_f[:, :, i],
                                        Xw_f, GAMw_f, nxw_f, Xw_r, GAMw_r, nxw_r, g.RCUT, LCUT)
            # Rear wing  
            Vncw_r[i, :] = n_vel_T_by_W(istep, nxt_r, XC_r[:, :, i], NC_r[:, :, i],
                                        Xw_f, GAMw_f, nxw_f, Xw_r, GAMw_r, nxw_r, g.RCUT, LCUT)

        if g.idebg:
            print(f"Vncw {istep + 1}")
            print(np.allclose(Vncw_f, data['Vncw_f'], atol=1e-16))
            print(np.allclose(Vncw_r, data['Vncw_r'], atol=1e-16))

        # Calculation of the time-dependent sub-matrices MVNs_ij (i~=j)
        MVNs_12 = cross_matrix(XC_f[:, :, 0], NC_f[:, :, 0], nxt_f, Xt_f[:, :, :, 1], nxt_f, g.RCUT)
        MVNs_13 = cross_matrix(XC_f[:, :, 0], NC_f[:, :, 0], nxt_f, Xt_r[:, :, :, 0], nxt_r, g.RCUT)
        MVNs_14 = cross_matrix(XC_f[:, :, 0], NC_f[:, :, 0], nxt_f, Xt_r[:, :, :, 1], nxt_r, g.RCUT)
        MVNs_21 = cross_matrix(XC_f[:, :, 1], NC_f[:, :, 1], nxt_f, Xt_f[:, :, :, 0], nxt_f, g.RCUT)
        MVNs_23 = cross_matrix(XC_f[:, :, 1], NC_f[:, :, 1], nxt_f, Xt_r[:, :, :, 0], nxt_r, g.RCUT)
        MVNs_24 = cross_matrix(XC_f[:, :, 1], NC_f[:, :, 1], nxt_f, Xt_r[:, :, :, 1], nxt_r, g.RCUT)
        MVNs_31 = cross_matrix(XC_r[:, :, 0], NC_r[:, :, 0], nxt_r, Xt_f[:, :, :, 0], nxt_f, g.RCUT)
        MVNs_32 = cross_matrix(XC_r[:, :, 0], NC_r[:, :, 0], nxt_r, Xt_f[:, :, :, 1], nxt_f, g.RCUT)
        MVNs_34 = cross_matrix(XC_r[:, :, 0], NC_r[:, :, 0], nxt_r, Xt_r[:, :, :, 1], nxt_r, g.RCUT)
        MVNs_41 = cross_matrix(XC_r[:, :, 1], NC_r[:, :, 1], nxt_r, Xt_f[:, :, :, 0], nxt_f, g.RCUT)
        MVNs_42 = cross_matrix(XC_r[:, :, 1], NC_r[:, :, 1], nxt_r, Xt_f[:, :, :, 1], nxt_f, g.RCUT)
        MVNs_43 = cross_matrix(XC_r[:, :, 1], NC_r[:, :, 1], nxt_r, Xt_r[:, :, :, 0], nxt_r, g.RCUT)

        # Assemble the total matrix using MVNs_f[:,:,1], MVNs_r[:,:,1], MVNs_ij[:,:]
        MVN = assemble_matrix(nxt_f, nxt_r, MVNs_f, MVNs_r,
                              MVNs_12, MVNs_13, MVNs_14,
                              MVNs_21, MVNs_23, MVNs_24,
                              MVNs_31, MVNs_32, MVNs_34,
                              MVNs_41, MVNs_42, MVNs_43)

        # Solve the system of equations
        GAMA = solution(nxt_f, nxt_r, MVN, Vnc_f, Vncw_f, Vnc_r, Vncw_r)

        if g.idebg:
            print(f"MVN and GAMA {istep + 1}:")
            print(np.allclose(MVN, data['MVN'], atol=1e-16))
            print(np.allclose(GAMA, data['GAMA'], atol=1e-16))

        # Split GAMA into 4 parts
        GAM_f = np.zeros((2, nxt_f))
        GAM_r = np.zeros((2, nxt_r))

        GAM_f[0, 0:nxt_f] = GAMA[0:nxt_f]  # Front right wing
        GAM_f[1, 0:nxt_f] = GAMA[nxt_f:(2 * nxt_f)]  # Front left  wing
        GAM_r[0, 0:nxt_r] = GAMA[(2 * nxt_f):(2 * nxt_f + nxt_r)]  # Rear right wing
        GAM_r[1, 0:nxt_r] = GAMA[(2 * nxt_f + nxt_r):(2 * nxt_f + 2 * nxt_r)]  # Rear left  wing

        # Plot GAMA at the collocation points of the elements
        # using the unit normal direction: positive up and negative down
        if g.gplot:
            for i in range(g.nwing):
                # Front wing

                iteration["GAM_front"] = {
                    "i": i,
                    "t": t,
                    "GAM_f": np.copy(GAM_f[i,:]),
                    "XC_f": np.copy(XC_f[:,:,i]),
                    "NC_f": np.copy(NC_f[:,:,i])
                }
                
                # plot_GAM(0, i, t, GAM_f[i,:], XC_f[:,:,i], NC_f[:,:,i])

                # Rear wing

                iteration["GAM_rear"] = {
                    "i": i,
                    "t": t,
                    "GAM_r": np.copy(GAM_r[i,:]),
                    "XC_r": np.copy(XC_r[:,:,i]),
                    "NC_r": np.copy(NC_r[:,:,i])
                }
                
                # plot_GAM(1, i, t, GAM_r[i,:], XC_r[:,:,i], NC_r[:,:,i])

        # Plot locations, Xb & Xw, of border & wake vortices (space-fixed sys)
        if g.wplot:
            iteration["wake"] = {
                "istep": istep,
                "nxb_f": nxb_f,
                "nxw_f": nxw_f,
                "Xb_f": np.copy(Xb_f),
                "Xw_f": np.copy(Xw_f),
                "nxb_r": nxb_r,
                "nxw_r": nxw_r,
                "Xb_r": np.copy(Xb_r),
                "Xw_r": np.copy(Xw_r)
            }
            # plot_WB(istep, g.nxb_f, nxw_f, Xb_f, Xw_f, g.nxb_r, nxw_r, Xb_r, Xw_r)

        if g.nstep > 3:  # At least 4 steps needed to calculate forces and moments
            # Calculate impulses in the body-translating system
            # Include all of the bound vortices and wake vortices
            # For istep=1, there are no wake vortices
            # Front wing
            limpa, aimpa, limpw, aimpw = \
                s_impulse_WT(istep, U, t, Xt_f, Xw_f, GAM_f, GAMw_f,
                             beta[0:2], phi[0:2], theta[0:2], a[0:2])
            for j in range(3):
                for w in range(g.nwing):
                    limpa_f[j, istep, w] = limpa[j, w]
                    aimpa_f[j, istep, w] = aimpa[j, w]
                    limpw_f[j, istep, w] = limpw[j, w]
                    aimpw_f[j, istep, w] = aimpw[j, w]
            # Rear wing
            limpa, aimpa, limpw, aimpw = \
                s_impulse_WT(istep, U, t, Xt_r, Xw_r, GAM_r, GAMw_r,
                             beta[2:4], phi[2:4], theta[2:4], a[2:4])
            for j in range(3):
                for w in range(g.nwing):
                    limpa_r[j, istep, w] = limpa[j, w]
                    aimpa_r[j, istep, w] = aimpa[j, w]
                    limpw_r[j, istep, w] = limpw[j, w]
                    aimpw_r[j, istep, w] = aimpw[j, w]

        if g.idebg:
            print(f"Impulse arrays {istep + 1}:")
            print(np.allclose(data['limpa_f'], limpa_f, atol=1e-16))
            print(np.allclose(data['aimpa_f'], aimpa_f, atol=1e-16))
            print(np.allclose(data['limpw_f'], limpw_f, atol=1e-16))
            print(np.allclose(data['aimpw_f'], aimpw_f, atol=1e-16))
            print(np.allclose(data['limpa_r'], limpa_r, atol=1e-16))
            print(np.allclose(data['aimpa_r'], aimpa_r, atol=1e-16))
            print(np.allclose(data['limpw_r'], limpw_r, atol=1e-16))
            print(np.allclose(data['aimpw_r'], aimpw_r, atol=1e-16))

        # Extract GAMAb (border & shed) from GAM
        GAMAb_f = GAM_f[:, :nxb_f].copy()
        GAMAb_r = GAM_r[:, :nxb_r].copy()

        if g.idebg:
            print(f"GAMAb {istep + 1}:")
            print(np.allclose(GAMAb_f, data['GAMAb_f'], atol=1e-16))
            print(np.allclose(GAMAb_r, data['GAMAb_r'], atol=1e-16))

        # Calculate velocity of border and wake vortices to be shed or convected
        # Influence coeff for the border elem vel due to the total wing elem
        # Self-influence coeff for each wing; calculated at each time step
        cVBT_f = b_vel_B_by_T_matrix(nxb_f, nxt_f, Xb_f, Xt_f, g.RCUT)
        cVBT_r = b_vel_B_by_T_matrix(nxb_r, nxt_r, Xb_r, Xt_r, g.RCUT)

        if g.idebg:
            print(f"cVBT {istep + 1}:")
            print(np.allclose(cVBT_f, data['cVBT_f'], atol=1e-16))
            print(np.allclose(cVBT_r, data['cVBT_r'], atol=1e-16))

        # Border element veocity due to the total wing elements: self-influence
        # VBTs_m(j,n,ixb,w);  vel on wing w due to total elem on wing w
        VBTs_f = vel_B_by_T(cVBT_f, GAM_f, nxt_f)
        VBTs_r = vel_B_by_T(cVBT_r, GAM_r, nxt_r)

        if g.idebg:
            print(f"VBTs {istep + 1}:")
            print(np.allclose(VBTs_f, data['VBTs_f'], atol=1e-16))
            print(np.allclose(VBTs_r, data['VBTs_r'], atol=1e-16))

        # Border element veocity due to the total wing elements: cross-influence
        VBTs_12 = cross_vel_B_by_T(Xb_f[..., 0], nxb_f, Xt_f[..., 1], GAM_f[1, :], nxt_f, g.RCUT, LCUT)
        VBTs_13 = cross_vel_B_by_T(Xb_f[..., 0], nxb_f, Xt_r[..., 0], GAM_r[0, :], nxt_r, g.RCUT, LCUT)
        VBTs_14 = cross_vel_B_by_T(Xb_f[..., 0], nxb_f, Xt_r[..., 1], GAM_r[1, :], nxt_r, g.RCUT, LCUT)
        VBTs_21 = cross_vel_B_by_T(Xb_f[..., 1], nxb_f, Xt_f[..., 0], GAM_f[0, :], nxt_f, g.RCUT, LCUT)
        VBTs_23 = cross_vel_B_by_T(Xb_f[..., 1], nxb_f, Xt_r[..., 0], GAM_r[0, :], nxt_r, g.RCUT, LCUT)
        VBTs_24 = cross_vel_B_by_T(Xb_f[..., 1], nxb_f, Xt_r[..., 1], GAM_r[1, :], nxt_r, g.RCUT, LCUT)
        VBTs_31 = cross_vel_B_by_T(Xb_r[..., 0], nxb_r, Xt_f[..., 0], GAM_f[0, :], nxt_f, g.RCUT, LCUT)
        VBTs_32 = cross_vel_B_by_T(Xb_r[..., 0], nxb_r, Xt_f[..., 1], GAM_f[1, :], nxt_f, g.RCUT, LCUT)
        VBTs_34 = cross_vel_B_by_T(Xb_r[..., 0], nxb_r, Xt_r[..., 1], GAM_r[1, :], nxt_r, g.RCUT, LCUT)
        VBTs_41 = cross_vel_B_by_T(Xb_r[..., 1], nxb_r, Xt_f[..., 0], GAM_f[0, :], nxt_f, g.RCUT, LCUT)
        VBTs_42 = cross_vel_B_by_T(Xb_r[..., 1], nxb_r, Xt_f[..., 1], GAM_f[1, :], nxt_f, g.RCUT, LCUT)
        VBTs_43 = cross_vel_B_by_T(Xb_r[..., 1], nxb_r, Xt_r[..., 0], GAM_r[0, :], nxt_r, g.RCUT, LCUT)

        # Assemble the total border element velocity due to two wings
        VBT_f, VBT_r = assemble_vel_B_by_T(nxb_f, VBTs_f, VBTs_12, VBTs_13, VBTs_14, VBTs_21, VBTs_23, VBTs_24,
                                           nxb_r, VBTs_r, VBTs_31, VBTs_32, VBTs_34, VBTs_41, VBTs_42, VBTs_43)

        if g.idebg:
            print(f"VBT {istep + 1}:")
            print(np.allclose(VBT_f, data['VBT_f'], atol=1e-16))
            print(np.allclose(VBT_r, data['VBT_r'], atol=1e-16))

        # Velocity from wake vortices
        if istep > 0:
            # Velocity of the border elements due to wake vortices
            for i in range(g.nwing):
                VBW_f[..., i] = vel_by(istep, Xb_f[..., i], nxb_f, Xw_f, GAMw_f, nxw_f, Xw_r,
                                                    GAMw_r, nxw_r, g.RCUT, LCUT)
                VBW_r[..., i] = vel_by(istep, Xb_r[..., i], nxb_r, Xw_f, GAMw_f, nxw_f, Xw_r,
                                                    GAMw_r, nxw_r, g.RCUT, LCUT)

            # Velocity of the wake elements due to total wing vortices
            for i in range(g.nwing):
                VWT_f[..., :istep * nxb_f, i] = vel_by(istep, Xw_f[..., i], nxw_f, Xt_f, GAM_f, nxt_f,
                                                              Xt_r, GAM_r, nxt_r, g.RCUT, LCUT)
                VWT_r[..., :istep * nxb_r, i] = vel_by(istep, Xw_r[:, :, :, i], nxw_r, Xt_f, GAM_f, nxt_f,
                                                              Xt_r, GAM_r, nxt_r, g.RCUT, LCUT)

            # Velocity of the wake elements due to wake elements
            for i in range(g.nwing):
                VWW_f[..., :istep * nxb_f, i] = vel_by(istep, Xw_f[..., i], nxw_f, Xw_f, GAMw_f, nxw_f,
                                                              Xw_r, GAMw_r, nxw_r, g.RCUT, LCUT)
                VWW_r[..., :istep * nxb_r, i] = vel_by(istep, Xw_r[..., i], nxw_r, Xw_f, GAMw_f, nxw_f,
                                                              Xw_r, GAMw_r, nxw_r, g.RCUT, LCUT)

        if g.idebg:
            print(f"VBW, VWT, VWW {istep + 1}:")
            print(np.allclose(VBW_f, data['VBW_f'], atol=1e-16))
            print(np.allclose(VBW_r, data['VBW_r'], atol=1e-16))
            if istep > 0:
                print(np.allclose(VWT_f[:, :, :nxw_f, :], data['VWT_f'], atol=1e-16))
                print(np.allclose(VWT_r[:, :, :nxw_r, :], data['VWT_r'], atol=1e-16))
                print(np.allclose(VWW_f[:, :, :nxw_f, :], data['VWW_f'], atol=1e-16))
                print(np.allclose(VWW_r[:, :, :nxw_r, :], data['VWW_r'], atol=1e-16))

        # Shed border vortex elements
        Xs_f = Xb_f + g.dt * (VBT_f + VBW_f)
        Xs_r = Xb_r + g.dt * (VBT_r + VBW_r)

        # Convect wake vortices
        if istep > 0:
            Xw_f = Xw_f + g.dt * (VWT_f + VWW_f)
            Xw_r = Xw_r + g.dt * (VWT_r + VWW_r)

        # Add shed vortices to wake vortex
        if istep == 0:
            # Front wings
            GAMw_f = GAMAb_f
            nxw_f = nxb_f
            Xw_f[:, :, :nxb_f, :] = Xs_f
            # Rear wings
            GAMw_r = GAMAb_r
            nxw_r = nxb_r
            Xw_r[:, :, :nxb_r, :] = Xs_r
        else:
            GAMw_f, nxw_f, Xw_f = add_wake(istep, nxb_f, GAMAb_f, Xs_f, GAMw_f, Xw_f)
            GAMw_r, nxw_r, Xw_r = add_wake(istep, nxb_r, GAMAb_r, Xs_r, GAMw_r, Xw_r)

        if g.idebg:
            print(f"End {istep + 1}:")
            print(np.allclose(GAMw_f, data['GAMw_f'], atol=1e-16))
            print(np.allclose(nxw_f, data['nxw_f'], atol=1e-16))
            print(np.allclose(Xw_f[:, :, :nxw_f, :], data['Xw_f'], atol=1e-16))
            print(np.allclose(GAMw_r, data['GAMw_r'], atol=1e-16))
            print(np.allclose(nxt_r, data['nxt_r'], atol=1e-16))
            print(np.allclose(Xw_r[:, :, :nxw_r, :], data['Xw_r'], atol=1e-16))

        g.iterations.append(iteration)
    # END TIME MARCH

    # Calculate the force and moment on the airfoil
    if g.nstep > 3:
        force_moment(g.rho_, v_[0], d_[0], g.nstep, g.dt, U,
                     limpa_f, limpa_r, aimpa_f, aimpa_r,
                     limpw_f, limpw_r, aimpw_f, aimpw_r)

    plot_graphs()


def create_directories():
    from pathlib import Path

    base_dir = Path(g.folder)
    if not base_dir.exists():
        base_dir.mkdir()

    mesh_dir = base_dir / Path("mesh")
    if not mesh_dir.exists():
        mesh_dir.mkdir()

    debug_dir = base_dir / Path("debug")
    if not debug_dir.exists():
        debug_dir.mkdir()

    wake_dir = base_dir / Path("wake")
    if not wake_dir.exists():
        wake_dir.mkdir()

    f_and_m = base_dir / Path("f&m")
    if not f_and_m.exists():
        f_and_m.mkdir()


if __name__ == "__main__":
    tombo()
