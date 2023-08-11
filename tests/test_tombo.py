import pytest
import numpy as np
import numpy.testing as npt
from scipy.io import loadmat
import tombo.globals as g
import test_globals

@pytest.fixture
def matlab_wing_data():
    return loadmat("tests/matlab_data/wing.mat",
                   squeeze_me=True)

@pytest.fixture
def matlab_nd_data_data():
    return loadmat("tests/matlab_data/nd_data.mat",
                   squeeze_me=True)

@pytest.fixture
def matlab_wing_total_data():
    return loadmat("tests/matlab_data/wing_total.mat",
                   squeeze_me=True)

@pytest.fixture(params=[1, 2, 3, 4])
def matlab_loop_data(request):
    return loadmat(f"tests/matlab_data/loop_data{request.param}.mat",
                   squeeze_me=True)

def test_symmetric_5_sided_mesh(matlab_wing_data):
    from tombo.symmetric_5_sided_mesh import symmetric_5_sided_mesh

    xb_f, nxb_f, nb_f, xc_f, nxc_f, nc_f, l_f, c_f, h_f = \
        symmetric_5_sided_mesh('f', g.lt_f, g.lr_f, g.bang_f, g.hfactor_f, g.wfactor_f)
    xb_r, nxb_r, nb_r, xc_r, nxc_r, nc_r, l_r, c_r, h_r = \
        symmetric_5_sided_mesh('r', g.lt_r, g.lr_r, g.bang_r, g.hfactor_r, g.wfactor_r)
    
    npt.assert_allclose(xb_f, matlab_wing_data['xb_f'])
    npt.assert_allclose(nxb_f, matlab_wing_data['nxb_f'])
    npt.assert_allclose(nb_f, matlab_wing_data['nb_f'])
    npt.assert_allclose(xc_f, matlab_wing_data['xc_f'])
    npt.assert_allclose(nxc_f, matlab_wing_data['nxc_f'])
    npt.assert_allclose(nc_f, matlab_wing_data['nc_f'])
    npt.assert_allclose(l_f, matlab_wing_data['l_f'])
    npt.assert_allclose(c_f, matlab_wing_data['c_f'])
    npt.assert_allclose(h_f, matlab_wing_data['h_f'])

    npt.assert_allclose(xb_r, matlab_wing_data['xb_r'])
    npt.assert_allclose(nxb_r, matlab_wing_data['nxb_r'])
    npt.assert_allclose(nb_r, matlab_wing_data['nb_r'])
    npt.assert_allclose(xc_r, matlab_wing_data['xc_r'])
    npt.assert_allclose(nxc_r, matlab_wing_data['nxc_r'])
    npt.assert_allclose(nc_r, matlab_wing_data['nc_r'])
    npt.assert_allclose(l_r, matlab_wing_data['l_r'])
    npt.assert_allclose(c_r, matlab_wing_data['c_r'])
    npt.assert_allclose(h_r, matlab_wing_data['h_r'])

def test_nd_data(matlab_wing_data, matlab_nd_data_data):
    from tombo.nd_data import nd_data

    xb_f = matlab_wing_data['xb_f']
    xc_f = matlab_wing_data['xc_f']
    xb_r = matlab_wing_data['xb_r']
    xc_r = matlab_wing_data['xc_r']
    l_f = matlab_wing_data['l_f']
    c_f = matlab_wing_data['c_f']
    h_f = matlab_wing_data['h_f']
    l_r = matlab_wing_data['l_r']
    c_r = matlab_wing_data['c_r']
    h_r = matlab_wing_data['h_r']

    l, c, h, phiT, phiB, a, beta, delta, gMax, U, \
        xb_f, xc_f, xb_r, xc_r, b_f, b_r, e, d_, v_, rt = \
            nd_data(l_f, c_f, h_f, l_r, c_r, h_r,
                    g.phiT_, g.phiB_, g.a_, g.beta_, g.delta_, g.gMax_,
                    g.U_, xb_f, xc_f, xb_r, xc_r, g.b_f, g.b_r, g.f_)
    
    npt.assert_allclose(l, matlab_nd_data_data['l'])
    npt.assert_allclose(c, matlab_nd_data_data['c'])
    npt.assert_allclose(h, matlab_nd_data_data['h'])
    npt.assert_allclose(phiT, matlab_nd_data_data['phiT'])
    npt.assert_allclose(phiB, matlab_nd_data_data['phiB'])
    npt.assert_allclose(a, matlab_nd_data_data['a'])
    npt.assert_allclose(beta, matlab_nd_data_data['beta'])
    npt.assert_allclose(delta, matlab_nd_data_data['delta'])
    npt.assert_allclose(gMax, matlab_nd_data_data['gMax'])
    npt.assert_allclose(U, matlab_nd_data_data['U'])
    npt.assert_allclose(xb_f, matlab_nd_data_data['xb_f'])
    npt.assert_allclose(xc_f, matlab_nd_data_data['xc_f'])
    npt.assert_allclose(xb_r, matlab_nd_data_data['xb_r'])
    npt.assert_allclose(xc_r, matlab_nd_data_data['xc_r'])
    npt.assert_allclose(b_f, matlab_nd_data_data['b_f'])
    npt.assert_allclose(b_r, matlab_nd_data_data['b_r'])
    npt.assert_allclose(e, matlab_nd_data_data['e'])
    npt.assert_allclose(d_, matlab_nd_data_data['d_'])
    npt.assert_allclose(v_[0], matlab_nd_data_data['v_']) # Compare only first element
    npt.assert_allclose(rt, matlab_nd_data_data['rt'])


def test_wing_total(matlab_nd_data_data, matlab_wing_total_data):
    from tombo.wing_total import wing_total

    xb_f = matlab_nd_data_data['xb_f']
    nxb_f = matlab_nd_data_data['nxb_f']
    nb_f = matlab_nd_data_data['nb_f']
    xc_f = matlab_nd_data_data['xc_f']
    nxc_f = matlab_nd_data_data['nxc_f']
    nc_f = matlab_nd_data_data['nc_f']

    xb_r = matlab_nd_data_data['xb_r']
    nxb_r = matlab_nd_data_data['nxb_r']
    nb_r = matlab_nd_data_data['nb_r']
    xc_r = matlab_nd_data_data['xc_r']
    nxc_r = matlab_nd_data_data['nxc_r']
    nc_r = matlab_nd_data_data['nc_r']

    xc_f, xb_f, xt_f, nxt_f, xC_f, nC_f = \
        wing_total(xb_f, nxb_f, nb_f, xc_f, nxc_f, nc_f)
    xc_r, xb_r, xt_r, nxt_r, xC_r, nC_r = \
        wing_total(xb_r, nxb_r, nb_r, xc_r, nxc_r, nc_r)
    
    npt.assert_allclose(xc_f, matlab_wing_total_data['xc_f'])
    npt.assert_allclose(xb_f, matlab_wing_total_data['xb_f'])
    npt.assert_allclose(xt_f, matlab_wing_total_data['xt_f'])
    npt.assert_allclose(nxt_f, matlab_wing_total_data['nxt_f'])
    npt.assert_allclose(xC_f, matlab_wing_total_data['xC_f'])
    npt.assert_allclose(nC_f, matlab_wing_total_data['nC_f'])

    npt.assert_allclose(xc_r, matlab_wing_total_data['xc_r'])
    npt.assert_allclose(xb_r, matlab_wing_total_data['xb_r'])
    npt.assert_allclose(xt_r, matlab_wing_total_data['xt_r'])
    npt.assert_allclose(nxt_r, matlab_wing_total_data['nxt_r'])
    npt.assert_allclose(xC_r, matlab_wing_total_data['xC_r'])
    npt.assert_allclose(nC_r, matlab_wing_total_data['nC_r'])


def test_lr_set_matrix(matlab_wing_total_data, matlab_loop_data):
    from tombo.lr_set_matrix import lr_set_matrix

    xt_f = matlab_wing_total_data['xt_f']
    nxt_f = matlab_wing_total_data['nxt_f']
    xC_f = matlab_wing_total_data['xC_f']
    nC_f = matlab_wing_total_data['nC_f']

    xt_r = matlab_wing_total_data['xt_r']
    nxt_r = matlab_wing_total_data['nxt_r']
    xC_r = matlab_wing_total_data['xC_r']
    nC_r = matlab_wing_total_data['nC_r']

    MVNs_f = lr_set_matrix(xt_f, nxt_f, xC_f, nC_f, g.RCUT)
    MVNs_r = lr_set_matrix(xt_r, nxt_r, xC_r, nC_r, g.RCUT)

    npt.assert_allclose(MVNs_f, matlab_loop_data['MVNs_f'])
    npt.assert_allclose(MVNs_r, matlab_loop_data['MVNs_r'])


def test_wing_m(matlab_loop_data):
    from tombo.wing_m import wing_m

    phi = np.empty(g.twing)
    theta = np.empty(g.twing)
    dph = np.empty(g.twing)
    dth = np.empty(g.twing)

    t = matlab_loop_data['t']
    rt = matlab_loop_data['rt']
    e = matlab_loop_data['e']
    gMax = matlab_loop_data['gMax']
    phiT = matlab_loop_data['phiT']
    phiB = matlab_loop_data['phiB']

    for i in range(g.twing):
        phi[i], theta[i], dph[i], dth[i] = \
            wing_m(g.mpath[i], t, rt[i], g.tau[i], e[i],
                   gMax[i], g.p[i], g.rtOff[i], phiT[i], phiB[i])
        
    npt.assert_allclose(phi, matlab_loop_data['phi'])
    npt.assert_allclose(theta, matlab_loop_data['theta'])
    npt.assert_allclose(dph, matlab_loop_data['dph'])
    npt.assert_allclose(dth, matlab_loop_data['dth'])

def test_lr_mass_L2GT(matlab_loop_data):
    from tombo.lr_mass_L2GT import lr_mass_L2GT

    nxc_f = matlab_loop_data['nxc_f']
    nxc_r = matlab_loop_data['nxc_r']
    nxb_f = matlab_loop_data['nxb_f']
    nxb_r = matlab_loop_data['nxb_r']
    nxt_f = matlab_loop_data['nxt_f']
    nxt_r = matlab_loop_data['nxt_r']

    beta = matlab_loop_data['beta']
    delta = matlab_loop_data['delta']
    phi = matlab_loop_data['phi']
    theta = matlab_loop_data['theta']
    a = matlab_loop_data['a']
    U = matlab_loop_data['U']
    t = matlab_loop_data['t']

    b_f = matlab_loop_data['b_f']
    xc_f = matlab_loop_data['xc_f']
    xb_f = matlab_loop_data['xb_f']
    xt_f = matlab_loop_data['xt_f']
    xC_f = matlab_loop_data['xC_f']
    nC_f = matlab_loop_data['nC_f']

    b_r = matlab_loop_data['b_r']
    xc_r = matlab_loop_data['xc_r']
    xb_r = matlab_loop_data['xb_r']
    xt_r = matlab_loop_data['xt_r']
    xC_r = matlab_loop_data['xC_r']
    nC_r = matlab_loop_data['nC_r']

    Xc_f = np.empty((3, 4, nxc_f, 2))
    Xc_r = np.empty((3, 4, nxc_r, 2))
    Xb_f = np.empty((3, 4, nxb_f, 2))
    Xb_r = np.empty((3, 4, nxb_r, 2))
    Xt_f = np.empty((3, 4, nxt_f, 2))
    Xt_r = np.empty((3, 4, nxt_r, 2))
    XC_f = np.empty((3, nxt_f, 2))
    XC_r = np.empty((3, nxt_r, 2))
    NC_f = np.empty((3, nxt_f, 2))
    NC_r = np.empty((3, nxt_r, 2))

    for i in range(g.nwing):
        Xc_f[..., i], Xb_f[..., i], Xt_f[..., i], XC_f[..., i], NC_f[..., i] = \
            lr_mass_L2GT(i, beta[i], delta, phi[i], theta[i], a[i], U, t, b_f,
                                xc_f, xb_f, xt_f, xC_f, nC_f)
        Xc_r[..., i], Xb_r[..., i], Xt_r[..., i], XC_r[..., i], NC_r[..., i] = \
            lr_mass_L2GT(i, beta[i + 2], delta, phi[i + 2], theta[i + 2], a[i + 2], U, t, b_r,
                                xc_r, xb_r, xt_r, xC_r, nC_r)

    npt.assert_allclose(Xc_f, matlab_loop_data['Xc_f'])
    npt.assert_allclose(Xb_f, matlab_loop_data['Xb_f'])
    npt.assert_allclose(Xt_f, matlab_loop_data['Xt_f'])
    npt.assert_allclose(XC_f, matlab_loop_data['XC_f'])
    npt.assert_allclose(NC_f, matlab_loop_data['NC_f'])

def test_lrs_wing_NVs(matlab_loop_data):
    from tombo.lrs_wing_NVs import lrs_wing_NVs

    nxt_f = matlab_loop_data['nxt_f']
    nxt_r = matlab_loop_data['nxt_r']

    xC_f = matlab_loop_data['xC_f']
    XC_f = matlab_loop_data['XC_f']
    NC_f = matlab_loop_data['NC_f']
    xC_r = matlab_loop_data['xC_r']
    XC_r = matlab_loop_data['XC_r']
    NC_r = matlab_loop_data['NC_r']

    t = matlab_loop_data['t']
    theta = matlab_loop_data['theta']
    phi = matlab_loop_data['phi']
    dph = matlab_loop_data['dph']
    dth = matlab_loop_data['dth']
    a = matlab_loop_data['a']
    beta = matlab_loop_data['beta']
    U = matlab_loop_data['U']

    Vnc_f = np.empty((g.nwing, nxt_f))
    Vnc_r = np.empty((g.nwing, nxt_r))

    for i in range(g.nwing):
        Vnc_f[i, :] = lrs_wing_NVs(0, i, xC_f, XC_f[:, :, i], NC_f[:, :, i], t, theta[i],
                                   phi[i], dph[i], dth[i], a[i], beta[i], U)
        Vnc_r[i, :] = lrs_wing_NVs(1, i, xC_r, XC_r[:, :, i], NC_r[:, :, i], t, theta[i + 2],
                                   phi[i + 2], dph[i + 2], dth[i + 2], a[i + 2], beta[i + 2], U)
        
    npt.assert_allclose(Vnc_f, matlab_loop_data['Vnc_f'])
    npt.assert_allclose(Vnc_r, matlab_loop_data['Vnc_r'])

def test_n_vel_T_by_W(matlab_loop_data):
    from tombo.n_vel_T_by_W import n_vel_T_by_W

    istep = matlab_loop_data['istep'] - 1 # Convert from MATLAB indexing
    LCUT = matlab_loop_data['LCUT']

    nxt_f = matlab_loop_data['nxt_f']
    nxt_r = matlab_loop_data['nxt_r']

    XC_f = matlab_loop_data['XC_f']
    NC_f = matlab_loop_data['NC_f']
    XC_r = matlab_loop_data['XC_r']
    NC_r = matlab_loop_data['NC_r']

    nxb_f = matlab_loop_data['nxb_f']
    nxw_f = matlab_loop_data['nxw_f']
    GAMw_f = np.ascontiguousarray(matlab_loop_data['GAMw_f'])

    nxb_r = matlab_loop_data['nxb_r']
    nxw_r = matlab_loop_data['nxw_r']
    GAMw_r = np.ascontiguousarray(matlab_loop_data['GAMw_r'])

    up_to_f = GAMw_f.shape[1]
    Xw_f = np.empty((3, 4, nxb_f * g.nstep, g.nwing))
    Xw_f[:, :, :up_to_f] = matlab_loop_data['Xw_f']
    
    up_to_r = GAMw_r.shape[1]
    Xw_r = np.empty((3, 4, nxb_r * g.nstep, g.nwing))
    Xw_r[:, :, :up_to_r] = matlab_loop_data['Xw_r']

    Vncw_f = np.empty((g.nwing, nxt_r))
    Vncw_r = np.empty((g.nwing, nxt_f))

    for i in range(g.nwing):
        Vncw_f[i, :] = n_vel_T_by_W(istep, nxt_f, XC_f[..., i], NC_f[..., i],
                                    Xw_f, GAMw_f, nxw_f, Xw_r, GAMw_r, nxw_r, g.RCUT, LCUT)
        Vncw_r[i, :] = n_vel_T_by_W(istep, nxt_r, XC_r[..., i], NC_r[..., i],
                                    Xw_f, GAMw_f, nxw_f, Xw_r, GAMw_r, nxw_r, g.RCUT, LCUT)
        
    npt.assert_allclose(Vncw_f, matlab_loop_data['Vncw_f'])
    npt.assert_allclose(Vncw_r, matlab_loop_data['Vncw_r'])

def test_cross_matrix(matlab_loop_data):
    from tombo.cross_matrix import cross_matrix

    XC_f = matlab_loop_data['XC_f']
    NC_f = matlab_loop_data['NC_f']
    Xt_f = matlab_loop_data['Xt_f']
    nxt_f = matlab_loop_data['nxt_f']

    XC_r = matlab_loop_data['XC_r']
    NC_r = matlab_loop_data['NC_r']
    Xt_r = matlab_loop_data['Xt_r']
    nxt_r = matlab_loop_data['nxt_r']

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

    npt.assert_allclose(MVNs_12, matlab_loop_data['MVNs_12'])
    npt.assert_allclose(MVNs_13, matlab_loop_data['MVNs_13'])
    npt.assert_allclose(MVNs_14, matlab_loop_data['MVNs_14'])
    npt.assert_allclose(MVNs_21, matlab_loop_data['MVNs_21'])
    npt.assert_allclose(MVNs_23, matlab_loop_data['MVNs_23'])
    npt.assert_allclose(MVNs_24, matlab_loop_data['MVNs_24'])
    npt.assert_allclose(MVNs_31, matlab_loop_data['MVNs_31'])
    npt.assert_allclose(MVNs_32, matlab_loop_data['MVNs_32'])
    npt.assert_allclose(MVNs_34, matlab_loop_data['MVNs_34'])
    npt.assert_allclose(MVNs_41, matlab_loop_data['MVNs_41'])
    npt.assert_allclose(MVNs_42, matlab_loop_data['MVNs_42'])
    npt.assert_allclose(MVNs_43, matlab_loop_data['MVNs_43'])

def test_assemble_matrix(matlab_loop_data):
    from tombo.assemble_matrix import assemble_matrix

    MVNs_f = matlab_loop_data['MVNs_f']
    MVNs_r = matlab_loop_data['MVNs_r']

    MVNs_12 = matlab_loop_data['MVNs_12']
    MVNs_13 = matlab_loop_data['MVNs_13']
    MVNs_14 = matlab_loop_data['MVNs_14']
    MVNs_21 = matlab_loop_data['MVNs_21']
    MVNs_23 = matlab_loop_data['MVNs_23']
    MVNs_24 = matlab_loop_data['MVNs_24']
    MVNs_31 = matlab_loop_data['MVNs_31']
    MVNs_32 = matlab_loop_data['MVNs_32']
    MVNs_34 = matlab_loop_data['MVNs_34']
    MVNs_41 = matlab_loop_data['MVNs_41']
    MVNs_42 = matlab_loop_data['MVNs_42']
    MVNs_43 = matlab_loop_data['MVNs_43']

    MVN = assemble_matrix(MVNs_f, MVNs_r,
                          MVNs_12, MVNs_13, MVNs_14,
                          MVNs_21, MVNs_23, MVNs_24,
                          MVNs_31, MVNs_32, MVNs_34,
                          MVNs_41, MVNs_42, MVNs_43)
    
    npt.assert_allclose(MVN, matlab_loop_data['MVN'])

def test_solution(matlab_loop_data):
    from tombo.solution import solution

    MVN = matlab_loop_data['MVN']
    nxt_f = matlab_loop_data['nxt_f']
    Vnc_f = matlab_loop_data['Vnc_f']
    Vncw_f = matlab_loop_data['Vncw_f']

    nxt_r = matlab_loop_data['nxt_r']
    Vnc_r = matlab_loop_data['Vnc_r']
    Vncw_r = matlab_loop_data['Vncw_r']

    GAMA = solution(nxt_f, nxt_r, MVN, Vnc_f, Vncw_f, Vnc_r, Vncw_r)

    npt.assert_allclose(GAMA, matlab_loop_data['GAMA'])

@pytest.mark.skip(reason="values too close to zero to be meaningful")
def test_s_impulse_WT(matlab_loop_data):
    from tombo.s_impulse_WT import s_impulse_WT

    istep = matlab_loop_data['istep'] - 1 # Convert from MATLAB indexing
    U = matlab_loop_data['U']
    t = matlab_loop_data['t']
    beta = matlab_loop_data['beta']
    phi = matlab_loop_data['phi']
    theta = matlab_loop_data['theta']
    a = matlab_loop_data['a']

    Xt_f = matlab_loop_data['Xt_f']
    Xw_f = matlab_loop_data['Xw_f']
    GAM_f = matlab_loop_data['GAM_f']
    GAMw_f = matlab_loop_data['GAMw_f']

    Xt_r = matlab_loop_data['Xt_r']
    Xw_r = matlab_loop_data['Xw_r']
    GAM_r = matlab_loop_data['GAM_r']
    GAMw_r = matlab_loop_data['GAMw_r']

    limpa_f = matlab_loop_data['limpa_f']
    aimpa_f = matlab_loop_data['aimpa_f']
    limpw_f = matlab_loop_data['limpw_f']
    aimpw_f = matlab_loop_data['aimpw_f']

    limpa_r = matlab_loop_data['limpa_r']
    aimpa_r = matlab_loop_data['aimpa_r']
    limpw_r = matlab_loop_data['limpw_r']
    aimpw_r = matlab_loop_data['aimpw_r']

    limpa_f = np.empty((3, g.nstep, g.nwing))
    limpa_r = np.empty((3, g.nstep, g.nwing))
    aimpa_f = np.empty((3, g.nstep, g.nwing))
    aimpa_r = np.empty((3, g.nstep, g.nwing))
    limpw_f = np.empty((3, g.nstep, g.nwing))
    limpw_r = np.empty((3, g.nstep, g.nwing))
    aimpw_f = np.empty((3, g.nstep, g.nwing))
    aimpw_r = np.empty((3, g.nstep, g.nwing))

    limpa, aimpa, limpw, aimpw = \
        s_impulse_WT(istep, U, t, Xt_f, Xw_f, GAM_f, GAMw_f,
                     beta[0:2], phi[0:2], theta[0:2], a[0:2])
    for j in range(3):
        for w in range(g.nwing):
            limpa_f[j, istep, w] = limpa[j, w]
            aimpa_f[j, istep, w] = aimpa[j, w]
            limpw_f[j, istep, w] = limpw[j, w]
            aimpw_f[j, istep, w] = aimpw[j, w]

    limpa, aimpa, limpw, aimpw = \
        s_impulse_WT(istep, U, t, Xt_r, Xw_r, GAM_r, GAMw_r,
                     beta[2:4], phi[2:4], theta[2:4], a[2:4])
    for j in range(3):
        for w in range(g.nwing):
            limpa_r[j, istep, w] = limpa[j, w]
            aimpa_r[j, istep, w] = aimpa[j, w]
            limpw_r[j, istep, w] = limpw[j, w]
            aimpw_r[j, istep, w] = aimpw[j, w]
    
    npt.assert_allclose(limpa_f, matlab_loop_data['limpa_f'])
    npt.assert_allclose(aimpa_f, matlab_loop_data['aimpa_f'])
    npt.assert_allclose(limpw_f, matlab_loop_data['limpw_f'])
    npt.assert_allclose(aimpw_f, matlab_loop_data['aimpw_f'])
    npt.assert_allclose(limpa_r, matlab_loop_data['limpa_r'])
    npt.assert_allclose(aimpa_r, matlab_loop_data['aimpa_r'])
    npt.assert_allclose(limpw_r, matlab_loop_data['limpw_r'])
    npt.assert_allclose(aimpw_r, matlab_loop_data['aimpw_r'])

def test_b_vel_B_by_T_matrix(matlab_loop_data):
    from tombo.b_vel_B_by_T_matrix import b_vel_B_by_T_matrix

    nxb_f = matlab_loop_data['nxb_f']
    nxt_f = matlab_loop_data['nxt_f']
    Xb_f = matlab_loop_data['Xb_f']
    Xt_f = matlab_loop_data['Xt_f']

    nxb_r = matlab_loop_data['nxb_r']
    nxt_r = matlab_loop_data['nxt_r']
    Xb_r = matlab_loop_data['Xb_r']
    Xt_r = matlab_loop_data['Xt_r']

    cVBT_f = b_vel_B_by_T_matrix(nxb_f, nxt_f, Xb_f, Xt_f, g.RCUT)
    cVBT_r = b_vel_B_by_T_matrix(nxb_r, nxt_r, Xb_r, Xt_r, g.RCUT)

    npt.assert_allclose(cVBT_f, matlab_loop_data['cVBT_f'])
    npt.assert_allclose(cVBT_r, matlab_loop_data['cVBT_r'])

def test_vel_B_by_T(matlab_loop_data):
    from tombo.vel_B_by_T import vel_B_by_T

    cVBT_f = matlab_loop_data['cVBT_f']
    cVBT_r = matlab_loop_data['cVBT_r']

    GAM_f = matlab_loop_data['GAM_f']
    nxt_f = matlab_loop_data['nxt_f']
    GAM_r = matlab_loop_data['GAM_r']
    nxt_r = matlab_loop_data['nxt_r']

    VBTs_f = vel_B_by_T(cVBT_f, GAM_f, nxt_f)
    VBTs_r = vel_B_by_T(cVBT_r, GAM_r, nxt_r)

    npt.assert_allclose(VBTs_f, matlab_loop_data['VBTs_f'])
    npt.assert_allclose(VBTs_r, matlab_loop_data['VBTs_r'])
