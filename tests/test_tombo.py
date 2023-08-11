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

@pytest.fixture(params=[1, 2, 3, 4])
def matlab_loop_data(request):
    return loadmat(f"tests/matlab_data/data{request.param}.mat",
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


def test_wing_total():
    pass

def test_lr_set_matrix():
    pass

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

