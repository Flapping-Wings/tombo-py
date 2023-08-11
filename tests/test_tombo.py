import pytest
import numpy as np
import numpy.testing as npt
from scipy.io import loadmat
import tombo.globals as g
import test_globals

@pytest.fixture
def matlab_wing_data():
    return loadmat(f"tests/matlab_data/wing.mat")

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

def test_nd_data():
    pass

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

