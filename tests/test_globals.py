"""Overwrite globals with test globals"""

import numpy as np
import tomli
import tombo.globals as g

with open('tests/test_config.toml', mode='rb') as file:
    config = tomli.load(file)

# General
# -------

g.solver = config['general']['solver']
g.idebg = config['general']['idebg']
g.save_data = config['general']['save_data']


# Plotting
# --------

g.output_folder = config['plotting']['output_folder']
g.data_folder = config['plotting']['data_folder']
g.plot_folder = config['plotting']['plot_folder']
g.labels = config['plotting']['labels']
g.plot_enabled = config['plotting']['plot_enabled']


# Time
# ----

g.dt = config['time']['dt']
g.nstep = config['time']['nstep']


# Body geometry
# -------------

g.twing = config['body_geometry']['twing']
g.nwing = config['body_geometry']['nwing']
g.delta_= config['body_geometry']['delta_']
g.b_f = config['body_geometry']['b_f']
g.b_r = config['body_geometry']['b_r']


# Wing geometry
# -------------

g.icamber = config['wing_geometry']['icamber']
g.acamber = config['wing_geometry']['acamber']

g.hfactor_f = config['wing_geometry']['hfactor_f']
g.wfactor_f = config['wing_geometry']['wfactor_f']
g.hfactor_r = config['wing_geometry']['hfactor_r']
g.wfactor_r = config['wing_geometry']['wfactor_r']

g.lt_f = config['wing_geometry']['lt_f']
g.lr_f = config['wing_geometry']['lr_f']
g.lt_r = config['wing_geometry']['lt_r']
g.lr_r = config['wing_geometry']['lr_r']

g.bang_f = config['wing_geometry']['bang_f']
g.bang_r = config['wing_geometry']['bang_r']

g.ielong = config['wing_geometry']['ielong']


# Wing motion
# -----------

g.phiT_ = np.array(config['wing_motion']['phiT_'])
g.phiB_ = np.array(config['wing_motion']['phiB_'])
g.a_    = np.array(config['wing_motion']['a_'])
g.beta_ = np.array(config['wing_motion']['beta_'])
g.f_    = np.array(config['wing_motion']['f_'])
g.gMax_ = np.array(config['wing_motion']['gMax_'])
g.p     = np.array(config['wing_motion']['p'])
g.rtOff = np.array(config['wing_motion']['rtOff'])
g.tau   = np.array(config['wing_motion']['tau'])
g.mpath = np.array(config['wing_motion']['mpath'])


# Fluid
# -----

g.rho_ = config['fluid']['rho_']
g.U_ = np.array(config['fluid']['U_'])


# Tolerance
# ---------------

g.RCUT = config['tolerance']['RCUT']
