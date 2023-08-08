import numpy as np
import tomli

with open('config.toml', mode='rb') as file:
    config = tomli.load(file)

# General
# -------

solver = config['general']['solver']
idebg = config['general']['idebg']


# Plotting
# --------

data_folder = config['plotting']['data_folder']
plot_folder = config['plotting']['plot_folder']
labels = config['plotting']['labels']
plot_enabled = config['plotting']['plot_enabled']


# Time
# ----

dt = config['time']['dt']
nstep = config['time']['nstep']


# Body geometry
# -------------

twing = config['body_geometry']['twing']
nwing = config['body_geometry']['nwing']
delta_= config['body_geometry']['delta_']
b_f = config['body_geometry']['b_f']
b_r = config['body_geometry']['b_r']


# Wing geometry
# -------------

icamber = config['wing_geometry']['icamber']
acamber = config['wing_geometry']['acamber']

hfactor_f = config['wing_geometry']['hfactor_f']
wfactor_f = config['wing_geometry']['wfactor_f']
hfactor_r = config['wing_geometry']['hfactor_r']
wfactor_r = config['wing_geometry']['wfactor_r']

lt_f = config['wing_geometry']['lt_f']
lr_f = config['wing_geometry']['lr_f']
lt_r = config['wing_geometry']['lt_r']
lr_r = config['wing_geometry']['lr_r']

bang_f = config['wing_geometry']['bang_f']
bang_r = config['wing_geometry']['bang_r']

ielong = config['wing_geometry']['ielong']


# Wing motion
# -----------

phiT_ = np.array(config['wing_motion']['phiT_'])
phiB_ = np.array(config['wing_motion']['phiB_'])
a_    = np.array(config['wing_motion']['a_'])
beta_ = np.array(config['wing_motion']['beta_'])
f_    = np.array(config['wing_motion']['f_'])
gMax_ = np.array(config['wing_motion']['gMax_'])
p     = np.array(config['wing_motion']['p'])
rtOff = np.array(config['wing_motion']['rtOff'])
tau   = np.array(config['wing_motion']['tau'])
mpath = np.array(config['wing_motion']['mpath'])


# Fluid
# -----

rho_ = config['fluid']['rho_']
U_ = np.array(config['fluid']['U_'])


# Tolerance
# ---------------

RCUT = config['tolerance']['RCUT']


"""Check config values"""

if np.any(p < 4):
    raise ValueError("p must >=4 for all wings")

if np.any(np.abs(rtOff) > 0.5):
    raise ValueError("-0.5 <= rtOff <= 0.5 must be satisfied for all wings")

if np.any((tau < 0) | (tau >= 2)):
    raise ValueError("0 <= tau < 2 must be satisfied for all wings")
