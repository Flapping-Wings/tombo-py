import numpy as np
import numpy.typing as npt

# General configuration
# ---------------------

# stime: bool = True
"""
Time specification:
- multiple times (False),
- time marching (True)
"""
solver: bool = False
"""
Linear equation solver: 
- forward elim/backward sub (False),
- backslash (True)
"""


# Plotting
# --------

data_folder: str = 'output/data'
"""Folder for data used for plotting"""
plot_folder: str = 'output/plots'
"""Folder for generated plots"""
labels: list = [['fr', 'fl'], ['rr', 'rl']]
"""
Labels used in filenames
- 'fr': front right
- 'fl': front left
- 'rr': rear right
- 'rl': rear left
"""
idebg: bool = True
"""Toggle debug prints"""


# Time
# ----

dt: float = 0.1
"""Time increment"""
nstep: int = 4
"""Number of time steps to iterate through"""


# Body geometry
# -------------

twing: int = 4
"""Total number of wings"""
nwing: int = twing // 2
"""Number of wings in fore/rear"""
delta_: float = 0.0
"""Body angle in degrees"""
b_f: float = -1.5
"""Fore wing location in cm"""
b_r: float = 1.5
"""Rear wing location in cm"""


# Wing geometry
# -------------

# General

icamber: int = 0
"""
Camber direction: 
- `0`: no camber
- `1`: x-direction
- `2`: y-direction
- `3`: both directions"""
acamber: float = 0.2
"""Camber amplitude"""

hfactor_f: float = 0.1
"""Ratio of border element height to wing chord length (front)"""
wfactor_f: float = 3
"""Ratio of border element width to border element height (front)"""
hfactor_r: float = 0.1
"""Ratio of border element height to wing chord length (rear)"""
wfactor_r: float = 3
"""Ratio of border element width to border element height (rear)"""

lt_f: float = 2
"""Length of tapered section of the wing in cm (front)"""
lr_f: float = 2
"""Length of straight section of the wing in cm (front)"""
lt_r: float = 1
"""Length of tapered section of the wing in cm (rear)"""
lr_r: float = 1
"""Length of straight section of the wing in cm (rear)"""

bang_f: float = 30
"""Base angle (angle between tapered edge and centerline) of the wing in degrees (front)"""
bang_r: float = 30
"""Base angle (angle between tapered edge and centerline) of the wing in degrees (rear)"""

ielong: bool  = False
"""Toggle fixed number of border elements"""


# Wing motion
# -----------

phiT_: npt.NDArray[np.floating] = np.array([80.0, 80.0, 80.0, 80.0])
"""Top stroke angle in degrees"""
phiB_: npt.NDArray[np.floating] = np.array([-45.0, -45.0, -45.0, -45.0])
"""Bottom stroke angle in degrees"""
a_: npt.NDArray[np.floating] = np.array([0.0, 0.0, 0.0, 0.0])
"""Rotation axis offset in cm"""
beta_: npt.NDArray[np.floating] = np.array([90.0, 90.0, 90.0, 90.0])
"""Stroke plane angle in degrees wrt the body axis"""
f_: npt.NDArray[np.floating] = np.array([30.0, 30.0, 30.0, 30.0])
"""Flapping frequency in Hz"""
gMax_: npt.NDArray[np.floating] = np.array([30.0, 30.0, 30.0, 30.0])
"""Max rotation amplitude in degrees; actual rotation is `2*gMax`"""
p: npt.NDArray[np.floating] = np.array([5.0, 5.0, 5.0, 5.0])
"""Rotation speed parameter (nondimensional); `p >= 4`"""
rtOff: npt.NDArray[np.floating] = np.array([0.0, 0.0, 0.0, 0.0])
"""
Rotation timing offset (nondimensional):
- `-0.5 < rtOff < 0.5`
- `rtOff < 0` (advanced), 
- `rtOff = 0` (symmetric),
- `rtOff > 0` (delayed);
"""
tau: npt.NDArray[np.floating] = np.array([0.0, 0.0, 0.0, 0.0])
"""
Phase shift for the time (nondimensional):
- `0 <= tau < 2`
- `tau = 0` (start from top, start with down stroke)
- `0 < tau < 1` (start in between, start with down stroke)
- `tau = 1` (start from bottom, start wtih up stroke)
- `1 < tau < 2` (start in between, start with up stroke)
"""
mpath: npt.NDArray[np.integer] = np.array([0, 0, 0, 0])
"""
Motion path parameter:
- `0`: noTail Standard flapping(sinusoidal) & rotation (smoothed step func)
- `1`: `td = 1` (DUTail), `nhp = 4` (2 periods)     
- `2`: `td = 2` (UDTail), `nhp = 4` (2 periods) 
- `3`: `td = 1` (DUDUTail), `nhp = 8` (4 periods) 
- `4`: `td = 2` (UDUDTail), `nhp = 8` (4 periods) 
- `5`: ctct(UDTailUTail), nhp = 4
"""


# Fluid
# -----

rho_: float = 0.001225
"""Air density"""
U_: npt.NDArray[np.floating] = np.array([100.0, 0.0, 0.0])
"""
Ambient velocity in (x, y, z); can be interpreted as the flight velocity
when the wind is calm
"""


# Cutoffs (error)
# ---------------

RCUT: float = 1.0e-10
"""Distance between source and observation points to be judged as zero"""


"""Check config values"""

if np.any(p < 4):
    raise ValueError("p must >=4 for all wings")

if np.any(np.abs(rtOff) > 0.5):
    raise ValueError("-0.5 <= rtOff <= 0.5 must be satisfied for all wings")

if np.any((tau < 0) | (tau >= 2)):
    raise ValueError("0 <= tau < 2 must be satisfied for all wings")
