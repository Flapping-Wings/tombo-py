# tombo-py

`tombo-py` is a simulation of flapping wings in 3D, for the purposes of studying insect flight. The original simulation was written by[ Dr. Mitsunori (Mitch) Denda](https://mae.rutgers.edu/mitsunori-mitch-denda) in MATLAB, and this repository is an optimized port of it to Python.

## Installation

Clone the repo and run the following commands in the directory (it is highly recommended that you use a Python virtual environment):

```shell
pip install -r requirements.txt
python3 -m build

# For users
pip install .
# For contributors
pip install --editable .[test]
```

If you do not have an interactive backend for Matplotlib installed, run
```
sudo apt install python3-tk
```

## Usage

`tombo-py` comes with a command line interface and is primarily used with subcommands. Run `tombo -h` to see the list of subcommands. The help option is available for each subcommand as well.

**IMPORTANT**: `tombo-py` relies on Numba-compiled functions that are cached on the first run for speed. The first simulation will therefore be much slower than subsequent simulations. We suggest first running a simulation with the `config.toml` that is shipped with `tombo-py`, which provides a minimal example, before running more intensive simulations.

### `sim`
Run the simulation, as specified in `config.toml`
```
tombo sim
```

### `plot`
Generates and saves plots using saved data from the simulation. If passed a path to a directory, it will generate plots from the data files in it. By default, it only plots those enabled in `config.toml`, but you can override this with the `--all` option.
```shell
# Generate plots specified in config.toml
tombo plot
# Generate plots only for force
tombo plot output/data/force
# Ignore settings in config.toml and generate all plots
tombo plot --all
```

### `simplot`
Convenience command to run simulation and generate plots at once. Simulation and plotting settings are taken from `config.toml`, but the `--all` option is supported for plotting. This subcommand does not take any arguments, as the plots are assumed to be for the simulation that was just run.
```shell
tombo simplot
```

### `view`
Opens a plot in an interactive viewer. It requires the path to the data file to be viewed to be passed as an argument.
```shell
tombo view output/data/wake/wake_0.npz
```

## Configuration
Settings for simulation and plotting can be configured in `config.toml`. Some of the user-relevant settings are described below.

### `plotting.plot_enabled`
Each of the boolean entries in `plotting.plot_enabled` acts as a toggle that enables or disables the plots of that category from being generated when calling the plotting comamnds.

### `wing_geometry.hfactor` and `wing_geometry.wfactor`
The overall resolution of the simulation is controlled through the resolution of the wing mesh. The mesh is made up of the border elements, which are user-configurable via the `wing_geometry.hfactor` and `wing_geometry.wfactor` settings, and the center elements, which are automatically determined by the border elements.

`hfactor` is the ratio of the height of each border element to the chord length of the wing. Smaller values lead to a higher resolution.

`wfactor` is the ratio of the width of each border element to its height. Therefore, with the value of `3` that `config.toml` ships with, each border element is a 3x1 rectangle. We recommend running the simulation once with the shipped settings to cache the compiled functions and then changing `wfactor` to `1` for a square mesh.

## Miscellaneous

Early development of `tombo-py` was done in [this repo](https://github.com/Flapping-Wings/Flapping-Wings). Refer to that repo if documentation of old pull requests or issues is needed.