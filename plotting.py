import os
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import globals as g

def plot_mesh_2D(wing, Xb, nXb, Xc, nXc, npoly=4, *, save=False):
    """
    Plot 2D view of wing mesh

    Parameters
    ----------
    wing: str
        'f' for front wing, or 'r' for rear wing
    Xb: ndarray[j, n, i]
        Border element coordinates
    nXb: int
        Number of border elements
    Xc: ndarray[j, n, i]
        Center element coordinates
    nXc: int
        Number of center elements
    npoly: int
        Order of drawn polygons TODO Ask about this
    folder: str
        Subfolder to save plot in
    save: bool
        If `True`, save plot as a png.
        If `False`, open plot in interactive viewer.
    """
    fig, ax = plt.subplots()
    ax.axis('equal')
    plot_mesh_2D_helper(ax, Xb, nXb, npoly, 'r')
    plot_mesh_2D_helper(ax, Xc, nXc, npoly, 'b')

    if save:
        fig.savefig(f'{g.plot_folder}/mesh2d/mesh2d_{wing}.png')
        plt.close(fig)
    else:
        plt.show()

def plot_mesh_2D_helper(ax, Xn, nXn, npoly, color):
    for i in range(nXn):
        x = np.zeros(npoly + 1)
        y = np.zeros(npoly + 1)

        x[:npoly] = Xn[0, :npoly, i]
        y[:npoly] = Xn[1, :npoly, i]
        x[npoly] = Xn[0, 0, i]
        y[npoly] = Xn[1, 0, i]

        cx = Xn[0, npoly, i]
        cy = Xn[1, npoly, i]

        ax.plot(x, y, color, linewidth=2)
        ax.plot(cx, cy, 'o')


def dummy():
    pass

# Relate each plot type to its corresponding function
plotting_funcs = {
    'mesh2d': plot_mesh_2D,
    'mesh3d': dummy,
    'airfoil_vel': dummy,
    'GAMA': dummy,
    'wake': dummy,
    'force': dummy,
    'moment': dummy
}

def create_directories(base_path):
    os.makedirs(base_path, exist_ok=True)
    base_dir = Path(base_path)

    for key in plotting_funcs.keys():
        dir = base_dir / Path(key)
        if not dir.exists():
            dir.mkdir()

# def plot():
#     create_directories(g.plot_folder)

#     with np.load(f'{g.data_folder}/mesh2d/mesh2d_f.npz') as data:
#         plotting_funcs['mesh2d'](*data.values(), save=False)


# if __name__ == '__main__':
#     plot()
