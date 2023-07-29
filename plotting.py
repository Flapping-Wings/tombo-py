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
    save: bool
        If `True`, save plot as a png.
        If `False`, open plot in interactive viewer.
    """
    fig, ax = plt.subplots()

    plot_mesh_2D_helper(ax, Xb, nXb, npoly, color='r')
    plot_mesh_2D_helper(ax, Xc, nXc, npoly, color='b')
    ax.axis('equal')

    if save:
        fig.savefig(f'{g.plot_folder}/mesh2d/mesh2d_{wing}.png')
        plt.close(fig)
    else:
        plt.show()

def plot_mesh_2D_helper(ax, X, nX, npoly, color):
    for i in range(nX):
        x = np.zeros(npoly + 1)
        y = np.zeros(npoly + 1)

        x[:npoly] = X[0, :npoly, i]
        y[:npoly] = X[1, :npoly, i]
        x[npoly] = X[0, 0, i]
        y[npoly] = X[1, 0, i]

        cx = X[0, npoly, i]
        cy = X[1, npoly, i]

        ax.plot(x, y, color, linewidth=2)
        ax.plot(cx, cy, 'o')


def plot_mesh_3D(wing, Xb, nXb, Nb, Xc, nXc, Nc, *, save=False):
    """
    Plot 3D view of wing mesh

    Parameters
    ----------
    wing: str
        'f' for front wing, or 'r' for rear wing
    Xb: ndarray[j, n, i]
        Border element coordinates
    nXb: int
        Number of border elements 
    Nb: ndarray[j, i]
        Unit normals to the border elements
    Xc: ndarray[j, n, i]
        Center element coordinates
    nXc: int
        Number of center elements
    Nc: ndarray[j, i]
        Unit normals to the center elements
    save: bool
        If `True`, save plot as a png.
        If `False`, open plot in interactive viewer.
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    plot_mesh_3D_helper(ax, Xb, nXb, Nb)
    plot_mesh_3D_helper(ax, Xc, nXc, Nc)
    ax.axis('equal')

    if save:
        fig.savefig(f'{g.plot_folder}/mesh3d/mesh3d_{wing}.png')
        plt.close()
    else:
        plt.show()

def plot_mesh_3D_helper(ax, X, nX, N):
    scale = 0.1
    Nline = np.zeros((2, 3))
    x = np.zeros(5)
    y = np.zeros(5)
    z = np.zeros(5)

    for i in range(nX):
        x[:4] = X[0, :4, i]
        y[:4] = X[1, :4, i]
        z[:4] = X[2, :4, i]
        
        x[4] = X[0, 0, i]
        y[4] = X[1, 0, i]
        z[4] = X[2, 0, i]

        ax.plot(x, y, z, 'k')
        
        Nline = np.array([X[:, 4, i], X[:, 4, i] + scale * N[:, i]])

        ax.plot(Nline[:, 0], Nline[:, 1], Nline[:, 2], color='r') 


def plot_airfoil_vel(XC, Vnc, NC, m, iwing, t, *, save=False):
    # End points for the normal velocity vector
    sf = 0.1    # Scale factor for the velocity plot
    xaif = XC[0,:]
    yaif = XC[1,:]
    zaif = XC[2,:]
    
    xtip = xaif + sf * Vnc * NC[0,:]
    ytip = yaif + sf * Vnc * NC[1,:]
    ztip = zaif + sf * Vnc * NC[2,:]
    
    # Plot normal velocity vectors at collocation points     
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    for i in range(len(xtip)):
      ax.plot([xaif[i], xtip[i]], [yaif[i], ytip[i]], [zaif[i], ztip[i]])

    ax.scatter(XC[0, :], XC[1, :], XC[2, :], marker='o')

    ax.set_title('Normal velocity vectors at collocation points')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1, 1, 1])
    ax.axis('equal')

    if save:
        labels = [['fr', 'fl'], ['rr', 'rl']]

        plt.savefig(f'{g.plot_folder}/airfoil_vel/airfoil_vel_{labels[m][iwing]}_{t:.4f}.png')
        plt.close(fig)
    else:
        plt.show()


def dummy():
    pass

# Relate each plot type to its corresponding function
plotting_funcs = {
    'mesh2d': plot_mesh_2D,
    'mesh3d': plot_mesh_3D,
    'airfoil_vel': plot_airfoil_vel,
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

def plot():
    create_directories(g.plot_folder)

    with np.load(f'{g.data_folder}/airfoil_vel/airfoil_vel_rl_0.0000.npz') as data:
        plotting_funcs['airfoil_vel'](*data.values(), save=False)


if __name__ == '__main__':
    plot()
