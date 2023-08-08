import os
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

import globals as g

def plot_mesh_2D(Xb, nXb, Xc, nXc, npoly=4, *, filename, save):
    """
    Plot 2D view of wing mesh

    Parameters
    ----------
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
    filename: str
        Filename (stem) of saved plot
    save: bool
        If `True`, save plot as a png.
        If `False`, open plot in interactive viewer.
    """
    fig, ax = plt.subplots()

    plot_mesh_2D_helper(ax, Xb, nXb, npoly, color='r')
    plot_mesh_2D_helper(ax, Xc, nXc, npoly, color='b')
    ax.axis('equal')

    if save:
        fig.savefig(f'{g.plot_folder}/mesh2d/{filename}.png')
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

def plot_mesh_3D(Xb, nXb, Nb, Xc, nXc, Nc, *, filename, save):
    """
    Plot 3D view of wing mesh

    Parameters
    ----------
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
    filename: str
        Filename (stem) of saved plot
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
        fig.savefig(f'{g.plot_folder}/mesh3d/{filename}.png')
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

def plot_airfoil_vel(Vnc, XC, NC, *, filename, save):
    """
    Plot normal velocity of the airfoil at collocation points

    Parameters
    ----------
    Vnc: ndarray[1, nXC]
        Normal velocity (global) 
    XC: ndarray[j, i]
        Collocation points (global)
    NC[j, i]
        Unit normals at collocation points (global)
    filename: str
        Filename (stem) of saved plot
    save: bool
        If `True`, save plot as a png.
        If `False`, open plot in interactive viewer.
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    scale_factor = 0.1
    plot_velocity(ax, scale_factor, Vnc, XC, NC)
    ax.set_title('Normal velocity vectors at collocation points')

    if save:
        plt.savefig(f'{g.plot_folder}/airfoil_vel/{filename}.png')
        plt.close(fig)
    else:
        plt.show()

def plot_GAMA(GAMA, XC, NC, *, filename, save):
    """
    Plot GAMA at the collocation points using the normal direction

    Parameters
    ----------
    GAMA: 
        TODO
    XC: ndarray[j, i]
        Collocation points (global)
    NC[j, i]
        Unit normals at collocation points (global)
    filename: str
        Filename (stem) of saved plot
    save: bool
        If `True`, save plot as a png.
        If `False`, open plot in interactive viewer.
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')
    
    scale_factor = 1.0
    plot_velocity(ax, scale_factor, GAMA, XC, NC)
    ax.set_title('GAMA at collocation points')

    if save:
        plt.savefig(f'{g.plot_folder}/GAMA/{filename}.png')
        plt.close(fig)
    else:
        plt.show()

def plot_velocity(ax, scale_factor, vel, XC, NC):
    """Helper for `plot_airfoil_vel` and `plot_GAMA`"""
    # End points for the normal velocity vector
    xaif = XC[0]
    yaif = XC[1]
    zaif = XC[2]
    
    xtip = xaif + scale_factor * vel * NC[0]
    ytip = yaif + scale_factor * vel * NC[1]
    ztip = zaif + scale_factor * vel * NC[2]

    for i in range(len(xtip)):
      ax.plot([xaif[i], xtip[i]], [yaif[i], ytip[i]], [zaif[i], ztip[i]])

    ax.scatter(xaif, yaif, zaif, marker='o')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_box_aspect([1, 1, 1])
    ax.axis('equal')

def plot_wake(nXb_f, nXw_f, Xb_f, Xw_f, nXb_r, nXw_r, Xb_r, Xw_r, *, filename, save):
    """
    Plot  wake elements along with the original border elements
    
    Parameters
    ----------
    nXb_f: int
        Number of border vortices (front)
    nXw_f: int
        Number of wake vortices (front)
    Xb_f: ndarray[j, inode, iXb, w]
        Coordinate j of node point for the border iXb (front)
    Xw_f: ndarray[j, inode, iXw, w]
        Coordinate j of node point for wake element iXb (front)
    nXb_r: int
        Number of border vortices (rear)
    nXw_r: int
        Number of wake vortices (rear)
    Xb_r: ndarray[j, inode, iXb, w]
        Coordinate j of node point for the border iXb (rear)
    Xw_r: ndarray[j, inode, iXw, w]
        Coordinate j of node point for wake element iXb (rear)
    filename: str
        Filename (stem) of saved plot
    save: bool
        If `True`, save plot as a png.
        If `False`, open plot in interactive viewer.
    """
    fig = plt.figure()
    ax = fig.add_subplot(projection='3d')

    plot_wing_set(ax, nXb_f, nXw_f, Xb_f, Xw_f)
    plot_wing_set(ax, nXb_r, nXw_r, Xb_r, Xw_r)
    ax.axis('equal')

    if save:      
        plt.savefig(f'{g.plot_folder}/wake/{filename}.png')
        plt.close(fig)
    else:
        plt.show()

def plot_wing_set(ax, nXb, nXw, Xb, Xw):
    """Helper for `plot_wake`"""
    # Original border elements
    for w in range(g.nwing):
        for i in range(nXb):
            x = [Xb[0,0,i,w], Xb[0,1,i,w], Xb[0,2,i,w], Xb[0,3,i,w], Xb[0,0,i,w]]
            y = [Xb[1,0,i,w], Xb[1,1,i,w], Xb[1,2,i,w], Xb[1,3,i,w], Xb[1,0,i,w]]
            z = [Xb[2,0,i,w], Xb[2,1,i,w], Xb[2,2,i,w], Xb[2,3,i,w], Xb[2,0,i,w]]

            ax.plot(x, y, z, 'r')
    
    # Wake elements
    for w in range(g.nwing):
        for i in range(nXw):
            x = [Xw[0,0,i,w], Xw[0,1,i,w], Xw[0,2,i,w], Xw[0,3,i,w], Xw[0,0,i,w]]
            y = [Xw[1,0,i,w], Xw[1,1,i,w], Xw[1,2,i,w], Xw[1,3,i,w], Xw[1,0,i,w]]
            z = [Xw[2,0,i,w], Xw[2,1,i,w], Xw[2,2,i,w], Xw[2,3,i,w], Xw[2,0,i,w]]

            if w == 0:
                ax.plot(x, y, z, 'k')
            else:
                ax.plot(x, y, z, 'b')

def plot_force(times, force, *, filename, save):
    """
    Plot force on the wings throughout the simulation

    Parameters
    ----------
    times: ndarray
        Points in time that were simulated
    force: ndarray
        Force on the wings at each corresponding point in time    
    """
    fig = plt.figure()
    plt.plot(times, force, 'x-k')
    plt.grid(True)

    if save:
        plt.savefig(f'{g.plot_folder}/force/{filename}.png')
        plt.close(fig)
    else:
        plt.show()

def plot_moment(times, moment, *, filename, save):
    """
    Plot moment on the wings throughout the simulation

    Parameters
    ----------
    times: ndarray
        Points in time that were simulated
    moment: ndarray
        Moment on the wings at each corresponding point in time    
    """
    fig = plt.figure()
    plt.plot(times, moment, 'o-r')
    plt.grid(True)

    if save:
        plt.savefig(f'{g.plot_folder}/moment/{filename}.png')
        plt.close(fig)
    else:
        plt.show()

# Relate each plot type to its corresponding function
plotting_funcs = {
    'mesh2d': plot_mesh_2D,
    'mesh3d': plot_mesh_3D,
    'airfoil_vel': plot_airfoil_vel,
    'GAMA': plot_GAMA,
    'wake': plot_wake,
    'force': plot_force,
    'moment': plot_moment
}

def create_directories(base_path):
    os.makedirs(base_path, exist_ok=True)
    base_dir = Path(base_path)

    for key in plotting_funcs.keys():
        dir = base_dir / Path(key)
        if not dir.exists():
            dir.mkdir()

def make_plot(full_path, save=True):
    path = Path(full_path)
    plot_type = path.parts[-2]

    with np.load(full_path) as data:
        plotting_funcs[plot_type](*data.values(),
                                  filename=path.stem,
                                  save=save)

def view_plot(path):
    """Generate single plot in interactive viewer"""
    make_plot(path, save=False)

def generate_plots(dir, chunksize=10):
    """Generate and save plot for each data file in dir"""
    # Construct list of all data files in directory
    data_files = [os.path.join(root, file)
                    for root, _, files in os.walk(dir)
                    for file in files]
    
    # Create plots
    if len(data_files) < chunksize:
        # Force map to make calls with empty loop
        for _ in map(make_plot, data_files):
            pass
    else:
        with Pool() as pool:
            # Force imap_unordered to make calls with empty loop
            for _ in pool.imap_unordered(make_plot, data_files, chunksize=chunksize):
                pass

def main(path):
    CHUNKSIZE = 10

    create_directories(g.plot_folder)

    if os.path.isfile(path):
        view_plot(path)
    elif os.path.isdir(path):
        generate_plots(path)
    else:
        raise ValueError("argument must be path to a directory or a file")

if __name__ == '__main__':
    main(sys.argv[1])
