import os
from pathlib import Path

def dummy():
    pass

# Relate each plot type to its corresponding function
plotting_funcs = {
    'mesh2d': dummy,
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
