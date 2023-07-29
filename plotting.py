import os
from pathlib import Path

def create_directories(base_path):
    os.makedirs(base_path, exist_ok=True)
    base_dir = Path(base_path)

    mesh2d_dir = base_dir / Path("mesh2d")
    if not mesh2d_dir.exists():
        mesh2d_dir.mkdir()

    mesh3d_dir = base_dir / Path("mesh3d")
    if not mesh3d_dir.exists():
        mesh3d_dir.mkdir()

    airfoil_vel_dir = base_dir / Path("airfoil_vel")
    if not airfoil_vel_dir.exists():
        airfoil_vel_dir.mkdir()

    GAMA_dir = base_dir / Path("GAMA")
    if not GAMA_dir.exists():
        GAMA_dir.mkdir()

    wake_dir = base_dir / Path("wake")
    if not wake_dir.exists():
        wake_dir.mkdir()

    force_dir = base_dir / Path("force")
    if not force_dir.exists():
        force_dir.mkdir()

    moment_dir = base_dir / Path("moment")
    if not moment_dir.exists():
        moment_dir.mkdir()