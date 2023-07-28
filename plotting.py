import os
from pathlib import Path

def create_directories(base_path):
    os.makedirs(base_path, exist_ok=True)
    base_dir = Path(base_path)

    mesh_dir = base_dir / Path("mesh")
    if not mesh_dir.exists():
        mesh_dir.mkdir()

    debug_dir = base_dir / Path("debug")
    if not debug_dir.exists():
        debug_dir.mkdir()

    wake_dir = base_dir / Path("wake")
    if not wake_dir.exists():
        wake_dir.mkdir()

    f_and_m = base_dir / Path("f&m")
    if not f_and_m.exists():
        f_and_m.mkdir()