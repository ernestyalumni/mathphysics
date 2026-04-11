from pathlib import Path
import sys

src_path = Path(__file__).resolve().parents[1] / "src"

if src_path.exists():
    if src_path not in sys.path:
        sys.path.append(str(src_path))
        print(f"Added {src_path} to sys.path")
else:
    raise FileNotFoundError(
        f"MHV amplitudes src path {src_path} does not exist")