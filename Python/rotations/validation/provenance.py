"""Bind tested native binaries, source trees and Python/CAS environment."""
import hashlib
import json
import pathlib
import subprocess
import sys
import sympy
import numpy
import scipy
repo,run,rust,cpp,*extra_binaries=map(pathlib.Path,sys.argv[1:])
paths=[]
for relative in ('Python/rotations','Rust/lie-groups/src','Rust/lie-groups/tests','Rust/lie-groups/examples','CppMath/LieGroups/include','CppMath/LieGroups/tests','CppMath/LieGroups/examples'):
    root=repo/relative
    if root.exists(): paths.extend(p for p in root.rglob('*') if p.is_file() and p.suffix in ('.py','.rs','.hpp','.cpp','.sh','.toml','.md'))
paths += [repo/'Rust/lie-groups/Cargo.toml',repo/'Rust/lie-groups/Cargo.lock',repo/'CppMath/LieGroups/CMakeLists.txt',rust,cpp,*extra_binaries]
print(json.dumps({'interpreter':sys.executable,'python':sys.version,'sympy':sympy.__version__,'numpy':numpy.__version__,'scipy':scipy.__version__,'git':subprocess.check_output(['git','-C',str(repo),'rev-parse','HEAD'],text=True).strip(),'hashes':{str(path):hashlib.sha256(path.read_bytes()).hexdigest() for path in sorted(set(paths))},'scope':'Source and binary hashes; numerical/CAS evidence is bounded, not all-input numerical proof'},indent=2))
