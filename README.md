# mathphysics
 Various notes and solutions on Math and Physics

## Contents and Table of Contents

### Notes and solutions only
| codename        | directory | Keywords | Description             | External links |
| --------------- | :------------------------------------- | :---------------------: | :------------------ | :-------------| 
| `the geometry of physics problems.tex` | `./LaTeX_and_pdfs/the geometry of physics problems`                     | Frankel, geometry, topology, physics,                           | Notes and solutions for Frankel's **The Geometry of Physics** | |
| `master.tex` | `./documents/mono-notes` | thermodynamics, statistical mechanics, Kittel, Kroemer | Mono notes, first part. Letter and widescreen from one source. | |

## Verified Rust QFT tools

The [exact Rust QFT kernel](Rust/qft-symbolics/README.md) starts the
Cadabra-inspired symbolic/Python migration: spinors, Clifford algebra,
a strict bounded typed-index kernel, selected Srednicki problems (including
formal-d gamma identities and evaluated chiral traces), and exact
half-collinear SMGA recurrence checks.
See its [corpus](Rust/qft-symbolics/corpus.json) and
[verification scope](Rust/qft-symbolics/docs/VERIFICATION.md).
[TheoPhysEval research-workflow lessons](agent/research/2026-10-16-theophyseval.md)
explain the verifier-first approach and the paper's limits.

## Frankel Lie groups — Rust and C++

[Scalar-generic Rust](Rust/lie-groups/README.md) and
[matching C++20](CppMath/LieGroups/README.md) implement numerical and exact symbolic
SO(3), so(3), SU(2), su(2), epsilon/cross products, exponentials, Pauli bases and
the double cover. Rust supports f32/f64 and software binary16; C++ supports
float/double/long double. Exact symbolic charts and polynomial Lie identities
are cross-checked with live SymPy, SageMath and Cadabra.
See the [source audit](Rust/lie-groups/docs/SOURCE-AUDIT.md) and
[verification contract](Rust/lie-groups/docs/VERIFICATION.md).
These are independently authored known-result tools, not complete CAS ports.

## Goldstein rotation libraries

The [Python rotation library and validation runner](Python/rotations/README.md),
[Rust API](Rust/lie-groups/docs/GOLDSTEIN-API.md), and
[C++ API](CppMath/LieGroups/README.md) implement the matrix, Euler-angle,
axis–angle, Hamilton-quaternion and Cayley–Klein dictionaries used in the
mono-notes. Independent numerical and exact symbolic checks cover signs,
composition order, singular branches, conjugation and Cayley–Hamilton.
See the [verification receipt](Rust/lie-groups/docs/GOLDSTEIN-VERIFICATION.md)
for executed checks, commands and domain limits.

## Interactive six-DOF rigid-body lab

Run `./run_rigid_body_lab.sh`, then open <http://127.0.0.1:8765>.
The Python backend advances fixed 1/240 s physical steps; the local WebGL
viewer provides spacecraft translation/rotation controls, free tumbling,
intermediate-axis instability, pause/step, and live conservation diagnostics.
See [the viewer and controls](JavaScript/rigid-body-viewer/README.md) and
[the dynamics model and validation](Python/rotations/DYNAMICS.md).
The model is an ideal constant-mass rigid body with explicit external wrenches;
measured spacecraft properties, gravity, flexibility, slosh, and actuator
hardware are separate model extensions.

# Creating and starting a virtual environment for Python 3

Create a directory for a virtual environment:

```
mathphysics]$ python3 -m venv ./venv/
```

Activate it:
```
mathphysics]$ source ./venv/bin/activate
```

Deactivate it:
```
deactivate
```

# Pip install (in the virtual environment) requirements.

Go to `Manifolds/` (you'll want the `requirements.txt` file accessible)

```
pip install -r requirements.txt
```

Running `pip freeze` **before** and after gives a good idea to the user of what `pip` packages have been installed.

## jupyter notebook in this virtual environment

You'll also want to do a 
```
pip install jupyter notebook
```
because you may be using a jupyter notebook for your local system; check this with
```
which jupyter
```

cf. https://www.codingforentrepreneurs.com/blog/install-jupyter-notebooks-virtualenv

If you want to use the virtual environments "version" or set of pip installed libraries, then from https://stackoverflow.com/questions/37891550/jupyter-notebook-running-kernel-in-different-env

```
pip install ipykernel

# and
# (but read further below because this might not be the right command for your setup)

python -m ipykernel install --user --name ENVNAME --display-name "Python (whatever you want to call it)"
```
e.g.

```
python -m ipykernel install --user --name venv --display-name "Python_venv"
```

Because I only have a **single Python 3 kernel**, the above command didn't work, but this did:

```
python -m ipykernel install --user
```
## Mathematical-physics master notes

The unified [collected notes and worked problems](documents/mono-notes/README.md)
combine algebra, geometry, mechanics, thermodynamics, and quantum/relativity
material in one source, with portrait and screen PDF editions.
