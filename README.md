# mathphysics
 Various notes and solutions on Math and Physics

## Contents and Table of Contents

### Notes and solutions only
| codename        | directory | Keywords | Description             | External links |
| --------------- | :------------------------------------- | :---------------------: | :------------------ | :-------------| 
| `the geometry of physics problems.tex` | `./LaTeX_and_pdfs/the geometry of physics problems`                     | Frankel, geometry, topology, physics,                           | Notes and solutions for Frankel's **The Geometry of Physics** | |

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