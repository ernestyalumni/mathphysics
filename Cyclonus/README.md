# Protocol Buffers

## Installation

cf. https://github.com/protocolbuffers/protobuf/blob/main/src/README.md

```
sudo apt-get install autoconf automake libtool curl make g++ unzip
```

```
git clone https://github.com/protocolbuffers/protobuf.git
cd protobuf
git submodule update --init --recursive
./autogen.sh
```

```
 ./configure
 make -j$(nproc) # $(nproc) ensures it uses all cores for compilation
 make check
 sudo make install
 sudo ldconfig # refresh shared library cache.
```

By default, the protobuf package will be installed to `/usr/local`.

### Compilation

#### Python

We had done this in the `Cyclonus/Cyclonus` subdirectory.

```
$ protoc -I=../ --python_out=./ ../ProtoFiles/PhysicalConstants/Fundamental.proto
```

## Running Jupyter notebooks for Protocol Buffers

**tl;dr**

Do:

```
# In mathphysics/Cyclonus/
python -m venv ./venv

# In mathphysics/Cyclonus/
jupyter notebook
```

When running jupyter notebooks (in the terminal, you'd launch or run it by typing `jupyter notebook`), I found that it matters from where (filepath) you launch it, and also where you've placed your virtual environment contents (e.g. the `venv` subdirectory).

For instance, I've created a virtual environment (Python) directory with contents of my virtual environment specifically in `mathphysics/Cyclonus/`, named `mathphysics/Cyclonus/venv`. When I launch Jupyter notebooks, I run, exactly in `mathphysics/Cyclonus/`

```
jupyter notebook
```
and only then will the server, and hence the GUI on your browser recognize where `mathphysics/Cyclonus/venv` is. Otherwise, if you run `jupyter notebook` from somewhere else, say `mathphysics/` the kernel won't be able to find the Python virtual environment to run your notebook in.
