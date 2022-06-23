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