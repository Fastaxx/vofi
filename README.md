# VOFI Library

The  VOFI library initializes the volume fraction scalar field in a computational mesh, given an analytic expression f(x,y,z) for the interface. The grid cells can be cuboids of variable size. The implicit function f(x,y,z) is specified by the user: the interface is represented by the zero level set, f(x,y,z)=0, and the reference phase is located where f(x,y,z)<0. Each routine in the directory 'src' contains a brief description of what it does and of the I/O variables. 

## Installation

```shell
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/where/to/install -DBUILD_SHARED_LIBS=ON
make
make install
```

