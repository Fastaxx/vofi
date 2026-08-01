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

Tests are built by default and run with `ctest` from the build directory
(configure with `-DBUILD_TESTING=OFF` to skip them):

* `regression_2d3d` compares every dimension against `tests/golden.dat`.
  `tests/gen_golden` regenerates that file, and a diff in it is a claim that
  a change in behaviour was intended.
* `hyper_4d` checks the 4D kernel against closed forms.
* `stress_4d` is a randomised self-consistency sweep -- each cell against its
  sixteen children -- over balls, ellipsoids, planes, two-sheeted slabs, tori
  and growing spheres, with cell aspect ratios spanning three decades. Run it
  with more cases and another seed to search harder: `stress_4d 5000 7`.

## Dimensions

`vofi_get_cc` and `vofi_get_cell_type` take `ndim0` = 1, 2, 3 or 4. See
`include/vofi.h` for how the caller's arrays must be sized in each case --
in particular `xex` grows to 5 reals and `npt` to 6 ints in 4D.

### 4D

A 4D cell is a hypercube; the motivating case is a space-time slab, three
space axes and one time axis. It is handled exactly as 2D and 3D are, with
one more level of the same nested height-function quadrature:

    2D:  s -> root along p
    3D:  t -> s -> root along p
    4D:  u -> t -> s -> root along p

All four directions are ordered once per cell by |grad f|, largest first, so
`pdir` carries the height and `udir` is the outer sweep; every hyperplane
`u = const` inherits that frame, the way every plane of a 3D cell inherits
`(p,s,t)`. The `u` axis is subdivided where the cross-section changes
topology: the crossings of the eight cell edges parallel to `u` (internal
limits) and the extrema of `u` on the interface inside the two 3-faces normal
to `pdir` (external limits). The latter are found the way 3D finds its
tertiary extrema -- by marching on the centre of the reference phase until
the section pinches off -- one dimension up, so the chord midpoint becomes
the centre of two chords.

Because the whole cell shares one frame, the interface is a single graph
`x_p = H(s,t,u)`, and its 3-volume is integrated on the same nodes as the
hypervolume rather than reconstructed from facets.

Accuracy, from `tests/test_4d.c`:

| quantity | test | error |
|---|---|---|
| fraction | hyperplane cut, 90 orientations/offsets | 3.3e-15 |
| fraction | 4-ball, 6^4 grid sum vs pi^2 R^4 / 2 | 4.4e-16 |
| fraction | space-time slab, growing ball | 9.4e-13 |
| centroid | axis-aligned cuts, against the exact box | 2.2e-16 |
| centroid | diagonal section of the unit hypercube (23/60) | 2.2e-16 |
| measure | hyperplane cut, relative | 1.7e-12 |
| measure | 3-sphere, 6^4 grid sum vs 2 pi^2 R^3 | 5.5e-14 |

Cost is about 1.9e4 function evaluations per cut cell, or 4.7e4 with the
interface measure requested. Most of that is the nested quadrature itself --
the limits and topology machinery is 2.3% of it -- so the way to trade
accuracy for speed is the `npt` array: capping every level at 10 points costs
2.9x less and still gives ~6e-12 on the 4-ball grid sum, and at 8 points,
4.4x less for ~4e-10.

As in 2D and 3D, an interface that does not reach the boundary of a cell is
not seen: the grid has to resolve it. And as in 3D, a single frame for the
whole cell is worse than a per-hyperplane one where the cross-section turns
sharply inside one cell -- that is the trade 2D and 3D already make.
