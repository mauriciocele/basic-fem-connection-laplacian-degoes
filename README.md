# basic-fem-connection-laplacian-degoes
Basic FEM to compute a vector field through parallel transporting some input vectors using connection laplacian, based off Vector Heat Method


## Geometric algebra version and the notes

`src/GAConnection.h` implements the method in geometric algebra (tangent bivectors, rotors between tangent
planes, the connection as a complex number of the face plane, hat gradients as a reciprocal frame), and
`src/ga_fem.cpp` uses it without graphics: it prints the worked examples of the notes, solves the cactus with the
real 2x2-block system of `FEM.cpp` and with the equivalent complex Hermitian system (they agree to 1e-15), and
validates the result on a sphere against exact parallel transport along great circles.

```sh
cd src
c++ -std=c++17 -O2 -I/path/to/eigen ga_fem.cpp HalfEdge/Mesh.cpp -o ga_fem
./ga_fem ../fem-basic-2d/data                                  # data for the figures
cd ../fem-basic-2d/figures && python3 make_figures.py ../data  # 3D figures as TikZ
```

The notes, `fem-basic-2d/fembasic2d.tex`, derive the method from the PDE to the code (Part I, matrices) and again
in geometric algebra with the intrinsic geometry made explicit (Part II); build them with pdflatex and bibtex from
`fem-basic-2d/` (the code listings are read from `../src`).
