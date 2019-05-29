# Here is something about the Version 5 of HanFem
Firstly, this package will still offer support for simple Rectanglar shape problem with P1/P2 method.

Secondly, the main idea about this FEM realization package is listed below:
1. Main Process: Generate Mesh, Define the FE space, Ensemble the matrices, Solve the problem;
2. The rotation information will be recorded as the space was generated.

Thirdly, in this version, things will be improved further, and here is some major improvements expected:
1. Rotation info will be sperated with the edge info
2. The function used in matrices ensembling process will be vectors of nodal values;
3. The cost of matrices ensembling will be reduced without complicate the structure too much;
4. The analysis of computation complexity will be attached;
5. Refine.m will be checked; Other simple shape problem may be studied.