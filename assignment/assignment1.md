# Assignment 1: Exterior Calculus

## Written

## Coding

### cotan

cotangent formula
$$
cot(A)=\frac{u\cdot v}{||u\times v||}
$$

```c++
double VertexPositionGeometry::cotan(Halfedge he) const {
    if (!he.isInterior()) return 0;

    Vector3 a = inputVertexPositions[he.next().tipVertex()],
            b = inputVertexPositions[he.tailVertex()],
            c = inputVertexPositions[he.tipVertex()];

    Vector3 u=b-a,
            v=c-a;

    return dot(u,v)/cross(u,v).norm();
}
```

### barycentricDualArea

```c++
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {
    double area = 0.;

    for (Face f : v.adjacentFaces()) area += faceArea(f) / 3;

    return area;
}
```

### buildHodgeStar0Form

```c++
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Tri;
```



```c++
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    std::vector<Tri> coefficients;

    for (Vertex v : mesh.vertices())
        coefficients.push_back(Tri(v.getIndex(), v.getIndex(), barycentricDualArea(v)));

    SpMat m(mesh.nVertices(), mesh.nVertices());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m;
}
```

### buildHodgeStar1Form

```c++
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    std::vector<Tri> coefficients;

    for (Edge e : mesh.edges()) {
        double area = (cotan(e.halfedge()) + cotan(e.halfedge().twin())) / 2;
        coefficients.push_back(Tri(e.getIndex(), e.getIndex(), area));
    }

    SpMat m(mesh.nEdges(), mesh.nEdges());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m; 
}
```

### buildHodgeStar2Form

```c++
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    std::vector<Tri> coefficients;

    for (Face f : mesh.faces())
        coefficients.push_back(Tri(f.getIndex(), f.getIndex(), 1/faceArea(f)));

    SpMat m(mesh.nFaces(), mesh.nFaces());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m;
}
```

### buildExteriorDerivative0Form

```C++
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
    std::vector<Tri> coefficients;

    for (Edge e : mesh.edges()) {
        coefficients.push_back(Tri(e.getIndex(), e.firstVertex().getIndex(), -1));
        coefficients.push_back(Tri(e.getIndex(), e.secondVertex().getIndex(), 1));
    }

    SpMat m(mesh.nEdges(), mesh.nVertices());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m;
}
```

### buildExteriorDerivative1Form

```c++
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    std::vector<Tri> coefficients;

    for (Face f : mesh.faces()) {
        for (Halfedge he : f.adjacentHalfedges()) {
            double orientation = he.getIndex() < he.twin().getIndex() ? 1 : - 1;
            coefficients.push_back(Tri(f.getIndex(), he.edge().getIndex(), orientation));
        }
    }


    SpMat m(mesh.nFaces(), mesh.nEdges());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m;
}
```

