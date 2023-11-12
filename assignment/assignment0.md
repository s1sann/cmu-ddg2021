# Assignment 0: Combinatorial Surfaces

## Written

## Coding

### buildVertexEdgeAdjacencyMatrix

```c++
typedef Eigen::Triplet<size_t> Tri;
typedef SparseMatrix<size_t> SpMat;
```



```c++
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();

    std::vector<Tri> cofficients;
    SpMat A(mesh->nEdges(), mesh->nVertices());

    for (Edge e : mesh->edges()) {
        size_t idx = geometry->edgeIndices[e];
        cofficients.push_back(Tri(idx, geometry->vertexIndices[e.firstVertex()], 1));
        cofficients.push_back(Tri(idx, geometry->vertexIndices[e.secondVertex()], 1));
    }

    A.setFromTriplets(cofficients.begin(), cofficients.end());

    return A;
}
```

### buildEdgeFaceAdjacencyMatrix

```c++
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {
    geometry->requireFaceIndices();
    geometry->requireEdgeIndices();

    std::vector<T> cofficients;
    SpMat A(mesh->nFaces(), mesh->nEdges());

    for (Face f : mesh->faces()) {
        size_t idx = geometry->faceIndices[f];
        for (Edge e : f.adjacentEdges()) cofficients.push_back(T(idx, geometry->edgeIndices[e], 1));
    }

    A.setFromTriplets(cofficients.begin(), cofficients.end());

    return A; 
}
```

### buildVertexVector

```c++
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {
    Vector<size_t> v = Vector<size_t>::Zero(mesh->nVertices());

    for (size_t idx : subset.vertices) v[idx] = 1;

    return v;
}
```

### buildEdgeVector

```c++
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {
    Vector<size_t> e = Vector<size_t>::Zero(mesh->nEdges());

    for (size_t idx : subset.edges) e[idx] = 1;

    return e;
}
```

### buildFaceVector

```c++
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {
    Vector<size_t> f = Vector<size_t>::Zero(mesh->nFaces());

    for (size_t idx : subset.edges) f[idx] = 1;

    return f;
}
```

### star

```c++
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {
    SpMat vertexEdgeSpMat = buildVertexEdgeAdjacencyMatrix();
    SpMat faceEdgeSpMat = buildFaceEdgeAdjacencyMatrix();

    MeshSubset star = subset.deepCopy();

    for (size_t idx : subset.vertices)
        for (int k = 0; k < vertexEdgeSpMat.outerSize(); k++)
            for (SpMat::InnerIterator it(vertexEdgeSpMat, k); it; ++it)
                if (it.col() == int(idx)) star.addEdge(it.row());

    for (size_t idx : star.edges)
        for (int k = 0; k < faceEdgeSpMat.outerSize(); k++)
            for (SpMat::InnerIterator it(faceEdgeSpMat, k); it; ++it)
                if (it.col() == int(idx)) star.addFace(it.row());

    return star;
}

```

### closure

```c++
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {
    SpMat vertexEdgeSpMat = buildVertexEdgeAdjacencyMatrix();
    SpMat faceEdgeSpMat = buildFaceEdgeAdjacencyMatrix();

    MeshSubset closure = subset.deepCopy();

    for (size_t idx : subset.faces)
        for (int k = 0; k < faceEdgeSpMat.outerSize(); k++)
            for (SpMat::InnerIterator it(faceEdgeSpMat, k); it; ++it)
                if (it.row() == int(idx)) closure.addEdge(it.col());

    for (size_t idx : closure.edges)
        for (int k = 0; k < vertexEdgeSpMat.outerSize(); k++)
            for (SpMat::InnerIterator it(vertexEdgeSpMat, k); it; ++it)
                if (it.row() == int(idx)) closure.addVertex(it.col());

    return closure;
}
```

### link

```c++
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {
    MeshSubset link;
    link = closure(star(subset));
    link.deleteSubset(star(closure(subset)));

    return link;
}
```

### isComplex

```c++
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {
    if (closure(subset).equals(subset)) return true;
    return false;
}
```

### isPureComplex

```c++
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {
    SpMat vertexEdgeSpMat = buildVertexEdgeAdjacencyMatrix();
    SpMat faceEdgeSpMat = buildFaceEdgeAdjacencyMatrix();   

    MeshSubset isPure;

    if (!isComplex(subset)) {
        return -1;
    } else if (!subset.faces.empty()) {
        isPure = subset;
        isPure.deleteEdges(subset.edges);
        isPure.deleteVertices(subset.vertices);
        isPure = closure(isPure);

        if (isPure.edges == subset.edges && isPure.vertices == subset.vertices)
            return 2;
        else
            return 1;
    } else if (!subset.edges.empty()) {
        isPure = subset;
        isPure.deleteVertices(subset.vertices);
        isPure = closure(isPure);

        if (isPure.vertices == subset.vertices)
            return 1;
        else
            return -1;
    } else if (!subset.vertices.empty())
        return 0;

    return -1;
}
```

### boundary

```c++
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {
    SpMat vertexEdgeSpMat = buildVertexEdgeAdjacencyMatrix();
    SpMat faceEdgeSpMat = buildFaceEdgeAdjacencyMatrix();

    MeshSubset boundary;

    if (isPureComplex(subset) == 2) {
        std::vector<size_t> allEdge;
        std::set<size_t> boundaryEdge;

        for (size_t idx : subset.faces)
            for (int k = 0; k < faceEdgeSpMat.outerSize(); k++)
                for (SpMat::InnerIterator it(faceEdgeSpMat, k); it; ++it)
                    if (it.row() == int(idx)) allEdge.push_back(it.col());

        for (size_t idx : subset.edges)
            if (std::count(allEdge.begin(), allEdge.end(), idx) == 1) boundaryEdge.insert(idx);

        boundary.addEdges(boundaryEdge);
        boundary = closure(boundary);
        return boundary;
    } else if (isPureComplex(subset) == 1) {
        std::vector<size_t> allVertex;
        std::set<size_t> boundaryVertex;

        for (size_t idx : subset.edges)
            for (int k = 0; k < vertexEdgeSpMat.outerSize(); k++)
                for (SpMat::InnerIterator it(vertexEdgeSpMat, k); it; ++it)
                    if (it.row() == int(idx)) allVertex.push_back(it.col());

        for (size_t idx : subset.edges)
            if (std::count(allVertex.begin(), allVertex.end(), idx) == 1) boundaryVertex.insert(idx);

        boundary.addEdges(boundaryVertex);
        boundary = closure(boundary);
        return boundary;
    }

    return boundary;
}
```

