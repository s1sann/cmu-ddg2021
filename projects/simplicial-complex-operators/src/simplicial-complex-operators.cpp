// Implement member functions for SimplicialComplexOperators class.
#include "simplicial-complex-operators.h"
#include <Eigen>
#include <iostream>
#include <vector>

typedef Eigen::Triplet<size_t> Tri;
typedef SparseMatrix<size_t> SpMat;

using namespace geometrycentral;
using namespace geometrycentral::surface;

/*
 * Assign a unique index to each vertex, edge, and face of a mesh.
 * All elements are 0-indexed.
 *
 * Input: None. Access geometry via the member variable <geometry>, and pointer to the mesh via <mesh>.
 * Returns: None.
 */
void SimplicialComplexOperators::assignElementIndices() {

    // Needed to access geometry->vertexIndices, etc. as cached quantities.
    // Not needed if you're just using v->getIndex(), etc.
    geometry->requireVertexIndices();
    geometry->requireEdgeIndices();
    geometry->requireFaceIndices();

    // You can set the index field of a vertex via geometry->vertexIndices[v], where v is a Vertex object (or an
    // integer). Similarly you can do edges and faces via geometry->edgeIndices, geometry->faceIndices, like so:
    size_t idx = 0;
    for (Vertex v : mesh->vertices()) {
        idx = geometry->vertexIndices[v];
    }

    for (Edge e : mesh->edges()) {
        idx = geometry->edgeIndices[e];
    }

    for (Face f : mesh->faces()) {
        idx = geometry->faceIndices[f];
    }

    // You can more easily get the indices of mesh elements using the function getIndex(), albeit less efficiently and
    // technically less safe (although you don't need to worry about it), like so:
    //
    //      v.getIndex()
    //
    // where v can be a Vertex, Edge, Face, Halfedge, etc. For example:

    for (Vertex v : mesh->vertices()) {
        idx = v.getIndex(); // == geometry->vertexIndices[v])
    }

    // Geometry Central already sets the indices for us, though, so this function is just here for demonstration.
    // You don't have to do anything :)
}

/*
 * Construct the unsigned vertex-edge adjacency matrix A0.
 *
 * Input:
 * Returns: The sparse vertex-edge adjacency matrix which gets stored in the global variable A0.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildVertexEdgeAdjacencyMatrix() const {

    // TODO
    // Note: You can build an Eigen sparse matrix from triplets, then return it as a Geometry Central SparseMatrix.
    // See <https://eigen.tuxfamily.org/dox/group__TutorialSparse.html> for documentation.
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

    return A; // placeholder
}

/*
 * Construct the unsigned face-edge adjacency matrix A1.
 *
 * Input:
 * Returns: The sparse face-edge adjacency matrix which gets stored in the global variable A1.
 */
SparseMatrix<size_t> SimplicialComplexOperators::buildFaceEdgeAdjacencyMatrix() const {

    // TODO
    geometry->requireFaceIndices();
    geometry->requireEdgeIndices();

    std::vector<T> cofficients;
    SpMat A(mesh->nFaces(), mesh->nEdges());

    for (Face f : mesh->faces()) {
        size_t idx = geometry->faceIndices[f];
        for (Edge e : f.adjacentEdges()) cofficients.push_back(T(idx, geometry->edgeIndices[e], 1));
    }

    A.setFromTriplets(cofficients.begin(), cofficients.end());

    return A; // placeholder
}

/*
 * Construct a vector encoding the vertices in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |V|, where |V| = # of vertices in the mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildVertexVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> v = Vector<size_t>::Zero(mesh->nVertices());

    for (size_t idx : subset.vertices) v[idx] = 1;

    return v;
}

/*
 * Construct a vector encoding the edges in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |E|, where |E| = # of edges in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildEdgeVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> e = Vector<size_t>::Zero(mesh->nEdges());

    for (size_t idx : subset.edges) e[idx] = 1;

    return e;
}

/*
 * Construct a vector encoding the faces in the selected subset of simplices.
 *
 * Input: Selected subset of simplices.
 * Returns: Vector of length |F|, where |F| = # of faces in mesh.
 */
Vector<size_t> SimplicialComplexOperators::buildFaceVector(const MeshSubset& subset) const {

    // TODO
    Vector<size_t> f = Vector<size_t>::Zero(mesh->nFaces());

    for (size_t idx : subset.edges) f[idx] = 1;

    return f;
}

/*
 * Compute the simplicial star St(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The star of the given subset.
 */
MeshSubset SimplicialComplexOperators::star(const MeshSubset& subset) const {

    // TODO
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

    return star; // placeholder
}


/*
 * Compute the closure Cl(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The closure of the given subset.
 */
MeshSubset SimplicialComplexOperators::closure(const MeshSubset& subset) const {

    // TODO
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

    return closure; // placeholder
}

/*
 * Compute the link Lk(S) of the selected subset of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The link of the given subset.
 */
MeshSubset SimplicialComplexOperators::link(const MeshSubset& subset) const {

    // TODO
    MeshSubset link;
    link = closure(star(subset));
    link.deleteSubset(star(closure(subset)));

    return link; // placeholder
}

/*
 * Return true if the selected subset is a simplicial complex, false otherwise.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: True if given subset is a simplicial complex, false otherwise.
 */
bool SimplicialComplexOperators::isComplex(const MeshSubset& subset) const {

    // TODO
    if (closure(subset).equals(subset)) return true;
    return false; // placeholder
}

/*
 * Check if the given subset S is a pure simplicial complex. If so, return the degree of the complex. Otherwise, return
 * -1.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: int representing the degree of the given complex (-1 if not pure)
 */
int SimplicialComplexOperators::isPureComplex(const MeshSubset& subset) const {

    // TODO
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

    return -1; // placeholder
}

/*
 * Compute the set of simplices contained in the boundary bd(S) of the selected subset S of simplices.
 *
 * Input: A MeshSubset object containing the indices of the currently active vertices, edges, and faces, respectively.
 * Returns: The boundary of the given subset.
 */
MeshSubset SimplicialComplexOperators::boundary(const MeshSubset& subset) const {

    // TODO
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

    return boundary; // placeholder
}