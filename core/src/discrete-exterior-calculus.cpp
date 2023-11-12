// PLEASE READ:
//
// This file additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because we are
// "inside" the class, we no longer have to call
//
//          geometry->inputVertexPositions[v], etc.
//
// We can just call
//
//          this->inputVertexPositions[v], etc.
//
// or simply
//
//          inputVertexPositions[v], etc.
//
// In addition, we no longer access the corresponding surface mesh via
//
//          mesh->vertices(), etc.
//
// but instead <mesh> is not a pointer anymore, so we use
//
//          mesh.vertices(), etc.
//
// Functions in this file can be called from other projects simply by using geometry->buildHodgeStar0Form(), etc. where
// "geometry" is a pointer to a VertexPositionGeometry. This avoids having to declare a GeometryRoutines object in every
// project, and also mimics the way that geometry routines are normally called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.

#include "geometrycentral/surface/vertex_position_geometry.h"

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Tri;

namespace geometrycentral {
namespace surface {


/*
 * Build Hodge operator on 0-forms.
 * By convention, the area of a vertex is 1.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar0Form() const {
    // TODO
    std::vector<Tri> coefficients;

    for (Vertex v : mesh.vertices()) coefficients.push_back(Tri(v.getIndex(), v.getIndex(), barycentricDualArea(v)));

    SpMat m(mesh.nVertices(), mesh.nVertices());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m; // placeholder
}

/*
 * Build Hodge operator on 1-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar1Form() const {
    std::vector<Tri> coefficients;

    for (Edge e : mesh.edges()) {
        double area = (cotan(e.halfedge()) + cotan(e.halfedge().twin())) / 2;
        coefficients.push_back(Tri(e.getIndex(), e.getIndex(), area));
    }

    SpMat m(mesh.nEdges(), mesh.nEdges());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m; // placeholder
}

/*
 * Build Hodge operator on 2-forms.
 *
 * Input:
 * Returns: A sparse diagonal matrix representing the Hodge operator that can be applied to discrete 2-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildHodgeStar2Form() const {
    // TODO
    std::vector<Tri> coefficients;

    for (Face f : mesh.faces()) coefficients.push_back(Tri(f.getIndex(), f.getIndex(), 1 / faceArea(f)));

    SpMat m(mesh.nFaces(), mesh.nFaces());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m; // placeholder
}

/*
 * Build exterior derivative on 0-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 0-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative0Form() const {
    // TODO
    std::vector<Tri> coefficients;

    for (Edge e : mesh.edges()) {
        coefficients.push_back(Tri(e.getIndex(), e.firstVertex().getIndex(), -1));
        coefficients.push_back(Tri(e.getIndex(), e.secondVertex().getIndex(), 1));
    }

    SpMat m(mesh.nEdges(), mesh.nVertices());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m; // placeholder
}

/*
 * Build exterior derivative on 1-forms.
 *
 * Input:
 * Returns: A sparse matrix representing the exterior derivative that can be applied to discrete 1-forms.
 */
SparseMatrix<double> VertexPositionGeometry::buildExteriorDerivative1Form() const {
    // TODO
    std::vector<Tri> coefficients;

    for (Face f : mesh.faces()) {
        for (Halfedge he : f.adjacentHalfedges()) {
            double orientation = he.getIndex() < he.twin().getIndex() ? 1 : - 1;
            coefficients.push_back(Tri(f.getIndex(), he.edge().getIndex(), orientation));
        }
    }

    SpMat m(mesh.nFaces(), mesh.nEdges());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m; // placeholder
}

} // namespace surface
} // namespace geometrycentral