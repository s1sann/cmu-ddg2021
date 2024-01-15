// PLEASE READ:
//
// This file implements additional geometry routines for the VertexPositionGeometry class in Geometry Central. Because
// we are "inside" the class, we no longer have to call
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
// Functions in this file can be called from other projects simply by using geometry->cotan(he),
// geometry->barycentricDualArea(v), etc. where "geometry" is a pointer to a VertexPositionGeometry. This avoids having
// to declare a GeometryRoutines object in every project, and also mimics the way that geometry routines are normally
// called in Geometry Central.
//
// Other notes: In this file, you can use the constant pi by using PI.


#include "geometrycentral/surface/vertex_position_geometry.h"
#include <complex>

typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Triplet<double> Tri;

namespace geometrycentral {
namespace surface {

/*
 * Compute the Euler characteristic of the mesh.
 */
int VertexPositionGeometry::eulerCharacteristic() const {
    return (int)mesh.nVertices() - (int)mesh.nEdges() + (int)mesh.nFaces();
}

/*
 * Compute the mean length of all the edges in the mesh.
 *
 * Input:
 * Returns: The mean edge length.
 */
double VertexPositionGeometry::meanEdgeLength() const {

    double total = 0.0;
    for (Edge e : mesh.edges()) {
        total += edgeLength(e);
    }
    return total / mesh.nEdges();
}

/*
 * Compute the total surface area of the mesh.
 *
 * Input:
 * Returns: The surface area of the mesh.
 */
double VertexPositionGeometry::totalArea() const {

    double total = 0.0;
    for (Face f : mesh.faces()) {
        total += faceArea(f);
    }
    return total;
}

/*
 * Computes the cotangent of the angle opposite to a halfedge. (Do NOT use built-in function for this)
 *
 * Input: The halfedge whose cotan weight is to be computed.
 * Returns: The cotan of the angle opposite the given halfedge.
 */
double VertexPositionGeometry::cotan(Halfedge he) const {

    // TODO
    if (!he.isInterior()) return 0;

    Vector3 a = inputVertexPositions[he.next().tipVertex()], b = inputVertexPositions[he.tailVertex()],
            c = inputVertexPositions[he.tipVertex()];

    Vector3 u = b - a, v = c - a;

    return dot(u, v) / cross(u, v).norm(); // placeholder
}

/*
 * Computes the barycentric dual area of a vertex.
 *
 * Input: The vertex whose barycentric dual area is to be computed.
 * Returns: The barycentric dual area of the given vertex.
 */
double VertexPositionGeometry::barycentricDualArea(Vertex v) const {

    // TODO
    double area = 0.;

    for (Face f : v.adjacentFaces()) area += faceArea(f) / 3;

    return area; // placeholder
}

/*
 * Computes the angle (in radians) at a given corner. (Do NOT use built-in function for this)
 *
 *
 * Input: The corner at which the angle needs to be computed.
 * Returns: The angle clamped between 0 and π.
 */
double VertexPositionGeometry::angle(Corner c) const {
    // TODO
    Vector3 i = inputVertexPositions[c.vertex()];
    Vector3 j = inputVertexPositions[c.halfedge().next().vertex()];
    Vector3 k = inputVertexPositions[c.halfedge().next().next().vertex()];

    Vector3 ij = j - i;
    Vector3 ik = k - i;

    return acos(dot(ij, ik) / norm(ij) / norm(ik)); // placeholder
}

/*
 * Computes the signed angle (in radians) between two adjacent faces. (Do NOT use built-in function for this)
 *
 * Input: The halfedge (shared by the two adjacent faces) on which the dihedral angle is computed.
 * Returns: The dihedral angle.
 */
double VertexPositionGeometry::dihedralAngle(Halfedge he) const {

    // TODO
    Vector3 ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];

    Vector3 n1 = faceNormal(he.face());
    Vector3 n2 = faceNormal(he.twin().face());

    return atan2(dot(ij / norm(ij), cross(n1, n2)), dot(n1, n2)); // placeholder
}

/*
 * Computes the normal at a vertex using the "equally weighted" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "equally weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalEquallyWeighted(Vertex v) const {
    // TODO
    Vector3 n = {0, 0, 0};
    for (Face f : v.adjacentFaces()) n += faceNormal(f);

    return n / norm(n); // placeholder
}

/*
 * Computes the normal at a vertex using the "tip angle weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "tip angle weights" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAngleWeighted(Vertex v) const {
    // TODO
    Vector3 n = {0, 0, 0};
    for (Corner c : v.adjacentCorners()) n += faceNormal(c.face()) * angle(c);

    return n / norm(n); // placeholder
}

/*
 * Computes the normal at a vertex using the "inscribed sphere" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "inscribed sphere" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalSphereInscribed(Vertex v) const {
    // TODO
    Vector3 n = {0, 0, 0};
    for (Corner c : v.adjacentCorners()) {
        Vector3 i = inputVertexPositions[c.vertex()];
        Vector3 j = inputVertexPositions[c.halfedge().next().vertex()];
        Vector3 k = inputVertexPositions[c.halfedge().next().next().vertex()];

        Vector3 ij = j - i;
        Vector3 ik = k - i;

        n += cross(ij, ik) / pow(norm(ij) * norm(ik),2);
    }

    return n / norm(n); // placeholder
}

/*
 * Computes the normal at a vertex using the "face area weights" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "face area weighted" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalAreaWeighted(Vertex v) const {
    // TODO
    Vector3 n = {0, 0, 0};
    for (Face f : v.adjacentFaces()) n += faceArea(f) * faceNormal(f);

    return n / norm(n); // placeholder
}

/*
 * Computes the normal at a vertex using the "Gauss curvature" method.
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "Gauss curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalGaussianCurvature(Vertex v) const {
    // TODO
    Vector3 n = {0, 0, 0};
    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];
        n += dihedralAngle(he) * ij / norm(ij);
    }

    return n / norm(n); // placeholder
}

/*
 * Computes the normal at a vertex using the "mean curvature" method (equivalent to the "area gradient" method).
 *
 * Input: The vertex on which the normal is to be computed.
 * Returns: The "mean curvature" normal vector.
 */
Vector3 VertexPositionGeometry::vertexNormalMeanCurvature(Vertex v) const {
    // TODO
    Vector3 n = {0, 0, 0};
    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];
        n += (cotan(he) + cotan(he.twin())) * ij;
    }

    return n / norm(n); // placeholder
}

/*
 * Computes the angle defect at a vertex.
 *
 * Input: The vertex whose angle defect is to be computed.
 * Returns: The angle defect of the given vertex.
 */
double VertexPositionGeometry::angleDefect(Vertex v) const {
    // TODO
    double k = 2 * PI;
    for (Corner c : v.adjacentCorners()) k -= angle(c);
    return k; // placeholder
}

/*
 * Computes the total angle defect of the mesh.
 *
 * Input:
 * Returns: The total angle defect
 */
double VertexPositionGeometry::totalAngleDefect() const {
    // TODO
    double k = 0.;
    for (Vertex v : mesh.vertices()) k += angleDefect(v);
    return k; // placeholder
}

/*
 * Computes the (integrated) scalar mean curvature at a vertex.
 *
 * Input: The vertex whose mean curvature is to be computed.
 * Returns: The mean curvature at the given vertex.
 */
double VertexPositionGeometry::scalarMeanCurvature(Vertex v) const {
    // TODO
    double h = 0.;
    for (Halfedge he : v.outgoingHalfedges()) {
        Vector3 ij = inputVertexPositions[he.tipVertex()] - inputVertexPositions[he.tailVertex()];
        h += dihedralAngle(he) * norm(ij);
    }

    return h * 0.5; // placeholder
}

/*
 * Computes the circumcentric dual area of a vertex.
 *
 * Input: The vertex whose circumcentric dual area is to be computed.
 * Returns: The circumcentric dual area of the given vertex.
 */
double VertexPositionGeometry::circumcentricDualArea(Vertex v) const {
    // TODO
    double a = 0.;
    for (Corner c : v.adjacentCorners()) {
        Vector3 i = inputVertexPositions[c.vertex()];
        Vector3 j = inputVertexPositions[c.halfedge().next().vertex()];
        Vector3 k = inputVertexPositions[c.halfedge().next().next().vertex()];

        Vector3 ij = j - i;
        Vector3 ik = k - i;

        a += norm2(ij) * cotan(c.halfedge()) + norm2(ik) * cotan(c.halfedge().next().next());
    }

    return a / 8.; // placeholder
}

/*
 * Computes the (pointwise) minimum and maximum principal curvature values at a vertex.
 *
 * Input: The vertex on which the principal curvatures need to be computed.
 * Returns: A std::pair containing the minimum and maximum principal curvature values at a vertex.
 */
std::pair<double, double> VertexPositionGeometry::principalCurvatures(Vertex v) const {
    // TODO
    double a = circumcentricDualArea(v);
    double g = angleDefect(v) / a;
    double m = 2 * scalarMeanCurvature(v) / a;

    double k2 = (m + sqrt(m * m - 4 * g)) / 2;
    double k1 = m - k2;
    return std::make_pair(std::min(k1, k2), std::max(k1, k2)); // placeholder
}


/*
 * Builds the sparse POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace matrix,
 * multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse positive definite Laplace matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::laplaceMatrix() const {
    // TODO
    std::vector<Tri> coefficients;

    for (Vertex v : mesh.vertices()) {
        double sum = 0.;
        for (Halfedge he : v.incomingHalfedges()) {
            double num = (cotan(he) + cotan(he.twin())) / 2;
            coefficients.push_back(Tri(v.getIndex(), he.vertex().getIndex(), -num));
            sum += num;
        }
        coefficients.push_back(Tri(v.getIndex(), v.getIndex(), sum+1e-8));
    }

    SpMat m(mesh.nVertices(), mesh.nVertices());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m; // placeholder
}

/*
 * Builds the sparse diagonal mass matrix containing the barycentric dual area of each vertex.
 *
 * Input:
 * Returns: Sparse mass matrix for the mesh.
 */
SparseMatrix<double> VertexPositionGeometry::massMatrix() const {
    // TODO
    std::vector<Tri> coefficients;

    for (Vertex v : mesh.vertices()) coefficients.push_back(Tri(v.getIndex(), v.getIndex(), barycentricDualArea(v)));

    SpMat m(mesh.nVertices(), mesh.nVertices());
    m.setFromTriplets(coefficients.begin(), coefficients.end());

    return m; // placeholder
}

/*
 * Builds the sparse complex POSITIVE DEFINITE Laplace matrix. Do this by building the negative semidefinite Laplace
 * matrix, multiplying by -1, and shifting the diagonal elements by a small constant (e.g. 1e-8).
 *
 * Input:
 * Returns: Sparse complex positive definite Laplace matrix for the mesh.
 */
SparseMatrix<std::complex<double>> VertexPositionGeometry::complexLaplaceMatrix() const {

    // TODO
    return identityMatrix<std::complex<double>>(1); // placeholder
}

/*
 * Compute the center of mass of a mesh.
 */
Vector3 VertexPositionGeometry::centerOfMass() const {

    // Compute center of mass.
    Vector3 center = {0.0, 0.0, 0.0};
    for (Vertex v : mesh.vertices()) {
        center += inputVertexPositions[v];
    }
    center /= mesh.nVertices();

    return center;
}

/*
 * Centers a mesh about the origin.
 * Also rescales the mesh to unit radius if <rescale> == true.
 */
void VertexPositionGeometry::normalize(const Vector3& origin, bool rescale) {

    // Compute center of mass.
    Vector3 center = centerOfMass();

    // Translate to origin [of original mesh].
    double radius = 0;
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] -= center;
        radius = std::max(radius, inputVertexPositions[v].norm());
    }

    // Rescale.
    if (rescale) {
        for (Vertex v : mesh.vertices()) {
            inputVertexPositions[v] /= radius;
        }
    }

    // Translate to origin [of original mesh].
    for (Vertex v : mesh.vertices()) {
        inputVertexPositions[v] += origin;
    }
}

} // namespace surface
} // namespace geometrycentral