// Implement member functions for MeanCurvatureFlow class.
#include "mean-curvature-flow.h"
#include "geometrycentral/numerical/linear_solvers.h"

/* Constructor
 * Input: The surface mesh <inputMesh> and geometry <inputGeo>.
 */
MeanCurvatureFlow::MeanCurvatureFlow(ManifoldSurfaceMesh* inputMesh, VertexPositionGeometry* inputGeo) {

    // Build member variables: mesh, geometry
    mesh = inputMesh;
    geometry = inputGeo;
}

/*
 * Build the mean curvature flow operator.
 *
 * Input: The mass matrix <M> of the mesh, and the timestep <h>.
 * Returns: A sparse matrix representing the mean curvature flow operator.
 */
SparseMatrix<double> MeanCurvatureFlow::buildFlowOperator(const SparseMatrix<double>& M, double h) const {
    // TODO
    SparseMatrix<double> m = geometry->laplaceMatrix();
    return M+(h*m); // placeholder
}

/*
 * Performs mean curvature flow.
 *
 * Input: The timestep <h>.
 * Returns:
 */
void MeanCurvatureFlow::integrate(double h) {

    // TODO
    // Note: Geometry Central has linear solvers: https://geometry-central.net/numerical/linear_solvers/
    // Note: Update positions via geometry->inputVertexPositions
    SparseMatrix<double> M = geometry->massMatrix();
    SparseMatrix<double> A = buildFlowOperator(M,h);

    std::vector<Vector<double>> positions(3, Vector<double>::Zero(mesh->nVertices()));
    for (Vertex v : mesh->vertices()) {
        positions[0][v.getIndex()] = geometry->inputVertexPositions[v].x; // placeholder
        positions[1][v.getIndex()] = geometry->inputVertexPositions[v].y; // placeholder
        positions[2][v.getIndex()] = geometry->inputVertexPositions[v].z; // placeholder
    }

    for(int i=0;i<3;i++) {
        positions[i] = M * positions[i];
        positions[i] = solvePositiveDefinite(A, positions[i]);
    }

    for (Vertex v : mesh->vertices()) {
        geometry->inputVertexPositions[v] = { 
            positions[0][v.getIndex()],
            positions[1][v.getIndex()],
            positions[2][v.getIndex()]
        };
    }
}