/*
	ga_fem.cpp: the method of FEM.cpp written with geometric algebra (GAConnection.h), as in Part II of
	fembasic2d.pdf, side by side with the matrix version of Part I, without any graphics.

	  1. prints the worked examples of the text;
	  2. solves the cactus of FEM.cpp twice, with 2x2 real blocks (Part I, the code of FEM.cpp) and with the
	     n x n complex Hermitian system (Part II), and compares the two;
	  3. solves the same problem on a sphere, where parallel transport along great circles is known exactly
	     (it is the smallest rotor from the source to the point), and measures the error;
	  4. writes CSV files for the figures of the text into the directory given as argument, and checks Part VI (end of main).

	Build, from src/:
	  c++ -std=c++17 -O2 -I<eigen> ga_fem.cpp HalfEdge/Mesh.cpp -o ga_fem && ./ga_fem ../fem-basic-2d/data
*/
#include <cstdio>
#include <fstream>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "GAConnection.h"
#include "HalfEdge/Mesh.h"

using Eigen::Vector3d;
using Eigen::Vector2d;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
typedef Eigen::Matrix<double, 3, 2> Matrix32d;
using ga::Complex;

// ---------------------------------------------------------------------------------------------------------------
// The geometry shared by both versions (Section 2): vertex and face normals and frames, as in FEM.cpp
// ---------------------------------------------------------------------------------------------------------------
struct Geometry {
	std::vector<Vector3d> faceNormal, vertexT, faceT;   // normals, first frame vectors t_i and t_f
	std::vector<double> faceArea;
	std::vector<Matrix32d> vertexFrame, faceFrame;      // T_i = [t_i  n_i x t_i], T_f = [t_f  n_f x t_f]
};

static void faceVertices(const Face& f, Vertex* v[3]) {
	v[0] = f.edge->vertex;
	v[1] = f.edge->next->vertex;
	v[2] = f.edge->next->next->vertex;
}

static Geometry computeGeometry(Mesh& mesh) {
	Geometry g;
	g.vertexT.resize(mesh.numVertices());
	g.vertexFrame.resize(mesh.numVertices());
	for (Vertex& vi : mesh.getVertices()) {   // computeVertexTangentSpace
		const Vector3d e = vi.edge->pair->vertex->p - vi.p;
		const Vector3d t = (e - e.dot(vi.n) * vi.n).normalized();
		g.vertexT[vi.ID] = t;
		g.vertexFrame[vi.ID].col(0) = t;
		g.vertexFrame[vi.ID].col(1) = vi.n.cross(t);
	}
	const int nf = mesh.numFaces();
	g.faceNormal.resize(nf); g.faceArea.resize(nf); g.faceT.resize(nf); g.faceFrame.resize(nf);
	for (Face& f : mesh.getFaces()) {         // computeFaceNormals, computeFaceTangentSpace
		Vertex* v[3];
		faceVertices(f, v);
		const Vector3d n = (v[1]->p - v[0]->p).cross(v[2]->p - v[0]->p);
		g.faceArea[f.ID] = 0.5 * n.norm();
		g.faceNormal[f.ID] = n.normalized();
		const Vector3d t = (v[1]->p - v[0]->p).normalized();
		g.faceT[f.ID] = t;
		g.faceFrame[f.ID].col(0) = t;
		g.faceFrame[f.ID].col(1) = g.faceNormal[f.ID].cross(t);
	}
	return g;
}

// ---------------------------------------------------------------------------------------------------------------
// Part I: the matrix version, the formulas of FEM.cpp
// ---------------------------------------------------------------------------------------------------------------
namespace matrixform {

/// alignVectors: the Rodrigues matrix of the smallest rotation taking nv to nf (eq. rodrigues)
static Matrix3d alignVectors(const Vector3d& nv, const Vector3d& nf) {
	const double c = nv.dot(nf);
	const Vector3d v = nv.cross(nf);
	Matrix3d skew;
	skew << 0, -v(2), v(1), v(2), 0, -v(0), -v(1), v(0), 0;
	return Matrix3d::Identity() + skew + skew * skew / (1.0 + c);
}

/// Rvf: R_i^f = T_f^T Q T_i (eq. connection)
static Matrix2d Rvf(const Geometry& g, const Vertex* v, const Face& f) {
	return g.faceFrame[f.ID].transpose() * alignVectors(v->n, g.faceNormal[f.ID]) * g.vertexFrame[v->ID];
}

/// gradN: B (B^T B)^-1 applied to the reference gradients (Proposition gradient)
static void gradients(Vertex* v[3], Vector3d gradN[3]) {
	Matrix32d B;
	B.col(0) = v[1]->p - v[0]->p;
	B.col(1) = v[2]->p - v[0]->p;
	const Matrix32d Binv = ((B.transpose() * B).inverse() * B.transpose()).transpose();
	gradN[0] = -Binv.col(0) - Binv.col(1);
	gradN[1] = Binv.col(0);
	gradN[2] = Binv.col(1);
}

/// AssembleMatrix2 + AssembleDiagonalMassMatrix + IplusMinvTimesA + PreFactor + solve, in Eigen
static Eigen::VectorXd solve(Mesh& mesh, const Geometry& g, double dt, int source, Vector2d value) {
	const int n = mesh.numVertices();
	std::vector<double> mass(n, 0.0);
	std::map<std::pair<int, int>, double> A;   // dt L, in the 2n x 2n real form
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		faceVertices(f, v);
		Vector3d gradN[3];
		gradients(v, gradN);
		const double area = g.faceArea[f.ID];
		for (int a = 0; a < 3; ++a) {
			mass[v[a]->ID] += area / 3.0;
			for (int b = 0; b < 3; ++b) {
				const Vector2d ga = g.faceFrame[f.ID].transpose() * gradN[a];
				const Vector2d gb = g.faceFrame[f.ID].transpose() * gradN[b];
				const Matrix2d block = a == b ? Matrix2d(Matrix2d::Identity())
				                              : Matrix2d(Rvf(g, v[a], f).transpose() * Rvf(g, v[b], f));
				const Matrix2d w = dt * area * ga.dot(gb) * block;
				for (int r = 0; r < 2; ++r)
					for (int c = 0; c < 2; ++c) A[{ 2 * v[a]->ID + r, 2 * v[b]->ID + c }] += w(r, c);
			}
		}
	}
	std::vector<Eigen::Triplet<double>> triplets;
	for (const auto& [ij, value_ij] : A) {
		const int i = ij.first, j = ij.second;
		if (i / 2 == source) continue;   // Dirichlet rows are replaced
		triplets.emplace_back(i, j, value_ij / mass[i / 2] + (i == j ? 1.0 : 0.0));   // I + M^-1 (dt L)
	}
	triplets.emplace_back(2 * source, 2 * source, 1.0);
	triplets.emplace_back(2 * source + 1, 2 * source + 1, 1.0);
	Eigen::SparseMatrix<double> S(2 * n, 2 * n);
	S.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::VectorXd rhs = Eigen::VectorXd::Zero(2 * n);
	rhs.segment<2>(2 * source) = value;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> lu(S);
	return lu.solve(rhs);
}

}   // namespace matrixform

// ---------------------------------------------------------------------------------------------------------------
// Part II: the geometric algebra version, an n x n complex Hermitian system
// ---------------------------------------------------------------------------------------------------------------
namespace gaform {

/// Section 16: L_ij += A_f (grad phi_i . grad phi_j) e^{T (rho_j - rho_i)}, M_ii += A_f / 3, then the same
/// scaling I + M^-1 (dt L) and the same Dirichlet row as the matrix version, now with complex numbers.
static Eigen::VectorXcd solve(Mesh& mesh, const Geometry& g, double dt, int source, Complex value) {
	const int n = mesh.numVertices();
	std::vector<double> mass(n, 0.0);
	std::map<std::pair<int, int>, Complex> A;
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		faceVertices(f, v);
		ga::Vec3 grad[3];
		ga::hatGradients(v[0]->p, v[1]->p, v[2]->p, grad);   // the reciprocal frame
		const double area = g.faceArea[f.ID];
		Complex rho[3];                                     // e^{T_f rho_a^f}
		for (int a = 0; a < 3; ++a)
			rho[a] = ga::connection(g.vertexT[v[a]->ID], v[a]->n, g.faceT[f.ID], g.faceNormal[f.ID]);
		for (int a = 0; a < 3; ++a) {
			mass[v[a]->ID] += area / 3.0;
			for (int b = 0; b < 3; ++b) {
				// transport from b to a through the face: e^{T (rho_b - rho_a)} = e^{T rho_b} conj(e^{T rho_a})
				const Complex transport = a == b ? Complex(1.0) : rho[b] * std::conj(rho[a]);
				A[{ v[a]->ID, v[b]->ID }] += dt * area * grad[a].dot(grad[b]) * transport;
			}
		}
	}
	std::vector<Eigen::Triplet<Complex>> triplets;
	for (const auto& [ij, value_ij] : A) {
		if (ij.first == source) continue;
		triplets.emplace_back(ij.first, ij.second, value_ij / mass[ij.first] + (ij.first == ij.second ? 1.0 : 0.0));
	}
	triplets.emplace_back(source, source, 1.0);
	Eigen::SparseMatrix<Complex> S(n, n);
	S.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(n);
	rhs(source) = value;
	Eigen::SparseLU<Eigen::SparseMatrix<Complex>> lu(S);
	return lu.solve(rhs);
}

}   // namespace gaform

// ---------------------------------------------------------------------------------------------------------------
// Part III: the intrinsic version, the same n x n complex system from the edge lengths alone
// ---------------------------------------------------------------------------------------------------------------
#include "IntrinsicConnection.h"

namespace intrinsicform {

/// Section 23: the assembly of gaform::solve, with the gradients, areas and connection of IntrinsicConnection.h.
static Eigen::VectorXcd solve(Mesh& mesh, const intrinsic::Geometry& G, double dt, int source, Complex value) {
	const int n = mesh.numVertices();
	std::vector<double> mass(n, 0.0);
	std::map<std::pair<int, int>, Complex> A;
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		const intrinsic::FaceGeometry& g = G.faces[f.ID];
		for (int a = 0; a < 3; ++a) {
			mass[v[a]->ID] += g.area / 3.0;
			for (int b = 0; b < 3; ++b) {
				const Complex transport = a == b ? Complex(1.0) : g.rho[b] * std::conj(g.rho[a]);
				A[{ v[a]->ID, v[b]->ID }] += dt * g.area * g.grad[a].dot(g.grad[b]) * transport;
			}
		}
	}
	std::vector<Eigen::Triplet<Complex>> triplets;
	for (const auto& [ij, value_ij] : A) {
		if (ij.first == source) continue;
		triplets.emplace_back(ij.first, ij.second, value_ij / mass[ij.first] + (ij.first == ij.second ? 1.0 : 0.0));
	}
	triplets.emplace_back(source, source, 1.0);
	Eigen::SparseMatrix<Complex> S(n, n);
	S.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(n);
	rhs(source) = value;
	Eigen::SparseLU<Eigen::SparseMatrix<Complex>> lu(S);
	return lu.solve(rhs);
}

}   // namespace intrinsicform

// ---------------------------------------------------------------------------------------------------------------
// Part IV: the edge connection, with its curvature inside the faces; Part V: discrete exterior calculus
// ---------------------------------------------------------------------------------------------------------------
#include "EdgeConnection.h"

/// The scaled system of Sections 16 and 23, I + M^-1 (dt L) with the Dirichlet row of the source, from the entries of
/// dt L and the lumped masses, solved by sparse LU.
static Eigen::VectorXcd solveScaled(int n, const std::map<std::pair<int, int>, Complex>& A, const std::vector<double>& mass,
                                    int source, Complex value) {
	std::vector<Eigen::Triplet<Complex>> triplets;
	for (const auto& [ij, value_ij] : A) {
		if (ij.first == source) continue;
		triplets.emplace_back(ij.first, ij.second, value_ij / mass[ij.first] + (ij.first == ij.second ? 1.0 : 0.0));
	}
	triplets.emplace_back(source, source, 1.0);
	Eigen::SparseMatrix<Complex> S(n, n);
	S.setFromTriplets(triplets.begin(), triplets.end());
	Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(n);
	rhs(source) = value;
	Eigen::SparseLU<Eigen::SparseMatrix<Complex>> lu(S);
	return lu.solve(rhs);
}

/// The half edge from v[b] to v[a] of a face, for a != b.
static const Edge* halfEdge(Edge* h[3], int b, int a) { return a == (b + 1) % 3 ? h[b] : h[a]->pair; }

namespace curvedform {

/// Section 27: assemble the curved elements (eq. curved_element) with the lumped mass A_f / 3.
static Eigen::VectorXcd solve(Mesh& mesh, const intrinsic::Geometry& G, const edgeconn::Geometry& E, double dt, int source,
                              Complex value, int quadrature = 8) {
	const int n = mesh.numVertices();
	std::vector<double> mass(n, 0.0);
	std::map<std::pair<int, int>, Complex> A;
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		Complex transportBA[3][3], L[3][3];   // e^{i rho_ba}: from the plane of v[b] to that of v[a]
		for (int a = 0; a < 3; ++a)
			for (int b = 0; b < 3; ++b) transportBA[a][b] = a == b ? Complex(1.0) : std::polar(1.0, E.transport.at(halfEdge(h, b, a)));
		edgeconn::curvedElement(G.faces[f.ID], E.holonomy[f.ID], transportBA, quadrature, L);
		for (int a = 0; a < 3; ++a) {
			mass[v[a]->ID] += G.faces[f.ID].area / 3.0;
			for (int b = 0; b < 3; ++b) A[{ v[a]->ID, v[b]->ID }] += dt * L[a][b];
		}
	}
	return solveScaled(n, A, mass, source, value);
}

}   // namespace curvedform

namespace decform {

/// Section 33: the covariant exterior derivative d^nabla (one row per edge a -> b: psi_b - r_ab psi_a), the diagonal Hodge
/// stars *1 (half the sum of the cotangents opposite to each edge) and *0 (a third of the area of the faces around each
/// vertex), and the connection Laplacian L = (d^nabla)^H *1 d^nabla.
static Eigen::SparseMatrix<Complex> laplacian(Mesh& mesh, const intrinsic::Geometry& G, const edgeconn::Geometry& E,
                                              std::vector<double>& star0) {
	std::vector<const Edge*> edges;
	for (const Edge* e : mesh.getEdges()) if (e->vertex->ID < e->pair->vertex->ID) edges.push_back(e);
	auto cotOpposite = [&](const Edge* e) {   // the cotangent of the angle opposite to e in the face to its left
		if (!e->face) return 0.0;
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(*e->face, v, h);
		int k = 0;
		while (h[k] != e) ++k;
		return 1.0 / std::tan(G.faces[e->face->ID].angle[(k + 2) % 3]);
	};
	std::vector<Eigen::Triplet<Complex>> d, star1;
	for (int k = 0; k < (int)edges.size(); ++k) {
		const Edge* e = edges[k];
		d.emplace_back(k, e->pair->vertex->ID, 1.0);                                       // + psi_b
		d.emplace_back(k, e->vertex->ID, -std::polar(1.0, E.transport.at(e)));            // - r_ab psi_a
		star1.emplace_back(k, k, 0.5 * (cotOpposite(e) + cotOpposite(e->pair)));
	}
	const int n = mesh.numVertices(), m = (int)edges.size();
	Eigen::SparseMatrix<Complex> D(m, n), S1(m, m);
	D.setFromTriplets(d.begin(), d.end());
	S1.setFromTriplets(star1.begin(), star1.end());
	star0.assign(n, 0.0);
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		for (int a = 0; a < 3; ++a) star0[v[a]->ID] += G.faces[f.ID].area / 3.0;
	}
	return Eigen::SparseMatrix<Complex>(D.adjoint() * S1 * D);
}

static Eigen::VectorXcd solve(Mesh& mesh, const intrinsic::Geometry& G, const edgeconn::Geometry& E, double dt, int source,
                              Complex value) {
	std::vector<double> star0;
	const Eigen::SparseMatrix<Complex> L = laplacian(mesh, G, E, star0);
	std::map<std::pair<int, int>, Complex> A;
	for (int k = 0; k < L.outerSize(); ++k)
		for (Eigen::SparseMatrix<Complex>::InnerIterator it(L, k); it; ++it) A[{ (int)it.row(), (int)it.col() }] += dt * it.value();
	return solveScaled(mesh.numVertices(), A, star0, source, value);
}

}   // namespace decform

// ---------------------------------------------------------------------------------------------------------------
// Meshes
// ---------------------------------------------------------------------------------------------------------------
static double meanEdgeLength(Mesh& mesh) {
	double sum = 0;
	int count = 0;
	for (auto& e : mesh.getEdges()) { sum += (e->pair->vertex->p - e->vertex->p).norm(); ++count; }
	return sum / count;
}

/// A geodesic sphere: an icosahedron subdivided `levels` times, projected on the unit sphere, written as OBJ.
static void writeIcosphere(const std::string& path, int levels) {
	const double p = (1.0 + std::sqrt(5.0)) / 2.0;
	std::vector<Vector3d> v = { {-1, p, 0}, {1, p, 0}, {-1, -p, 0}, {1, -p, 0}, {0, -1, p}, {0, 1, p},
	                            {0, -1, -p}, {0, 1, -p}, {p, 0, -1}, {p, 0, 1}, {-p, 0, -1}, {-p, 0, 1} };
	for (auto& x : v) x.normalize();
	std::vector<std::array<int, 3>> f = { {0,11,5}, {0,5,1}, {0,1,7}, {0,7,10}, {0,10,11}, {1,5,9}, {5,11,4},
	    {11,10,2}, {10,7,6}, {7,1,8}, {3,9,4}, {3,4,2}, {3,2,6}, {3,6,8}, {3,8,9}, {4,9,5}, {2,4,11}, {6,2,10},
	    {8,6,7}, {9,8,1} };
	for (int l = 0; l < levels; ++l) {
		std::map<std::pair<int, int>, int> mid;
		auto midpoint = [&](int a, int b) {
			auto key = std::minmax(a, b);
			auto it = mid.find(key);
			if (it != mid.end()) return it->second;
			v.push_back((v[a] + v[b]).normalized());
			return mid[key] = (int)v.size() - 1;
		};
		std::vector<std::array<int, 3>> g;
		for (auto& t : f) {
			const int a = midpoint(t[0], t[1]), b = midpoint(t[1], t[2]), c = midpoint(t[2], t[0]);
			g.push_back({ t[0], a, c }); g.push_back({ t[1], b, a }); g.push_back({ t[2], c, b }); g.push_back({ a, b, c });
		}
		f = g;
	}
	std::ofstream out(path);
	for (auto& x : v) out << "v " << x.x() << " " << x.y() << " " << x.z() << "\n";
	for (auto& t : f) out << "f " << t[0] + 1 << " " << t[1] + 1 << " " << t[2] + 1 << "\n";
}

// ---------------------------------------------------------------------------------------------------------------
// The worked examples of the text
// ---------------------------------------------------------------------------------------------------------------
static void workedExamples() {
	printf("== Worked example 1: one triangle in space\n");
	const Vector3d x1(0, 0, 0), x2(2, 0, 0), x3(0, 1, 1);
	Matrix32d B;
	B.col(0) = x2 - x1;
	B.col(1) = x3 - x1;
	const Matrix32d Binv = B * (B.transpose() * B).inverse();
	const double area = 0.5 * (x2 - x1).cross(x3 - x1).norm();
	ga::Vec3 grad[3];
	ga::hatGradients(x1, x2, x3, grad);
	printf("  B^T B = [%g %g; %g %g], area = %.6f (sqrt 2 = %.6f)\n", (B.transpose() * B)(0, 0), (B.transpose() * B)(0, 1),
	       (B.transpose() * B)(1, 0), (B.transpose() * B)(1, 1), area, std::sqrt(2.0));
	for (int a = 0; a < 3; ++a) printf("  grad phi_%d = (%g, %g, %g)\n", a + 1, grad[a].x(), grad[a].y(), grad[a].z());
	printf("  columns of Binv: (%g %g %g), (%g %g %g)\n", Binv(0, 0), Binv(1, 0), Binv(2, 0), Binv(0, 1), Binv(1, 1), Binv(2, 1));
	const Vector3d n = (x2 - x1).cross(x3 - x1).normalized();
	const ga::Vec3 g2 = ga::hatGradient(x3, x1, ga::tangentBivector(n), area);
	printf("  grad phi_2 from the opposite edge turned a quarter turn: (%g, %g, %g)\n", g2.x(), g2.y(), g2.z());
	printf("  stiffness A grad phi_a . grad phi_b:\n");
	for (int a = 0; a < 3; ++a)
		printf("    %9.6f %9.6f %9.6f\n", area * grad[a].dot(grad[0]), area * grad[a].dot(grad[1]), area * grad[a].dot(grad[2]));
	const Vector3d u = x1 - x3, w = x2 - x3;
	printf("  -cot(theta_3)/2 = %.6f\n", -0.5 * u.dot(w) / u.cross(w).norm());

	printf("== Worked example 2: transport from a vertex plane into the face plane\n");
	const Vector3d ni(0, 0, 1), ti(0, 1, 0), tf = (x2 - x1).normalized();
	const ga::Rotor Q = ga::rotorBetween(ni, n);
	const Vector3d c = Q.rotate(ti);
	const Complex rho = ga::connection(ti, ni, tf, n);
	printf("  n_f = (%g, %g, %g), Q = %.6f + I (%.6f, %.6f, %.6f)\n", n.x(), n.y(), n.z(), Q.w, Q.b.x(), Q.b.y(), Q.b.z());
	printf("  Q t_i ~Q = (%.6f, %.6f, %.6f), e^{T rho} = %.6f + %.6f T, rho = %.1f degrees\n", c.x(), c.y(), c.z(),
	       rho.real(), rho.imag(), std::arg(rho) * 180.0 / M_PI);
	Matrix32d Ti, Tf;
	Ti.col(0) = ti; Ti.col(1) = ni.cross(ti);
	Tf.col(0) = tf; Tf.col(1) = n.cross(tf);
	const Matrix2d R = Tf.transpose() * matrixform::alignVectors(ni, n) * Ti;
	printf("  R_i^f = T_f^T Q T_i = [%.6f %.6f; %.6f %.6f]\n", R(0, 0), R(0, 1), R(1, 0), R(1, 1));
}


// ---------------------------------------------------------------------------------------------------------------
// Part III: comparing the extrinsic and the intrinsic connection
// ---------------------------------------------------------------------------------------------------------------
static double angleDeg(const Vector3d& a, const Vector3d& b) { return std::atan2(a.cross(b).norm(), a.dot(b)) * 180.0 / M_PI; }

/// The unit vector of face f at the given angle from its frame vector t_f, counterclockwise about n_f.
static Vector3d faceVector(const Geometry& g, int f, double angle) {
	return std::cos(angle) * g.faceT[f] + std::sin(angle) * g.faceNormal[f].cross(g.faceT[f]);
}

/// What each method puts in face f for the spinor coordinate z of its corner a: the vector of angle arg z + rho_a^f
/// (eq. ga_transport), which is what the finite element method uses in that face. Comparing these vectors face by face
/// needs no choice of a plane at the vertex (Section 24).
struct Methods {
	Mesh& mesh;
	const Geometry& g;
	const intrinsic::Geometry& G;
	Vector3d extrinsic(int f, int a, Complex z) const {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(mesh.faceAt(f), v, h);
		const int i = v[a]->ID;
		const Complex rho = ga::connection(g.vertexT[i], mesh.vertexAt(i).n, g.faceT[f], g.faceNormal[f]);
		return faceVector(g, f, std::arg(z) + std::arg(rho));
	}
	Vector3d intrinsic(int f, int a, Complex z) const { return faceVector(g, f, std::arg(z) + std::arg(G.faces[f].rho[a])); }
	/// Parts IV and V transport along edges: the vector keeps its angle to the edge h[a] leaving the vertex.
	Vector3d edge(int f, int a, Complex z) const {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(mesh.faceAt(f), v, h);
		return faceVector(g, f, G.faces[f].edgeAngle[a] + std::arg(z) - edgeconn::vertexAngle(G, h[a]));
	}
	/// Method k = 0, 1, 2, 3: Part II, III, IV, V.
	Vector3d of(int k, int f, int a, Complex z) const { return k == 0 ? extrinsic(f, a, z) : k == 1 ? intrinsic(f, a, z) : edge(f, a, z); }
};

static const char* methodName[4] = { "extrinsic (II)", "intrinsic (III)", "curved (IV)", "DEC (V)" };

/// The vector heat method with the four connections.
static std::vector<Eigen::VectorXcd> solveAll(Mesh& mesh, const Geometry& g, const intrinsic::Geometry& G,
                                              const edgeconn::Geometry& E, double dt, int source, Complex value) {
	return { gaform::solve(mesh, g, dt, source, value), intrinsicform::solve(mesh, G, dt, source, value),
	         curvedform::solve(mesh, G, E, dt, source, value), decform::solve(mesh, G, E, dt, source, value) };
}

/// The face of the reference edge of vertex i, and the corner of i in it.
static void firstCorner(Mesh& mesh, int i, int& f, int& a) {
	const Edge* e = mesh.vertexAt(i).edge;
	f = e->face->ID;
	Vertex* v[3];
	Edge* h[3];
	intrinsic::faceCorners(*e->face, v, h);
	for (a = 0; h[a] != e; ++a) {}
}

/// The discrete Gauss-Bonnet sum of Section 21: sum of 2 pi - Theta_i over interior vertices, pi - Theta_i on the boundary.
static double totalCurvature(const intrinsic::Geometry& G) {
	double sum = 0;
	for (const auto& v : G.vertices) sum += (v.boundary ? M_PI : 2.0 * M_PI) - v.totalAngle;
	return sum;
}

static double largestDefect(const intrinsic::Geometry& G) {
	double d = 0;
	for (const auto& v : G.vertices) if (!v.boundary) d = std::max(d, std::abs(2.0 * M_PI - v.totalAngle));
	return d * 180.0 / M_PI;
}

struct Errors { double mean = 0, max = 0; int count = 0;
	void add(double e) { mean += e; max = std::max(max, e); ++count; }
	double average() const { return count ? mean / count : 0.0; } };

/// Exact parallel transport on an intrinsically flat mesh: carry a vector of face f0 to every face by unfolding across
/// the edges, the rotation about the shared edge (Section 14). On a flat mesh the result does not depend on the path.
static std::vector<Vector3d> unfoldFrom(Mesh& mesh, const Geometry& g, int f0, const Vector3d& u0) {
	std::vector<Vector3d> u(mesh.numFaces(), Vector3d::Zero());
	std::vector<bool> seen(mesh.numFaces(), false);
	std::vector<int> queue{ f0 };
	u[f0] = u0;
	seen[f0] = true;
	for (size_t k = 0; k < queue.size(); ++k) {
		const int f = queue[k];
		Edge* e = mesh.faceAt(f).edge;
		for (int m = 0; m < 3; ++m, e = e->next) {
			const Face* other = e->pair->face;
			if (!other || seen[other->ID]) continue;
			seen[other->ID] = true;
			u[other->ID] = ga::rotorBetween(g.faceNormal[f], g.faceNormal[other->ID]).rotate(u[f]);
			queue.push_back(other->ID);
		}
	}
	return u;
}

/// Section 24: the vector heat method on a sphere, with the four connections, corner by corner against exact transport
/// along great circles, carried from the tangent plane of the smooth sphere into each face by the smallest rotation.
static void sphereAll(const std::string& dir, int levels, Errors err[4], double& h, int& nv, double& defect, double& face) {
	const std::string obj = dir + "/icosphere_" + std::to_string(levels) + ".obj";
	writeIcosphere(obj, levels);
	Mesh mesh;
	mesh.readOBJ(obj.c_str());
	std::remove(obj.c_str());
	mesh.computeNormals();
	const Geometry g = computeGeometry(mesh);
	const intrinsic::Geometry G = intrinsic::compute(mesh);
	const edgeconn::Geometry E = edgeconn::compute(mesh, G);
	const Methods M{ mesh, g, G };
	h = meanEdgeLength(mesh);
	nv = mesh.numVertices();
	defect = largestDefect(G);
	face = 0;
	for (double o : E.holonomy) face = std::max(face, std::abs(o) * 180 / M_PI);
	int source = 0;
	for (Vertex& v : mesh.getVertices()) if (v.p.z() > mesh.vertexAt(source).p.z()) source = v.ID;
	const Vector3d s = mesh.vertexAt(source).p.normalized();
	const Complex value = std::polar(1.0, 2.0 * M_PI / 3.0);
	const std::vector<Eigen::VectorXcd> u = solveAll(mesh, g, G, E, h * h, source, value);
	// the source vector: what each method puts in the first face of the source, rolled into the plane of the sphere
	int f0, a0;
	firstCorner(mesh, source, f0, a0);
	Vector3d src[4];
	for (int k = 0; k < 4; ++k) src[k] = ga::rotorBetween(g.faceNormal[f0], s).rotate(M.of(k, f0, a0, value));
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* hh[3];
		intrinsic::faceCorners(f, v, hh);
		for (int a = 0; a < 3; ++a) {
			const int i = v[a]->ID;
			const Vector3d x = v[a]->p.normalized();
			if (std::acos(std::max(-1.0, std::min(1.0, s.dot(x)))) > 0.75 * M_PI) continue;
			const ga::Rotor along = ga::rotorBetween(s, x), into = ga::rotorBetween(x, g.faceNormal[f.ID]);
			for (int k = 0; k < 4; ++k) err[k].add(angleDeg(M.of(k, f.ID, a, u[k](i)), into.rotate(along.rotate(src[k]))));
		}
	}
}

/// Section 24: a flat sheet of nx x ny squares, folded along every stripColumns-th column of vertices: parallel folds.
static void writeSheet(const std::string& path, int nx, int ny, int stripColumns, double foldDeg) {
	const double dx = 1.0 / ny;   // square cells, the sheet is nx / ny wide and 1 high
	std::vector<Vector3d> v;
	const double turn = foldDeg * M_PI / 180.0;
	for (int j = 0; j <= ny; ++j)
		for (int i = 0; i <= nx; ++i) {
			// walk along the polygon: every stripColumns cells the direction turns by `turn`
			Vector3d p(0, 0, j * dx);
			double heading = 0;
			for (int k = 0; k < i; ++k) {
				if (k > 0 && k % stripColumns == 0) heading += turn;
				p += dx * Vector3d(std::cos(heading), std::sin(heading), 0);
			}
			v.push_back(p);
		}
	std::ofstream out(path);
	out.precision(17);
	for (auto& x : v) out << "v " << x.x() << " " << x.y() << " " << x.z() << "\n";
	auto id = [&](int i, int j) { return j * (nx + 1) + i + 1; };
	for (int j = 0; j < ny; ++j)
		for (int i = 0; i < nx; ++i) {
			out << "f " << id(i, j) << " " << id(i, j + 1) << " " << id(i + 1, j + 1) << "\n";
			out << "f " << id(i, j) << " " << id(i + 1, j + 1) << " " << id(i + 1, j) << "\n";
		}
}

/// Section 24: a Miura-ori fold, X(i, j) = P_i + Q_j with the zigzags P_i = (i u, 0, (i mod 2) w) and
/// Q_j = ((j mod 2) r, j t, 0). Every face is a parallelogram, cut into two triangles, and the four angles at every
/// vertex add up to 2 pi exactly (acos x + acos(-x) = pi twice), so the folded sheet is isometric to a flat one.
static void writeMiura(const std::string& path, int nx, int ny, double u, double w, double r, double t) {
	std::ofstream out(path);
	out.precision(17);
	for (int j = 0; j <= ny; ++j)
		for (int i = 0; i <= nx; ++i)
			out << "v " << i * u + (j % 2) * r << " " << j * t << " " << (i % 2) * w << "\n";
	auto id = [&](int i, int j) { return j * (nx + 1) + i + 1; };
	for (int j = 0; j < ny; ++j)
		for (int i = 0; i < nx; ++i) {
			out << "f " << id(i, j) << " " << id(i + 1, j) << " " << id(i + 1, j + 1) << "\n";
			out << "f " << id(i, j) << " " << id(i + 1, j + 1) << " " << id(i, j + 1) << "\n";
		}
}

/// Section 24: the vector heat method on an intrinsically flat folded sheet, with the four connections, corner by corner
/// against the exact field obtained by unfolding. Writes the faces with their largest error for Part II.
static void foldExperiment(const char* name, const std::string& obj, int source, Complex value, const std::string& csv) {
	Mesh m;
	m.readOBJ(obj.c_str());
	std::remove(obj.c_str());
	m.computeNormals();
	const Geometry g = computeGeometry(m);
	const intrinsic::Geometry G = intrinsic::compute(m);
	const edgeconn::Geometry E = edgeconn::compute(m, G);
	const Methods M{ m, g, G };
	const double h = meanEdgeLength(m);
	const std::vector<Eigen::VectorXcd> u = solveAll(m, g, G, E, h * h, source, value);
	// the source vector, as each method puts it in the first face of the source, carried exactly to every face
	int f0, a0;
	firstCorner(m, source, f0, a0);
	std::vector<Vector3d> exact[4];
	for (int k = 0; k < 4; ++k) exact[k] = unfoldFrom(m, g, f0, M.of(k, f0, a0, value));
	Errors err[4];
	std::ofstream out;
	if (!csv.empty()) { out.open(csv); out << "x1,y1,z1,x2,y2,z2,x3,y3,z3,ux,uy,uz,error_extrinsic\n"; }
	for (Face& f : m.getFaces()) {
		Vertex* v[3];
		Edge* hh[3];
		intrinsic::faceCorners(f, v, hh);
		double worst = 0;
		for (int a = 0; a < 3; ++a)
			for (int k = 0; k < 4; ++k) {
				const double d = angleDeg(M.of(k, f.ID, a, u[k](v[a]->ID)), exact[k][f.ID]);
				err[k].add(d);
				if (k == 0) worst = std::max(worst, d);
			}
		if (!csv.empty()) {
			for (int a = 0; a < 3; ++a) out << v[a]->p.x() << "," << v[a]->p.y() << "," << v[a]->p.z() << ",";
			out << exact[1][f.ID].x() << "," << exact[1][f.ID].y() << "," << exact[1][f.ID].z() << "," << worst << "\n";
		}
	}
	printf("== %s (%d vertices, largest interior angle defect %.0e degrees):", name, m.numVertices(), largestDefect(G));
	for (int k = 0; k < 4; ++k) printf(" %s mean %.3g max %.3g;", methodName[k], err[k].average(), err[k].max);
	printf(" degrees\n");
}

/// The worked examples of Part III.
static void workedExamplesIntrinsic(const std::string& dir) {
	printf("== Worked example 3: the triangle of Worked example 1 from its three lengths\n");
	{
		intrinsic::FaceGeometry g;
		g.length[0] = 2.0; g.length[1] = std::sqrt(6.0); g.length[2] = std::sqrt(2.0);   // l12, l23, l31
		const double l12 = g.length[0], l23 = g.length[1], l13 = g.length[2];
		printf("  Gram matrix from lengths: [%g %g; %g %g]\n", l12 * l12, 0.5 * (l12 * l12 + l13 * l13 - l23 * l23),
		       0.5 * (l12 * l12 + l13 * l13 - l23 * l23), l13 * l13);
		intrinsic::layOut(g);
		printf("  Heron: area = %.6f; angles %.4f %.4f %.4f degrees\n", g.area, g.angle[0] * 180 / M_PI,
		       g.angle[1] * 180 / M_PI, g.angle[2] * 180 / M_PI);
		printf("  layout: (%g, %g), (%g, %g), (%.6f, %.6f)\n", g.x[0].x(), g.x[0].y(), g.x[1].x(), g.x[1].y(), g.x[2].x(), g.x[2].y());
		for (int a = 0; a < 3; ++a) printf("  grad phi_%d = (%.6f, %.6f), |grad|^2 = %.6f\n", a + 1, g.grad[a].x(), g.grad[a].y(),
		                                   g.grad[a].squaredNorm());
		printf("  stiffness A grad phi_a . grad phi_b:\n");
		for (int a = 0; a < 3; ++a)
			printf("    %9.6f %9.6f %9.6f\n", g.area * g.grad[a].dot(g.grad[0]), g.area * g.grad[a].dot(g.grad[1]),
			       g.area * g.grad[a].dot(g.grad[2]));
	}
	// a vertex with five neighbours: the apex of a pyramid over a regular pentagon, upright and then leaning
	for (int pass = 0; pass < 2; ++pass) {
		const double lean = pass == 0 ? 0.0 : 0.4;
		const std::string obj = dir + "/pyramid.obj";
		{
			std::ofstream out(obj);
			out.precision(17);
			const double h = std::sqrt(std::pow(2.0 * std::sin(M_PI / 5.0), 2) - 1.0);   // equilateral faces when upright
			out << "v " << lean << " 0 " << h << "\n";
			for (int k = 0; k < 5; ++k) out << "v " << std::cos(2 * M_PI * k / 5) << " " << std::sin(2 * M_PI * k / 5) << " 0\n";
			for (int k = 0; k < 5; ++k) out << "f 1 " << k + 2 << " " << (k + 1) % 5 + 2 << "\n";
		}
		Mesh mesh;
		mesh.readOBJ(obj.c_str());
		std::remove(obj.c_str());
		mesh.computeNormals();
		const Geometry g = computeGeometry(mesh);
		const intrinsic::Geometry G = intrinsic::compute(mesh);
		const Vertex& apex = mesh.vertexAt(0);
		const intrinsic::VertexGeometry& va = G.vertices[0];
		printf("== Worked example %d: the apex of a pyramid over a pentagon, %s\n", 4 + pass, pass == 0 ? "upright" : "leaning");
		printf("  angle sum %.4f degrees, defect %.4f degrees, s = %.6f, vertex normal (%.4f, %.4f, %.4f)\n",
		       va.totalAngle * 180 / M_PI, 360 - va.totalAngle * 180 / M_PI, va.scale, apex.n.x(), apex.n.y(), apex.n.z());
		// the corners in counterclockwise order from the reference edge, where each method puts them in the vertex plane
		struct Corner { double alpha, cone, rolled; };
		std::vector<Corner> corners;
		for (Face& f : mesh.getFaces()) {
			Vertex* v[3];
			Edge* h[3];
			intrinsic::faceCorners(f, v, h);   // the apex is the first vertex of every face
			const double alpha = G.faces[f.ID].angle[0];
			const double cone = std::fmod(G.outAngle.at(h[0]) + va.totalAngle, va.totalAngle);
			// Part II: roll the face into the plane of the vertex normal and read the direction of its first edge
			const Vector3d d = ga::rotorBetween(g.faceNormal[f.ID], apex.n).rotate((v[1]->p - v[0]->p).normalized());
			double rolled = std::atan2(g.vertexT[0].cross(d).dot(apex.n), g.vertexT[0].dot(d));
			if (rolled < -1e-9) rolled += 2 * M_PI;
			corners.push_back({ alpha, cone, rolled });
		}
		std::sort(corners.begin(), corners.end(), [](const Corner& a, const Corner& b) { return a.cone < b.cone; });
		double rhoDifference = 0;
		for (Face& f : mesh.getFaces()) {
			const double rhoI = std::arg(G.faces[f.ID].rho[0]);
			const double rhoE = std::arg(ga::connection(g.vertexT[0], apex.n, g.faceT[f.ID], g.faceNormal[f.ID]));
			rhoDifference = std::max(rhoDifference, std::abs(std::arg(std::polar(1.0, rhoI - rhoE))) * 180 / M_PI);
		}
		printf("  largest difference of the two connections rho at the apex: %.3e degrees\n", rhoDifference);
		printf("  corner  alpha    | intrinsic: from   width  | extrinsic: from   width   gap to the next (degrees)\n");
		const double deg = 180 / M_PI;
		for (size_t k = 0; k < corners.size(); ++k) {
			const Corner& c = corners[k];
			const Corner& n = corners[(k + 1) % corners.size()];
			double gap = n.rolled - (c.rolled + c.alpha);
			while (gap > M_PI) gap -= 2 * M_PI;
			while (gap < -M_PI) gap += 2 * M_PI;
			printf("  %zu      %7.3f   |  %8.3f  %7.3f     |  %8.3f  %7.3f   %7.3f\n", k, c.alpha * deg, va.scale * c.cone * deg,
			       va.scale * c.alpha * deg, c.rolled * deg, c.alpha * deg, gap * deg);
		}
	}
}

/// The worked examples of Parts IV and V.
static void workedExamplesEdge(const std::string& dir) {
	// the upright pyramid of Worked example 4: edge transports and face holonomies
	const std::string obj = dir + "/pyramid.obj";
	{
		std::ofstream out(obj);
		out.precision(17);
		const double h = std::sqrt(std::pow(2.0 * std::sin(M_PI / 5.0), 2) - 1.0);
		out << "v 0 0 " << h << "\n";
		for (int k = 0; k < 5; ++k) out << "v " << std::cos(2 * M_PI * k / 5) << " " << std::sin(2 * M_PI * k / 5) << " 0\n";
		for (int k = 0; k < 5; ++k) out << "f 1 " << k + 2 << " " << (k + 1) % 5 + 2 << "\n";
	}
	Mesh mesh;
	mesh.readOBJ(obj.c_str());
	std::remove(obj.c_str());
	mesh.computeNormals();
	const intrinsic::Geometry G = intrinsic::compute(mesh);
	const edgeconn::Geometry E = edgeconn::compute(mesh, G);
	const double deg = 180 / M_PI;
	printf("== Worked example 6: the edge connection of the upright pyramid\n");
	{
		Face& f = mesh.faceAt(0);
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		auto turn = [](double x) { return std::fmod(std::fmod(x, 360.0) + 360.0, 360.0); };   // in [0, 360)
		auto half = [&](double x) { const double t = turn(x); return t > 180 ? t - 360 : t; };   // in (-180, 180]
		double sum = 0;
		for (int a = 0; a < 3; ++a) {
			const double rho = half(E.transport.at(h[a]) * deg);
			sum += rho;
			printf("  edge %d -> %d: phi_ab = %7.3f, phi_ba = %7.3f, rho_ab = %7.3f degrees\n", v[a]->ID, v[(a + 1) % 3]->ID,
			       turn(edgeconn::vertexAngle(G, h[a]) * deg), turn(edgeconn::vertexAngle(G, h[a]->pair) * deg), rho);
		}
		printf("  sum of the three transports %.3f degrees; holonomy Omega_f = sum (s - 1) alpha = %.3f degrees\n", sum,
		       E.holonomy[0] * deg);
		double total = 0;
		for (double o : E.holonomy) total += o;
		printf("  the five face holonomies add up to %.3f degrees, the angle defect of the apex\n", total * deg);
		// the curved element of the face, in the gauge of the face (no vertex frames), against the flat one
		Complex one[3][3], L[3][3], L0[3][3];
		for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) one[a][b] = 1.0;
		edgeconn::curvedElement(G.faces[0], E.holonomy[0], one, 8, L);
		edgeconn::curvedElement(G.faces[0], 0.0, one, 8, L0);
		printf("  kappa = Omega_f / A_f = %.6f; the integrals of the curved element (left) and the flat one (right):\n",
		       E.holonomy[0] / G.faces[0].area);
		for (int a = 0; a < 3; ++a)
			printf("    %8.5f%+8.5fi %8.5f%+8.5fi %8.5f%+8.5fi   | %8.5f %8.5f %8.5f\n", L[a][0].real(), L[a][0].imag(),
			       L[a][1].real(), L[a][1].imag(), L[a][2].real(), L[a][2].imag(), L0[a][0].real(), L0[a][1].real(), L0[a][2].real());
	}
	printf("== Worked example 7: discrete exterior calculus on the triangle of Worked example 1\n");
	{
		intrinsic::FaceGeometry g;
		g.length[0] = 2.0; g.length[1] = std::sqrt(6.0); g.length[2] = std::sqrt(2.0);
		intrinsic::layOut(g);
		// edges 1->2, 2->3, 3->1 opposite to the vertices 3, 1, 2
		const char* name[3] = { "12", "23", "31" };
		for (int a = 0; a < 3; ++a)
			printf("  *1 on edge %s = cot(%.4f deg) / 2 = %.6f\n", name[a], g.angle[(a + 2) % 3] * deg, 0.5 / std::tan(g.angle[(a + 2) % 3]));
		Eigen::Matrix3d d0;   // rows: edges 12, 23, 31; columns: vertices 1, 2, 3
		d0 << -1, 1, 0, 0, -1, 1, 1, 0, -1;
		Eigen::Vector3d star1;
		for (int a = 0; a < 3; ++a) star1(a) = 0.5 / std::tan(g.angle[(a + 2) % 3]);
		const Eigen::Matrix3d L = d0.transpose() * star1.asDiagonal() * d0;
		printf("  d0^T *1 d0 =\n");
		for (int a = 0; a < 3; ++a) printf("    %9.6f %9.6f %9.6f\n", L(a, 0), L(a, 1), L(a, 2));
		// the circumcenter of the layout, and the distances to the edge midpoints: the dual edge pieces
		const Eigen::Vector2d A = g.x[0], B = g.x[1], C = g.x[2];
		const double D = 2 * (A.x() * (B.y() - C.y()) + B.x() * (C.y() - A.y()) + C.x() * (A.y() - B.y()));
		const Eigen::Vector2d cc((A.squaredNorm() * (B.y() - C.y()) + B.squaredNorm() * (C.y() - A.y()) + C.squaredNorm() * (A.y() - B.y())) / D,
		                         (A.squaredNorm() * (C.x() - B.x()) + B.squaredNorm() * (A.x() - C.x()) + C.squaredNorm() * (B.x() - A.x())) / D);
		printf("  circumcenter of the layout (%.6f, %.6f); distance to the midpoints of 12, 23, 31: %.6f %.6f %.6f\n", cc.x(), cc.y(),
		       (cc - 0.5 * (A + B)).norm(), (cc - 0.5 * (B + C)).norm(), (cc - 0.5 * (C + A)).norm());
	}
	printf("== Worked example 8: d^nabla d^nabla on a face of the upright pyramid\n");
	{
		const Complex H = std::polar(1.0, E.holonomy[0]);
		printf("  1 - e^{i Omega_f} = %.6f %+.6fi, |1 - e^{i Omega_f}| = %.6f = 2 sin(Omega_f / 2)\n", (1.0 - H).real(), (1.0 - H).imag(),
		       std::abs(1.0 - H));
		// check with the matrices: d^nabla on 0-forms, then on 1-forms, around the face 0 of the pyramid
		Face& f = mesh.faceAt(0);
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		Complex psi[3] = { Complex(0.3, -0.7), Complex(1.1, 0.2), Complex(-0.4, 0.5) }, eta[3];
		for (int a = 0; a < 3; ++a) eta[a] = psi[(a + 1) % 3] - std::polar(1.0, E.transport.at(h[a])) * psi[a];   // in the frame of the end
		const Complex r0 = std::polar(1.0, E.transport.at(h[0])), r1 = std::polar(1.0, E.transport.at(h[1])),
		              r2 = std::polar(1.0, E.transport.at(h[2]));
		const Complex dd = eta[2] + r2 * eta[1] + r2 * r1 * eta[0];   // all carried to the first vertex
		printf("  for psi = (%.1f%+.1fi, %.1f%+.1fi, %.1f%+.1fi): (d d psi)_f = %.6f%+.6fi, (1 - e^{i Omega_f}) psi_1 = %.6f%+.6fi\n",
		       psi[0].real(), psi[0].imag(), psi[1].real(), psi[1].imag(), psi[2].real(), psi[2].imag(), dd.real(), dd.imag(),
		       ((1.0 - r2 * r1 * r0) * psi[0]).real(), ((1.0 - r2 * r1 * r0) * psi[0]).imag());
	}
}

// ---------------------------------------------------------------------------------------------------------------
// Part VI: how the pieces fit. The leaning pyramid traced through the four methods, the corner transport of Part III
// written as the edge transport of Part IV plus a correction, and the energies of a rotation field on the sphere.
// ---------------------------------------------------------------------------------------------------------------

/// The pyramid of Worked examples 4 and 5: the apex over a regular pentagon of radius 1, moved by `lean` along x.
static void loadPyramid(Mesh& mesh, const std::string& dir, double lean) {
	const std::string obj = dir + "/pyramid.obj";
	{
		std::ofstream out(obj);
		out.precision(17);
		const double h = std::sqrt(std::pow(2.0 * std::sin(M_PI / 5.0), 2) - 1.0);   // equilateral faces when upright
		out << "v " << lean << " 0 " << h << "\n";
		for (int k = 0; k < 5; ++k) out << "v " << std::cos(2 * M_PI * k / 5) << " " << std::sin(2 * M_PI * k / 5) << " 0\n";
		for (int k = 0; k < 5; ++k) out << "f 1 " << k + 2 << " " << (k + 1) % 5 + 2 << "\n";
	}
	mesh.readOBJ(obj.c_str());
	std::remove(obj.c_str());
	mesh.computeNormals();
}

/// An angle in degrees, wrapped to (-180, 180].
static double wrapDeg(double radians) { return std::arg(std::polar(1.0, radians)) * 180.0 / M_PI; }

/// The corner connections e^{T rho_a^f} of a face, with the extrinsic rule (k = 0, Part II) or the intrinsic one (k = 1,
/// Part III).
static void cornerConnections(Mesh& mesh, const Geometry& g, const intrinsic::Geometry& G, int k, const Face& f, Complex rho[3]) {
	Vertex* v[3];
	Edge* h[3];
	intrinsic::faceCorners(f, v, h);
	for (int a = 0; a < 3; ++a)
		rho[a] = k == 0 ? ga::connection(g.vertexT[v[a]->ID], v[a]->n, g.faceT[f.ID], g.faceNormal[f.ID]) : G.faces[f.ID].rho[a];
}

/// Section 40: the corner transport from v[a] to v[a+1] through the face of h[a], divided by the edge transport of the
/// same half edge. Returns its angle and the prediction of Proposition corner_edge, -[(s_a - 1) alpha_a + (s_b - 1) alpha_b]/2.
static void cornerOverEdge(const intrinsic::Geometry& G, const edgeconn::Geometry& E, const Edge* e, const Complex rho[3],
                           double& measured, double& predicted) {
	Vertex* v[3];
	Edge* h[3];
	intrinsic::faceCorners(*e->face, v, h);
	int a = 0;
	while (h[a] != e) ++a;
	const int b = (a + 1) % 3;
	const Complex corner = rho[a] * std::conj(rho[b]);   // from v[a] into the face, then out to v[b]
	measured = std::arg(corner * std::polar(1.0, -E.transport.at(e)));
	const intrinsic::FaceGeometry& fg = G.faces[e->face->ID];
	predicted = -0.5 * ((G.vertices[v[a]->ID].scale - 1) * fg.angle[a] + (G.vertices[v[b]->ID].scale - 1) * fg.angle[b]);
}

/// Section 46, the running example: the leaning pyramid through every step of the pipeline, with the four methods.
static void workedExamplesPieces(const std::string& dir) {
	const double deg = 180 / M_PI;
	for (int pass = 0; pass < 2; ++pass) {
		const double lean = pass == 0 ? 0.0 : 0.4;
		Mesh mesh;
		loadPyramid(mesh, dir, lean);
		const Geometry g = computeGeometry(mesh);
		const intrinsic::Geometry G = intrinsic::compute(mesh);
		const edgeconn::Geometry E = edgeconn::compute(mesh, G);
		const intrinsic::VertexGeometry& apex = G.vertices[0];
		printf("== Worked example %d: the %s pyramid through the pipeline\n", 9 + pass, pass == 0 ? "upright" : "leaning");
		printf("  geometry: Theta = %.3f, defect kappa = %.3f degrees, s = %.6f\n", apex.totalAngle * deg,
		       (2 * M_PI - apex.totalAngle) * deg, apex.scale);
		// step 2 and 3: the corner connections of the apex (II, III) and the edge transports leaving the apex (IV, V)
		printf("  face  alpha_apex | rho II   rho III  | edge 0->k: rho_ab (IV, V) | Omega_f (IV, V) | (s-1) alpha (III)\n");
		double sumOmega = 0;
		for (Face& f : mesh.getFaces()) {
			Vertex* v[3];
			Edge* h[3];
			intrinsic::faceCorners(f, v, h);   // the apex is v[0] of every face
			Complex rhoII[3], rhoIII[3];
			cornerConnections(mesh, g, G, 0, f, rhoII);
			cornerConnections(mesh, g, G, 1, f, rhoIII);
			sumOmega += E.holonomy[f.ID];
			printf("  %d     %7.3f      | %7.3f  %7.3f  | 0->%d: %8.3f               | %7.3f         | %7.3f\n", f.ID,
			       G.faces[f.ID].angle[0] * deg, wrapDeg(std::arg(rhoII[0])), wrapDeg(std::arg(rhoIII[0])), v[1]->ID,
			       wrapDeg(E.transport.at(h[0])), E.holonomy[f.ID] * deg, (apex.scale - 1) * G.faces[f.ID].angle[0] * deg);
		}
		printf("  the face holonomies add up to %.3f degrees, the defect of the apex\n", sumOmega * deg);
		// step 3, the transports between two vertices: along each edge apex -> base k, through the face to its left (f)
		// and to its right (g), with the corner connections of II and III, against the edge transport of IV and V
		printf("  edge 0->k | III through f, g minus edge transport (predicted)        | loop f,g: II      III\n");
		for (Face& f : mesh.getFaces()) {
			Vertex* v[3];
			Edge* h[3];
			intrinsic::faceCorners(f, v, h);
			const Edge* e = h[0];             // apex -> v[1], face f to its left
			const Edge* back = e->pair;       // v[1] -> apex, face g to its left
			Complex rf[2][3], rg[2][3];
			for (int k = 0; k < 2; ++k) {
				cornerConnections(mesh, g, G, k, f, rf[k]);
				cornerConnections(mesh, g, G, k, *back->face, rg[k]);
			}
			double mf, pf, mg, pg;
			cornerOverEdge(G, E, e, rf[1], mf, pf);
			cornerOverEdge(G, E, back, rg[1], mg, pg);   // from v[1] to the apex; reversed, its angle changes sign
			// the loop apex -> f -> v[1] -> g -> apex, with the corner connections of each method
			double loop[2];
			for (int k = 0; k < 2; ++k) {
				Vertex* w[3];
				Edge* hg[3];
				intrinsic::faceCorners(*back->face, w, hg);
				int a = 0;
				while (hg[a] != back) ++a;   // w[a] = v[1], w[a+1] = apex
				const Complex viaF = rf[k][0] * std::conj(rf[k][1]);                    // apex -> v[1] through f
				const Complex viaG = rg[k][a] * std::conj(rg[k][(a + 1) % 3]);          // v[1] -> apex through g
				loop[k] = std::arg(viaF * viaG);
			}
			printf("  0->%d      | f: %8.4f (%8.4f)   g: %8.4f (%8.4f)              | %8.4f  %8.4f\n", v[1]->ID, mf * deg, pf * deg,
			       -mg * deg, -pg * deg, loop[0] * deg, loop[1] * deg);
		}
		// step 5, the element: the entry (apex, v[1]) of the element matrix of face 0 with the four methods
		{
			Face& f = mesh.faceAt(0);
			Vertex* v[3];
			Edge* h[3];
			intrinsic::faceCorners(f, v, h);
			const intrinsic::FaceGeometry& fg = G.faces[0];
			Complex rho[2][3], entry[4], transportBA[3][3], L[3][3];
			for (int k = 0; k < 2; ++k) {
				cornerConnections(mesh, g, G, k, f, rho[k]);
				entry[k] = fg.area * fg.grad[0].dot(fg.grad[1]) * rho[k][1] * std::conj(rho[k][0]);
			}
			for (int a = 0; a < 3; ++a)
				for (int b = 0; b < 3; ++b) transportBA[a][b] = a == b ? Complex(1.0) : std::polar(1.0, E.transport.at(halfEdge(h, b, a)));
			edgeconn::curvedElement(fg, E.holonomy[0], transportBA, 8, L);
			entry[2] = L[0][1];
			entry[3] = fg.area * fg.grad[0].dot(fg.grad[1]) * transportBA[0][1];
			printf("  element of face 0, entry (apex, %d):", v[1]->ID);
			for (int k = 0; k < 4; ++k) printf(" %s %.5f at %.3f deg;", methodName[k], std::abs(entry[k]), wrapDeg(std::arg(entry[k])));
			printf(" diagonal IV %.5f vs flat %.5f\n", L[0][0].real(), fg.area * fg.grad[0].squaredNorm());
		}
		// step 7, the solution: source at the apex, all four methods, compared corner by corner in the faces
		{
			const Methods M{ mesh, g, G };
			const double h = meanEdgeLength(mesh);
			const Complex value = std::polar(1.0, 2.0 * M_PI / 3.0);
			const std::vector<Eigen::VectorXcd> u = solveAll(mesh, g, G, E, h * h, 0, value);
			double worst[4] = { 0, 0, 0, 0 };
			for (Face& f : mesh.getFaces()) {
				Vertex* v[3];
				Edge* hh[3];
				intrinsic::faceCorners(f, v, hh);
				for (int a = 1; a < 3; ++a)   // the base corners
					for (int k = 0; k < 4; ++k)
						worst[k] = std::max(worst[k], angleDeg(M.of(k, f.ID, a, u[k](v[a]->ID)), M.of(1, f.ID, a, u[1](v[a]->ID))));
			}
			printf("  solution (source at the apex): largest angle to III at the base corners:");
			for (int k = 0; k < 4; ++k) printf(" %s %.3f;", methodName[k], worst[k]);
			printf(" degrees; |u| at base vertex 1: %.4f\n", std::abs(u[1](1)));
		}
	}
	// Proposition corner_edge on a large mesh: the corner transport of III is the edge transport of IV turned by half the
	// stretching of the two corners at the ends of the edge
	{
		Mesh mesh;
		mesh.readOBJ("cactus1.obj");
		mesh.CenterAndNormalize();
		mesh.computeNormals();
		const Geometry g = computeGeometry(mesh);
		const intrinsic::Geometry G = intrinsic::compute(mesh);
		const edgeconn::Geometry E = edgeconn::compute(mesh, G);
		double worst = 0, largest = 0, loopSum = 0;
		for (Face& f : mesh.getFaces()) {
			Vertex* v[3];
			Edge* h[3];
			intrinsic::faceCorners(f, v, h);
			Complex rho[3];
			cornerConnections(mesh, g, G, 1, f, rho);
			for (int a = 0; a < 3; ++a) {
				double m, p;
				cornerOverEdge(G, E, h[a], rho, m, p);
				worst = std::max(worst, std::abs(wrapDeg(m - p)));
				largest = std::max(largest, std::abs(m) * deg);
				loopSum += m;   // each edge twice, once per face: the sum of all the edge loops of III
			}
		}
		printf("== Cactus: corner transport (III) = edge transport (IV) turned by -[(s_a-1)alpha_a + (s_b-1)alpha_b]/2: largest "
		       "deviation %.2e degrees (largest correction %.2f degrees); sum of all edge loops %.6f = -2 pi chi = %.6f\n", worst,
		       largest, loopSum, -2.0 * M_PI * (mesh.numVertices() - mesh.numEdges() + mesh.numFaces()));
	}
}

/// Section 43: the rotation field U(x) = e_z x x on the unit sphere. Its Dirichlet energy int |grad U|^2 = 8 pi / 3, its
/// Hodge energy int (div U)^2 + (curl U)^2 = 16 pi / 3, and their difference int K |U|^2 = 8 pi / 3 (Weitzenboeck). The
/// first is computed with the four connection Laplacians, the second with scalar DEC on 1-forms, the third with the
/// angle defects.
static void killingEnergies(const std::string& dir) {
	std::ofstream out(dir + "/killing_energy.csv");
	out << "vertices,extrinsic,intrinsic,curved,dec,hodge,divergence,curvature\n";
	const Vector3d z(0, 0, 1);
	for (int levels = 1; levels <= 5; ++levels) {
		const std::string obj = dir + "/icosphere_k.obj";
		writeIcosphere(obj, levels);
		Mesh mesh;
		mesh.readOBJ(obj.c_str());
		std::remove(obj.c_str());
		mesh.computeNormals();
		const Geometry g = computeGeometry(mesh);
		const intrinsic::Geometry G = intrinsic::compute(mesh);
		const edgeconn::Geometry E = edgeconn::compute(mesh, G);
		const int n = mesh.numVertices();
		auto U = [&](int i) { const Vector3d x = mesh.vertexAt(i).p.normalized(); return Vector3d(z.cross(x)); };
		// the spinor coordinates of U in each method's vertex planes. Part II: the vertex frame. Parts III to V: the
		// direction of U in each face, carried back to the vertex plane by the method's own rule, averaged over the corners.
		std::vector<Eigen::VectorXcd> psi(4, Eigen::VectorXcd::Zero(n));
		for (int i = 0; i < n; ++i) {
			const Vector3d u = U(i);
			psi[0](i) = Complex(g.vertexFrame[i].col(0).dot(u), g.vertexFrame[i].col(1).dot(u));
		}
		for (Face& f : mesh.getFaces()) {
			Vertex* v[3];
			Edge* h[3];
			intrinsic::faceCorners(f, v, h);
			for (int a = 0; a < 3; ++a) {
				const int i = v[a]->ID;
				const Vector3d u = U(i);
				const double beta = std::atan2(g.faceNormal[f.ID].cross(g.faceT[f.ID]).dot(u), g.faceT[f.ID].dot(u));
				const double w = G.faces[f.ID].angle[a];
				psi[1](i) += w * std::polar(1.0, beta - std::arg(G.faces[f.ID].rho[a]));
				psi[2](i) += w * std::polar(1.0, beta - G.faces[f.ID].edgeAngle[a] + edgeconn::vertexAngle(G, h[a]));
			}
		}
		for (int i = 0; i < n; ++i) {
			const double m = U(i).norm();
			for (int k = 1; k <= 2; ++k) psi[k](i) = m * std::polar(1.0, std::arg(psi[k](i)));
		}
		psi[3] = psi[2];
		// the four stiffness matrices, without dt and without boundary rows
		std::vector<std::map<std::pair<int, int>, Complex>> L(3);
		for (Face& f : mesh.getFaces()) {
			Vertex* v[3];
			Edge* h[3];
			intrinsic::faceCorners(f, v, h);
			const intrinsic::FaceGeometry& fg = G.faces[f.ID];
			Complex rho[2][3], transportBA[3][3], C[3][3];
			for (int k = 0; k < 2; ++k) cornerConnections(mesh, g, G, k, f, rho[k]);
			for (int a = 0; a < 3; ++a)
				for (int b = 0; b < 3; ++b) transportBA[a][b] = a == b ? Complex(1.0) : std::polar(1.0, E.transport.at(halfEdge(h, b, a)));
			edgeconn::curvedElement(fg, E.holonomy[f.ID], transportBA, 8, C);
			for (int a = 0; a < 3; ++a)
				for (int b = 0; b < 3; ++b) {
					const double s = fg.area * fg.grad[a].dot(fg.grad[b]);
					for (int k = 0; k < 2; ++k) L[k][{ v[a]->ID, v[b]->ID }] += s * (a == b ? Complex(1.0) : rho[k][b] * std::conj(rho[k][a]));
					L[2][{ v[a]->ID, v[b]->ID }] += C[a][b];
				}
		}
		double energy[4];
		for (int k = 0; k < 3; ++k) {
			Complex e = 0;
			for (const auto& [ij, value] : L[k]) e += std::conj(psi[k](ij.first)) * value * psi[k](ij.second);
			energy[k] = e.real();
		}
		std::vector<double> star0;
		const Eigen::SparseMatrix<Complex> Ldec = decform::laplacian(mesh, G, E, star0);
		energy[3] = (psi[3].adjoint() * (Ldec * psi[3]))(0).real();
		// the Hodge energy of the 1-form U^flat: eta_ab = int_ab U . dx = z . (x_a x x_b), exact along the chord
		std::map<const Edge*, double> eta;
		for (const Edge* e : mesh.getEdges()) eta[e] = z.dot(e->vertex->p.cross(e->pair->vertex->p));
		double curl = 0, divergence = 0, curvature = 0;
		for (Face& f : mesh.getFaces()) {
			Vertex* v[3];
			Edge* h[3];
			intrinsic::faceCorners(f, v, h);
			const double circulation = eta[h[0]] + eta[h[1]] + eta[h[2]];   // d1 eta, the integral of curl U over the face
			curl += circulation * circulation / G.faces[f.ID].area;         // *2 = 1 / A_f
		}
		std::vector<double> codiff(n, 0.0);   // d0^T *1 eta: minus the integral of div U over the dual cell
		auto cotOpposite = [&](const Edge* e) {
			Vertex* v[3];
			Edge* h[3];
			intrinsic::faceCorners(*e->face, v, h);
			int k = 0;
			while (h[k] != e) ++k;
			return 1.0 / std::tan(G.faces[e->face->ID].angle[(k + 2) % 3]);
		};
		for (const Edge* e : mesh.getEdges()) {
			if (e->vertex->ID > e->pair->vertex->ID) continue;
			const double w = 0.5 * (cotOpposite(e) + cotOpposite(e->pair)) * eta[e];
			codiff[e->pair->vertex->ID] += w;
			codiff[e->vertex->ID] -= w;
		}
		for (int i = 0; i < n; ++i) {
			divergence += codiff[i] * codiff[i] / star0[i];
			curvature += (2.0 * M_PI - G.vertices[i].totalAngle) * U(i).squaredNorm();
		}
		printf("== Rotation field on the sphere, %5d vertices: int |grad U|^2: II %.5f, III %.5f, IV %.5f, V %.5f (8 pi/3 = %.5f); "
		       "Hodge %.5f + %.1e (16 pi/3 = %.5f); sum kappa |U|^2 = %.5f\n", n, energy[0], energy[1], energy[2], energy[3],
		       8 * M_PI / 3, curl, divergence, 16 * M_PI / 3, curvature);
		out << n << "," << energy[0] << "," << energy[1] << "," << energy[2] << "," << energy[3] << "," << curl + divergence << ","
		    << divergence << "," << curvature << "\n";
	}
}

int main(int argc, char** argv) {
	const std::string dir = argc > 1 ? argv[1] : "../fem-basic-2d/data";
	workedExamples();

	// --- the cactus of FEM.cpp, both versions -----------------------------------------------------------------
	{
		Mesh mesh;
		mesh.readOBJ("cactus1.obj");
		mesh.CenterAndNormalize();
		mesh.computeNormals();
		const Geometry g = computeGeometry(mesh);
		const double h = meanEdgeLength(mesh);
		int source = 0;
		for (Vertex& v : mesh.getVertices())
			if (v.p.norm() < 2.5e-2 || v.p.z() > 0.83) { source = v.ID; break; }
		const Complex value = std::polar(1.0, 2.0 * M_PI / 3.0);
		const Eigen::VectorXd ur = matrixform::solve(mesh, g, h * h, source, Vector2d(value.real(), value.imag()));
		const Eigen::VectorXcd uc = gaform::solve(mesh, g, h * h, source, value);
		double diff = 0;
		for (int i = 0; i < mesh.numVertices(); ++i)
			diff = std::max(diff, std::abs(Complex(ur(2 * i), ur(2 * i + 1)) - uc(i)));
		printf("== Cactus (%d vertices, source %d): max |real form - complex form| = %.2e\n", mesh.numVertices(), source, diff);
		std::ofstream out(dir + "/cactus_vectors.csv");
		out << "x,y,z,ux,uy,uz\n";
		for (Vertex& v : mesh.getVertices()) {
			const Complex z = uc(v.ID) / std::abs(uc(v.ID));
			const Vector3d u = g.vertexT[v.ID] * z.real() + ga::quarterTurn(g.vertexT[v.ID], ga::tangentBivector(v.n)) * z.imag();
			out << v.p.x() << "," << v.p.y() << "," << v.p.z() << "," << u.x() << "," << u.y() << "," << u.z() << "\n";
		}
		std::ofstream faces(dir + "/cactus_faces.csv");
		faces << "x,y,z\n";
		for (Face& f : mesh.getFaces()) {
			Vertex* v[3];
			faceVertices(f, v);
			for (int a = 0; a < 3; ++a) faces << v[a]->p.x() << "," << v[a]->p.y() << "," << v[a]->p.z() << "\n";
		}
	}

	// --- the sphere: exact parallel transport along great circles ---------------------------------------------
	{
		const std::string obj = dir + "/icosphere.obj";
		writeIcosphere(obj, 4);
		Mesh mesh;
		mesh.readOBJ(obj.c_str());
		mesh.computeNormals();
		const Geometry g = computeGeometry(mesh);
		const double h = meanEdgeLength(mesh);
		int source = 0;
		for (Vertex& v : mesh.getVertices()) if (v.p.z() > mesh.vertexAt(source).p.z()) source = v.ID;
		const Vector3d s = mesh.vertexAt(source).p.normalized();
		const Complex value = std::polar(1.0, 2.0 * M_PI / 3.0);
		const Eigen::VectorXcd uc = gaform::solve(mesh, g, h * h, source, value);
		const Vector3d us = g.vertexT[source] * value.real() + mesh.vertexAt(source).n.cross(g.vertexT[source]) * value.imag();
		std::ofstream out(dir + "/sphere_vectors.csv");
		out << "x,y,z,ux,uy,uz,ex,ey,ez,distance,error_deg,magnitude\n";
		double worst = 0, sum = 0;
		int count = 0;
		for (Vertex& v : mesh.getVertices()) {
			const Vector3d x = v.p.normalized();
			const double d = std::acos(std::max(-1.0, std::min(1.0, s.dot(x))));
			const Complex z = uc(v.ID) / std::abs(uc(v.ID));
			const Vector3d u = g.vertexT[v.ID] * z.real() + v.n.cross(g.vertexT[v.ID]) * z.imag();
			// exact: the source vector carried along the great circle from s to x, by the rotor from s to x
			Vector3d e = d < M_PI - 1e-6 ? ga::rotorBetween(s, x).rotate(us) : Vector3d(us);
			e = (e - e.dot(v.n) * v.n).normalized();
			const double err = std::atan2(u.cross(e).norm(), u.dot(e)) * 180.0 / M_PI;
			if (d < 0.75 * M_PI) { worst = std::max(worst, err); sum += err; ++count; }
			out << x.x() << "," << x.y() << "," << x.z() << "," << u.x() << "," << u.y() << "," << u.z() << ","
			    << e.x() << "," << e.y() << "," << e.z() << "," << d << "," << err << "," << std::abs(uc(v.ID)) << "\n";
		}
		printf("== Sphere (%d vertices, mean edge %.4f): angle to exact transport within 3/4 of the way to the antipode: "
		       "mean %.3f, max %.3f degrees\n", mesh.numVertices(), h, sum / count, worst);
	}

	// --- Part III: the intrinsic connection ---------------------------------------------------------------------
	workedExamplesIntrinsic(dir);
	workedExamplesEdge(dir);
	{   // the sphere under refinement
		std::ofstream out(dir + "/sphere_convergence.csv");
		out << "levels,vertices,h,defect,holonomy,mean_extrinsic,max_extrinsic,mean_intrinsic,max_intrinsic,mean_curved,max_curved,"
		       "mean_dec,max_dec\n";
		for (int levels = 1; levels <= 5; ++levels) {
			Errors err[4];
			double h, defect, face;
			int nv;
			sphereAll(dir, levels, err, h, nv, defect, face);
			printf("== Sphere, %5d vertices, h = %.4f, largest defect %5.2f, largest face holonomy %5.2f:", nv, h, defect, face);
			out << levels << "," << nv << "," << h << "," << defect << "," << face;
			for (int k = 0; k < 4; ++k) {
				printf(" %s %.3f/%.3f", methodName[k], err[k].average(), err[k].max);
				out << "," << err[k].average() << "," << err[k].max;
			}
			printf(" degrees (mean/max)\n");
			out << "\n";
		}
	}
	{   // the cactus: how far apart are the connections?
		Mesh mesh;
		mesh.readOBJ("cactus1.obj");
		mesh.CenterAndNormalize();
		mesh.computeNormals();
		const Geometry g = computeGeometry(mesh);
		const intrinsic::Geometry G = intrinsic::compute(mesh);
		const edgeconn::Geometry E = edgeconn::compute(mesh, G);
		const Methods M{ mesh, g, G };
		const double h = meanEdgeLength(mesh);
		int source = 0;
		for (Vertex& v : mesh.getVertices())
			if (v.p.norm() < 2.5e-2 || v.p.z() > 0.83) { source = v.ID; break; }
		// all coordinates are angles from the same reference edge, so the same value is given to all
		const Complex value = std::polar(1.0, 2.0 * M_PI / 3.0);
		const std::vector<Eigen::VectorXcd> u = solveAll(mesh, g, G, E, h * h, source, value);
		double faceHolonomy = 0, totalHolonomy = 0;
		for (double o : E.holonomy) { faceHolonomy = std::max(faceHolonomy, std::abs(o) * 180 / M_PI); totalHolonomy += o; }
		printf("== Cactus: total curvature %.6f, 2 pi chi = %.6f, sum of face holonomies %.6f, largest defect %.1f, largest "
		       "face holonomy %.1f degrees\n", totalCurvature(G), 2.0 * M_PI * (mesh.numVertices() - mesh.numEdges() + mesh.numFaces()),
		       totalHolonomy, largestDefect(G), faceHolonomy);
		const int pairs[3][2] = { { 0, 1 }, { 1, 2 }, { 2, 3 } };
		std::ofstream out(dir + "/cactus_intrinsic.csv");
		out << "difference_deg,magnitude,difference_curved,difference_dec\n";
		std::vector<double> perVertex[3];
		for (int p = 0; p < 3; ++p) {
			std::vector<double> corners;
			perVertex[p].assign(mesh.numVertices(), 0.0);
			for (Face& f : mesh.getFaces()) {
				Vertex* v[3];
				Edge* hh[3];
				intrinsic::faceCorners(f, v, hh);
				for (int a = 0; a < 3; ++a) {
					const int i = v[a]->ID, k = pairs[p][0], l = pairs[p][1];
					corners.push_back(angleDeg(M.of(k, f.ID, a, u[k](i)), M.of(l, f.ID, a, u[l](i))));
					perVertex[p][i] = std::max(perVertex[p][i], corners.back());
				}
			}
			std::sort(corners.begin(), corners.end());
			printf("   %s vs %s, corner by corner: median %.2f, 90%% %.2f, max %.2f degrees\n", methodName[pairs[p][0]],
			       methodName[pairs[p][1]], corners[corners.size() / 2], corners[corners.size() * 9 / 10], corners.back());
		}
		for (int i = 0; i < mesh.numVertices(); ++i)
			out << perVertex[0][i] << "," << std::abs(u[1](i)) << "," << perVertex[1][i] << "," << perVertex[2][i] << "\n";
	}
	{   // folded sheets: exact bendings of a flat sheet
		const int nx = 48, ny = 16, strip = 4;
		const std::string obj = dir + "/sheet.obj";
		writeSheet(obj, nx, ny, strip, 0.0);
		foldExperiment("Flat sheet        ", obj, (ny / 2) * (nx + 1) + nx / 4, std::polar(1.0, M_PI / 3.0), "");
		writeSheet(obj, nx, ny, strip, 30.0);
		foldExperiment("Parallel folds    ", obj, (ny / 2) * (nx + 1) + nx / 4, std::polar(1.0, M_PI / 3.0), "");
		const int mx = 24, my = 16;
		writeMiura(obj, mx, my, 1.0 / 12, 0.6 / 12, 0.6 / 12, 1.0 / 12);
		foldExperiment("Miura-ori fold    ", obj, (my / 2) * (mx + 1) + mx / 2, std::polar(1.0, M_PI / 3.0),
		               dir + "/miura_faces.csv");
	}

	// --- Part VI: how the pieces fit ------------------------------------------------------------------------------
	workedExamplesPieces(dir);
	killingEnergies(dir);
	return 0;
}
