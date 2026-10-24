/*
	ga_fem.cpp: the method of FEM.cpp written with geometric algebra (GAConnection.h), as in Part II of
	fembasic2d.pdf, side by side with the matrix version of Part I, without any graphics.

	  1. prints the worked examples of the text;
	  2. solves the cactus of FEM.cpp twice, with 2x2 real blocks (Part I, the code of FEM.cpp) and with the
	     n x n complex Hermitian system (Part II), and compares the two;
	  3. solves the same problem on a sphere, where parallel transport along great circles is known exactly
	     (it is the smallest rotor from the source to the point), and measures the error;
	  4. writes CSV files for the figures of the text into the directory given as argument.

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

/// Section 15: L_ij += A_f (grad phi_i . grad phi_j) e^{T (rho_j - rho_i)}, M_ii += A_f / 3, then the same
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
	return 0;
}
