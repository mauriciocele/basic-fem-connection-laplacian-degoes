#pragma once
/*
	The edge connection of Parts IV and V of fembasic2d.pdf: parallel transport along the edges of the mesh, one angle
	per edge, computed from the intrinsic angles of IntrinsicConnection.h, with its curvature, the holonomy of each face.

	Part IV uses it in a finite element whose connection is curved inside each face (curvedElement); Part V uses it in
	the connection Laplacian of discrete exterior calculus (Section 33).

	Every function names the section and equation of the text it implements.
*/
#include <cmath>
#include <complex>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include "IntrinsicConnection.h"

namespace edgeconn {

using Complex = std::complex<double>;
using Vec2 = Eigen::Vector2d;

/// The angle of a half edge leaving vertex i, in the plane of the vertex: phi = s_i theta (eq. scale).
inline double vertexAngle(const intrinsic::Geometry& G, const Edge* e) {
	return G.vertices[e->vertex->ID].scale * G.outAngle.at(e);
}

struct Geometry {
	std::map<const Edge*, double> transport;   // rho_ab for the half edge a -> b (eq. edge_transport)
	std::vector<double> holonomy;              // Omega_f, the holonomy of each face (eq. face_holonomy)
};

/// Section 26: the transport along every edge and the holonomy of every face.
inline Geometry compute(Mesh& mesh, const intrinsic::Geometry& G) {
	Geometry E;
	for (const Edge* e : mesh.getEdges()) {
		// carry a vector from a to b along the edge: keep its angle to the edge, which points along phi_ab at a and
		// along phi_ba + pi at b
		E.transport[e] = vertexAngle(G, e->pair) + M_PI - vertexAngle(G, e);
	}
	E.holonomy.resize(mesh.numFaces());
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		double omega = 0;
		for (int a = 0; a < 3; ++a) omega += (G.vertices[v[a]->ID].scale - 1.0) * G.faces[f.ID].angle[a];
		E.holonomy[f.ID] = omega;
	}
	return E;
}

/// Gauss-Legendre nodes and weights on [0, 1].
inline void gaussLegendre(int n, std::vector<double>& x, std::vector<double>& w) {
	x.resize(n);
	w.resize(n);
	for (int k = 0; k < n; ++k) {
		double t = std::cos(M_PI * (k + 0.75) / (n + 0.5)), dp = 0;
		for (int it = 0; it < 100; ++it) {
			double p0 = 1, p1 = t;
			for (int m = 2; m <= n; ++m) { const double p2 = ((2 * m - 1) * t * p1 - (m - 1) * p0) / m; p0 = p1; p1 = p2; }
			dp = n * (t * p1 - p0) / (t * t - 1);
			const double dt = p1 / dp;
			t -= dt;
			if (std::abs(dt) < 1e-16) break;
		}
		x[k] = 0.5 * (1 - t);
		w[k] = 1.0 / ((1 - t * t) * dp * dp);   // 2 / ((1 - t^2) P'^2), halved for [0, 1]
	}
}

/// Section 27, eq. curved_element: the 3 x 3 Hermitian element matrix of a face whose connection has the constant
/// curvature kappa = Omega_f / A_f, for coordinates in the planes of the three vertices:
///   L_ab = e^{i rho_ba} \int_f e^{i kappa A(x_a, x_b, x)} (grad phi_a - i phi_a w_a) . (grad phi_b + i phi_b w_b) dA,
///   w_a(x) = (kappa / 2) J (x_a - x),
/// with A(x_a, x_b, x) the signed area of the triangle x_a x_b x. Integrated with an n x n Gauss rule on the triangle,
/// seen as a collapsed square.
inline void curvedElement(const intrinsic::FaceGeometry& g, double omega, const Complex transportBA[3][3], int n,
                          Complex L[3][3]) {
	const double kappa = omega / g.area;
	std::vector<double> x, w;
	gaussLegendre(n, x, w);
	auto J = [](const Vec2& v) { return Vec2(-v.y(), v.x()); };
	auto cross = [](const Vec2& p, const Vec2& q) { return p.x() * q.y() - p.y() * q.x(); };
	for (int a = 0; a < 3; ++a) for (int b = 0; b < 3; ++b) L[a][b] = 0;
	for (int i = 0; i < n; ++i)
		for (int j = 0; j < n; ++j) {
			// the collapsed square: lambda = (1 - s, s (1 - t), s t), dA = 2 A s ds dt
			const double s = x[i], t = x[j];
			const double phi[3] = { 1 - s, s * (1 - t), s * t };
			const double weight = w[i] * w[j] * 2.0 * g.area * s;
			const Vec2 p = phi[0] * g.x[0] + phi[1] * g.x[1] + phi[2] * g.x[2];
			Eigen::Vector2cd D[3];   // grad phi_a + i phi_a w_a
			for (int a = 0; a < 3; ++a) {
				const Vec2 wa = 0.5 * kappa * J(g.x[a] - p);
				D[a] = g.grad[a].cast<Complex>() + Complex(0, 1) * phi[a] * wa.cast<Complex>();
			}
			for (int a = 0; a < 3; ++a)
				for (int b = 0; b < 3; ++b) {
					const double area = 0.5 * cross(g.x[b] - g.x[a], p - g.x[a]);
					L[a][b] += weight * std::polar(1.0, kappa * area) * D[a].dot(D[b]);   // dot conjugates D[a]
				}
		}
	for (int a = 0; a < 3; ++a)
		for (int b = 0; b < 3; ++b) L[a][b] *= transportBA[a][b];
}

}   // namespace edgeconn
