#pragma once
/*
	The intrinsic connection of Part III of fembasic2d.pdf: everything the finite element method needs, computed from
	the edge lengths of the mesh alone. No normal, no rotor of space, no position in R^3 enters the computation; the
	positions are read once, to measure the lengths (Section 20), and once more at the end, to draw the result.

	Conventions are those of HalfEdge/: a half edge points away from its vertex (e->vertex is its origin), the face
	of a half edge lies to its left, and e->pair->next is the next outgoing half edge clockwise around the origin.

	Every function names the section and equation of the text it implements.
*/
#include <cmath>
#include <complex>
#include <map>
#include <vector>
#include <Eigen/Dense>
#include "HalfEdge/Mesh.h"

namespace intrinsic {

using Complex = std::complex<double>;
using Vec2 = Eigen::Vector2d;

/// The three vertices of a face, counterclockwise, and the three half edges leaving them inside the face:
/// h[a] goes from v[a] to v[a+1].
inline void faceCorners(const Face& f, Vertex* v[3], Edge* h[3]) {
	h[0] = f.edge; h[1] = f.edge->next; h[2] = f.edge->next->next;
	for (int a = 0; a < 3; ++a) v[a] = h[a]->vertex;
}

/// Heron's formula (eq. heron): the area of a triangle from its three edge lengths.
inline double heronArea(double l0, double l1, double l2) {
	const double s = 0.5 * (l0 + l1 + l2);
	return std::sqrt(std::max(0.0, s * (s - l0) * (s - l1) * (s - l2)));
}

/// The law of cosines (eq. law_of_cosines): the angle opposite to the edge of length c in a triangle with sides a, b, c.
inline double cornerAngle(double a, double b, double c) {
	return std::acos(std::max(-1.0, std::min(1.0, (a * a + b * b - c * c) / (2.0 * a * b))));
}

/// What the method needs of one face, from its lengths (Section 20).
struct FaceGeometry {
	double length[3];   // length[a] = |h[a]|, the edge from v[a] to v[a+1]
	double angle[3];    // the corner angle at v[a]
	double area;
	Vec2 x[3];          // the layout: an isometric copy of the triangle in the plane, v[0] at 0, v[1] on the x axis
	Vec2 grad[3];       // gradients of the hat functions in the layout (eq. intrinsic_gradient)
	double edgeAngle[3];// direction of h[a] in the layout, measured from the x axis (the face frame t_f)
	Complex rho[3];     // the connection from the plane of v[a] into the face, e^{T rho_a^f} (eq. intrinsic_rho)
};

/// What the method needs of one vertex (Section 21).
struct VertexGeometry {
	double totalAngle;  // Theta_i, the sum of the corner angles at the vertex
	double scale;       // s_i = 2 pi / Theta_i for an interior vertex, 1 at the boundary (eq. scale)
	bool boundary;
};

struct Geometry {
	std::vector<FaceGeometry> faces;
	std::vector<VertexGeometry> vertices;
	std::map<const Edge*, double> outAngle;   // the cone angle theta of every outgoing half edge (eq. cone_angle)
};

/// Section 20: lay a face out in the plane from its lengths alone and compute its area, angles and gradients.
inline void layOut(FaceGeometry& g) {
	const double l01 = g.length[0], l12 = g.length[1], l20 = g.length[2];
	g.angle[0] = cornerAngle(l01, l20, l12);
	g.angle[1] = cornerAngle(l01, l12, l20);
	g.angle[2] = M_PI - g.angle[0] - g.angle[1];
	g.area = heronArea(l01, l12, l20);
	g.x[0] = Vec2(0, 0);
	g.x[1] = Vec2(l01, 0);
	g.x[2] = l20 * Vec2(std::cos(g.angle[0]), std::sin(g.angle[0]));
	for (int a = 0; a < 3; ++a) {
		const Vec2 e = g.x[(a + 1) % 3] - g.x[a];
		g.edgeAngle[a] = std::atan2(e.y(), e.x());
		// the opposite edge, from v[a+1] to v[a+2], turned a quarter turn counterclockwise, over twice the area
		const Vec2 opposite = g.x[(a + 2) % 3] - g.x[(a + 1) % 3];
		g.grad[a] = Vec2(-opposite.y(), opposite.x()) / (2.0 * g.area);
	}
}

/// Sections 20 to 22: the whole intrinsic geometry. The positions are used only to measure the edge lengths.
inline Geometry compute(Mesh& mesh) {
	Geometry G;
	G.faces.resize(mesh.numFaces());
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		faceCorners(f, v, h);
		for (int a = 0; a < 3; ++a) G.faces[f.ID].length[a] = (v[(a + 1) % 3]->p - v[a]->p).norm();   // the only use of p
		layOut(G.faces[f.ID]);
	}
	// the corner angle of the face to the left of an outgoing half edge, at its origin
	auto angleAtOrigin = [&](const Edge* e) {
		Vertex* v[3];
		Edge* h[3];
		faceCorners(*e->face, v, h);
		for (int a = 0; a < 3; ++a) if (h[a] == e) return G.faces[e->face->ID].angle[a];
		return 0.0;
	};
	// Section 21: walk clockwise around each vertex from its reference edge, accumulating the cone angle
	G.vertices.resize(mesh.numVertices());
	for (Vertex& vi : mesh.getVertices()) {
		const Edge* e = vi.edge;
		double theta = 0.0;
		G.outAngle[e] = theta;
		bool boundary = false;
		while (true) {
			const Edge* pair = e->pair;
			if (!pair->face) { boundary = true; break; }   // past the last face of a boundary vertex
			const Edge* next = pair->next;                  // the next outgoing half edge, clockwise
			theta -= angleAtOrigin(next);                    // the face between e and next lies to the left of next
			if (next == vi.edge) break;
			G.outAngle[next] = theta;
			e = next;
		}
		VertexGeometry& g = G.vertices[vi.ID];
		g.totalAngle = -theta;
		g.boundary = boundary || !vi.edge->face;
		g.scale = g.boundary ? 1.0 : 2.0 * M_PI / g.totalAngle;
	}
	// eq. intrinsic_rho: compare the bisector of each corner in the face and in the plane of the vertex
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		faceCorners(f, v, h);
		FaceGeometry& g = G.faces[f.ID];
		for (int a = 0; a < 3; ++a) {
			const double s = G.vertices[v[a]->ID].scale;
			const double theta = G.outAngle.at(h[a]);   // the corner spans [theta, theta + angle] in the cone
			const double inFace = g.edgeAngle[a] + 0.5 * g.angle[a];
			const double inVertex = s * (theta + 0.5 * g.angle[a]);
			g.rho[a] = std::polar(1.0, inFace - inVertex);
		}
	}
	return G;
}

/// Section 23, drawing the result: the unit vector of R^3 that the spinor coordinate z at vertex i points to. The
/// direction of z in the plane of the vertex is the cone angle arg(z) / s_i; it is drawn in the face whose corner
/// contains that angle, turned from the face's first edge by the remaining angle.
inline Eigen::Vector3d lift(Mesh& mesh, const Geometry& G, int i, Complex z) {
	Vertex& vi = mesh.vertexAt(i);
	const VertexGeometry& gv = G.vertices[i];
	double c = std::arg(z) / gv.scale;   // the cone angle, in (-Theta, Theta) for an interior vertex
	if (!gv.boundary && c > 0) c -= gv.totalAngle;   // the corners cover [-Theta, 0]
	const Edge* best = nullptr;
	double bestDistance = 1e300;
	for (Vertex::EdgeAroundIterator it = vi.iterator(); !it.end(); it++) {
		const Edge* e = it.edge_out();
		if (!e->face) continue;
		Vertex* v[3];
		Edge* h[3];
		faceCorners(*e->face, v, h);
		int a = 0;
		while (h[a] != e) ++a;
		const double from = G.outAngle.at(e), to = from + G.faces[e->face->ID].angle[a];
		const double distance = c < from ? from - c : c > to ? c - to : 0.0;
		if (distance < bestDistance) { bestDistance = distance; best = e; }
	}
	const Eigen::Vector3d p0 = best->vertex->p, p1 = best->pair->vertex->p, p2 = best->next->next->vertex->p;
	const Eigen::Vector3d n = (p1 - p0).cross(p2 - p0).normalized();
	const Eigen::Vector3d edge = (p1 - p0).normalized();
	const double turn = c - G.outAngle.at(best);
	return std::cos(turn) * edge + std::sin(turn) * n.cross(edge);
}

}   // namespace intrinsic
