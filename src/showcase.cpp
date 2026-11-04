/*
	showcase.cpp: two demonstrations of the parallel transport of fembasic2d.pdf, rendered to images without any
	graphics library.

	  1. A decal on the cactus, pasted with the logarithmic map: geodesic polar coordinates around one point. The
	     distance r comes from the heat method (one scalar heat step, then a Poisson solve), and the angle theta from
	     the vector heat method with the intrinsic connection of Part III: theta is the angle between the direction in
	     which r grows and the transported source direction. The texture coordinates are r (cos theta, sin theta).
	  2. The Miura-ori fold, with the same logarithmic map computed with the connection of Part II (vertex normals) and
	     with that of Part III (edge lengths). The sheet is intrinsically flat, so the exact logarithmic map is the
	     unfolding of the paper. The distances are taken from that unfolding, so the two pictures differ only in the
	     connection.

	Build, from src/:
	  c++ -std=c++17 -O2 -I<eigen> showcase.cpp HalfEdge/Mesh.cpp -o showcase && ./showcase <output directory>
	It writes showcase_cactus.ppm, showcase_miura.ppm (extrinsic left, intrinsic right) and prints the errors.
*/
#include <algorithm>
#include <array>
#include <cstdio>
#include <fstream>
#include <functional>
#include <string>
#include <vector>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include "GAConnection.h"
#include "HalfEdge/Mesh.h"
#include "IntrinsicConnection.h"

using Eigen::Vector2d;
using Eigen::Vector3d;
using Complex = std::complex<double>;
using Corners = std::vector<std::array<Complex, 3>>;

// ---------------------------------------------------------------------------------------------------------------
// The connections and the solvers
// ---------------------------------------------------------------------------------------------------------------

/// The corner connections e^{T rho_a^f} of every face, in the frame of the face whose first axis is its first edge.
/// That frame is both the layout of Part III and the face frame t_f of Part II. extrinsic: Part II; else Part III.
static Corners cornerConnections(Mesh& mesh, const intrinsic::Geometry& G, bool extrinsic) {
	Corners rho(mesh.numFaces());
	std::vector<Vector3d> t(mesh.numVertices());   // Part II: the reference edge projected into the vertex plane
	for (Vertex& v : mesh.getVertices()) {
		const Vector3d e = v.edge->pair->vertex->p - v.p;
		t[v.ID] = (e - e.dot(v.n) * v.n).normalized();
	}
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		const Vector3d nf = (v[1]->p - v[0]->p).cross(v[2]->p - v[0]->p).normalized();
		const Vector3d tf = (v[1]->p - v[0]->p).normalized();
		for (int a = 0; a < 3; ++a)
			rho[f.ID][a] = extrinsic ? ga::connection(t[v[a]->ID], v[a]->n, tf, nf) : G.faces[f.ID].rho[a];
	}
	return rho;
}

/// One backward Euler step of the vector heat equation, (M + dt L) psi = 0 away from the source and psi = 1 at the
/// source: the reference direction of the source, carried to every vertex (Stages 6 to 8 of the text).
static Eigen::VectorXcd vectorHeat(Mesh& mesh, const intrinsic::Geometry& G, const Corners& rho, double dt, int source) {
	const int n = mesh.numVertices();
	std::vector<Eigen::Triplet<Complex>> T;
	std::vector<double> mass(n, 0.0);
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		const intrinsic::FaceGeometry& g = G.faces[f.ID];
		for (int a = 0; a < 3; ++a) {
			mass[v[a]->ID] += g.area / 3.0;
			if (v[a]->ID == source) continue;
			for (int b = 0; b < 3; ++b) {
				const Complex transport = a == b ? Complex(1.0) : rho[f.ID][b] * std::conj(rho[f.ID][a]);
				T.emplace_back(v[a]->ID, v[b]->ID, dt * g.area * g.grad[a].dot(g.grad[b]) * transport);
			}
		}
	}
	for (int i = 0; i < n; ++i) T.emplace_back(i, i, i == source ? 1.0 : mass[i]);
	Eigen::SparseMatrix<Complex> A(n, n);
	A.setFromTriplets(T.begin(), T.end());
	Eigen::VectorXcd rhs = Eigen::VectorXcd::Zero(n);
	rhs(source) = 1.0;
	Eigen::SparseLU<Eigen::SparseMatrix<Complex>> lu(A);
	return lu.solve(rhs);
}

/// The scalar stiffness K_ab = sum_f A_f grad phi_a . grad phi_b, with the rows of the pinned vertex replaced.
static Eigen::VectorXd scalarSolve(Mesh& mesh, const intrinsic::Geometry& G, double dt, bool withMass, int pinned,
                                   double pinnedValue, const Eigen::VectorXd& rhs) {
	const int n = mesh.numVertices();
	std::vector<Eigen::Triplet<double>> T;
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		const intrinsic::FaceGeometry& g = G.faces[f.ID];
		for (int a = 0; a < 3; ++a) {
			if (v[a]->ID == pinned) continue;
			if (withMass) T.emplace_back(v[a]->ID, v[a]->ID, g.area / 3.0);
			for (int b = 0; b < 3; ++b) T.emplace_back(v[a]->ID, v[b]->ID, dt * g.area * g.grad[a].dot(g.grad[b]));
		}
	}
	T.emplace_back(pinned, pinned, 1.0);
	Eigen::SparseMatrix<double> A(n, n);
	A.setFromTriplets(T.begin(), T.end());
	Eigen::VectorXd b = rhs;
	b(pinned) = pinnedValue;
	Eigen::SparseLU<Eigen::SparseMatrix<double>> lu(A);
	return lu.solve(b);
}

/// The gradient of a piecewise linear function on face f, in the layout (face) frame.
static Vector2d faceGradient(Mesh& mesh, const intrinsic::Geometry& G, const Face& f, const Eigen::VectorXd& u) {
	Vertex* v[3];
	Edge* h[3];
	intrinsic::faceCorners(f, v, h);
	const intrinsic::FaceGeometry& g = G.faces[f.ID];
	return u(v[0]->ID) * g.grad[0] + u(v[1]->ID) * g.grad[1] + u(v[2]->ID) * g.grad[2];
}

/// The heat method for geodesic distance: one scalar heat step from the source, the unit vector field X = -grad u /
/// |grad u| on the faces, and the function whose gradient is closest to X, K phi = b with b_a = sum_f A_f grad phi_a . X.
static Eigen::VectorXd heatDistance(Mesh& mesh, const intrinsic::Geometry& G, double dt, int source) {
	const int n = mesh.numVertices();
	const Eigen::VectorXd u = scalarSolve(mesh, G, dt, true, source, 1.0, Eigen::VectorXd::Zero(n));
	Eigen::VectorXd b = Eigen::VectorXd::Zero(n);
	for (Face& f : mesh.getFaces()) {
		const Vector2d gu = faceGradient(mesh, G, f, u);
		const Vector2d X = -gu / std::max(gu.norm(), 1e-300);
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		for (int a = 0; a < 3; ++a) b(v[a]->ID) += G.faces[f.ID].area * G.faces[f.ID].grad[a].dot(X);
	}
	return scalarSolve(mesh, G, 1.0, false, source, 0.0, b);
}

/// The logarithmic map: at every vertex, r e^{T theta}, where theta is the angle from the transported source direction
/// Y to the direction in which r grows, R. Both are compared inside each face around the vertex (Y carried into the
/// face with the corner connection), and the corners are averaged with their angles as weights.
static std::vector<Vector2d> logMap(Mesh& mesh, const intrinsic::Geometry& G, const Corners& rho, const Eigen::VectorXcd& Y,
                                    const Eigen::VectorXd& r, const std::function<Complex(const Face&, int)>& radial) {
	const int n = mesh.numVertices();
	std::vector<Complex> sum(n, 0.0);
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		for (int a = 0; a < 3; ++a) {
			const Complex R = radial(f, a), y = rho[f.ID][a] * Y(v[a]->ID);
			if (std::abs(R) < 1e-14 || std::abs(y) < 1e-300) continue;
			sum[v[a]->ID] += G.faces[f.ID].angle[a] * (R / std::abs(R)) * std::conj(y / std::abs(y));
		}
	}
	std::vector<Vector2d> uv(n);
	for (int i = 0; i < n; ++i) {
		const double theta = std::abs(sum[i]) > 0 ? std::arg(sum[i]) : 0.0;
		uv[i] = r(i) * Vector2d(std::cos(theta), std::sin(theta));
	}
	return uv;
}

/// An intrinsically flat mesh laid out in the plane, face after face across the edges: the exact unfolding. Returns
/// the plane position of every vertex, and for every face the angle of the rotation from its layout into the plane.
static std::vector<Vector2d> unfold(Mesh& mesh, const intrinsic::Geometry& G, std::vector<double>& faceRotation) {
	std::vector<Vector2d> p(mesh.numVertices(), Vector2d::Zero());
	faceRotation.assign(mesh.numFaces(), 0.0);
	std::vector<bool> done(mesh.numFaces(), false);
	std::vector<int> queue{ 0 };
	{
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(mesh.faceAt(0), v, h);
		for (int a = 0; a < 3; ++a) p[v[a]->ID] = G.faces[0].x[a];
	}
	done[0] = true;
	for (size_t k = 0; k < queue.size(); ++k) {
		Face& f = mesh.faceAt(queue[k]);
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		for (int a = 0; a < 3; ++a) {
			const Face* g = h[a]->pair->face;
			if (!g || done[g->ID]) continue;
			Vertex* w[3];
			Edge* hg[3];
			intrinsic::faceCorners(*g, w, hg);
			int i = 0;
			while (hg[i] != h[a]->pair) ++i;   // hg[i] goes from w[i] = v[a+1] to w[i+1] = v[a]
			const intrinsic::FaceGeometry& lg = G.faces[g->ID];
			const Vector2d A = p[w[i]->ID], B = p[w[(i + 1) % 3]->ID];
			const Vector2d dl = lg.x[(i + 1) % 3] - lg.x[i], dp = B - A;
			const double gamma = std::atan2(dp.y(), dp.x()) - std::atan2(dl.y(), dl.x());
			const Eigen::Rotation2Dd Rg(gamma);
			p[w[(i + 2) % 3]->ID] = A + Rg * (lg.x[(i + 2) % 3] - lg.x[i]);
			faceRotation[g->ID] = gamma;
			done[g->ID] = true;
			queue.push_back(g->ID);
		}
	}
	return p;
}

// ---------------------------------------------------------------------------------------------------------------
// A small z-buffered renderer with texture coordinates and diffuse light
// ---------------------------------------------------------------------------------------------------------------
struct Image {
	int w, h;
	std::vector<Vector3d> c;
	std::vector<double> z;
	Image(int w_, int h_) : w(w_), h(h_), c(w_ * h_, Vector3d(1, 1, 1)), z(w_ * h_, -1e300) {}
};

struct Camera {
	Vector3d right, up, toward;   // toward points from the scene to the eye
	Camera(double elevDeg, double azimDeg) {
		const double e = elevDeg * M_PI / 180, a = azimDeg * M_PI / 180;
		right = Vector3d(-std::sin(a), std::cos(a), 0);
		up = Vector3d(-std::sin(e) * std::cos(a), -std::sin(e) * std::sin(a), std::cos(e));
		toward = Vector3d(std::cos(e) * std::cos(a), std::cos(e) * std::sin(a), std::sin(e));
	}
};

using Texture = std::function<Vector3d(const Vector2d&)>;

/// The colour of a pixel, from the face, its three vertices and the barycentric coordinates of the pixel in it.
using ColorAt = std::function<Vector3d(const Face&, Vertex* const v[3], const double l[3])>;

/// Rasterize the mesh into the columns [x0, x0 + width), with a z-buffer and diffuse light. The view is fitted to the
/// vertices accepted by `fit` (all by default).
static void raster(Mesh& mesh, const ColorAt& colorAt, const Camera& cam, Image& img, int x0, int width, bool smooth,
                   double ambient, const std::function<bool(int)>& fit = [](int) { return true; }) {
	double minx = 1e300, maxx = -1e300, miny = 1e300, maxy = -1e300;
	for (Vertex& v : mesh.getVertices()) {
		if (!fit(v.ID)) continue;
		const double x = v.p.dot(cam.right), y = v.p.dot(cam.up);
		minx = std::min(minx, x); maxx = std::max(maxx, x); miny = std::min(miny, y); maxy = std::max(maxy, y);
	}
	const double scale = 0.92 * std::min(width / (maxx - minx), img.h / (maxy - miny));
	const double cx = x0 + 0.5 * width, cy = 0.5 * img.h, mx = 0.5 * (minx + maxx), my = 0.5 * (miny + maxy);
	auto project = [&](const Vector3d& p) {
		return Vector3d(cx + scale * (p.dot(cam.right) - mx), cy - scale * (p.dot(cam.up) - my), p.dot(cam.toward));
	};
	const Vector3d light = (cam.toward + 0.6 * cam.up - 0.4 * cam.right).normalized();
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		const Vector3d s[3] = { project(v[0]->p), project(v[1]->p), project(v[2]->p) };
		const Vector3d nf = (v[1]->p - v[0]->p).cross(v[2]->p - v[0]->p).normalized();
		const double area = (s[1].x() - s[0].x()) * (s[2].y() - s[0].y()) - (s[2].x() - s[0].x()) * (s[1].y() - s[0].y());
		if (std::abs(area) < 1e-12) continue;
		const int xa = std::max(x0, (int)std::floor(std::min({ s[0].x(), s[1].x(), s[2].x() })));
		const int xb = std::min(x0 + width - 1, (int)std::ceil(std::max({ s[0].x(), s[1].x(), s[2].x() })));
		const int ya = std::max(0, (int)std::floor(std::min({ s[0].y(), s[1].y(), s[2].y() })));
		const int yb = std::min(img.h - 1, (int)std::ceil(std::max({ s[0].y(), s[1].y(), s[2].y() })));
		for (int py = ya; py <= yb; ++py)
			for (int px = xa; px <= xb; ++px) {
				const double qx = px + 0.5, qy = py + 0.5;
				double l[3];
				for (int k = 0; k < 3; ++k) {   // barycentric coordinates from the signed areas
					const Vector3d& p = s[(k + 1) % 3];
					const Vector3d& q = s[(k + 2) % 3];
					l[k] = ((p.x() - qx) * (q.y() - qy) - (q.x() - qx) * (p.y() - qy)) / area;
				}
				if (l[0] < 0 || l[1] < 0 || l[2] < 0) continue;
				const double depth = l[0] * s[0].z() + l[1] * s[1].z() + l[2] * s[2].z();
				const int idx = py * img.w + px;
				if (depth <= img.z[idx]) continue;
				img.z[idx] = depth;
				const Vector3d n = smooth ? Vector3d((l[0] * v[0]->n + l[1] * v[1]->n + l[2] * v[2]->n).normalized()) : nf;
				img.c[idx] = (ambient + (1 - ambient) * std::abs(n.dot(light))) * colorAt(f, v, l);
			}
	}
}

/// Texture coordinates of one face, with the guard: where the map is stretched more than 4 times in area (past the
/// cut locus of the logarithmic map) the face gets no texture, signalled by coordinates far outside any decal.
static Vector2d texCoord(const std::vector<Vector2d>& uv, Vertex* const v[3], const double l[3], bool guard) {
	if (guard) {
		const Vector3d cr = (v[1]->p - v[0]->p).cross(v[2]->p - v[0]->p);
		const Vector2d du = uv[v[1]->ID] - uv[v[0]->ID], dv = uv[v[2]->ID] - uv[v[0]->ID];
		if (std::abs(du.x() * dv.y() - du.y() * dv.x()) > 4.0 * cr.norm()) return Vector2d(1e9, 1e9);
	}
	return l[0] * uv[v[0]->ID] + l[1] * uv[v[1]->ID] + l[2] * uv[v[2]->ID];
}

/// Render the mesh with one texture map.
static void render(Mesh& mesh, const std::vector<Vector2d>& uv, const Texture& tex, const Camera& cam, Image& img,
                   int x0, int width, bool smooth, double ambient, bool guard = false,
                   const std::function<bool(int)>& fit = [](int) { return true; }) {
	raster(mesh, [&](const Face&, Vertex* const v[3], const double l[3]) { return tex(texCoord(uv, v, l, guard)); }, cam, img,
	       x0, width, smooth, ambient, fit);
}

/// A soft vertical gradient behind the columns [x0, x0 + width).
static void background(Image& img, int x0, int width, const Vector3d& top, const Vector3d& bottom) {
	for (int y = 0; y < img.h; ++y)
		for (int x = x0; x < x0 + width; ++x) img.c[y * img.w + x] = top + (bottom - top) * (double)y / (img.h - 1);
}

/// Write the image downsampled by `ss` in each direction (box filter: supersampled antialiasing), as binary PPM.
static void writePPM(const Image& img, int ss, const std::string& path) {
	const int w = img.w / ss, h = img.h / ss;
	std::ofstream out(path, std::ios::binary);
	out << "P6\n" << w << " " << h << "\n255\n";
	for (int y = 0; y < h; ++y)
		for (int x = 0; x < w; ++x) {
			Vector3d c = Vector3d::Zero();
			for (int j = 0; j < ss; ++j)
				for (int i = 0; i < ss; ++i) c += img.c[(y * ss + j) * img.w + x * ss + i];
			c /= ss * ss;
			for (int k = 0; k < 3; ++k) out.put((char)(unsigned char)std::max(0.0, std::min(255.0, 255.0 * c(k))));
		}
}

/// A decal: a checkerboard of cell `cell` inside the radius `radius`, with an arrow along the source direction (red)
/// and a shorter one a quarter turn to its left (blue), so that both orientation and handedness are visible.
static Texture decal(double radius, double cell, const Vector3d& outside, double arrow) {
	return [=](const Vector2d& t) -> Vector3d {
		const double r = t.norm();
		if (r > radius) return outside;
		if (r > radius - 0.012 * arrow) return Vector3d(0.1, 0.1, 0.1);
		const double w = 0.045 * arrow;
		if (t.x() > 0 && t.x() < 0.62 * arrow && std::abs(t.y()) < w) return Vector3d(0.85, 0.1, 0.1);
		if (t.x() >= 0.62 * arrow && t.x() < 0.8 * arrow && std::abs(t.y()) < 2.4 * w * (0.8 * arrow - t.x()) / (0.18 * arrow))
			return Vector3d(0.85, 0.1, 0.1);
		if (t.y() > 0 && t.y() < 0.32 * arrow && std::abs(t.x()) < 0.8 * w) return Vector3d(0.1, 0.25, 0.85);
		if (t.y() >= 0.32 * arrow && t.y() < 0.44 * arrow && std::abs(t.x()) < 2.0 * w * (0.44 * arrow - t.y()) / (0.12 * arrow))
			return Vector3d(0.1, 0.25, 0.85);
		const int k = (int)std::floor(t.x() / cell) + (int)std::floor(t.y() / cell);
		return (k & 1) ? Vector3d(0.98, 0.83, 0.45) : Vector3d(1.0, 0.98, 0.92);
	};
}

/// A polar decal: rings of constant geodesic distance every `cell` and spokes every 30 degrees, on alternating bands.
static Texture polar(double radius, double cell, const Vector3d& light, const Vector3d& dark, const Vector3d& line,
                     const Vector3d& outside) {
	return [=](const Vector2d& t) -> Vector3d {
		const double r = t.norm();
		if (r > radius) return outside;
		if (r > radius - 0.012 * radius) return Vector3d(0.1, 0.1, 0.1);
		if (r < 0.05 * radius) return line;
		const double ring = r / cell - std::floor(r / cell + 0.5);
		if (std::abs(ring) * cell < 0.012 * radius) return line;
		const double a = std::atan2(t.y(), t.x()) / (M_PI / 6);
		if (std::abs(a - std::floor(a + 0.5)) * (M_PI / 6) * r < 0.01 * radius) return line;
		return ((int)std::floor(r / cell) & 1) ? dark : light;
	};
}

/// Graph paper: thin lines every `cell` and thicker ones every 4 cells, with the two arrows at the origin. Straight
/// lines on the surface mean that the map is the unfolding of the paper.
static Texture gridPaper(double cell, double arrow) {
	const Texture arrows = decal(1e300, 1e300, Vector3d(1, 1, 1), arrow);
	return [=](const Vector2d& t) -> Vector3d {
		const Vector3d a = arrows(t);
		if ((a - Vector3d(1.0, 0.98, 0.92)).norm() > 1e-9 && (a - Vector3d(0.98, 0.83, 0.45)).norm() > 1e-9) return a;
		auto dist = [&](double x, double c) { const double f = x / c - std::floor(x / c + 0.5); return std::abs(f) * c; };
		if (dist(t.x(), 4 * cell) < 0.006 || dist(t.y(), 4 * cell) < 0.006) return Vector3d(0.1, 0.2, 0.45);
		if (dist(t.x(), cell) < 0.0025 || dist(t.y(), cell) < 0.0025) return Vector3d(0.35, 0.5, 0.75);
		return Vector3d(1.0, 0.99, 0.95);
	};
}

// ---------------------------------------------------------------------------------------------------------------
// The two demonstrations
// ---------------------------------------------------------------------------------------------------------------
static double meanEdgeLength(Mesh& mesh) {
	double sum = 0;
	int count = 0;
	for (auto& e : mesh.getEdges()) { sum += (e->pair->vertex->p - e->vertex->p).norm(); ++count; }
	return sum / count;
}

static void writeMiura(const std::string& path, int nx, int ny, double u, double w, double r, double t) {
	std::ofstream out(path);
	out.precision(17);
	for (int j = 0; j <= ny; ++j)
		for (int i = 0; i <= nx; ++i) out << "v " << i * u + (j % 2) * r << " " << j * t << " " << (i % 2) * w << "\n";
	auto id = [&](int i, int j) { return j * (nx + 1) + i + 1; };
	for (int j = 0; j < ny; ++j)
		for (int i = 0; i < nx; ++i) {
			out << "f " << id(i, j) << " " << id(i + 1, j) << " " << id(i + 1, j + 1) << "\n";
			out << "f " << id(i, j) << " " << id(i + 1, j + 1) << " " << id(i, j + 1) << "\n";
		}
}

/// The logarithmic map around `source` with the intrinsic connection, turned so that the source direction points up
/// (along e_z projected into the surface). Returns the texture coordinates, and the distances in r.
static std::vector<Vector2d> uprightLogMap(Mesh& mesh, const intrinsic::Geometry& G, const Corners& rho, double h,
                                           int source, Eigen::VectorXd& r) {
	const Eigen::VectorXcd Y = vectorHeat(mesh, G, rho, h * h, source);
	r = heatDistance(mesh, G, h * h, source);
	auto radial = [&](const Face& f, int) {
		const Vector2d g = faceGradient(mesh, G, f, r);
		return Complex(g.x(), g.y());
	};
	std::vector<Vector2d> uv = logMap(mesh, G, rho, Y, r, radial);
	// turn the decal upright: rotate the texture coordinates so that the source direction (the red arrow) points up
	// the trunk. In the first face of the source, the source direction is the face vector at the angle arg(rho), and
	// "up" is e_z projected into the face.
	{
		const Edge* e = mesh.vertexAt(source).edge;
		const Face& f = *e->face;
		Vertex* v[3];
		Edge* hh[3];
		intrinsic::faceCorners(f, v, hh);
		int a = 0;
		while (hh[a] != e) ++a;
		const Vector3d nf = (v[1]->p - v[0]->p).cross(v[2]->p - v[0]->p).normalized();
		const Vector3d tf = (v[1]->p - v[0]->p).normalized(), bf = nf.cross(tf);
		const Vector3d up = Vector3d(0, 0, 1) - nf.z() * nf;
		const double alpha = std::atan2(up.dot(bf), up.dot(tf)) - std::arg(rho[f.ID][a]);
		const Eigen::Rotation2Dd turn(-alpha);
		for (Vector2d& t : uv) t = turn * t;
	}
	return uv;
}

static void cactusDemo(const std::string& dir, int source, double radius, double elev, double azim) {
	Mesh mesh;
	mesh.readOBJ("cactus2.obj");   // the high-resolution cactus (5261 vertices)
	mesh.CenterAndNormalize();
	mesh.computeNormals();
	const intrinsic::Geometry G = intrinsic::compute(mesh);
	const double h = meanEdgeLength(mesh);
	const Corners rho = cornerConnections(mesh, G, false);
	Eigen::VectorXd r;
	const std::vector<Vector2d> uv = uprightLogMap(mesh, G, rho, h, source, r);
	const Camera cam(elev, azim);
	Image img(2 * 900, 2 * 900);
	const Texture tex = decal(radius, radius / 5, Vector3d(0.55, 0.78, 0.5), radius);
	render(mesh, uv, tex, cam, img, 0, img.w, true, 0.35, true);
	writePPM(img, 2, dir + "/showcase_cactus.ppm");
	Image zoom(2 * 900, 2 * 900);   // a close-up of the decal
	render(mesh, uv, tex, cam, zoom, 0, zoom.w, true, 0.35, true, [&](int i) { return r(i) < 1.6 * radius; });
	writePPM(zoom, 2, dir + "/showcase_cactus_zoom.ppm");
	printf("cactus: %d vertices, source %d at (%.3f, %.3f, %.3f), mean edge %.4f, largest distance %.3f\n",
	       mesh.numVertices(), source, mesh.vertexAt(source).p.x(), mesh.vertexAt(source).p.y(), mesh.vertexAt(source).p.z(), h,
	       r.maxCoeff());
}

static void miuraDemo(const std::string& dir) {
	const std::string obj = dir + "/showcase_miura.obj";
	const int nx = 24, ny = 16;
	writeMiura(obj, nx, ny, 1.0 / 12, 0.6 / 12, 0.6 / 12, 1.0 / 12);
	Mesh mesh;
	mesh.readOBJ(obj.c_str());
	std::remove(obj.c_str());
	mesh.computeNormals();
	const intrinsic::Geometry G = intrinsic::compute(mesh);
	const double h = meanEdgeLength(mesh);
	const int source = (ny / 2) * (nx + 1) + nx / 2;
	std::vector<double> rotation;
	const std::vector<Vector2d> p = unfold(mesh, G, rotation);
	Eigen::VectorXd r(mesh.numVertices());
	for (int i = 0; i < mesh.numVertices(); ++i) r(i) = (p[i] - p[source]).norm();
	// the exact direction in which r grows at the corner (f, a), in the layout frame of f
	auto radial = [&](const Face& f, int a) {
		Vertex* v[3];
		Edge* hh[3];
		intrinsic::faceCorners(f, v, hh);
		const Vector2d d = Eigen::Rotation2Dd(-rotation[f.ID]) * (p[v[a]->ID] - p[source]);
		return Complex(d.x(), d.y());
	};
	const Camera cam(60, -62);
	Image img(2 * 1600, 2 * 760);
	const char* name[2] = { "extrinsic (Part II)", "intrinsic (Part III)" };
	for (int k = 0; k < 2; ++k) {
		const Corners rho = cornerConnections(mesh, G, k == 0);
		const Eigen::VectorXcd Y = vectorHeat(mesh, G, rho, h * h, source);
		const std::vector<Vector2d> uv = logMap(mesh, G, rho, Y, r, radial);
		// the exact map is a rotation of the unfolding: uv = e^{-T alpha} (p - p_s); find alpha from a neighbour of the
		// source and measure the angle of every vertex against it
		Complex ref = 0;
		for (int i = 0; i < mesh.numVertices(); ++i)
			if (i != source && r(i) < 3 * h) {
				const Vector2d d = p[i] - p[source];
				ref += Complex(d.x(), d.y()) * std::conj(Complex(uv[i].x(), uv[i].y()));
			}
		const Complex rot = ref / std::abs(ref);
		double mean = 0, worst = 0;
		int count = 0;
		for (int i = 0; i < mesh.numVertices(); ++i) {
			if (i == source) continue;
			const Vector2d d = p[i] - p[source];
			const double err = std::abs(std::arg(Complex(d.x(), d.y()) * std::conj(rot * Complex(uv[i].x(), uv[i].y())))) * 180 / M_PI;
			mean += err; worst = std::max(worst, err); ++count;
		}
		printf("miura, %s: angle of the log map against the unfolded paper: mean %.3g, max %.3g degrees\n", name[k],
		       mean / count, worst);
		render(mesh, uv, gridPaper(1.0 / 16, 0.35), cam, img, k * img.w / 2, img.w / 2, false, 0.8);
	}
	writePPM(img, 2, dir + "/showcase_miura.ppm");
}

/// The picture below the abstract: the cactus with three decals (a checkerboard with the two arrows, and two polar
/// grids whose rings are geodesic circles), next to the folded Miura sheet carrying straight graph paper.
static void teaser(const std::string& dir) {
	Mesh mesh;
	mesh.readOBJ("cactus2.obj");
	mesh.CenterAndNormalize();
	mesh.computeNormals();
	const intrinsic::Geometry G = intrinsic::compute(mesh);
	const double h = meanEdgeLength(mesh);
	const Corners rho = cornerConnections(mesh, G, false);
	const Camera cam(12, -64);
	// a vertex facing the camera, nearest to the point (fx, fy) of the projected bounding box (fractions, y upwards)
	double minx = 1e300, maxx = -1e300, miny = 1e300, maxy = -1e300;
	for (Vertex& v : mesh.getVertices()) {
		minx = std::min(minx, v.p.dot(cam.right)); maxx = std::max(maxx, v.p.dot(cam.right));
		miny = std::min(miny, v.p.dot(cam.up)); maxy = std::max(maxy, v.p.dot(cam.up));
	}
	auto facing = [&](double fx, double fy) {
		const Vector2d target(minx + fx * (maxx - minx), miny + fy * (maxy - miny));
		int best = 0;
		double d = 1e300;
		for (Vertex& v : mesh.getVertices()) {
			if (v.n.dot(cam.toward) < 0.75) continue;
			const double e = (Vector2d(v.p.dot(cam.right), v.p.dot(cam.up)) - target).norm();
			if (e < d) { d = e; best = v.ID; }
		}
		return best;
	};
	struct Decal { int source; double radius; Texture tex; std::vector<Vector2d> uv; };
	const Vector3d green(0.55, 0.78, 0.5);
	std::vector<Decal> decals = {
		{ facing(0.52, 0.27), 0.15, decal(0.15, 0.03, green, 0.15), {} },
		{ facing(0.52, 0.80), 0.12, polar(0.12, 0.024, Vector3d(0.80, 0.88, 1.0), Vector3d(0.62, 0.76, 0.98),
		                                  Vector3d(0.08, 0.2, 0.5), green), {} },
		{ facing(0.13, 0.66), 0.07, polar(0.07, 0.0175, Vector3d(1.0, 0.86, 0.9), Vector3d(0.96, 0.66, 0.76),
		                                  Vector3d(0.5, 0.08, 0.25), green), {} },
	};
	for (Decal& d : decals) {
		Eigen::VectorXd r;
		d.uv = uprightLogMap(mesh, G, rho, h, d.source, r);
	}
	Image img(2 * 1800, 2 * 1000);
	background(img, 0, img.w, Vector3d(0.97, 0.98, 1.0), Vector3d(0.84, 0.88, 0.94));
	raster(mesh, [&](const Face&, Vertex* const v[3], const double l[3]) -> Vector3d {
		for (const Decal& d : decals) {
			const Vector2d t = texCoord(d.uv, v, l, true);
			if (t.norm() <= d.radius) return d.tex(t);
		}
		return green;
	}, cam, img, 0, 2 * 760, true, 0.35);
	// the Miura sheet with the intrinsic connection: graph paper, straight across every fold
	const std::string obj = dir + "/teaser_miura.obj";
	const int nx = 24, ny = 16;
	writeMiura(obj, nx, ny, 1.0 / 12, 0.6 / 12, 0.6 / 12, 1.0 / 12);
	Mesh miura;
	miura.readOBJ(obj.c_str());
	std::remove(obj.c_str());
	miura.computeNormals();
	const intrinsic::Geometry Gm = intrinsic::compute(miura);
	const int source = (ny / 2) * (nx + 1) + nx / 2;
	std::vector<double> rotation;
	const std::vector<Vector2d> p = unfold(miura, Gm, rotation);
	Eigen::VectorXd r(miura.numVertices());
	for (int i = 0; i < miura.numVertices(); ++i) r(i) = (p[i] - p[source]).norm();
	const Corners rm = cornerConnections(miura, Gm, false);
	const Eigen::VectorXcd Y = vectorHeat(miura, Gm, rm, meanEdgeLength(miura) * meanEdgeLength(miura), source);
	const std::vector<Vector2d> uv = logMap(miura, Gm, rm, Y, r, [&](const Face& f, int a) {
		Vertex* v[3];
		Edge* hh[3];
		intrinsic::faceCorners(f, v, hh);
		const Vector2d d = Eigen::Rotation2Dd(-rotation[f.ID]) * (p[v[a]->ID] - p[source]);
		return Complex(d.x(), d.y());
	});
	render(miura, uv, gridPaper(1.0 / 16, 0.35), Camera(52, -62), img, 2 * 760, 2 * 1040, false, 0.78);
	writePPM(img, 2, dir + "/teaser.ppm");
	printf("teaser: decal sources %d, %d, %d\n", decals[0].source, decals[1].source, decals[2].source);
}

int main(int argc, char** argv) {
	const std::string dir = argc > 1 ? argv[1] : ".";
	const int source = argc > 2 ? std::atoi(argv[2]) : -1;
	const double radius = argc > 3 ? std::atof(argv[3]) : 0.17;
	const double elev = argc > 4 ? std::atof(argv[4]) : 15, azim = argc > 5 ? std::atof(argv[5]) : -70;
	const double height = argc > 6 ? std::atof(argv[6]) : -0.35;   // the height of the default source on the trunk
	if (source >= 0) cactusDemo(dir, source, radius, elev, azim);
	else {
		// default: the vertex facing the camera most directly, half way up the trunk
		Mesh mesh;
		mesh.readOBJ("cactus2.obj");
		mesh.CenterAndNormalize();
		const Camera cam(elev, azim);
		// the axis of the trunk at that height: the centroid of the vertices in a thin slab
		Vector3d axis = Vector3d::Zero();
		int count = 0;
		for (Vertex& v : mesh.getVertices())
			if (std::abs(v.p.z() - height) < 0.04) { axis += v.p; ++count; }
		axis /= std::max(count, 1);
		const Vector3d facing = Vector3d(cam.toward.x(), cam.toward.y(), 0).normalized();
		int best = 0;
		double score = -1e300;
		for (Vertex& v : mesh.getVertices()) {
			if (std::abs(v.p.z() - height) > 0.04) continue;
			const double s = (v.p - axis).dot(facing);
			if (s > score) { score = s; best = v.ID; }
		}
		cactusDemo(dir, best, radius, elev, azim);
	}
	miuraDemo(dir);
	teaser(dir);
	return 0;
}
