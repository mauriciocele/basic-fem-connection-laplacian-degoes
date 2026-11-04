/*
	showcase_render.cpp: the pictures of showcase.cpp, path traced (PathTracer.h) instead of rasterized: soft light
	from a large area light and the sky, glossy coats, contact shadows on a seamless studio floor, and the cactus in
	a terracotta pot. The mathematics is that of showcase.cpp, which this file includes unchanged: the logarithmic
	map from the heat method and the vector heat method with the intrinsic connection, and the Miura sheet with the
	extrinsic and intrinsic connections.

	Build, from src/:
	  c++ -std=c++17 -O3 -I<eigen> showcase_render.cpp HalfEdge/Mesh.cpp -o showcase_render
	  ./showcase_render <output directory> [samples per pixel] [scale]
	It writes teaser.ppm, showcase_cactus.ppm, showcase_cactus_zoom.ppm, showcase_miura.ppm and showcase_whole.ppm
	(the checkerboard of one logarithmic map over the whole cactus, front and back), and showcase_atlas.ppm (the same
	with an atlas of eight stitched charts).
*/
#define main showcase_main
#include "showcase.cpp"
#undef main
#include "PathTracer.h"

using pt::V3;

namespace {

/// A rigid motion with uniform scale, x -> s R x + t.
struct Placement {
	Eigen::Matrix3d R = Eigen::Matrix3d::Identity();
	double s = 1;
	V3 t = V3::Zero();
	V3 point(const V3& x) const { return s * (R * x) + t; }
	V3 normal(const V3& n) const { return (R * n).normalized(); }
};

/// Add the faces of a mesh to the scene, each triangle remembering its face, with the corners in the order of
/// intrinsic::faceCorners (so that the barycentric coordinates match those of the textures).
void addMesh(pt::Scene& scene, Mesh& mesh, const Placement& P, int material) {
	for (Face& f : mesh.getFaces()) {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(f, v, h);
		pt::Triangle T;
		for (int a = 0; a < 3; ++a) { T.p[a] = P.point(v[a]->p); T.n[a] = P.normal(v[a]->n); }
		T.material = material;
		T.face = f.ID;
		scene.tris.push_back(T);
	}
}

/// A surface of revolution about the vertical axis through `c`, from a profile of (radius, height) points; each
/// segment of the profile is smooth around the axis and sharp at its ends.
void addLathe(pt::Scene& scene, const V3& c, const std::vector<Eigen::Vector2d>& profile, int material, int segments = 160) {
	for (size_t k = 0; k + 1 < profile.size(); ++k) {
		const Eigen::Vector2d a = profile[k], b = profile[k + 1], d = b - a;
		const Eigen::Vector2d nrm = Eigen::Vector2d(d.y(), -d.x()).normalized();   // outwards for an upward profile
		for (int i = 0; i < segments; ++i) {
			const double t0 = 2 * M_PI * i / segments, t1 = 2 * M_PI * (i + 1) / segments;
			auto P = [&](const Eigen::Vector2d& q, double t) -> V3 { return c + V3(q.x() * std::cos(t), q.x() * std::sin(t), q.y()); };
			auto N = [&](double t) -> V3 { return V3(nrm.x() * std::cos(t), nrm.x() * std::sin(t), nrm.y()).normalized(); };
			pt::Triangle T1{ { P(a, t0), P(b, t0), P(b, t1) }, { N(t0), N(t0), N(t1) }, material, 0 };
			pt::Triangle T2{ { P(a, t0), P(b, t1), P(a, t1) }, { N(t0), N(t1), N(t1) }, material, 0 };
			scene.tris.push_back(T1);
			scene.tris.push_back(T2);
		}
	}
}

/// A terracotta pot whose soil surface is at height `soil`, centred under the trunk.
void addPot(pt::Scene& scene, const V3& c, double bottom, double soil, double rTop, int clay, int earth) {
	const double rim = 0.07, rBottom = 0.78 * rTop;
	addLathe(scene, c, {
		{ 0.0, bottom }, { rBottom - 0.01, bottom }, { rBottom, bottom + 0.01 },   // the base
		{ rTop, soil + 0.02 - rim },                                                 // the tapered wall
		{ rTop + 0.025, soil + 0.02 - rim }, { rTop + 0.03, soil + 0.025 - rim },    // the rim
		{ rTop + 0.03, soil + 0.025 }, { rTop + 0.02, soil + 0.035 },
		{ rTop - 0.01, soil + 0.035 }, { rTop - 0.02, soil + 0.025 },
		{ rTop - 0.03, soil },                                                       // the inner wall, down to the soil
	}, clay);
	addLathe(scene, c, { { rTop - 0.03, soil }, { 0.0, soil + 0.004 } }, earth);
}

/// The studio backdrop: a soft radial gradient, lighter behind the subject.
pt::V3 backdrop(double x, double y) {
	const double r = std::hypot((x - 0.42) * 1.2, (y - 0.38));
	const double t = std::min(1.0, r / 0.9);
	const V3 inner(0.985, 0.98, 0.97), outer(0.80, 0.84, 0.90);
	return inner + (outer - inner) * (t * t * (3 - 2 * t));
}

/// A soft box: an area light of the given size centred at `center`, facing `target`.
void softbox(pt::Scene& scene, const V3& center, const V3& target, double size, double radiance) {
	const V3 dir = (center - target).normalized();
	const V3 u = dir.cross(V3(0, 0, 1)).normalized();
	const V3 v = dir.cross(u);
	scene.lightU = size * u;
	scene.lightV = size * v;
	scene.lightCorner = center - 0.5 * scene.lightU - 0.5 * scene.lightV;
	scene.lightRadiance = V3(1.0, 0.97, 0.92) * (1.5 * radiance);
}

void setSky(pt::Scene& scene, double strength) {
	scene.sky = [strength](const V3& d) {
		const double z = std::max(-1.0, std::min(1.0, d.z()));
		const V3 zenith(0.62, 0.72, 0.92), horizon(0.95, 0.93, 0.90), ground(0.55, 0.52, 0.50);
		const V3 c = z > 0 ? V3(horizon + (zenith - horizon) * std::sqrt(z)) : V3(horizon + (ground - horizon) * std::sqrt(-z));
		return V3(strength * c);
	};
}

void save(const std::vector<V3>& img, int w, int h, const std::string& path) {
	pt::writePPM(img, w, h, path);
	printf("wrote %s (%d x %d)\n", path.c_str(), w, h);
}

/// The cactus of showcase.cpp, normalized, with its intrinsic geometry and connection.
struct Cactus {
	Mesh mesh;
	intrinsic::Geometry G;
	Corners rho;
	double h;
	Cactus() {
		mesh.readOBJ("cactus2.obj");
		mesh.CenterAndNormalize();
		mesh.computeNormals();
		G = intrinsic::compute(mesh);
		rho = cornerConnections(mesh, G, false);
		h = meanEdgeLength(mesh);
	}
};

/// One decal: its source, radius, texture and texture coordinates (the upright logarithmic map).
struct Decal {
	int source;
	double radius;
	Texture tex;
	std::vector<Vector2d> uv;
	Eigen::VectorXd r;
};

const V3 kGreen(0.47, 0.72, 0.44);
const V3 kClay(0.80, 0.45, 0.30), kEarth(0.30, 0.22, 0.16);
const double kExposure = 1.9;

/// The cactus, its decals and its pot, added to the scene. The decal textures are sRGB colours.
void addCactus(pt::Scene& scene, Cactus& C, const std::vector<Decal>& decals, const Placement& P, bool atlas = false,
               double seam = 0.0) {
	const Cactus* cp = &C;
	const std::vector<Decal>* dp = &decals;
	pt::Material skin;
	skin.albedo = [cp, dp, atlas, seam](int face, const double l[3]) -> V3 {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(const_cast<Mesh&>(cp->mesh).faceAt(face), v, h);
		if (atlas) {   // an atlas of charts: each point takes the chart of its nearest source (geodesically)
			int best = -1;
			double r1 = 1e300, r2 = 1e300;
			for (size_t k = 0; k < dp->size(); ++k) {
				const Eigen::VectorXd& r = (*dp)[k].r;
				const double rk = l[0] * r(v[0]->ID) + l[1] * r(v[1]->ID) + l[2] * r(v[2]->ID);
				if (rk < r1) { r2 = r1; r1 = rk; best = (int)k; } else if (rk < r2) r2 = rk;
			}
			if (r2 - r1 < seam) return V3(0.30, 0.45, 0.30);   // a thin stitch where two charts meet
			const Decal& d = (*dp)[best];
			return d.tex(texCoord(d.uv, v, l, false));
		}
		for (const Decal& d : *dp) {
			const Vector2d t = texCoord(d.uv, v, l, true);
			if (t.norm() <= d.radius) return d.tex(t);
		}
		return kGreen;
	};
	skin.f0 = 0.05;
	skin.roughness = 0.22;
	scene.materials.push_back(skin);
	addMesh(scene, C.mesh, P, (int)scene.materials.size() - 1);
	pt::Material clay;
	clay.albedo = [](int, const double*) { return kClay; };
	clay.f0 = 0.03;
	clay.roughness = 0.5;
	scene.materials.push_back(clay);
	pt::Material earth;
	earth.albedo = [](int, const double*) { return kEarth; };
	earth.f0 = 0.0;
	scene.materials.push_back(earth);
	const int nm = (int)scene.materials.size();
	// the pot is placed for the cactus as it is loaded (P is the identity in all the scenes)
	addPot(scene, V3(-0.004, 0.04, 0), -1.10, -0.74, 0.27, nm - 2, nm - 1);
	scene.floorZ = -1.10;
}

/// The Miura sheet with the logarithmic map from its centre: with the extrinsic connection (Part II) or the
/// intrinsic one (Part III), and the exact distances of the unfolding.
struct Miura {
	Mesh mesh;
	std::vector<Vector2d> uv;
	explicit Miura(const std::string& dir, bool extrinsic) {
		const std::string obj = dir + "/render_miura.obj";
		const int nx = 24, ny = 16;
		writeMiura(obj, nx, ny, 1.0 / 12, 0.6 / 12, 0.6 / 12, 1.0 / 12);
		mesh.readOBJ(obj.c_str());
		std::remove(obj.c_str());
		mesh.computeNormals();
		const intrinsic::Geometry G = intrinsic::compute(mesh);
		const int source = (ny / 2) * (nx + 1) + nx / 2;
		std::vector<double> rotation;
		const std::vector<Vector2d> p = unfold(mesh, G, rotation);
		Eigen::VectorXd r(mesh.numVertices());
		for (int i = 0; i < mesh.numVertices(); ++i) r(i) = (p[i] - p[source]).norm();
		const Corners rho = cornerConnections(mesh, G, extrinsic);
		const double h = meanEdgeLength(mesh);
		const Eigen::VectorXcd Y = vectorHeat(mesh, G, rho, h * h, source);
		uv = logMap(mesh, G, rho, Y, r, [&](const Face& f, int a) {
			Vertex* v[3];
			Edge* hh[3];
			intrinsic::faceCorners(f, v, hh);
			const Vector2d d = Eigen::Rotation2Dd(-rotation[f.ID]) * (p[v[a]->ID] - p[source]);
			return Complex(d.x(), d.y());
		});
		// centre the sheet at the origin, in its own plane
		Vector3d c = Vector3d::Zero();
		for (Vertex& v : mesh.getVertices()) c += v.p;
		c /= mesh.numVertices();
		for (Vertex& v : mesh.getVertices()) v.p -= c;
	}
};

void addMiura(pt::Scene& scene, Miura& M, const Placement& P) {
	const Miura* mp = &M;
	const Texture paper = gridPaper(1.0 / 16, 0.35);
	pt::Material m;
	m.albedo = [mp, paper](int face, const double l[3]) -> V3 {
		Vertex* v[3];
		Edge* h[3];
		intrinsic::faceCorners(const_cast<Mesh&>(mp->mesh).faceAt(face), v, h);
		return paper(texCoord(mp->uv, v, l, false));
	};
	m.f0 = 0.02;
	m.roughness = 0.6;
	m.smooth = false;
	scene.materials.push_back(m);
	addMesh(scene, M.mesh, P, (int)scene.materials.size() - 1);
}

/// The rotation whose columns are the images of the sheet's axes: x along `along`, z along `normal`.
Eigen::Matrix3d frameFrom(const V3& along, const V3& normal) {
	const V3 n = normal.normalized();
	const V3 x = (along - along.dot(n) * n).normalized();
	Eigen::Matrix3d R;
	R.col(0) = x;
	R.col(1) = n.cross(x);
	R.col(2) = n;
	return R;
}

/// The rotation that shows an object to `cam` as showcase.cpp's camera `old` showed it: the old camera's right, up
/// and towards-the-eye axes go to those of the new one.
Eigen::Matrix3d seenAs(const ::Camera& old, const pt::Camera& cam) {
	Eigen::Matrix3d A, B;
	A << old.right, old.up, old.toward;
	B << cam.right, cam.up, -cam.forward;
	return B * A.transpose();
}

V3 viewDirection(double elevDeg, double azimDeg) {   // from the scene towards the eye, as in showcase.cpp
	const double e = elevDeg * M_PI / 180, a = azimDeg * M_PI / 180;
	return V3(std::cos(e) * std::cos(a), std::cos(e) * std::sin(a), std::sin(e));
}

std::vector<Decal> teaserDecals(Cactus& C) {
	const Vector3d green = kGreen;
	// the sources chosen by teaser() of showcase.cpp: vertices facing its camera at fractions of the picture
	std::vector<Decal> decals = {
		{ 1369, 0.15, decal(0.15, 0.03, green, 0.15), {}, {} },
		{ 979, 0.12, polar(0.12, 0.024, Vector3d(0.80, 0.88, 1.0), Vector3d(0.62, 0.76, 0.98), Vector3d(0.08, 0.2, 0.5), green), {}, {} },
		{ 2133, 0.07, polar(0.07, 0.0175, Vector3d(1.0, 0.86, 0.9), Vector3d(0.96, 0.66, 0.76), Vector3d(0.5, 0.08, 0.25), green), {}, {} },
	};
	for (Decal& d : decals) d.uv = uprightLogMap(C.mesh, C.G, C.rho, C.h, d.source, d.r);
	return decals;
}

void teaserScene(const std::string& dir, int spp, double scale) {
	Cactus C;
	const std::vector<Decal> decals = teaserDecals(C);
	Miura M(dir, false);
	pt::Scene scene;
	const V3 view = viewDirection(12, -64);
	const V3 right = V3(0, 0, 1).cross(view).normalized();   // the camera's right, horizontal
	addCactus(scene, C, decals, Placement());
	const V3 target = 0.62 * right + V3(0, 0, -0.12);
	const pt::Camera cam(target + 6.2 * (view + V3(0, 0, 0.05)).normalized(), target, 21);
	Placement sheet;   // the sheet as the old teaser showed it, from 52 degrees above its plane
	sheet.s = 0.70;
	sheet.R = seenAs(::Camera(52, -62), cam);
	sheet.t = 1.32 * right + V3(0, 0, -0.10) - 0.1 * view;
	addMiura(scene, M, sheet);
	scene.build();
	setSky(scene, 0.55);
	scene.floorAlbedo = V3(0.75, 0.75, 0.76);
	softbox(scene, target + 4.0 * (0.3 * view - 0.85 * right + V3(0, 0, 1.0)).normalized(), target, 2.6, 2.8);
	const int w = (int)(1800 * scale), h = (int)(1000 * scale);
	save(pt::render(scene, cam, w, h, spp, backdrop, kExposure), w, h, dir + "/teaser.ppm");
}

void cactusScene(const std::string& dir, int spp, double scale) {
	Cactus C;
	std::vector<Decal> decals = { { 1369, 0.17, decal(0.17, 0.17 / 5, kGreen, 0.17), {}, {} } };
	decals[0].uv = uprightLogMap(C.mesh, C.G, C.rho, C.h, decals[0].source, decals[0].r);
	pt::Scene scene;
	addCactus(scene, C, decals, Placement());
	scene.build();
	setSky(scene, 0.55);
	const V3 view = viewDirection(15, -70);
	const V3 right = V3(0, 0, 1).cross(view).normalized();
	const V3 target(0, 0, -0.12);
	softbox(scene, target + 4.0 * (0.5 * view - 0.7 * right + V3(0, 0, 0.9)).normalized(), target, 2.4, 2.6);
	const int w = (int)(900 * scale), h = (int)(900 * scale);
	save(pt::render(scene, pt::Camera(target + 6.5 * view, target, 20), w, h, spp, backdrop, kExposure), w, h, dir + "/showcase_cactus.ppm");
	// the close-up: looking at the source along its normal, a little from above
	const V3 s = C.mesh.vertexAt(1369).p, n = C.mesh.vertexAt(1369).n;
	const V3 eye = s + 1.6 * (n + V3(0, 0, 0.25)).normalized();
	save(pt::render(scene, pt::Camera(eye, s, 17), w, h, spp, backdrop, kExposure), w, h, dir + "/showcase_cactus_zoom.ppm");
}

void miuraScene(const std::string& dir, int spp, double scale) {
	Miura E(dir, true), I(dir, false);
	pt::Scene scene;
	const V3 view = viewDirection(30, -90);
	const V3 right = V3(0, 0, 1).cross(view).normalized();
	const V3 target(0, 0, 0.45);
	const pt::Camera cam(target + 8.0 * view, target, 20);
	for (int k = 0; k < 2; ++k) {   // each sheet as the rasterized figure showed it, from 60 degrees above its plane
		Placement P;
		P.R = seenAs(::Camera(60, -62), cam);
		P.t = target + (k == 0 ? -1.22 : 1.22) * right;
		addMiura(scene, k == 0 ? E : I, P);
	}
	scene.build();
	scene.floorZ = -0.35;
	setSky(scene, 0.55);
	scene.floorAlbedo = V3(0.75, 0.75, 0.76);
	softbox(scene, target + 5.0 * (0.6 * view - 0.5 * right + V3(0, 0, 1.0)).normalized(), target, 3.0, 2.8);
	const int w = (int)(1600 * scale), h = (int)(760 * scale);
	save(pt::render(scene, cam, w, h, spp, backdrop, kExposure), w, h, dir + "/showcase_miura.ppm");
}

/// The checkerboard of the logarithmic map from one source, over the whole cactus, seen from the front (the source)
/// and from the back (the cut locus, where the faces stretched past the guard of texCoord stay green). The two
/// views are written side by side into one picture.
void wholeScene(const std::string& dir, int spp, double scale, double cell) {
	Cactus C;
	const double azimuth = -25;   // the front view: the arms spread out, and the source faces the camera
	int source = 0;
	double best = -1e300;
	for (Vertex& v : C.mesh.getVertices()) {   // on the trunk, a third of the way up, facing the front view
		if (std::abs(v.p.z() + 0.30) > 0.03) continue;
		const double score = v.n.dot(viewDirection(0, azimuth));
		if (score > best) { best = score; source = v.ID; }
	}
	printf("whole: source %d\n", source);
	std::vector<Decal> decals = { { source, 1e8, decal(1e8, cell, kGreen, 0.17), {}, {} } };
	decals[0].uv = uprightLogMap(C.mesh, C.G, C.rho, C.h, decals[0].source, decals[0].r);
	pt::Scene scene;
	addCactus(scene, C, decals, Placement());
	scene.build();
	setSky(scene, 0.55);
	const int w = (int)(800 * scale), h = (int)(1000 * scale);
	std::vector<V3> both(2 * w * h);
	for (int k = 0; k < 2; ++k) {
		const V3 view = viewDirection(14, azimuth + 180 * k);
		const V3 right = V3(0, 0, 1).cross(view).normalized();
		const V3 target(0, 0, -0.14);
		softbox(scene, target + 4.0 * (0.5 * view - 0.7 * right + V3(0, 0, 0.9)).normalized(), target, 2.4, 2.6);
		const std::vector<V3> img = pt::render(scene, pt::Camera(target + 6.5 * view, target, 22), w, h, spp,
		                                       [k](double x, double y) { return backdrop(0.5 * (x + k), y); }, kExposure);
		for (int y = 0; y < h; ++y)
			for (int x = 0; x < w; ++x) both[y * 2 * w + k * w + x] = img[y * w + x];
	}
	save(both, 2 * w, h, dir + "/showcase_whole.ppm");
}

/// The checkerboard over the whole cactus with an atlas of charts: each chart is the logarithmic map from its own
/// source, and each point uses the chart whose source is geodesically nearest. With a source on the front and one
/// on the back of each cylinder-like part, every chart stops at the sides of its part, well before its cut locus.
void atlasScene(const std::string& dir, int spp, double scale, double cell, int charts, double seam) {
	Cactus C;
	const double azimuth = -25;
	const V3 front = viewDirection(0, azimuth);
	// the nearest vertex to a point on an axis of the cactus, on the side facing `dir`
	auto onSide = [&](const V3& axis, const V3& dir) {
		int best = 0;
		double score = -1e300;
		for (Vertex& v : C.mesh.getVertices()) {
			if ((v.p - axis).norm() > 0.2 || v.n.dot(dir) < 0.5) continue;
			const double s = (v.p - axis).dot(dir) - 2 * (v.p - axis - (v.p - axis).dot(dir) * dir).norm();
			if (s > score) { score = s; best = v.ID; }
		}
		return best;
	};
	std::vector<V3> axes = { V3(-0.005, 0.03, -0.35) };   // the trunk; then its top, and the two arms
	if (charts > 2) { axes.push_back(V3(-0.005, 0.0, 0.50)); axes.push_back(V3(0, -0.45, 0.45)); axes.push_back(V3(0, 0.45, 0.55)); }
	std::vector<Decal> decals;
	for (const V3& a : axes)
		for (int side = 0; side < 2; ++side)
			decals.push_back({ onSide(a, side == 0 ? front : V3(-front)), 1e8,
			                   decal(1e8, cell, kGreen, decals.empty() ? 0.17 : 0.0), {}, {} });   // arrows on the first
	for (Decal& d : decals) {
		d.uv = uprightLogMap(C.mesh, C.G, C.rho, C.h, d.source, d.r);
		printf("atlas: chart at vertex %d\n", d.source);
	}
	// stitch the charts: every chart is already turned upright (its red arrow along e_z), so only the phase of its
	// checkerboard is free; translate each chart to match the charts placed before it, on average over the vertices
	// near their common boundary (where its distance is within a band of the nearest placed one)
	const int n = C.mesh.numVertices();
	for (size_t k = 1; k < decals.size(); ++k) {
		std::vector<Vector2d> a, b;
		for (int i = 0; i < n; ++i) {
			int owner = -1;
			double rOwner = 1e300;
			for (size_t j = 0; j < k; ++j) if (decals[j].r(i) < rOwner) { rOwner = decals[j].r(i); owner = (int)j; }
			bool nearestOfAll = true;   // the vertex must be nearer to chart k or to its owner than to any later chart
			for (size_t j = k + 1; j < decals.size(); ++j) if (decals[j].r(i) < std::min(rOwner, decals[k].r(i))) nearestOfAll = false;
			if (nearestOfAll && std::abs(decals[k].r(i) - rOwner) < 3 * C.h) { a.push_back(decals[k].uv[i]); b.push_back(decals[owner].uv[i]); }
		}
		if (a.empty()) continue;
		Vector2d shift = Vector2d::Zero();
		for (size_t m = 0; m < a.size(); ++m) shift += b[m] - a[m];
		shift /= a.size();
		for (Vector2d& t : decals[k].uv) t += shift;
		printf("atlas: chart %zu stitched on %zu vertices\n", k, a.size());
	}
	pt::Scene scene;
	addCactus(scene, C, decals, Placement(), true, seam);
	scene.build();
	setSky(scene, 0.55);
	const int w = (int)(800 * scale), h = (int)(1000 * scale);
	std::vector<V3> both(2 * w * h);
	for (int k = 0; k < 2; ++k) {
		const V3 view = viewDirection(14, azimuth + 180 * k);
		const V3 right = V3(0, 0, 1).cross(view).normalized();
		const V3 target(0, 0, -0.14);
		softbox(scene, target + 4.0 * (0.5 * view - 0.7 * right + V3(0, 0, 0.9)).normalized(), target, 2.4, 2.6);
		const std::vector<V3> img = pt::render(scene, pt::Camera(target + 6.5 * view, target, 22), w, h, spp,
		                                       [k](double x, double y) { return backdrop(0.5 * (x + k), y); }, kExposure);
		for (int y = 0; y < h; ++y)
			for (int x = 0; x < w; ++x) both[y * 2 * w + k * w + x] = img[y * w + x];
	}
	save(both, 2 * w, h, dir + (charts == 8 ? std::string("/showcase_atlas.ppm") : "/showcase_atlas" + std::to_string(charts) + ".ppm"));
}

}   // namespace

int main(int argc, char** argv) {
	const std::string dir = argc > 1 ? argv[1] : ".";
	const int spp = argc > 2 ? std::atoi(argv[2]) : 256;
	const double scale = argc > 3 ? std::atof(argv[3]) : 1.0;
	const std::string only = argc > 4 ? argv[4] : "all";
	if (only == "all" || only == "teaser") teaserScene(dir, spp, scale);
	if (only == "all" || only == "cactus") cactusScene(dir, spp, scale);
	if (only == "all" || only == "miura") miuraScene(dir, spp, scale);
	if (only == "all" || only == "atlas" || only == "atlas2")
		atlasScene(dir, spp, scale, 0.045, only == "atlas2" ? 2 : 8, argc > 5 ? std::atof(argv[5]) : 0.0);
	if (only == "all" || only == "whole") wholeScene(dir, spp, scale, argc > 5 ? std::atof(argv[5]) : 0.045);
	return 0;
}
