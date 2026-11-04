#pragma once
/*
	PathTracer.h: a small physically based renderer for the pictures of showcase_render.cpp, with no dependency beyond
	Eigen and the standard library.

	  - triangles in a bounding volume hierarchy, with smooth or flat normals and a colour callback per object;
	  - one rectangular area light (sampled directly) and a sky that lights the scene from every direction;
	  - a diffuse base with a glossy clear coat whose strength follows Schlick's Fresnel term;
	  - a "shadow catcher" floor: seen directly it shows the backdrop, darkened only by the shadows and the occlusion
	    of the objects, so the objects stand on a seamless studio background; seen indirectly it is a real surface;
	  - antialiasing by jittered samples, a filmic tone curve and sRGB output, on all cores.
*/
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <string>
#include <thread>
#include <vector>
#include <Eigen/Dense>

namespace pt {

using V3 = Eigen::Vector3d;

inline double srgbToLinear(double c) { return std::pow(std::max(0.0, c), 2.2); }
inline V3 srgbToLinear(const V3& c) { return V3(srgbToLinear(c.x()), srgbToLinear(c.y()), srgbToLinear(c.z())); }

/// One object: its colour at a point of a face (given by the face's user index and barycentric coordinates, in sRGB),
/// and its clear coat.
struct Material {
	std::function<V3(int face, const double l[3])> albedo;
	double f0 = 0.04;          // reflectance of the coat at normal incidence
	double roughness = 0.25;   // 0 is a mirror
	bool smooth = true;        // interpolate the vertex normals
};

struct Triangle {
	V3 p[3], n[3];
	int material, face;
};

/// The pseudo-random numbers of one thread (xorshift64*).
struct Rng {
	uint64_t s;
	explicit Rng(uint64_t seed) : s(seed * 0x9E3779B97F4A7C15ull + 1) {}
	double operator()() {
		s ^= s >> 12; s ^= s << 25; s ^= s >> 27;
		return ((s * 0x2545F4914F6CDD1Dull) >> 11) * (1.0 / 9007199254740992.0);
	}
};

/// An orthonormal basis around a unit vector n.
inline void basis(const V3& n, V3& t, V3& b) {
	t = (std::abs(n.x()) > 0.5 ? V3(0, 1, 0) : V3(1, 0, 0)).cross(n).normalized();
	b = n.cross(t);
}

inline V3 cosineSample(const V3& n, Rng& rng) {
	const double u = rng(), v = rng(), r = std::sqrt(u), phi = 2 * M_PI * v;
	V3 t, b;
	basis(n, t, b);
	return (r * std::cos(phi) * t + r * std::sin(phi) * b + std::sqrt(std::max(0.0, 1 - u)) * n).normalized();
}

struct Scene {
	std::vector<Triangle> tris;
	std::vector<Material> materials;
	double floorZ = -1e300;                   // the shadow-catcher floor, the plane z = floorZ
	V3 floorAlbedo = V3(0.8, 0.8, 0.8);       // what the objects see of it (linear)
	V3 lightCorner, lightU, lightV;            // the area light: corner + s U + t V, emitting towards -(U x V)
	V3 lightRadiance = V3(6, 6, 6);
	std::function<V3(const V3&)> sky = [](const V3& d) -> V3 { return V3(0.6, 0.65, 0.75) * (0.6 + 0.4 * std::max(0.0, d.z())); };

	// ------------------------------------------------------------------ the bounding volume hierarchy
	struct Node { V3 lo, hi; int left = -1, right = -1, first = 0, count = 0; };
	std::vector<Node> nodes;
	std::vector<int> order;

	void build() {
		order.resize(tris.size());
		for (size_t i = 0; i < tris.size(); ++i) order[i] = (int)i;
		nodes.clear();
		nodes.reserve(2 * tris.size());
		buildNode(0, (int)tris.size());
	}

	int buildNode(int first, int count) {
		const int id = (int)nodes.size();
		nodes.emplace_back();
		V3 lo = V3::Constant(1e300), hi = V3::Constant(-1e300), clo = lo, chi = hi;
		for (int k = first; k < first + count; ++k)
			for (const V3& p : tris[order[k]].p) { lo = lo.cwiseMin(p); hi = hi.cwiseMax(p); }
		for (int k = first; k < first + count; ++k) {
			const V3 c = centroid(order[k]);
			clo = clo.cwiseMin(c); chi = chi.cwiseMax(c);
		}
		nodes[id].lo = lo; nodes[id].hi = hi;
		if (count <= 4) { nodes[id].first = first; nodes[id].count = count; return id; }
		int axis;
		(chi - clo).maxCoeff(&axis);
		const int mid = first + count / 2;
		std::nth_element(order.begin() + first, order.begin() + mid, order.begin() + first + count,
		                 [&](int a, int b) { return centroid(a)(axis) < centroid(b)(axis); });
		const int l = buildNode(first, mid - first);
		const int r = buildNode(mid, first + count - mid);
		nodes[id].left = l; nodes[id].right = r;
		return id;
	}

	V3 centroid(int t) const { return (tris[t].p[0] + tris[t].p[1] + tris[t].p[2]) / 3.0; }

	static bool box(const Node& n, const V3& o, const V3& inv, double tmax) {
		double t0 = 0, t1 = tmax;
		for (int a = 0; a < 3; ++a) {
			double ta = (n.lo(a) - o(a)) * inv(a), tb = (n.hi(a) - o(a)) * inv(a);
			if (ta > tb) std::swap(ta, tb);
			t0 = std::max(t0, ta); t1 = std::min(t1, tb);
			if (t0 > t1) return false;
		}
		return true;
	}

	/// The nearest triangle hit along o + t d, t in (eps, tmax): returns its index (or -1), t and barycentrics.
	int intersect(const V3& o, const V3& d, double tmax, double& tHit, double& b1, double& b2) const {
		const V3 inv(1 / d.x(), 1 / d.y(), 1 / d.z());
		int stack[64], top = 0, best = -1;
		stack[top++] = 0;
		tHit = tmax;
		while (top) {
			const Node& n = nodes[stack[--top]];
			if (!box(n, o, inv, tHit)) continue;
			if (n.left < 0) {
				for (int k = n.first; k < n.first + n.count; ++k) {
					const Triangle& T = tris[order[k]];
					const V3 e1 = T.p[1] - T.p[0], e2 = T.p[2] - T.p[0], pv = d.cross(e2);
					const double det = e1.dot(pv);
					if (std::abs(det) < 1e-14) continue;
					const double id = 1 / det;
					const V3 tv = o - T.p[0];
					const double u = tv.dot(pv) * id;
					if (u < 0 || u > 1) continue;
					const V3 qv = tv.cross(e1);
					const double v = d.dot(qv) * id;
					if (v < 0 || u + v > 1) continue;
					const double t = e2.dot(qv) * id;
					if (t > 1e-7 && t < tHit) { tHit = t; b1 = u; b2 = v; best = order[k]; }
				}
			} else {
				stack[top++] = n.left;
				stack[top++] = n.right;
			}
		}
		return best;
	}

	bool occluded(const V3& o, const V3& d, double tmax) const {
		double t, u, v;
		return intersect(o, d, tmax, t, u, v) >= 0;
	}

	/// The distance along o + t d to the area light, or a negative number.
	double hitLight(const V3& o, const V3& d) const {
		const V3 n = lightU.cross(lightV);
		const double den = d.dot(n);
		if (std::abs(den) < 1e-12) return -1;
		const double t = (lightCorner - o).dot(n) / den;
		if (t <= 0) return -1;
		const V3 q = o + t * d - lightCorner;
		const double s = q.dot(lightU) / lightU.squaredNorm(), r = q.dot(lightV) / lightV.squaredNorm();
		return (s >= 0 && s <= 1 && r >= 0 && r <= 1) ? t : -1;
	}

	/// Direct light from the area light arriving at p with normal n, times the cosine at p (no BRDF factor).
	V3 directLight(const V3& p, const V3& n, Rng& rng, bool testOccluders = true) const {
		const V3 q = lightCorner + rng() * lightU + rng() * lightV;
		V3 w = q - p;
		const double dist2 = w.squaredNorm(), dist = std::sqrt(dist2);
		w /= dist;
		const V3 ln = lightU.cross(lightV);
		const double area = ln.norm();
		const double cosL = -w.dot(ln) / area, cosP = w.dot(n);
		if (cosL <= 0 || cosP <= 0) return V3::Zero();
		if (testOccluders && occluded(p + 1e-6 * n, w, dist)) return V3::Zero();
		return lightRadiance * (cosP * cosL * area / dist2);
	}

	/// The radiance arriving along the ray o + t d (t > 0), by path tracing. `primaryFloor` is set when the ray is a
	/// camera ray that meets no object (the floor or nothing): then nothing is returned, and the caller composites the
	/// backdrop, darkened by the ratio shadow / unshadowed of the light the floor receives with and without the objects.
	V3 radiance(V3 o, V3 d, Rng& rng, bool* primaryFloor = nullptr, double* shadow = nullptr, double* unshadowed = nullptr) const {
		V3 L = V3::Zero(), beta = V3::Ones();
		bool specular = true;   // the camera ray, or a glossy bounce: it may see the light directly
		for (int depth = 0; depth < 6; ++depth) {
			double t, b1 = 0, b2 = 0;
			const int k = intersect(o, d, 1e300, t, b1, b2);
			const double tf = d.z() < -1e-12 ? (floorZ - o.z()) / d.z() : -1;
			const bool floorFirst = tf > 1e-7 && (k < 0 || tf < t);
			if (k < 0 && !floorFirst) {
				if (depth == 0 && primaryFloor) { *primaryFloor = true; *shadow = *unshadowed = 1; return V3::Zero(); }
				if (specular && depth > 0) { const double tl = hitLight(o, d); if (tl > 0) return L + beta.cwiseProduct(lightRadiance); }
				return L + beta.cwiseProduct(sky(d));
			}
			if (floorFirst) {
				const V3 p = o + tf * d, n(0, 0, 1);
				if (depth == 0 && primaryFloor) {
					// shadow catcher: the light received with the objects, against the light received without them
					*primaryFloor = true;
					Rng r2 = rng;
					const V3 e = directLight(p, n, rng, true), e0 = directLight(p, n, r2, false);
					const V3 w = cosineSample(n, rng);
					const double sk = sky(w).mean(), vis = occluded(p + 1e-6 * n, w, 1e300) ? 0.0 : 1.0;
					*shadow = e.mean() / M_PI + sk * vis;
					*unshadowed = e0.mean() / M_PI + sk;
					return V3::Zero();
				}
				L += beta.cwiseProduct(floorAlbedo).cwiseProduct(directLight(p, n, rng)) / M_PI;
				beta = beta.cwiseProduct(floorAlbedo);
				o = p + 1e-6 * n;
				d = cosineSample(n, rng);
				specular = false;
			} else {
				const Triangle& T = tris[k];
				const Material& m = materials[T.material];
				const double l[3] = { 1 - b1 - b2, b1, b2 };
				const V3 p = o + t * d;
				V3 ng = (T.p[1] - T.p[0]).cross(T.p[2] - T.p[0]).normalized();
				V3 ns = m.smooth ? V3((l[0] * T.n[0] + l[1] * T.n[1] + l[2] * T.n[2]).normalized()) : ng;
				if (ng.dot(d) > 0) ng = -ng;
				if (ns.dot(ng) < 0) ns = -ns;
				const V3 albedo = srgbToLinear(m.albedo(T.face, l));
				const double c = std::max(0.0, -d.dot(ns)), F = m.f0 + (1 - m.f0) * std::pow(1 - c, 5);
				if (rng() < F) {   // the clear coat: a glossy reflection
					const V3 mirror = d - 2 * d.dot(ns) * ns;
					V3 w;
					int tries = 0;
					do {
						V3 j(rng() - 0.5, rng() - 0.5, rng() - 0.5);
						w = (mirror + 2 * m.roughness * j).normalized();
					} while (w.dot(ng) <= 0 && ++tries < 8);
					if (w.dot(ng) <= 0) break;
					o = p + 1e-6 * ng;
					d = w;
					specular = true;
				} else {           // the diffuse base, with the area light sampled directly
					L += beta.cwiseProduct(albedo).cwiseProduct(directLight(p + 1e-6 * ng, ns, rng)) / M_PI;
					beta = beta.cwiseProduct(albedo);
					o = p + 1e-6 * ng;
					d = cosineSample(ns, rng);
					if (d.dot(ng) <= 0) break;
					specular = false;
				}
			}
			if (depth >= 2) {   // Russian roulette
				const double q = std::min(0.95, beta.maxCoeff());
				if (rng() > q) break;
				beta /= q;
			}
		}
		return L;
	}
};

/// A pinhole camera looking from `eye` at `target`, with vertical field of view `fovDeg`.
struct Camera {
	V3 eye, forward, right, up;
	double tanHalf;
	Camera(const V3& eye_, const V3& target, double fovDeg) : eye(eye_) {
		forward = (target - eye).normalized();
		right = forward.cross(V3(0, 0, 1)).normalized();
		up = right.cross(forward);
		tanHalf = std::tan(0.5 * fovDeg * M_PI / 180);
	}
	V3 ray(double x, double y, int w, int h) const {   // x, y in pixels from the top left
		const double aspect = (double)w / h;
		const double sx = (2 * x / w - 1) * tanHalf * aspect, sy = (1 - 2 * y / h) * tanHalf;
		return (forward + sx * right + sy * up).normalized();
	}
};

/// A filmic tone curve (the fit of the ACES curve by Narkowicz), on linear values.
inline double filmic(double x) {
	x *= 0.6;
	return std::max(0.0, std::min(1.0, (x * (2.51 * x + 0.03)) / (x * (2.43 * x + 0.59) + 0.14)));
}

/// Render to an sRGB image (row-major, values in [0, 1]). `backdrop(x, y)` is the sRGB colour of the studio
/// background at a pixel (x, y in [0, 1]); the floor shows it, shaded by the shadows of the objects.
inline std::vector<V3> render(const Scene& scene, const Camera& cam, int w, int h, int spp,
                              const std::function<V3(double, double)>& backdrop, double exposure = 1.0,
                              double shadowStrength = 0.85) {
	std::vector<V3> out(w * h);
	std::atomic<int> nextRow(0);
	auto work = [&]() {
		for (int y; (y = nextRow++) < h;) {
			Rng rng(1234567 + 7919 * y);
			for (int x = 0; x < w; ++x) {
				V3 sumObj = V3::Zero();
				double sumShadow = 0, sumUnshadowed = 0;
				int floorSamples = 0;
				for (int s = 0; s < spp; ++s) {
					const V3 d = cam.ray(x + rng(), y + rng(), w, h);
					bool primaryFloor = false;
					double sh = 0, un = 0;
					const V3 L = scene.radiance(cam.eye, d, rng, &primaryFloor, &sh, &un);
					if (primaryFloor) { ++floorSamples; sumShadow += sh; sumUnshadowed += un; }
					else {
						const V3 e = exposure * L;
						sumObj += V3(filmic(e.x()), filmic(e.y()), filmic(e.z()));
					}
				}
				const V3 bg = srgbToLinear(backdrop((x + 0.5) / w, (y + 0.5) / h));
				const double ratio = sumUnshadowed > 0 ? sumShadow / sumUnshadowed : 1.0;
				const double shade = 1 - shadowStrength * (1 - std::min(1.0, ratio));
				// pixels that see no object at all are pure backdrop (missed rays look at the backdrop too)
				const V3 c = (sumObj + floorSamples * shade * bg) / spp;
				for (int k = 0; k < 3; ++k) out[y * w + x](k) = std::pow(std::max(0.0, std::min(1.0, c(k))), 1 / 2.2);
			}
		}
	};
	std::vector<std::thread> pool;
	const unsigned n = std::max(1u, std::thread::hardware_concurrency());
	for (unsigned i = 0; i < n; ++i) pool.emplace_back(work);
	for (auto& t : pool) t.join();
	return out;
}

inline void writePPM(const std::vector<V3>& img, int w, int h, const std::string& path) {
	std::ofstream out(path, std::ios::binary);
	out << "P6\n" << w << " " << h << "\n255\n";
	for (const V3& c : img)
		for (int k = 0; k < 3; ++k) out.put((char)(unsigned char)std::lround(255 * std::max(0.0, std::min(1.0, c(k)))));
}

}   // namespace pt
