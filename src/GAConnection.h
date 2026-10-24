#pragma once
/*
	The connection Laplacian of FEM.cpp in the language of geometric algebra (Part II of fembasic2d.pdf).

	A minimal geometric algebra of R^3: only the grades the method needs.
	  - vectors:   Eigen::Vector3d
	  - bivectors: Bivector, stored by its dual vector b, B = I b (so e23, e31, e12 have coefficients b.x, b.y, b.z)
	  - rotors:    Rotor = w + B, an even element with w^2 + |b|^2 = 1
	  - complex numbers of a plane: std::complex<double>, alpha + beta T for the unit bivector T of that plane

	Every function names the section and equation of the text it implements.
*/
#include <cmath>
#include <complex>
#include <Eigen/Dense>

namespace ga {

using Vec3 = Eigen::Vector3d;
using Complex = std::complex<double>;

/// A bivector B = I b, stored by its dual vector b. Unit bivectors square to -1.
struct Bivector {
	Vec3 b;
	Bivector operator*(double s) const { return { b * s }; }
	double norm() const { return b.norm(); }
};

/// The outer product a ^ c = I (a x c): the oriented plane of a and c.
inline Bivector wedge(const Vec3& a, const Vec3& c) { return { a.cross(c) }; }

/// The inner product of a vector and a bivector, a . B, a vector of the plane of B, orthogonal to a.
/// For B = I b: a . (I b) = I (a ^ b) = I I (a x b) = b x a.
inline Vec3 dot(const Vec3& a, const Bivector& B) { return B.b.cross(a); }

/// B^-1 = -B / |B|^2, since B^2 = -|B|^2.
inline Bivector inverse(const Bivector& B) { return B * (-1.0 / B.b.squaredNorm()); }

/// A rotor R = w + I b (Section 12). R a ~R rotates the vector a.
struct Rotor {
	double w;
	Vec3 b;
	/// R a ~R. With R = w + I b, the rotation is the one of the unit quaternion (w, -b):
	/// a' = a + 2 w (u x a) + 2 u x (u x a), u = -b.
	Vec3 rotate(const Vec3& a) const {
		const Vec3 u = -b;
		const Vec3 t = 2.0 * u.cross(a);
		return a + w * t + u.cross(t);
	}
};

/// The geometric product of two vectors, c a = c . a + c ^ a, a rotor up to normalization.
inline Rotor product(const Vec3& c, const Vec3& a) { return { c.dot(a), c.cross(a) }; }

/// The tangent bivector T = I n of a unit normal n (eq. ga_T): the tangent plane as an oriented plane element.
inline Bivector tangentBivector(const Vec3& n) { return { n }; }

/// a T for a tangent vector a: the quarter turn n x a (eq. ga_quarter_turn).
/// For a in the plane of T the geometric product equals the inner product, a T = a . T.
inline Vec3 quarterTurn(const Vec3& a, const Bivector& T) { return dot(a, T); }

/// The smallest rotor taking the unit vector a to the unit vector c (Proposition ga_rotor):
/// Q = (1 + c a) / |1 + c a|. Undefined for c = -a, where 1 + c a = 0.
inline Rotor rotorBetween(const Vec3& a, const Vec3& c) {
	const Rotor ca = product(c, a);
	const double w = 1.0 + ca.w;
	const double norm = std::sqrt(w * w + ca.b.squaredNorm());   // |1 + c a| = sqrt(2 (1 + a.c))
	return { w / norm, ca.b / norm };
}

/// The complex number of the plane with normal n given by the product of two unit vectors of that plane:
/// t c = t . c + t ^ c = cos(rho) + T sin(rho), with t ^ c = (t x c . n) T (eq. ga_connection).
inline Complex planeProduct(const Vec3& t, const Vec3& c, const Vec3& n) {
	return { t.dot(c), t.cross(c).dot(n) };
}

/// The discrete connection from vertex i to face f, as a unit complex number (eq. ga_connection):
///   e^{T_f rho_i^f} = t_f (Q t_i ~Q),  Q = rotorBetween(n_i, n_f).
/// Transporting the spinor coordinate psi of vertex i into the face is multiplication by it (eq. ga_transport).
inline Complex connection(const Vec3& t_i, const Vec3& n_i, const Vec3& t_f, const Vec3& n_f) {
	const Rotor Q = rotorBetween(n_i, n_f);
	return planeProduct(t_f, Q.rotate(t_i), n_f);
}

/// The gradients of the three hat functions of the triangle x1, x2, x3 (Proposition ga_reciprocal):
/// the reciprocal frame e^1 = e2 . E^-1, e^2 = -e1 . E^-1 of the edges e1 = x2 - x1, e2 = x3 - x1, with
/// E = e1 ^ e2; then grad phi_2 = e^1, grad phi_3 = e^2, grad phi_1 = -e^1 - e^2.
inline void hatGradients(const Vec3& x1, const Vec3& x2, const Vec3& x3, Vec3 grad[3]) {
	const Vec3 e1 = x2 - x1, e2 = x3 - x1;
	const Bivector Einv = inverse(wedge(e1, e2));
	const Vec3 r1 = dot(e2, Einv);
	const Vec3 r2 = -dot(e1, Einv);
	grad[0] = -r1 - r2;
	grad[1] = r1;
	grad[2] = r2;
}

/// The same gradient, one vertex at a time (eq. ga_grad_edge): the opposite edge, from b to c in the
/// counterclockwise order a, b, c, turned a quarter turn in the plane of the face, divided by twice the area.
inline Vec3 hatGradient(const Vec3& xb, const Vec3& xc, const Bivector& Tf, double area) {
	return quarterTurn(xc - xb, Tf) / (2.0 * area);
}

}   // namespace ga
