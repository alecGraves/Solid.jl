// mucad.cpp  ----------------------------------------------------
//    C wrappers for the OCCT-based MUCAD library.
//-----------------------------------------------------------------

#include "mucad.h"   // the c header - it defines Shape, etc.
#include "mucad.hpp" // the c++ header - it defines c++ functions

#include <cstring>
#include <stddef.h> // size_t
#include <string>

//------------------------------------------------------------------
// Wrappers are wrapped in an extern "C" block so that the
// generated symbols are unmangled and can be linked from C-abi
// -----------------------------------------------------------------
extern "C" {

//----------  Memory management ------------------------------------
void mucad_free(Shape s) {
	mucadcpp_free(s);
}

//----------  2-D primitives ---------------------------------------
Shape mucad_circle(double cx, double cy, double radius) {
	return mucadcpp_circle(cx, cy, radius);
}

Shape mucad_rectangle(double x1, double y1, double x2, double y2) {
	return mucadcpp_rectangle(x1, y1, x2, y2);
}

Shape mucad_triangle(double x1, double y1,
		     double x2, double y2,
		     double x3, double y3) {
	return mucadcpp_triangle(x1, y1, x2, y2, x3, y3);
}

Shape mucad_polygon(const double *pts, size_t npts) {
	return mucadcpp_polygon(pts, npts);
}

//----------  3-D primitives ---------------------------------------
Shape mucad_sphere(double cx, double cy, double cz, double radius) {
	return mucadcpp_sphere(cx, cy, cz, radius);
}

Shape mucad_box(double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return mucadcpp_box(x1, y1, z1, x2, y2, z2);
}

//----------  Transformations --------------------------------------
Shape mucad_translate(Shape s, double dx, double dy, double dz) {
	return mucadcpp_translate(s, dx, dy, dz);
}

Shape mucad_rotate(Shape s,
		   double ax, double ay, double az,
		   double ux, double uy, double uz,
		   double angle) {
	return mucadcpp_rotate(s, ax, ay, az, ux, uy, uz, angle);
}

Shape mucad_transform(Shape s, const double mat[16]) {
	return mucadcpp_transform(s, mat);
}

Shape mucad_scale(Shape s, double sx, double sy, double sz) {
	return mucadcpp_scale(s, sx, sy, sz);
}

Shape mucad_mirror(Shape s,
		   double ax, double ay, double az,
		   double ux, double uy, double uz) {
	return mucadcpp_mirror(s, ax, ay, az, ux, uy, uz);
}

//----------  Sweep / Extrude / Revolve ----------------------------
Shape mucad_extrude(Shape shape2d, double ux, double uy, double uz) {
	return mucadcpp_extrude(shape2d, ux, uy, uz);
}

Shape mucad_sweep(Shape profile, Shape path) {
	return mucadcpp_sweep(profile, path);
}

Shape mucad_revolve(Shape s,
		    double ax, double ay, double az,
		    double ux, double uy, double uz,
		    double angle) {
	return mucadcpp_revolve(s, ax, ay, az, ux, uy, uz, angle);
}

//----------  Loft -------------------------------------------------
Shape mucad_loft(const Shape *shapes, size_t n,
		 int solid, int ruled) {
	// C callers pass 0/1 for booleans - cast to bool here.--
	return mucadcpp_loft(shapes, n, solid != 0, ruled != 0);
}

//----------  Boolean operations -----------------------------------
Shape mucad_union(Shape a, Shape b) {
	return mucadcpp_union(a, b);
}

Shape mucad_difference(Shape a, Shape b) {
	return mucadcpp_difference(a, b);
}

//----------  Splines & Advanced -----------------------------------
Shape mucad_quadratic_spline_wire(const double *pts, size_t npts) {
	return mucadcpp_quadratic_spline_wire(pts, npts);
}

Shape mucad_face_from_wire(Shape wire) {
	return mucadcpp_face_from_wire(wire);
}

//----------  Analysis ---------------------------------------------
int mucad_bounding_sphere(Shape s,
			  double *x, double *y, double *z, double *r) {
	return mucadcpp_bounding_sphere(s, x, y, z, r);
}

double mucad_volume(Shape s) {
	return mucadcpp_volume(s);
}

double mucad_surface_area(Shape s) {
	return mucadcpp_surface_area(s);
}

int mucad_centroid(Shape s,
		   double *cx, double *cy, double *cz) {
	return mucadcpp_centroid(s, cx, cy, cz);
}

int mucad_inertia(Shape s, double density,
		  double *Ixx, double *Iyy, double *Izz) {
	return mucadcpp_inertia(s, density, Ixx, Iyy, Izz);
}

int mucad_shape_type(Shape s) {
	return mucadcpp_shape_type(s);
}

//----------  Simplify ----------------------------------------------
Shape mucad_simplify(const Shape s) {
	return mucadcpp_simplify(s);
}

Shape mucad_simplify_to_solid(const Shape s) {
	return mucadcpp_simplify_to_solid(s);
}


//----------  I/O ---------------------------------------------------
static inline std::string make_string(const char *s, int len) {
	if(len < 0) { // negative -> use C-string length
		len = static_cast<int>(std::strlen(s));
	}
	return std::string(s, len);
}

int mucad_write_step(Shape s, const char *filename, int nlen) {
	std::string fname = make_string(filename, nlen);
	return mucadcpp_write_step(s, fname.c_str());
}

int mucad_write_stl(Shape s, const char *filename, int nlen, float resolution) {
	std::string fname = make_string(filename, nlen);
	return mucadcpp_write_stl(s, fname.c_str(), resolution);
}

int mucad_write_iges(Shape s, const char *filename, int nlen) {
	std::string fname = make_string(filename, nlen);
	return mucadcpp_write_iges(s, fname.c_str());
}

int mucad_write_obj(Shape s, const char *filename, int nlen, float resolution) {
	std::string fname = make_string(filename, nlen);
	return mucadcpp_write_obj(s, fname.c_str(), resolution);
}

} //extern "C"--
