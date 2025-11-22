// ccad.cpp  ----------------------------------------------------
//    C wrappers for the OCCT-based CCAD library.
//-----------------------------------------------------------------

#include "ccad.h"   // the c header - it defines Shape, etc.
#include "ccad.hpp" // the c++ header - it defines c++ functions

#include <cstring>
#include <stddef.h> // size_t
#include <string>

//------------------------------------------------------------------
// Wrappers are wrapped in an extern "C" block so that the
// generated symbols are unmangled and can be linked from C-abi
// -----------------------------------------------------------------
extern "C" {

//----------  Memory management ------------------------------------
void ccad_free(Shape s) {
	ccadcpp_free(s);
}

//----------  2-D primitives ---------------------------------------
Shape ccad_circle(double cx, double cy, double radius) {
	return ccadcpp_circle(cx, cy, radius);
}

Shape ccad_rectangle(double x1, double y1, double x2, double y2) {
	return ccadcpp_rectangle(x1, y1, x2, y2);
}

Shape ccad_triangle(double x1, double y1,
		     double x2, double y2,
		     double x3, double y3) {
	return ccadcpp_triangle(x1, y1, x2, y2, x3, y3);
}

Shape ccad_polygon(const double *pts, size_t npts) {
	return ccadcpp_polygon(pts, npts);
}

//----------  3-D primitives ---------------------------------------
Shape ccad_sphere(double cx, double cy, double cz, double radius) {
	return ccadcpp_sphere(cx, cy, cz, radius);
}

Shape ccad_box(double x1, double y1, double z1,
		double x2, double y2, double z2) {
	return ccadcpp_box(x1, y1, z1, x2, y2, z2);
}

//----------  Transformations --------------------------------------
Shape ccad_translate(Shape s, double dx, double dy, double dz) {
	return ccadcpp_translate(s, dx, dy, dz);
}

Shape ccad_rotate(Shape s,
		   double ax, double ay, double az,
		   double ux, double uy, double uz,
		   double angle) {
	return ccadcpp_rotate(s, ax, ay, az, ux, uy, uz, angle);
}

Shape ccad_transform(Shape s, const double mat[16]) {
	return ccadcpp_transform(s, mat);
}

Shape ccad_scale(Shape s, double sx, double sy, double sz) {
	return ccadcpp_scale(s, sx, sy, sz);
}

Shape ccad_mirror(Shape s,
		   double ax, double ay, double az,
		   double ux, double uy, double uz) {
	return ccadcpp_mirror(s, ax, ay, az, ux, uy, uz);
}

//----------  Sweep / Extrude / Revolve ----------------------------
Shape ccad_extrude(Shape shape2d, double ux, double uy, double uz) {
	return ccadcpp_extrude(shape2d, ux, uy, uz);
}

Shape ccad_sweep(Shape profile, Shape path) {
	return ccadcpp_sweep(profile, path);
}

Shape ccad_revolve(Shape s,
		    double ax, double ay, double az,
		    double ux, double uy, double uz,
		    double angle) {
	return ccadcpp_revolve(s, ax, ay, az, ux, uy, uz, angle);
}

//----------  Loft -------------------------------------------------
Shape ccad_loft(const Shape *shapes, size_t n,
		 int solid, int ruled) {
	// C callers pass 0/1 for booleans - cast to bool here.--
	return ccadcpp_loft(shapes, n, solid != 0, ruled != 0);
}

//----------  Boolean operations -----------------------------------
Shape ccad_union(Shape a, Shape b) {
	return ccadcpp_union(a, b);
}

Shape ccad_difference(Shape a, Shape b) {
	return ccadcpp_difference(a, b);
}

//----------  Splines & Advanced -----------------------------------
Shape ccad_quadratic_spline_wire(const double *pts, size_t npts) {
	return ccadcpp_quadratic_spline_wire(pts, npts);
}

Shape ccad_face_from_wire(Shape wire) {
	return ccadcpp_face_from_wire(wire);
}

//----------  Analysis ---------------------------------------------
int ccad_bounding_sphere(Shape s,
			  double *x, double *y, double *z, double *r) {
	return ccadcpp_bounding_sphere(s, x, y, z, r);
}

double ccad_volume(Shape s) {
	return ccadcpp_volume(s);
}

double ccad_surface_area(Shape s) {
	return ccadcpp_surface_area(s);
}

int ccad_centroid(Shape s,
		   double *cx, double *cy, double *cz) {
	return ccadcpp_centroid(s, cx, cy, cz);
}

int ccad_inertia(Shape s, double density,
		  double *Ixx, double *Iyy, double *Izz) {
	return ccadcpp_inertia(s, density, Ixx, Iyy, Izz);
}

int ccad_shape_type(Shape s) {
	return ccadcpp_shape_type(s);
}

//----------  Simplify ----------------------------------------------
Shape ccad_simplify(const Shape s) {
	return ccadcpp_simplify(s);
}

Shape ccad_simplify_to_solid(const Shape s) {
	return ccadcpp_simplify_to_solid(s);
}


//----------  I/O ---------------------------------------------------
static inline std::string make_string(const char *s, int len) {
	if(len < 0) { // negative -> use C-string length
		len = static_cast<int>(std::strlen(s));
	}
	return std::string(s, len);
}

int ccad_write_step(Shape s, const char *filename, int nlen) {
	std::string fname = make_string(filename, nlen);
	return ccadcpp_write_step(s, fname.c_str());
}

int ccad_write_stl(Shape s, const char *filename, int nlen, float resolution) {
	std::string fname = make_string(filename, nlen);
	return ccadcpp_write_stl(s, fname.c_str(), resolution);
}

int ccad_write_iges(Shape s, const char *filename, int nlen) {
	std::string fname = make_string(filename, nlen);
	return ccadcpp_write_iges(s, fname.c_str());
}

int ccad_write_obj(Shape s, const char *filename, int nlen, float resolution) {
	std::string fname = make_string(filename, nlen);
	return ccadcpp_write_obj(s, fname.c_str(), resolution);
}

} //extern "C"--
