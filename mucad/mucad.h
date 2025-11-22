#ifndef MUCAD_H
#define MUCAD_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/* Opaque handle - the implementation keeps it as a TopoDS_Shape* */
typedef void *Shape;

/* ---- shape creation -------------------------------------------------- */
Shape mucad_circle(double cx, double cy, double radius);
Shape mucad_rectangle(double x1, double y1, double x2, double y2);
Shape mucad_triangle(double x1, double y1, double x2, double y2, double x3, double y3);
Shape mucad_polygon(const double *pts, size_t npts);

Shape mucad_sphere(double cx, double cy, double cz, double radius);
Shape mucad_box(double x1, double y1, double z1, double x2, double y2, double z2);
Shape mucad_cylinder(double cx, double cy, double cz, double r, double h);
Shape mucad_cone(double cx, double cy, double cz, double r1, double r2, double h);
Shape mucad_torus(double cx, double cy, double cz, double R, double r);
Shape mucad_arc(double x1, double y1, double x2, double y2, double x3, double y3);

/* ---- transformations ------------------------------------------------- */
Shape mucad_translate(Shape s, double dx, double dy, double dz);
Shape mucad_rotate(Shape s, double ax, double ay, double az,
		   double ux, double uy, double uz, double angle);
Shape mucad_scale(Shape s, double sx, double sy, double sz);
Shape mucad_mirror(Shape s, double ax, double ay, double az,
		   double ux, double uy, double uz);
Shape mucad_transform(Shape s, const double mat[16]);

/* ---- boolean / set ops ---------------------------------------------- */
Shape mucad_union(Shape a, Shape b);
Shape mucad_difference(Shape a, Shape b);

/* ---- extrude / sweep / revolve / loft -------------------------------- */
Shape mucad_extrude(Shape shape2d, double ux, double uy, double uz);
Shape mucad_sweep(Shape profile, Shape path);
Shape mucad_revolve(Shape s, double ax, double ay, double az,
		    double ux, double uy, double uz, double angle);
Shape mucad_loft(const Shape *shapes, size_t n, int solid, int ruled);

/* ---- advanced primitives --------------------------------------------- */
Shape mucad_quadratic_spline_wire(const double *pts, size_t npts);
Shape mucad_face_from_wire(Shape wire);

/* ---- mesh / export --------------------------------------------------- */
int mucad_write_step(Shape s, const char *filename, int nlen);
int mucad_write_stl(Shape s, const char *filename, int nlen, float resolution);
int mucad_write_iges(Shape s, const char *filename, int nlen);
int mucad_write_obj(Shape s, const char *filename, int nlen, float resolution);

/* ---- shape cleanup --------------------------------------------------- */
Shape mucad_simplify(const Shape s);
Shape mucad_simplify_to_solid(const Shape s);

/* ---- analysis / diagnostics ------------------------------------------- */
int mucad_bounding_sphere(Shape s, double *x, double *y, double *z,
			  double *r);
double mucad_volume(Shape s);
double mucad_surface_area(Shape s);
int mucad_centroid(Shape s, double *cx, double *cy, double *cz);
int mucad_inertia(Shape s, double density,
		  double *Ixx, double *Iyy, double *Izz);

// Enumerator
// 1 Compound
// 2 Compsolid
// 3 Solid
// 4 Shell
// 5 Face
// 6 Wire
// 7 Edge
// 8 Vertex
// 9 Shape
int mucad_shape_type(Shape s);

/* ---- cleanup -------------------------------------------------------- */
void mucad_free(Shape s);

#ifdef __cplusplus
}
#endif
#endif /* MUCAD_H */
