#ifndef CCAD_H
#define CCAD_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/* Opaque handle - the implementation keeps it as a TopoDS_Shape* */
typedef void *Shape;

/* ---- shape creation -------------------------------------------------- */
Shape ccad_circle(double cx, double cy, double radius);
Shape ccad_rectangle(double x1, double y1, double x2, double y2);
Shape ccad_triangle(double x1, double y1, double x2, double y2, double x3, double y3);
Shape ccad_polygon(const double *pts, size_t npts);

Shape ccad_sphere(double cx, double cy, double cz, double radius);
Shape ccad_box(double x1, double y1, double z1, double x2, double y2, double z2);
Shape ccad_cylinder(double cx, double cy, double cz, double r, double h);
Shape ccad_cone(double cx, double cy, double cz, double r1, double r2, double h);
Shape ccad_torus(double cx, double cy, double cz, double R, double r);
Shape ccad_arc(double x1, double y1, double x2, double y2, double x3, double y3);

/* ---- transformations ------------------------------------------------- */
Shape ccad_translate(Shape s, double dx, double dy, double dz);
Shape ccad_rotate(Shape s, double ax, double ay, double az,
		   double ux, double uy, double uz, double angle);
Shape ccad_scale(Shape s, double sx, double sy, double sz);
Shape ccad_mirror(Shape s, double ax, double ay, double az,
		   double ux, double uy, double uz);
Shape ccad_transform(Shape s, const double mat[16]);

/* ---- boolean / set ops ---------------------------------------------- */
Shape ccad_union(Shape a, Shape b);
Shape ccad_difference(Shape a, Shape b);

/* ---- extrude / sweep / revolve / loft -------------------------------- */
Shape ccad_extrude(Shape shape2d, double ux, double uy, double uz);
Shape ccad_sweep(Shape profile, Shape path);
Shape ccad_revolve(Shape s, double ax, double ay, double az,
		    double ux, double uy, double uz, double angle);
Shape ccad_loft(const Shape *shapes, size_t n, int solid, int ruled);

/* ---- advanced primitives --------------------------------------------- */
Shape ccad_quadratic_spline_wire(const double *pts, size_t npts);
Shape ccad_face_from_wire(Shape wire);

/* ---- mesh / export --------------------------------------------------- */
int ccad_write_step(Shape s, const char *filename, int nlen);
int ccad_write_stl(Shape s, const char *filename, int nlen, float resolution);
int ccad_write_iges(Shape s, const char *filename, int nlen);
int ccad_write_obj(Shape s, const char *filename, int nlen, float resolution);

/* ---- shape cleanup --------------------------------------------------- */
Shape ccad_simplify(const Shape s);
Shape ccad_simplify_to_solid(const Shape s);

/* ---- analysis / diagnostics ------------------------------------------- */
int ccad_bounding_sphere(Shape s, double *x, double *y, double *z,
			  double *r);
double ccad_volume(Shape s);
double ccad_surface_area(Shape s);
int ccad_centroid(Shape s, double *cx, double *cy, double *cz);
int ccad_inertia(Shape s, double density,
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
int ccad_shape_type(Shape s);

/* ---- cleanup -------------------------------------------------------- */
void ccad_free(Shape s);

#ifdef __cplusplus
}
#endif
#endif /* CCAD_H */
