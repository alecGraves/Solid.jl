#ifndef CCAD_HPP
#define CCAD_HPP

#include <stddef.h>

/* Opaque handle - the implementation keeps it as a TopoDS_Shape* */
typedef void *Shape;

/* ---- shape creation -------------------------------------------------- */
Shape ccadcpp_circle(double cx, double cy, double radius);
Shape ccadcpp_rectangle(double x1, double y1, double x2, double y2);
Shape ccadcpp_triangle(double x1, double y1, double x2, double y2, double x3, double y3);
Shape ccadcpp_polygon(const double *pts, size_t npts);

Shape ccadcpp_sphere(double cx, double cy, double cz, double radius);
Shape ccadcpp_box(double x1, double y1, double z1, double x2, double y2, double z2);
Shape ccadcpp_cylinder(double cx, double cy, double cz, double r, double h);
Shape ccadcpp_cone(double cx, double cy, double cz, double r1, double r2, double h);
Shape ccadcpp_torus(double cx, double cy, double cz, double R, double r);
Shape ccadcpp_arc(double x1, double y1, double x2, double y2, double x3, double y3);

/* ---- transformations ------------------------------------------------- */
Shape ccadcpp_translate(Shape s, double dx, double dy, double dz);
Shape ccadcpp_rotate(Shape s, double ax, double ay, double az,
		      double ux, double uy, double uz, double angle);
Shape ccadcpp_scale(Shape s, double sx, double sy, double sz);
Shape ccadcpp_mirror(Shape s, double ax, double ay, double az,
		      double ux, double uy, double uz);
Shape ccadcpp_transform(Shape s, const double mat[16]);

/* ---- boolean / set ops ---------------------------------------------- */
Shape ccadcpp_union(Shape a, Shape b);
Shape ccadcpp_difference(Shape a, Shape b);

/* ---- extrude / sweep / revolve / loft -------------------------------- */
Shape ccadcpp_extrude(Shape shape2d, double ux, double uy, double uz);
Shape ccadcpp_sweep(Shape profile, Shape path);
Shape ccadcpp_revolve(Shape s, double ax, double ay, double az,
		       double ux, double uy, double uz, double angle);
Shape ccadcpp_loft(const Shape *shapes, size_t n, bool solid, bool ruled);

/* ---- advanced primitives --------------------------------------------- */
Shape ccadcpp_quadratic_spline_wire(const double *pts, size_t npts);
Shape ccadcpp_face_from_wire(Shape wire);

/* ---- mesh / export --------------------------------------------------- */
int ccadcpp_write_step(Shape s, const char *filename);
int ccadcpp_write_stl(Shape s, const char *filename, float resolution);
int ccadcpp_write_iges(Shape s, const char *filename);
int ccadcpp_write_obj(Shape s, const char *filename, float resolution);

/* ---- shape cleanup --------------------------------------------------- */
Shape ccadcpp_simplify(const Shape s);
Shape ccadcpp_simplify_to_solid(const Shape s);

/* ---- analysis / diagnostics ------------------------------------------- */
int ccadcpp_bounding_sphere(Shape s, double *x, double *y, double *z,
			     double *r);
double ccadcpp_volume(Shape s);
double ccadcpp_surface_area(Shape s);
int ccadcpp_centroid(Shape s, double *cx, double *cy, double *cz);
int ccadcpp_inertia(Shape s, double density,
		     double *Ixx, double *Iyy, double *Izz);

// Enumerator
// 0 = Edge, 1 = Wire, 2 = Face, 3 = Solid, 4 = Other
int ccadcpp_shape_type(Shape s);

/* ---- cleanup -------------------------------------------------------- */
void ccadcpp_free(Shape s);

#endif /* CCAD_HPP */
