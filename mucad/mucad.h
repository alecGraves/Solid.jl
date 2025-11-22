#ifndef MUCAD_H
#define MUCAD_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stddef.h>

/* Opaque handle - the implementation keeps it as a TopoDS_Shape* */
typedef void* Shape;

/* ---- shape creation -------------------------------------------------- */
Shape mucad_circle(double cx, double cy, double radius);
Shape mucad_rectangle(double x1,double y1,double x2,double y2);
Shape mucad_triangle(double x1,double y1,double x2,double y2,double x3,double y3);
Shape mucad_polygon(const double* pts,size_t npts);

Shape mucad_sphere(double cx,double cy,double cz,double radius);
Shape mucad_box(double x1,double y1,double z1,double x2,double y2,double z2);
Shape mucad_cylinder(double cx,double cy,double cz,double r,double h);
Shape mucad_cone(double cx,double cy,double cz,double r1,double r2,double h);
Shape mucad_torus(double cx,double cy,double cz,double R,double r);
Shape mucad_arc(double x1,double y1,double x2,double y2,double x3,double y3);

/* ---- transformations ------------------------------------------------- */
Shape mucad_translate(Shape s,double dx,double dy,double dz);
Shape mucad_rotate(Shape s,double ax,double ay,double az,
                   double ux,double uy,double uz,double angle);
Shape mucad_scale(Shape s,double sx,double sy,double sz);
Shape mucad_mirror(Shape s,double ax,double ay,double az,
                   double ux,double uy,double uz);
Shape mucad_transform(Shape s,const double mat[16]);

/* ---- boolean / set ops ---------------------------------------------- */
Shape mucad_union(Shape a,Shape b);
Shape mucad_difference(Shape a,Shape b);

/* ---- extrude / sweep / revolve / loft -------------------------------- */
Shape mucad_extrude(Shape shape2d,double ux,double uy,double uz);
Shape mucad_sweep(Shape profile,double x0,double y0,double z0,
                  double x1,double y1,double z1);
Shape mucad_revolve(Shape s,double ax,double ay,double az,
                    double ux,double uy,double uz,double angle);
Shape mucad_loft(const Shape* shapes,size_t n, int solid, int ruled);

/* ---- advanced primitives --------------------------------------------- */
Shape mucad_quadratic_spline_wire(const double* pts,size_t npts);
Shape mucad_face_from_wire(Shape wire);

/* ---- mesh / export --------------------------------------------------- */
int mucad_write_step(Shape s,const char* filename);
int mucad_write_stl(Shape s,const char* filename, float resolution);
int mucad_write_iges(Shape s,const char* filename);
int mucad_write_obj(Shape s,const char* filename);

/* ---- analysis / diagnostics ------------------------------------------- */
int mucad_bounding_sphere(Shape s, double* x, double* y, double* z,
        double* r);
double mucad_volume(Shape s);
double mucad_surface_area(Shape s);
int mucad_centroid(Shape s,double* cx,double* cy,double* cz);
int mucad_inertia(Shape s,double density,
                  double* Ixx,double* Iyy,double* Izz);

// Enumerator
// 1 TopAbs_COMPOUND     
// 2 TopAbs_COMPSOLID    
// 3 TopAbs_SOLID    
// 4 TopAbs_SHELL    
// 5 TopAbs_FACE     
// 6 TopAbs_WIRE     
// 7 TopAbs_EDGE     
// 8 TopAbs_VERTEX   
// 9 TopAbs_SHAPE 
int mucad_shape_type(Shape s); // enum: 0=WIRE,1=FACE,2=SHELL,3=SOLID,4=COMPOUND

/* ---- cleanup -------------------------------------------------------- */
void mucad_free(Shape s);

/* ---- optional helpers ------------------------------------------------ */
Shape mucad_clone(Shape s);
int mucad_set_tolerance(double tol);
double mucad_get_tolerance(void);

#ifdef __cplusplus
}
#endif
#endif /* MUCAD_H */
