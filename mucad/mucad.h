#ifndef MUCAD_H
#define MUCAD_H

#ifdef __cplusplus
extern "C" {
#endif

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
Shape mucad_transform(Shape s,const double mat[12]);

/* ---- boolean / set ops ---------------------------------------------- */
Shape mucad_union(Shape a,Shape b);
Shape mucad_difference(Shape a,Shape b);
Shape mucad_intersection(Shape a,Shape b);
Shape mucad_symdiff(Shape a,Shape b);
Shape mucad_union_many(const Shape* shapes,size_t n);
Shape mucad_difference_many(Shape base,const Shape* sub,size_t n);
Shape mucad_intersection_many(const Shape* shapes,size_t n);

/* ---- extrude / sweep / revolve / loft -------------------------------- */
Shape mucad_extrude(Shape shape2d,double ux,double uy,double uz);
Shape mucad_sweep(Shape profile,double x0,double y0,double z0,
                  double x1,double y1,double z1);
Shape mucad_revolve(Shape s,double ax,double ay,double az,
                    double ux,double uy,double uz,double angle);
Shape mucad_loft(const Shape* shapes,size_t n,bool solid,bool ruled);

/* ---- advanced primitives --------------------------------------------- */
Shape mucad_catmull_rom_wire(const double* pts,size_t npts,
                            double tension,bool closed,bool lineClose);
Shape mucad_quadratic_spline_wire(const double* pts,size_t npts);
Shape mucad_face_from_wire(Shape wire);

/* ---- mesh / export --------------------------------------------------- */
int mucad_mesh(Shape s,double deflection);
int mucad_write_step(Shape s,const char* filename);
int mucad_write_stl(Shape s,const char* filename);
int mucad_write_iges(Shape s,const char* filename);
int mucad_write_obj(Shape s,const char* filename);

/* ---- analysis / diagnostics ------------------------------------------- */
int mucad_bounding_box(Shape s,double* minx,double* miny,double* minz,
                       double* maxx,double* maxy,double* maxz);
double mucad_volume(Shape s);
double mucad_surface_area(Shape s);
int mucad_centroid(Shape s,double* cx,double* cy,double* cz);
int mucad_inertia(Shape s,double density,
                  double* Ixx,double* Iyy,double* Izz);
int mucad_shape_type(Shape s); // enum: 0=WIRE,1=FACE,2=SHELL,3=SOLID,4=COMPOUND
const char* mucad_last_error(void);

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
