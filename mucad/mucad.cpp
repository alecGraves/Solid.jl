// Build instructions (example with g++ on Linux):
//      g++ -c mucad.cpp -I. -o mucad.o
//      g++ main.c mucad.o  -L<OCCT-lib>  -lTKernel -lTKMath -lTKPrim 
//          -lTKGeomBase -lTKGeomAlgo -lTKBRep -lTKTopAlgo 
//          -lTKOffset -lTKBool -lTKShHealing -lTKFix 
//          -lTKMesh -lTKSTEP -lTKSTEPAttr -lTKSTL -lTKIGES -o demo

#include "mucad.h"

// OCCT Includes
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepOffsetAPI_MakePipeShell.hxx>
#include <BRepOffsetAPI_ThruSections.hxx> // <--- Needed for Loft
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepAlgoAPI_Fuse.hxx>
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <gp_GTrsf.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <GeomAPI_Interpolate.hxx>        // <--- Replacement for manual Spline
#include <TColgp_HArray1OfPnt.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Wire.hxx>
#include <TopoDS_Face.hxx>
#include <TopExp_Explorer.hxx>

#include <TColgp_Array1OfPnt.hxx>  // quadratic spline
#include <TColStd_Array1OfReal.hxx>
#include <TColStd_Array1OfInteger.hxx>

#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_GTrsf.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <STEPControl_Writer.hxx>
#include <IGESControl_Writer.hxx>
#include <StlAPI_Writer.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <BRep_Tool.hxx>
#include <Poly_Triangulation.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <fstream>
#include <iomanip>

#include <GCPnts_UniformAbscissa.hxx>
#include <BRepAdaptor_CompCurve.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepOffsetAPI_MakeFilling.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <ShapeFix_Solid.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <vector>
#include <cmath>
#include <algorithm>

#include <vector>

// ----------------------------------------------------------------
// Additional includes for missing types
// ----------------------------------------------------------------
#include <BRepBndLib.hxx>               // BRepBndLib::Add
#include <Bnd_Box.hxx>                  // Bnd_Box
#include <GProp_GProps.hxx>             // GProp_GProps
#include <BRepGProp.hxx>                // BRepGProp
#include <Precision.hxx>                // Precision::Confusion, Precision::Infinite
#include <ShapeAnalysis_FreeBounds.hxx> // ShapeAnalysis_FreeBounds::ConnectEdgesToWires
#include <ShapeFix_Solid.hxx>           // ShapeFix_Solid
#include <BRepTools_WireExplorer.hxx>   // BRepTools_WireExplorer
#include <Standard_Boolean.hxx>         // Standard_False, Standard_True

// Add new includes
#include <BRep_Builder.hxx>
#include <TopTools_HSequenceOfShape.hxx>


// ----------------------------------------------------------------
// Memory management configuration
// ----------------------------------------------------------------
// For environment variables
#ifdef _WIN32
    #define SETENV(name, value) _putenv_s(name, value)
#else
    #include <stdlib.h>
    #define SETENV(name, value) setenv(name, value, 1)
#endif


struct OCCTConfigurator {
    OCCTConfigurator() {
        // 1. Enable Optimized Pool Allocator
        SETENV("MMGT_OPT", "1");

        // 2. Disable Zeroing (Performance)
        SETENV("MMGT_CLEAR", "0");

        // 3. Increase Small Block Pool limit to 256 (4 cache lines)
        SETENV("MMGT_CELLSIZE", "256");

        // 4. Increase Recycle Threshold to 64KB
        SETENV("MMGT_THRESHOLD", "65536");
        
        // 5. Keep MMGT_MMAP default (1)
    }
};

static OCCTConfigurator auto_config = {};

// ----------------------------------------------------------------
// Helper: Shape Handle Conversion
// ----------------------------------------------------------------
static inline Shape toShape(const TopoDS_Shape& s) {
    return new TopoDS_Shape(s);
}

static inline const TopoDS_Shape& fromShape(Shape s) {
    return *static_cast<const TopoDS_Shape*>(s);
}

// ----------------------------------------------------------------
// Cleanup (all returns)
// ----------------------------------------------------------------
void mucadcpp_free(Shape s) {
    if (s) delete static_cast<TopoDS_Shape*>(s);
}

// ----------------------------------------------------------------
// 2-D Primitives
// ----------------------------------------------------------------
Shape mucadcpp_circle(double cx, double cy, double radius) {
    try {
        // 1. Define the Axis (Position and Orientation)
        gp_Pnt center(cx, cy, 0.0);
        gp_Dir dir(0.0, 0.0, 1.0); // Z-up
        gp_Ax2 axis(center, dir);

        // 2. Create the Geometric Circle
        gp_Circ circle_geom(axis, radius);

        // 3. Build the Topological Edge
        BRepBuilderAPI_MakeEdge mkEdge(circle_geom);
        
        // Check if construction failed (optional but recommended)
        if (!mkEdge.IsDone()) return nullptr;

        TopoDS_Edge edge = mkEdge.Edge();

        // 4. Wrap in a Wire (to match your original intent)
        BRepBuilderAPI_MakeWire mkWire(edge);
        
        if (!mkWire.IsDone()) return nullptr;

        TopoDS_Wire w = mkWire.Wire();

        return toShape(w);
    } catch (...) { return nullptr; }
}

Shape mucadcpp_rectangle(double x1, double y1, double x2, double y2) {
    try {
        gp_Pnt p1(x1, y1, 0.0);
        gp_Pnt p2(x2, y1, 0.0);
        gp_Pnt p3(x2, y2, 0.0);
        gp_Pnt p4(x1, y2, 0.0);

        BRepBuilderAPI_MakeWire mkWire;
        mkWire.Add(BRepBuilderAPI_MakeEdge(p1, p2).Edge());
        mkWire.Add(BRepBuilderAPI_MakeEdge(p2, p3).Edge());
        mkWire.Add(BRepBuilderAPI_MakeEdge(p3, p4).Edge());
        mkWire.Add(BRepBuilderAPI_MakeEdge(p4, p1).Edge());
        
        return toShape(mkWire.Wire());
    } catch (...) { return nullptr; }
}

Shape mucadcpp_triangle(double x1, double y1, double x2, double y2, double x3, double y3) {
    try {
        gp_Pnt p1(x1,y1,0), p2(x2,y2,0), p3(x3,y3,0);
        BRepBuilderAPI_MakeWire mkWire;
        mkWire.Add(BRepBuilderAPI_MakeEdge(p1,p2).Edge());
        mkWire.Add(BRepBuilderAPI_MakeEdge(p2,p3).Edge());
        mkWire.Add(BRepBuilderAPI_MakeEdge(p3,p1).Edge());
        return toShape(mkWire.Wire());
    } catch (...) { return nullptr; }
}

Shape mucadcpp_polygon(const double* pts, size_t npts) {
    if (npts < 3) return nullptr;
    try {
        BRepBuilderAPI_MakeWire mkWire;
        for (size_t i=0; i<npts; ++i) {
            gp_Pnt p1(pts[2*i], pts[2*i+1], 0.0);
            gp_Pnt p2(pts[2*((i+1)%npts)], pts[2*((i+1)%npts)+1], 0.0);
            mkWire.Add(BRepBuilderAPI_MakeEdge(p1, p2).Edge());
        }
        return toShape(mkWire.Wire());
    } catch (...) { return nullptr; }
}

// ----------------------------------------------------------------
// 3-D Primitives
// ----------------------------------------------------------------
Shape mucadcpp_sphere(double cx, double cy, double cz, double radius) {
    try {
        gp_Pnt center(cx, cy, cz);
        return toShape(BRepPrimAPI_MakeSphere(center, radius).Shape());
    } catch (...) { return nullptr; }
}

Shape mucadcpp_box(double x1, double y1, double z1, double x2, double y2, double z2) {
    try {
        gp_Pnt p1(x1, y1, z1);
        gp_Pnt p2(x2, y2, z2);
        return toShape(BRepPrimAPI_MakeBox(p1, p2).Shape());
    } catch (...) { return nullptr; }
}

// ----------------------------------------------------------------
// Transformations
// ----------------------------------------------------------------
Shape mucadcpp_translate(Shape s, double dx, double dy, double dz) {
    if (!s) return nullptr;
    try {
        gp_Trsf trsf; 
        trsf.SetTranslation(gp_Vec(dx, dy, dz));
        BRepBuilderAPI_Transform transf(fromShape(s), trsf);
        return toShape(transf.Shape());
    } catch (...) { return nullptr; }
}

Shape mucadcpp_rotate(Shape s, double ax, double ay, double az, 
                  double ux, double uy, double uz, double angle) {
    if (!s) return nullptr;
    try {
        gp_Ax1 axis(gp_Pnt(ax, ay, az), gp_Dir(ux, uy, uz));
        gp_Trsf trsf; 
        trsf.SetRotation(axis, angle);
        BRepBuilderAPI_Transform transf(fromShape(s), trsf);
        return toShape(transf.Shape());
    } catch (...) { return nullptr; }
}

Shape mucadcpp_transform(Shape s, const double mat[16])
{
    if (!s || !mat) return nullptr;
    try {
        gp_GTrsf gtrsf;

        // OCCT matrices are 1-based.
        // We map the input row-major 4x4 matrix to the Affine gp_GTrsf.
        // [ R  T ]
        // [ 0  1 ]
        
        // Row 1
        gtrsf.SetValue(1, 1, mat[0]);
        gtrsf.SetValue(1, 2, mat[1]);
        gtrsf.SetValue(1, 3, mat[2]);
        gtrsf.SetValue(1, 4, mat[3]); // Tx

        // Row 2
        gtrsf.SetValue(2, 1, mat[4]);
        gtrsf.SetValue(2, 2, mat[5]);
        gtrsf.SetValue(2, 3, mat[6]);
        gtrsf.SetValue(2, 4, mat[7]); // Ty

        // Row 3
        gtrsf.SetValue(3, 1, mat[8]);
        gtrsf.SetValue(3, 2, mat[9]);
        gtrsf.SetValue(3, 3, mat[10]);
        gtrsf.SetValue(3, 4, mat[11]); // Tz

        // BRepBuilderAPI_GTransform handles the heavy lifting.
        // It will convert analytic surfaces to BSplines if the matrix contains 
        // shear or non-uniform scaling.
        // copy = Standard_True (safe default)
        BRepBuilderAPI_GTransform applier(fromShape(s), gtrsf, Standard_True);
        
        return toShape(applier.Shape());

    } catch (...) {
        return nullptr;
    }
}

Shape mucadcpp_scale(Shape s, double sx, double sy, double sz)
{
    if (!s) return nullptr;
    try {
        // We use GTrsf (General Transform) because standard gp_Trsf 
        // only supports Uniform scaling (sx=sy=sz).
        gp_GTrsf gtrsf;
        
        // Set diagonal elements for scaling along global axes
        gtrsf.SetValue(1, 1, sx);
        gtrsf.SetValue(2, 2, sy);
        gtrsf.SetValue(3, 3, sz);
        
        // Use GTransform to allow geometry deformation (e.g., Sphere -> Ellipsoid)
        BRepBuilderAPI_GTransform applier(fromShape(s), gtrsf, Standard_True);
        
        return toShape(applier.Shape());
    } catch (...) {
        return nullptr;
    }
}

Shape mucadcpp_mirror(Shape s, double ax, double ay, double az, 
                            double ux, double uy, double uz)
{
    if (!s) return nullptr;
    try {
        gp_Pnt p(ax, ay, az);
        gp_Dir d(ux, uy, uz);
        gp_Ax2 plane(p, d); // Defines a coordinate system (Main Axis = Normal)

        gp_Trsf trsf;
        trsf.SetMirror(plane); // Mirror across the XY plane of the given axis

        // Standard Transform is sufficient (Mirror is a rigid body transform)
        BRepBuilderAPI_Transform applier(fromShape(s), trsf);
        
        return toShape(applier.Shape());
    } catch (...) { return nullptr; }
}

// ----------------------------------------------------------------
// Operations: Extrude, Sweep, Revolve, Loft
// ----------------------------------------------------------------
Shape mucadcpp_extrude(Shape shape2d, double ux, double uy, double uz) {
    if (!shape2d) return nullptr;
    try {
        const TopoDS_Shape& w = fromShape(shape2d);
        TopoDS_Face f;
        
        if (w.ShapeType() == TopAbs_FACE) {
            f = TopoDS::Face(w);
        } else if (w.ShapeType() == TopAbs_WIRE) {
            // BRepBuilderAPI_MakeFace assumes wire is planar
            f = BRepBuilderAPI_MakeFace(TopoDS::Wire(w)).Face();
        } else {
            return nullptr; 
        }
        
        gp_Vec dir(ux, uy, uz);
        return toShape(BRepPrimAPI_MakePrism(f, dir).Shape());
    } catch (...) { return nullptr; }
}

Shape mucadcpp_sweep(Shape profile, double x0, double y0, double z0, double x1, double y1, double z1) {
    if (!profile) return nullptr;
    try {
        // 1. Create the Path (Spine)
        gp_Pnt p0(x0,y0,z0);
        gp_Pnt p1(x1,y1,z1);
        TopoDS_Edge pathEdge = BRepBuilderAPI_MakeEdge(p0, p1).Edge();
        TopoDS_Wire pathWire = BRepBuilderAPI_MakeWire(pathEdge).Wire();

        // 2. Create the Pipe maker using the Spine
        BRepOffsetAPI_MakePipeShell pipe(pathWire);

        // 3. Add the profile (the cross section)
        // Note: The profile should ideally be perpendicular to the start of the path
        pipe.Add(fromShape(profile));
        
        pipe.Build();
        if (!pipe.IsDone()) return nullptr;
        return toShape(pipe.Shape());
    } catch (...) { return nullptr; }
}

Shape mucadcpp_revolve(Shape s, double ax, double ay, double az, 
                   double ux, double uy, double uz, double angle) {
    if (!s) return nullptr;
    try {
        gp_Pnt p(ax, ay, az);
        gp_Dir d(ux, uy, uz);
        gp_Ax1 axis(p, d);
        BRepPrimAPI_MakeRevol revol(fromShape(s), axis, angle);
        return toShape(revol.Shape());
    } catch (...) { return nullptr; }
}


// ----------------------------------------------------------------
// Loft Implementation
// ----------------------------------------------------------------
#include <GCPnts_QuasiUniformAbscissa.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <GCPnts_AbscissaPoint.hxx>

// #include <TopoDS.hxx>
// #include <TopoDS_Shape.hxx>
// #include <TopoDS_Wire.hxx>
// #include <TopExp_Explorer.hxx>

// // OCCT Modeling
// #include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepOffsetAPI_ThruSections.hxx>
#include <BRepFill_CompatibleWires.hxx>
// #include <ShapeFix_Solid.hxx>
// #include <TopTools_SequenceOfShape.hxx>

/* ----------------------------------------------------------------
 * Main loft routine
 * * Logic changed: Instead of manual point sampling (which smooths corners),
 * we now use BRepFill_CompatibleWires. This OCCT tool automatically 
 * cuts edges (e.g., splits a circle into 4 arcs) to match the topology 
 * of the target shape (e.g., a square) while preserving sharp vertices.
 * ---------------------------------------------------------------- */
Shape mucadcpp_loft(const Shape* shapes,
                    size_t n,
                    bool solid,
                    bool ruled)
{
    if (!shapes || n < 2) return nullptr;

    try {
        /* 1. Collect input Wires */
        TopTools_SequenceOfShape inputWires;

        for (size_t i = 0; i < n; ++i) {
            TopoDS_Shape s = fromShape(shapes[i]);
            if (s.IsNull()) return nullptr;

            // Extract the Wire from the Shape
            if (s.ShapeType() != TopAbs_WIRE) {
                TopExp_Explorer exp(s, TopAbs_WIRE);
                if (exp.More()) {
                    s = exp.Current();
                } else {
                    // If no wire found (e.g. passing an Edge), try to make one
                    TopExp_Explorer expE(s, TopAbs_EDGE);
                    if (expE.More()) {
                        BRepBuilderAPI_MakeWire mkWire;
                        mkWire.Add(TopoDS::Edge(expE.Current()));
                        s = mkWire.Wire();
                    } else {
                        return nullptr; // Invalid input
                    }
                }
            }
            inputWires.Append(TopoDS::Wire(s));
        }

        /* 2. Make Wires Compatible 
           This logic aligns the wires and splits edges so that a Circle (1 edge)
           can loft to a Square (4 edges) without smoothing the Square's corners. */
        BRepFill_CompatibleWires fixer(inputWires);
        fixer.SetPercent(0.1); // Precision for matching vertices
        fixer.Perform();

        TopTools_SequenceOfShape correctedWires;
        bool compatibilityMode = Standard_True;

        if (fixer.IsDone()) {
            // Use the processed wires (Circle is now split into arcs)
            correctedWires = fixer.Shape();
        } else {
            // Fallback: Use original wires if compatibility alg fails.
            // We disable compatibility check in ThruSections so it doesn't crash,
            // though the mesh might be twisted.
            correctedWires = inputWires;
            compatibilityMode = Standard_False;
        }

        /* 3. Build the Loft */
        // Param 3 (1e-6) is precision
        BRepOffsetAPI_ThruSections generator(solid, ruled, 1e-6);
        
        // If wires were successfully made compatible, we turn on the check.
        // This creates a cleaner mesh. If not, we turn it off to force creation.
        generator.CheckCompatibility(compatibilityMode);

        for (int i = 1; i <= correctedWires.Length(); i++) {
            generator.AddWire(TopoDS::Wire(correctedWires.Value(i)));
        }

        generator.Build();
        if (!generator.IsDone()) return nullptr;

        TopoDS_Shape result = generator.Shape();

        /* 4. Fix Solid topology if requested */
        if (solid && result.ShapeType() == TopAbs_SOLID) {
            ShapeFix_Solid fix;
            fix.Init(TopoDS::Solid(result));
            fix.Perform();
            result = fix.Solid();
        }

        return toShape(result);

    } catch (...) {
        return nullptr;
    }
}

// ----------------------------------------------------------------
// Booleans
// ----------------------------------------------------------------
Shape mucadcpp_union(Shape a, Shape b) {
    if (!a || !b) return nullptr;
    try {
        return toShape(BRepAlgoAPI_Fuse(fromShape(a), fromShape(b)).Shape());
    } catch (...) { return nullptr; }
}

Shape mucadcpp_difference(Shape a, Shape b) {
    if (!a || !b) return nullptr;
    try {
        return toShape(BRepAlgoAPI_Cut(fromShape(a), fromShape(b)).Shape());
    } catch (...) { return nullptr; }
}

// ----------------------------------------------------------------
// Splines & Advanced
// ----------------------------------------------------------------


/* --------------------------------------------------------------------
   Create a wire from 3-node quadratic segments (FEM style).
   
   The array is interpreted as:
   Segment 1: P[0] (Start), P[1] (Mid), P[2] (End)
   Segment 2: P[2] (Start), P[3] (Mid), P[4] (End)
   ...
   
   pts  : Flat array of coords {x0,y0,z0, x1,y1,z1...}
   npts : Total points. Must be odd (2*N_segments + 1).
          E.g., 3 points for 1 segment, 5 points for 2 segments.

Example:
double nodes[] = {
    0.0, 0.0, 0.0,   // Start
    0.5, 0.2, 0.0,   // Mid 1 (Curve goes through here)
    1.0, 0.0, 0.0,   // End 1 / Start 2
    1.5, -0.2, 0.0,  // Mid 2 (Curve goes through here)
    2.0, 0.0, 0.0    // End 2
};

// 5 points = 2 segments
Shape fem_wire = occt_quadratic_spline_wire(nodes, 5);
-------------------------------------------------------------------- */
Shape mucadcpp_quadratic_spline_wire(const double* pts, size_t npts)
{
    // We need at least 3 points (1 segment) and an odd number of points
    // to maintain the Start-Mid-End chain continuity.
    if (!pts || npts < 3 || (npts % 2 == 0)) return nullptr;

    try {
        BRepBuilderAPI_MakeWire mkWire;
        
        // We construct a Quadratic B-Spline (Degree 2) for each triplet.
        // A degree 2 curve through 3 points P0, P1, P2 is mathematically defined.
        // 
        // To do this explicitly in OCCT without an approximator (faster):
        // We calculate the geometric "Poles" (Control Points) required to make
        // the curve pass through the Mid point.
        //
        // For a rational B-Spline passing through P0(t=0), P1(t=0.5), P2(t=1):
        // The middle Control Pole (C1) is NOT P1. 
        // C1 = 2*P1 - 0.5*P0 - 0.5*P2
        
        for (size_t i = 0; i < npts - 2; i += 2) {
            gp_Pnt P0(pts[3*i],     pts[3*i+1],     pts[3*i+2]);
            gp_Pnt Pmid(pts[3*(i+1)], pts[3*(i+1)+1], pts[3*(i+1)+2]);
            gp_Pnt P2(pts[3*(i+2)], pts[3*(i+2)+1], pts[3*(i+2)+2]);

            gp_Pnt C1;
            C1.SetX( 2.0 * Pmid.X() - 0.5 * P0.X() - 0.5 * P2.X() );
            C1.SetY( 2.0 * Pmid.Y() - 0.5 * P0.Y() - 0.5 * P2.Y() );
            C1.SetZ( 2.0 * Pmid.Z() - 0.5 * P0.Z() - 0.5 * P2.Z() );

            TColgp_Array1OfPnt poles(1,3);
            poles(1) = P0;  poles(2) = C1;  poles(3) = P2;

            TColStd_Array1OfReal knots(1,2);
            knots(1) = 0.0; knots(2) = 1.0;

            TColStd_Array1OfInteger mults(1,2);
            mults(1) = 3; mults(2) = 3;

            Handle(Geom_Curve) curve = new Geom_BSplineCurve(poles, knots, mults, 2);

            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve).Edge();
            mkWire.Add(edge);
        }

        if (!mkWire.IsDone()) return nullptr;
        return toShape(mkWire.Wire());

    } catch (...) {
        return nullptr;
    }
}

Shape mucadcpp_face_from_wire(Shape wire) {
    if (!wire) return nullptr;
    try {
        return toShape(BRepBuilderAPI_MakeFace(TopoDS::Wire(fromShape(wire))).Face());
    } catch (...) { return nullptr; }
}


// ----------------------------------------------------------------
// Analysis
// ----------------------------------------------------------------
int mucadcpp_bounding_sphere(Shape s, double* x, double* y, double* z, double* r)
{
    if (!s || !x || !y || !z || !r) return -1;
    try {
        // 1. Compute Axis Aligned Bounding Box (AABB)
        Bnd_Box box;
        // Use triangulation if available, otherwise geometry controls
        // true = use triangulation (more precise if meshed), false = geometry
        BRepBndLib::Add(fromShape(s), box, false); 
        
        if (box.IsVoid()) return -1;

        // 2. Get Min/Max
        double xmin, ymin, zmin, xmax, ymax, zmax;
        box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

        // 3. Compute Center
        *x = (xmin + xmax) * 0.5;
        *y = (ymin + ymax) * 0.5;
        *z = (zmin + zmax) * 0.5;

        // 4. Compute Radius (Half Diagonal)
        // This is a "loose" bounding sphere. Minimal enclosing sphere 
        // is much harder to compute, but this covers the object safely.
        double dx = xmax - xmin;
        double dy = ymax - ymin;
        double dz = zmax - zmin;
        
        *r = 0.5 * std::sqrt(dx*dx + dy*dy + dz*dz);
        
        return 0;
    } catch (...) { return -1; }
}

double mucadcpp_volume(Shape s)
{
    if (!s) return 0.0;
    try {
        GProp_GProps props;
        BRepGProp::VolumeProperties(fromShape(s), props);
        return props.Mass(); // Mass = Volume when density is 1 (default)
    } catch (...) { return 0.0; }
}

double mucadcpp_surface_area(Shape s)
{
    if (!s) return 0.0;
    try {
        GProp_GProps props;
        BRepGProp::SurfaceProperties(fromShape(s), props);
        return props.Mass(); // Mass = Area for surface props
    } catch (...) { return 0.0; }
}

int mucadcpp_centroid(Shape s, double* cx, double* cy, double* cz)
{
    if (!s || !cx || !cy || !cz) return -1;
    try {
        const TopoDS_Shape& shape = fromShape(s);
        GProp_GProps props;

        // Determine best property type based on shape dimension
        // 1. Try Volume (Solid)
        BRepGProp::VolumeProperties(shape, props);
        
        // If volume is negligible (e.g. it's a Face or Wire), try Surface
        if (props.Mass() < Precision::Confusion()) {
            BRepGProp::SurfaceProperties(shape, props);
        }
        
        // If Area is negligible (e.g. it's a Wire/Edge), try Linear
        if (props.Mass() < Precision::Confusion()) {
             BRepGProp::LinearProperties(shape, props);
        }

        // If still zero, we can't compute a centroid
        if (props.Mass() < Precision::Confusion()) return -1;

        gp_Pnt center = props.CentreOfMass();
        *cx = center.X();
        *cy = center.Y();
        *cz = center.Z();
        
        return 0;
    } catch (...) { return -1; }
}

int mucadcpp_inertia(Shape s, double density, 
                  double* Ixx, double* Iyy, double* Izz)
{
    if (!s || !Ixx || !Iyy || !Izz) return -1;
    try {
        GProp_GProps props;
        
        // We assume the user wants Inertia of a SOLID volume.
        // (Using VolumeProperties calculates geometric inertia based on Volume)
        BRepGProp::VolumeProperties(fromShape(s), props);

        // Get the Matrix of Inertia relative to the GLOBAL Origin.
        // Note: This is the geometric inertia (integral of position^2 * dV).
        gp_Mat matrix = props.MatrixOfInertia();

        // Extract diagonal elements (Moments about global X, Y, Z)
        // Scale by density (Mass = Volume * Density)
        *Ixx = matrix(1, 1) * density;
        *Iyy = matrix(2, 2) * density;
        *Izz = matrix(3, 3) * density;

        return 0;
    } catch (...) {
        return -1;
    }
}

int mucadcpp_shape_type(Shape s)
{
    TopoDS_Shape shape = fromShape(s);
    if (shape.IsNull())
        return -1;  // invalid or null shape

    // Return the OCCT TopAbs shape type as an integer.
    return static_cast<int>(shape.ShapeType());
}

// ----------------------------------------------------------------
// I/O and Cleanup
// ----------------------------------------------------------------

int mucadcpp_write_step(Shape s, const char* filename) {
    if (!s || !filename) return -1;
    try {
        STEPControl_Writer writer;
        writer.Transfer(fromShape(s), STEPControl_AsIs);
        return (writer.Write(filename) == IFSelect_RetDone) ? 0 : -1;
    } catch (...) { return -1; }
}

int mucadcpp_write_stl(Shape s, const char* filename, float resolution) {
    if (!s || !filename) return -1;
    try {
        const TopoDS_Shape& shape = fromShape(s);

        // Force triangulation with a reasonable deflection
        BRepMesh_IncrementalMesh mesh(shape, resolution);
        if (!mesh.IsDone()) { // lazy evaluates.
            return -1; // Meshing failed
        }

        StlAPI_Writer writer;
        // Optionally set to ASCII for debugging (use Standard_False)
        writer.ASCIIMode() = false; // true = ASCII, false = Binary

        Standard_Boolean success = writer.Write(shape, filename);
        return success ? 0 : -1;
    } catch (...) {
        return -1;
    }
}

int mucadcpp_write_iges(Shape s, const char* filename)
{
    if (!s || !filename) return -1;
    try {
        IGESControl_Writer writer;
        // 0 = Metric (Millimeters), 1 = Inches. OCCT default is usually MM.
        // const int unit_mm = 0; 
        // writer.ChangeGlobalSection().SetUnitFlag(unit_mm);
        
        Standard_Boolean ok = writer.AddShape(fromShape(s));
        if (!ok) return -1;
        
        writer.ComputeModel();
        Standard_Boolean success = writer.Write(filename);
        return success ? 0 : -1;
    } catch (...) { return -1; }
}

int mucadcpp_write_obj(Shape s, const char* filename)
{
    if (!s || !filename) return -1;
    try {
        const TopoDS_Shape& shape = fromShape(s);

        // 1. Ensure the shape is triangulated (Meshed)
        //    Deflection controls quality (lower = smoother).
        BRepMesh_IncrementalMesh meshGen(shape, 0.01); 
        if (!meshGen.IsDone()) return -1;

        std::ofstream file(filename);
        if (!file.is_open()) return -1;

        file << "# Generated by mucadcpp/OCCT Wrapper\n";
        file << std::fixed << std::setprecision(6);

        // OBJ indices are global and 1-based. We must track offsets per face.
        int global_vertex_offset = 1;

        // 2. Iterate over every face to extract its triangulation
        TopExp_Explorer expl(shape, TopAbs_FACE);
        for (; expl.More(); expl.Next()) {
            const TopoDS_Face& face = TopoDS::Face(expl.Current());
            TopLoc_Location loc;
            Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);

            if (tri.IsNull()) continue;

            // Get nodes (vertices)
            // OCCT Poly_Triangulation nodes are 1-based in the array
            Standard_Integer nbNodes = tri->NbNodes();
            for (Standard_Integer i = 1; i <= nbNodes; ++i) {
                gp_Pnt p = tri->Node(i);
                // Apply location transformation if the face is instanced/moved
                p.Transform(loc.Transformation());
                file << "v " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
            }

            // Get triangles (faces)
            Standard_Integer nbTriangles = tri->NbTriangles();
            for (Standard_Integer i = 1; i <= nbTriangles; ++i) {
                const Poly_Triangle& t = tri->Triangle(i);
                Standard_Integer n1, n2, n3;
                t.Get(n1, n2, n3);
                
                // Check face orientation. If reversed, flip normal by swapping indices
                if (face.Orientation() == TopAbs_REVERSED) {
                     file << "f " << (n1 + global_vertex_offset - 1) << " " 
                                  << (n3 + global_vertex_offset - 1) << " " 
                                  << (n2 + global_vertex_offset - 1) << "\n";
                } else {
                     file << "f " << (n1 + global_vertex_offset - 1) << " " 
                                  << (n2 + global_vertex_offset - 1) << " " 
                                  << (n3 + global_vertex_offset - 1) << "\n";
                }
            }

            // Increment offset for the next face
            global_vertex_offset += nbNodes;
        }

        file.close();
        return 0;
    } catch (...) { return -1; }
}
