// Build instructions (example with g++ on Linux):
//		g++ -c mucad.cpp -I. -o mucad.o
// 		g++ main.c mucad.o  -L<OCCT-lib>  -lTKernel -lTKMath -lTKPrim 
// 			-lTKGeomBase -lTKGeomAlgo -lTKBRep -lTKTopAlgo 
// 			-lTKOffset -lTKBool -lTKShHealing -lTKFix 
// 			-lTKMesh -lTKSTEP -lTKSTEPAttr -lTKSTL -lTKIGES -o demo

#include "mucad.h"

// OCCT Includes
#include <BRepPrimAPI_MakeCircle.hxx>
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
void mucad_free(Shape s) {
    if (s) delete static_cast<TopoDS_Shape*>(s);
}

// ----------------------------------------------------------------
// 2-D Primitives
// ----------------------------------------------------------------
Shape mucad_circle(double cx, double cy, double radius) {
    try {
        gp_Pnt center(cx, cy, 0.0);
        gp_Dir dir(0.0, 0.0, 1.0);
        BRepPrimAPI_MakeCircle mk(center, dir, radius);
        TopoDS_Wire w = BRepBuilderAPI_MakeWire(mk.Edge()).Wire();
        return toShape(w);
    } catch (...) { return nullptr; }
}

Shape mucad_rectangle(double x1, double y1, double x2, double y2) {
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

Shape mucad_triangle(double x1, double y1, double x2, double y2, double x3, double y3) {
    try {
        gp_Pnt p1(x1,y1,0), p2(x2,y2,0), p3(x3,y3,0);
        BRepBuilderAPI_MakeWire mkWire;
        mkWire.Add(BRepBuilderAPI_MakeEdge(p1,p2).Edge());
        mkWire.Add(BRepBuilderAPI_MakeEdge(p2,p3).Edge());
        mkWire.Add(BRepBuilderAPI_MakeEdge(p3,p1).Edge());
        return toShape(mkWire.Wire());
    } catch (...) { return nullptr; }
}

Shape mucad_polygon(const double* pts, size_t npts) {
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
Shape mucad_sphere(double cx, double cy, double cz, double radius) {
    try {
        gp_Pnt center(cx, cy, cz);
        return toShape(BRepPrimAPI_MakeSphere(center, radius).Shape());
    } catch (...) { return nullptr; }
}

Shape mucad_box(double x1, double y1, double z1, double x2, double y2, double z2) {
    try {
        gp_Pnt p1(x1, y1, z1);
        gp_Pnt p2(x2, y2, z2);
        return toShape(BRepPrimAPI_MakeBox(p1, p2).Shape());
    } catch (...) { return nullptr; }
}

// ----------------------------------------------------------------
// Transformations
// ----------------------------------------------------------------
Shape mucad_translate(Shape s, double dx, double dy, double dz) {
    if (!s) return nullptr;
    try {
        gp_Trsf trsf; 
        trsf.SetTranslation(gp_Vec(dx, dy, dz));
        BRepBuilderAPI_Transform transf(fromShape(s), trsf);
        return toShape(transf.Shape());
    } catch (...) { return nullptr; }
}

Shape mucad_rotate(Shape s, double ax, double ay, double az, 
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

Shape mucad_transform(Shape s, const double mat[16])
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

Shape mucad_scale(Shape s, double sx, double sy, double sz)
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

Shape mucad_mirror(Shape s, double ax, double ay, double az, 
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
    } catch (...) {
        return nullptr;
    }
}

// ----------------------------------------------------------------
// Operations: Extrude, Sweep, Revolve, Loft
// ----------------------------------------------------------------
Shape mucad_extrude(Shape shape2d, double ux, double uy, double uz) {
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

Shape mucad_sweep(Shape profile, double x0, double y0, double z0, double x1, double y1, double z1) {
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
        
        if (!pipe.Build()) return nullptr;
        return toShape(pipe.Shape());
    } catch (...) { return nullptr; }
}

Shape mucad_revolve(Shape s, double ax, double ay, double az, 
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

// // 
// Shape mucad_loft(const Shape* shapes, size_t n, bool solid, bool ruled) {
//     if (!shapes || n < 2) return nullptr;
//     try {
//         // BRepOffsetAPI_ThruSections(isSolid, isRuled, precision)
//         BRepOffsetAPI_ThruSections generator(solid, ruled, 1.0e-06);

//         for (size_t i = 0; i < n; ++i) {
//             const TopoDS_Shape& s = fromShape(shapes[i]);
            
//             if (s.ShapeType() == TopAbs_WIRE) {
//                 generator.AddWire(TopoDS::Wire(s));
//             } 
//             else if (s.ShapeType() == TopAbs_VERTEX) {
//                 generator.AddVertex(TopoDS::Vertex(s));
//             }
//             else {
//                 // Try to find a wire inside (e.g. if a Face was passed)
//                 TopExp_Explorer exp(s, TopAbs_WIRE);
//                 if (exp.More()) {
//                     generator.AddWire(TopoDS::Wire(exp.Current()));
//                 }
//             }
//         }

//         generator.CheckCompatibility(Standard_False); // Optional: Helps if wires have different edge counts
//         generator.Build();
        
//         if (!generator.IsDone()) return nullptr;
//         return toShape(generator.Shape());
//     } catch (...) { return nullptr; }
// }

// Shape mucad_loft(const Shape* shapes, size_t n, bool solid, bool ruled)
// {
//     if (!shapes || n < 2) return nullptr;

//     try {
//         // ---------------------------------------------------------
//         // Step 1: Determine Global Resolution (Target Points)
//         // ---------------------------------------------------------
//         Standard_Integer maxVerts = 0;
//         std::vector<TopoDS_Wire> inputWires;
//         inputWires.reserve(n);

//         for (size_t i = 0; i < n; ++i) {
//             const TopoDS_Shape& s = fromShape(shapes[i]);
//             if (s.ShapeType() == TopAbs_WIRE) {
//                 inputWires.push_back(TopoDS::Wire(s));
//             } else {
//                 // Try to extract wire
//                 TopExp_Explorer exp(s, TopAbs_WIRE);
//                 if (exp.More()) inputWires.push_back(TopoDS::Wire(exp.Current()));
//                 else return nullptr; // Invalid input
//             }

//             // Count vertices
//             TopExp_Explorer vExp(inputWires.back(), TopAbs_VERTEX);
//             Standard_Integer vCount = 0;
//             for (; vExp.More(); vExp.Next()) vCount++;
//             if (vCount > maxVerts) maxVerts = vCount;
//         }

//         // Heuristic: 3x the max complexity, but at least 12 points for a circle
//         Standard_Integer N = std::max(12, maxVerts * 3);

//         // ---------------------------------------------------------
//         // Step 2: Resample All Wires into Point Clouds
//         // ---------------------------------------------------------
//         std::vector<std::vector<gp_Pnt>> pointGrid(n);

//         for (size_t i = 0; i < n; ++i) {
//             BRepAdaptor_CompCurve curveAdaptor(inputWires[i]);
            
//             // GCPnts_UniformAbscissa handles complex wires efficiently
//             GCPnts_UniformAbscissa uniform(curveAdaptor, N, Precision::Confusion());
            
//             if (!uniform.IsDone()) return nullptr; // Sampling failed
            
//             pointGrid[i].reserve(N);
//             for (Standard_Integer p = 1; p <= N; ++p) {
//                 // Note: GCPnts is 1-based
//                 Standard_Real param = uniform.Parameter(p);
//                 pointGrid[i].push_back(curveAdaptor.Value(param));
//             }
//         }

//         // ---------------------------------------------------------
//         // Step 3: Align Clouds (Minimize Twist)
//         // ---------------------------------------------------------
//         // Rotate grid[i] to match grid[i-1]
//         for (size_t i = 1; i < n; ++i) {
//             const auto& prev = pointGrid[i-1];
//             auto& curr = pointGrid[i];

//             Standard_Real minDist = Precision::Infinite();
//             Standard_Integer bestShift = 0;

//             // Check all N possible rotations
//             for (Standard_Integer shift = 0; shift < N; ++shift) {
//                 Standard_Real distSum = 0.0;
//                 // Sample points for speed (check every 2nd point)
//                 for (Standard_Integer k = 0; k < N; k += 2) {
//                     const gp_Pnt& pA = prev[k];
//                     const gp_Pnt& pB = curr[(k + shift) % N];
//                     distSum += pA.SquareDistance(pB);
//                 }

//                 if (distSum < minDist) {
//                     minDist = distSum;
//                     bestShift = shift;
//                 }
//             }

//             // Apply Shift
//             if (bestShift != 0) {
//                 std::vector<gp_Pnt> temp = curr;
//                 for (Standard_Integer k = 0; k < N; ++k) {
//                     curr[k] = temp[(k + bestShift) % N];
//                 }
//             }
//         }

//         // ---------------------------------------------------------
//         // Step 4: Generate Topology (Vertices)
//         // ---------------------------------------------------------
//         // We need explicit TopoDS_Vertices to ensure the patch edges connect
//         std::vector<std::vector<TopoDS_Vertex>> vertexGrid(n);
//         for (size_t i = 0; i < n; ++i) {
//             vertexGrid[i].reserve(N);
//             for (Standard_Integer j = 0; j < N; ++j) {
//                 vertexGrid[i].push_back(BRepBuilderAPI_MakeVertex(pointGrid[i][j]));
//             }
//         }

//         // ---------------------------------------------------------
//         // Step 5: Build Surfaces (Filling / Patching)
//         // ---------------------------------------------------------
//         BRepBuilderAPI_Sewing sewer(1.0e-6);

//         // Loop through slices
//         for (size_t i = 0; i < n - 1; ++i) {
//             // Loop through segments around the ring
//             for (Standard_Integer j = 0; j < N; ++j) {
//                 Standard_Integer nextJ = (j + 1) % N;

//                 // The 4 corners of the patch
//                 TopoDS_Vertex v00 = vertexGrid[i][j];
//                 TopoDS_Vertex v01 = vertexGrid[i][nextJ];
//                 TopoDS_Vertex v10 = vertexGrid[i+1][j];
//                 TopoDS_Vertex v11 = vertexGrid[i+1][nextJ];

//                 // Build 4 Edges sharing these exact vertices
//                 TopoDS_Edge eBottom = BRepBuilderAPI_MakeEdge(v00, v01);
//                 TopoDS_Edge eTop    = BRepBuilderAPI_MakeEdge(v10, v11);
//                 TopoDS_Edge eLeft   = BRepBuilderAPI_MakeEdge(v00, v10);
//                 TopoDS_Edge eRight  = BRepBuilderAPI_MakeEdge(v01, v11);

//                 if (ruled) {
//                     // Simple ruled face (linear interpolation)
//                     BRepBuilderAPI_MakeWire mkW(eBottom, eRight, eTop, eLeft); // Order matters?
//                     // Usually better to build a Filling even for ruled to ensure closure
//                     BRepOffsetAPI_MakeFilling fill(1, 4, 2, Standard_False, 1e-5, 1e-4, 0.1, 0.01);
//                     fill.Add(eBottom, GeomAbs_C0);
//                     fill.Add(eTop,    GeomAbs_C0);
//                     fill.Add(eLeft,   GeomAbs_C0);
//                     fill.Add(eRight,  GeomAbs_C0);
//                     fill.Build();
//                     if (fill.IsDone()) sewer.Add(fill.Shape());
//                 } else {
//                     // Smooth Filling (Coons Patch / Plate)
//                     // Degree 2 (Quadratic) is robust
//                     BRepOffsetAPI_MakeFilling fill(2, 10, 2, Standard_False, 1e-5, 1e-4, 0.1, 0.01);
//                     fill.Add(eBottom, GeomAbs_C0);
//                     fill.Add(eTop,    GeomAbs_C0);
//                     fill.Add(eLeft,   GeomAbs_C0);
//                     fill.Add(eRight,  GeomAbs_C0);
//                     fill.Build();
//                     if (fill.IsDone()) sewer.Add(fill.Shape());
//                 }
//             }
//         }

//         // ---------------------------------------------------------
//         // Step 6: Cap Ends (if Solid)
//         // ---------------------------------------------------------
//         if (solid) {
//             // Cap Start
//             BRepBuilderAPI_MakeWire startWire;
//             for (Standard_Integer j = 0; j < N; ++j) {
//                 startWire.Add(BRepBuilderAPI_MakeEdge(vertexGrid[0][j], vertexGrid[0][(j+1)%N]));
//             }
//             // Try planar first, then filling
//             BRepBuilderAPI_MakeFace mkStart(startWire.Wire(), true);
//             if (mkStart.IsDone()) sewer.Add(mkStart.Face());
//             else sewer.Add(occt_fill_face(toShape(startWire.Wire()))); // Use our existing fill helper

//             // Cap End
//             BRepBuilderAPI_MakeWire endWire;
//             for (Standard_Integer j = 0; j < N; ++j) {
//                 endWire.Add(BRepBuilderAPI_MakeEdge(vertexGrid[n-1][j], vertexGrid[n-1][(j+1)%N]));
//             }
//             BRepBuilderAPI_MakeFace mkEnd(endWire.Wire(), true);
//             if (mkEnd.IsDone()) sewer.Add(mkEnd.Face());
//             else sewer.Add(occt_fill_face(toShape(endWire.Wire())));
//         }

//         // ---------------------------------------------------------
//         // Step 7: Final Assembly & Simplification
//         // ---------------------------------------------------------
//         sewer.Perform();
//         TopoDS_Shape sewed = sewer.SewedShape();

//         TopoDS_Shape result = sewed;

//         if (solid) {
//             // Convert Shell to Solid
//             TopExp_Explorer shellExp(sewed, TopAbs_SHELL);
//             if (shellExp.More()) {
//                 TopoDS_Shell sh = TopoDS::Shell(shellExp.Current());
//                 BRepBuilderAPI_MakeSolid mkSolid(sh);
//                 if (mkSolid.IsDone()) result = mkSolid.Solid();
//             }
            
//             // Fix Orientation
//             ShapeFix_Solid fixer;
//             if (result.ShapeType() == TopAbs_SOLID) {
//                  fixer.Init(TopoDS::Solid(result));
//                  fixer.Perform();
//                  result = fixer.Solid();
//             }
//         }

//         // Optimization: Remove "unnecessary points" by merging faces
//         // This collapses the NxM grid faces back into single smooth faces where possible
//         ShapeUpgrade_UnifySameDomain simplifier(result, Standard_True, Standard_True, Standard_True);
//         simplifier.Build();
        
//         return toShape(simplifier.Shape());

//     } catch (...) {
//         return nullptr;
//     }
// }


// ----------------------------------------------------------------
// Helper: Create a Quadratic (Degree 2) B-Spline Edge from 3 points
// This ensures curvature matching as requested.
// ----------------------------------------------------------------
TopoDS_Edge make_quadratic_edge(const gp_Pnt& p1, const gp_Pnt& pMid, const gp_Pnt& p2)
{
    // To pass EXACTLY through pMid with a rational quadratic spline,
    // we compute the middle Control Pole (C1).
    // Formula: C1 = 2*Pmid - 0.5*Pstart - 0.5*Pend
    gp_Pnt C1;
    C1.SetX(2.0 * pMid.X() - 0.5 * p1.X() - 0.5 * p2.X());
    C1.SetY(2.0 * pMid.Y() - 0.5 * p1.Y() - 0.5 * p2.Y());
    C1.SetZ(2.0 * pMid.Z() - 0.5 * p1.Z() - 0.5 * p2.Z());

    TColgp_Array1OfPnt poles(1, 3);
    poles(1) = p1; poles(2) = C1; poles(3) = p2;

    TColStd_Array1OfReal knots(1, 2);
    knots(1) = 0.0; knots(2) = 1.0;

    TColStd_Array1OfInteger mults(1, 2);
    mults(1) = 3; mults(2) = 3;

    Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, mults, 2);
    return BRepBuilderAPI_MakeEdge(curve);
}

// ----------------------------------------------------------------
// Helper: Structure to hold a sampled shape
// ----------------------------------------------------------------
struct SampledWire {
    TopoDS_Wire wire;
    std::vector<gp_Pnt> cloud;      // Dense points for matching
    std::vector<double> params;     // Parameters on the CompCurve
    BRepAdaptor_CompCurve adaptor;  // For evaluating points later
};

// ----------------------------------------------------------------
// Loft Implementation
// ----------------------------------------------------------------
Shape mucad_loft_robust(const Shape* shapes, size_t n, bool solid, bool ruled)
{
    if (!shapes || n < 2) return nullptr;

    try {
        std::vector<SampledWire> wires(n);
        Standard_Integer maxVerts = 0;

        // -------------------------------------------------
        // 1. Pre-process: Build Dense Clouds (Vertices + Uniform)
        // -------------------------------------------------
        for (size_t i = 0; i < n; ++i) {
            const TopoDS_Shape& s = fromShape(shapes[i]);
            if (s.ShapeType() == TopAbs_WIRE) wires[i].wire = TopoDS::Wire(s);
            else {
                TopExp_Explorer exp(s, TopAbs_WIRE);
                if (exp.More()) wires[i].wire = TopoDS::Wire(exp.Current());
                else return nullptr;
            }

            wires[i].adaptor.Initialize(wires[i].wire);
            
            // Count Vertices to determine density
            int vCount = 0;
            TopExp_Explorer exp(wires[i].wire, TopAbs_VERTEX);
            for(; exp.More(); exp.Next()) vCount++;
            if (vCount > maxVerts) maxVerts = vCount;
        }

        // Target resolution: High enough to catch curvature, but used only for finding matches.
        // 10x max vertices ensures we have plenty of resolution to find a "nearest" point.
        int cloudSize = std::max(50, maxVerts * 10);

        for (size_t i = 0; i < n; ++i) {
            SampledWire& w = wires[i];
            w.cloud.reserve(cloudSize);
            w.params.reserve(cloudSize);

            // Uniform Abscissa on the whole wire
            GCPnts_UniformAbscissa uniform(w.adaptor, cloudSize, Precision::Confusion());
            
            if (uniform.IsDone()) {
                for (int p = 1; p <= cloudSize; ++p) {
                    double u = uniform.Parameter(p);
                    w.params.push_back(u);
                    w.cloud.push_back(w.adaptor.Value(u));
                }
            } else {
                // Fallback (should not happen often)
                double first = w.adaptor.FirstParameter();
                double last  = w.adaptor.LastParameter();
                for (int p = 0; p < cloudSize; ++p) {
                    double u = first + (last - first) * ((double)p/(cloudSize-1));
                    w.params.push_back(u);
                    w.cloud.push_back(w.adaptor.Value(u));
                }
            }
        }

        BRepBuilderAPI_Sewing sewer(1.0e-6);

        // -------------------------------------------------
        // 2. Lofting Loop
        // -------------------------------------------------
        for (size_t i = 0; i < n - 1; ++i) {
            SampledWire& w1 = wires[i];
            SampledWire& w2 = wires[i+1];

            // A. Extract Vertices from W1 (The Features)
            //    We strictly want to connect these to W2.
            std::vector<gp_Pnt> v1_pnts;
            std::vector<double> v1_params_on_w1;

            // Use WireExplorer to get vertices in order
            BRepTools_WireExplorer exp(w1.wire);
            for (; exp.More(); exp.Next()) {
                // We map the START of every edge
                TopoDS_Vertex v = exp.CurrentVertex();
                gp_Pnt p = BRep_Tool::Pnt(v);
                v1_pnts.push_back(p);
                
                // Find param on w1 adaptor (approximate via projection for safety)
                double u = 0.0;
                // (Optimization: In a real CompCurve, vertices are at knot intervals, 
                // but projection is safer if orientation is tricky)
                // Simple scan of our own cloud to find param:
                double minDist = Precision::Infinite();
                for (size_t k=0; k<w1.cloud.size(); ++k) {
                    double d = p.SquareDistance(w1.cloud[k]);
                    if(d < minDist) { minDist = d; u = w1.params[k]; }
                }
                v1_params_on_w1.push_back(u);
            }

            // B. Match Vertices of W1 to Cloud of W2
            std::vector<int> w2_indices;
            std::vector<double> w2_params;
            
            // Anti-twist logic:
            // We calculate the "Shift" that minimizes total distance of Vertices(W1) -> Cloud(W2)
            int bestShift = 0;
            double minTotalDist = Precision::Infinite();
            int N1 = (int)v1_pnts.size();
            int N2 = (int)w2.cloud.size();

            // We assume W1 vertices are sparse, W2 cloud is dense.
            // We check shifts in the W2 cloud indices.
            // Since N1 << N2, we step through N2.
            
            // To save time, check every (N2/N1) points roughly
            int step = std::max(1, N2 / N1);
            
            for (int shift = 0; shift < N2; shift += step) {
                double dist = 0;
                for (int k = 0; k < N1; ++k) {
                    // Where does vertex k map to?
                    // We distribute k evenly over the circle: index ~ (k * N2 / N1)
                    int targetIdx = (shift + (k * N2 / N1)) % N2;
                    dist += v1_pnts[k].SquareDistance(w2.cloud[targetIdx]);
                }
                if (dist < minTotalDist) {
                    minTotalDist = dist;
                    bestShift = shift;
                }
            }

            // Now Perform the actual Nearest Neighbor search constrained by that general alignment.
            // This prevents a vertex "jumping" backward.
            
            int currentSearchIdx = bestShift;
            for (int k = 0; k < N1; ++k) {
                // Search local neighborhood in W2 cloud starting from currentSearchIdx
                // We wrap around N2.
                double locMin = Precision::Infinite();
                int bestLocIdx = -1;
                
                // Search forward up to (N2/N1 * 1.5) points to allow flex but not folding
                int searchRange = (N2 / N1) * 2; 
                
                for (int s = 0; s < searchRange; ++s) {
                    int idx = (currentSearchIdx + s) % N2;
                    double d = v1_pnts[k].SquareDistance(w2.cloud[idx]);
                    if (d < locMin) {
                        locMin = d;
                        bestLocIdx = idx;
                    }
                }
                w2_indices.push_back(bestLocIdx);
                w2_params.push_back(w2.params[bestLocIdx]);
                
                // Advance search start so next vertex must be 'after' this one
                currentSearchIdx = bestLocIdx; 
            }

            // C. Build Patches Segment by Segment
            //    Segment k is between Vertex k and Vertex k+1
            for (int k = 0; k < N1; ++k) {
                int nextK = (k + 1) % N1;

                // 1. Get 4 Corners
                gp_Pnt p00 = v1_pnts[k];            // W1 Start
                gp_Pnt p10 = v1_pnts[nextK];        // W1 End
                
                gp_Pnt p01 = w2.cloud[w2_indices[k]];     // W2 Start
                gp_Pnt p11 = w2.cloud[w2_indices[nextK]]; // W2 End

                // 2. Construct Boundaries
                //    Edge 1 (Bottom): Actual edge of W1 (or approximation)
                //    To be robust, we rebuild it as a Quadratic Spline using the Midpoint from W1 parameters
                double u1_start = v1_params_on_w1[k];
                double u1_end   = v1_params_on_w1[nextK];
                if (u1_end < u1_start) u1_end += (w1.adaptor.LastParameter() - w1.adaptor.FirstParameter()); // Wrap
                double u1_mid = (u1_start + u1_end) * 0.5;
                // Wrap mid param back to range if needed for evaluation
                // (BRepAdaptor handles period, but let's be safe)
                gp_Pnt pMid_Bottom = w1.adaptor.Value(u1_mid);

                //    Edge 2 (Top): Segment of W2
                double u2_start = w2_params[k];
                double u2_end   = w2_params[nextK];
                if (u2_end < u2_start) u2_end += (w2.adaptor.LastParameter() - w2.adaptor.FirstParameter());
                double u2_mid = (u2_start + u2_end) * 0.5;
                gp_Pnt pMid_Top = w2.adaptor.Value(u2_mid);

                //    Edge 3/4 (Sides): Straight lines connecting W1->W2
                gp_Pnt pMid_Left  = p00.Translated(gp_Vec(p00, p01) * 0.5);
                gp_Pnt pMid_Right = p10.Translated(gp_Vec(p10, p11) * 0.5);

                TopoDS_Edge eBottom = make_quadratic_edge(p00, pMid_Bottom, p10);
                TopoDS_Edge eTop    = make_quadratic_edge(p01, pMid_Top,    p11);
                TopoDS_Edge eLeft   = ruled ? BRepBuilderAPI_MakeEdge(p00, p01) : make_quadratic_edge(p00, pMid_Left, p01);
                TopoDS_Edge eRight  = ruled ? BRepBuilderAPI_MakeEdge(p10, p11) : make_quadratic_edge(p10, pMid_Right, p11);

                // 3. Create Surface
                BRepOffsetAPI_MakeFilling fill(2, 10, 2, Standard_False, 1e-5, 1e-4, 0.1, 0.01);
                fill.Add(eBottom, GeomAbs_C0);
                fill.Add(eTop,    GeomAbs_C0);
                fill.Add(eLeft,   GeomAbs_C0);
                fill.Add(eRight,  GeomAbs_C0);
                
                fill.Build();
                if (fill.IsDone()) sewer.Add(fill.Shape());
            }
        }

        // -------------------------------------------------
        // 3. Capping (if Solid)
        // -------------------------------------------------
        if (solid) {
            // We reuse the wire/filling helper logic
            TopoDS_Wire startW = wires[0].wire;
            TopoDS_Wire endW   = wires[n-1].wire;
            
            BRepBuilderAPI_MakeFace mkStart(startW, true);
            if(mkStart.IsDone()) sewer.Add(mkStart.Face());
            else sewer.Add(occt_fill_face(toShape(startW)));

            BRepBuilderAPI_MakeFace mkEnd(endW, true);
            if(mkEnd.IsDone()) sewer.Add(mkEnd.Face());
            else sewer.Add(occt_fill_face(toShape(endW)));
        }

        // -------------------------------------------------
        // 4. Finalize
        // -------------------------------------------------
        sewer.Perform();
        TopoDS_Shape sewed = sewer.SewedShape();
        
        TopoDS_Shape result = sewed;
        if (solid) {
            TopExp_Explorer shExp(sewed, TopAbs_SHELL);
            if(shExp.More()) {
                 BRepBuilderAPI_MakeSolid mkS(TopoDS::Shell(shExp.Current()));
                 if(mkS.IsDone()) result = mkS.Solid();
            }
            ShapeFix_Solid fix; fix.Init(TopoDS::Solid(result)); fix.Perform();
            result = fix.Solid();
        }

        // "Remove all unnecessary points" -> Unify Same Domain
        ShapeUpgrade_UnifySameDomain unifier(result, Standard_True, Standard_True, Standard_True);
        unifier.Build();
        
        return toShape(unifier.Shape());

    } catch (...) {
        return nullptr;
    }
}

// ----------------------------------------------------------------
// Booleans
// ----------------------------------------------------------------
Shape mucad_union(Shape a, Shape b) {
    if (!a || !b) return nullptr;
    try {
        return toShape(BRepAlgoAPI_Fuse(fromShape(a), fromShape(b)).Shape());
    } catch (...) { return nullptr; }
}

Shape mucad_difference(Shape a, Shape b) {
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
Shape mucad_quadratic_spline_wire(const double* pts, size_t npts)
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

            // Calculate the intermediate control pole for exact interpolation
            // Based on the property of a parabola at t=0.5
            gp_Pnt C1;
            C1.SetX(2.0 * Pmid.X() - 0.5 * P0.X() - 0.5 * P2.X());
            C1.SetY(2.0 * Pmid.Y() - 0.5 * P0.Y() - 0.5 * P2.Y());
            C1.SetZ(2.0 * Pmid.Z() - 0.5 * P0.Z() - 0.5 * P2.Z());

            // Setup B-Spline Arrays
            TColgp_Array1OfPnt poles(1, 3);
            poles(1) = P0;
            poles(2) = C1; // The calculated pole, not the mid point!
            poles(3) = P2;

            TColStd_Array1OfReal knots(1, 2);
            knots(1) = 0.0;
            knots(2) = 1.0;

            TColStd_Array1OfInteger mults(1, 2);
            mults(1) = 3; // Multiplicity = Degree + 1 at ends
            mults(2) = 3;

            Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(
                poles, knots, mults, 2 // Degree 2 (Quadratic)
            );

            TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve).Edge();
            mkWire.Add(edge);
        }

        if (!mkWire.IsDone()) return nullptr;
        return toShape(mkWire.Wire());

    } catch (...) {
        return nullptr;
    }
}

Shape mucad_face_from_wire(Shape wire) {
    if (!wire) return nullptr;
    try {
        return toShape(BRepBuilderAPI_MakeFace(TopoDS::Wire(fromShape(wire))).Face());
    } catch (...) { return nullptr; }
}

// ----------------------------------------------------------------
// Analysis
// ----------------------------------------------------------------
int mucad_bounding_sphere(Shape s, double* x, double* y, double* z, double* r)
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

double mucad_volume(Shape s)
{
    if (!s) return 0.0;
    try {
        GProp_GProps props;
        BRepGProp::VolumeProperties(fromShape(s), props);
        return props.Mass(); // Mass = Volume when density is 1 (default)
    } catch (...) { return 0.0; }
}

double mucad_surface_area(Shape s)
{
    if (!s) return 0.0;
    try {
        GProp_GProps props;
        BRepGProp::SurfaceProperties(fromShape(s), props);
        return props.Mass(); // Mass = Area for surface props
    } catch (...) { return 0.0; }
}

int mucad_centroid(Shape s, double* cx, double* cy, double* cz)
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

int mucad_inertia(Shape s, double density, 
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
        // I_mass = I_geometric * density
        *Ixx = matrix(1, 1) * density;
        *Iyy = matrix(2, 2) * density;
        *Izz = matrix(3, 3) * density;

        return 0;
    } catch (...) { return -1; }
}

// ----------------------------------------------------------------
// I/O and Cleanup
// ----------------------------------------------------------------

int mucad_write_step(Shape s, const char* filename) {
    if (!s || !filename) return -1;
    try {
        STEPControl_Writer writer;
        writer.Transfer(fromShape(s), STEPControl_AsIs);
        return (writer.Write(filename) == IFSelect_RetDone) ? 0 : -1;
    } catch (...) { return -1; }
}

int mucad_write_stl(Shape s, const char* filename) {
    if (!s || !filename) return -1;
    try {
        StlAPI_Writer writer;
        writer.Write(fromShape(s), filename);
        return 0;
    } catch (...) { return -1; }
}

// Simplify functions
Shape mucad_simplify(Shape s) 
{
    if (!s) return nullptr;
    try {
        // ShapeUpgrade_UnifySameDomain is the standard tool to remove 
        // "seam" edges on planar faces or merge co-planar faces.
        // Arguments: Shape, UnifyFaces, UnifyEdges, ConcatBSplines
        ShapeUpgrade_UnifySameDomain unifier(fromShape(s), Standard_True, Standard_True, Standard_True);
        unifier.Build();
        return toShape(unifier.Shape());
    } catch (...) {
        return nullptr;
    }
}

/* --------------------------------------------------------------------
   Simplify to Solid (Sewing + Capping + Healing)
   
   This function attempts to turn a loose collection of shapes into a 
   single watertight solid.
   
   1. It sews faces together within the given tolerance.
   2. It detects open boundaries (holes).
   3. It attempts to cap planar holes (e.g., closing a pipe).
   4. It simplifies the geometry (merging coplanar faces).
   
   shapes : array of input shapes (faces, shells, etc.)
   n      : number of shapes
   tol    : tolerance for sewing (e.g., 1.0e-6 to 1.0e-3)
--------------------------------------------------------------------- */
Shape mucad_simplify_to_solid(const Shape* shapes, size_t n, double tol) 
{
    if (!shapes || n == 0) return nullptr;

    try {
        // ---------------------------------------------------------
        // Step 1: Sew the input faces together
        // ---------------------------------------------------------
        BRepBuilderAPI_Sewing sewing(tol);
        
        for (size_t i = 0; i < n; ++i) {
            const TopoDS_Shape& s = fromShape(shapes[i]);
            sewing.Add(s);
        }
        
        sewing.Perform();
        TopoDS_Shape sewedShape = sewing.SewedShape();

        // ---------------------------------------------------------
        // Step 2: Detect and Cap Holes (Planar Capping)
        // ---------------------------------------------------------
        // ShapeAnalysis_FreeBounds extracts the open wires (boundaries)
        TopoDS_Compound closedWires, openWires;
        ShapeAnalysis_FreeBounds::ConnectEdgesToWires(sewedShape, tol, Standard_False, closedWires, openWires);
        
        TopTools_ListOfShape caps;
        TopExp_Explorer wireExp(closedWires, TopAbs_WIRE);
        
        for (; wireExp.More(); wireExp.Next()) {
            const TopoDS_Wire& holeWire = TopoDS::Wire(wireExp.Current());
            
            // Try to make a planar face to fill the hole
            BRepBuilderAPI_MakeFace mkFace(holeWire, true); // true = only planar
            if (mkFace.IsDone()) {
                caps.Append(mkFace.Face());
            }
        }

        // If we found caps, we need to sew them onto the main body
        if (!caps.IsEmpty()) {
            BRepBuilderAPI_Sewing resewing(tol);
            resewing.Add(sewedShape);
            
            TopTools_ListIteratorOfListOfShape it(caps);
            for (; it.More(); it.Next()) {
                resewing.Add(it.Value());
            }
            
            resewing.Perform();
            sewedShape = resewing.SewedShape();
        }

        // ---------------------------------------------------------
        // Step 3: Convert to Solid
        // ---------------------------------------------------------
        // The sewed shape is likely a Shell or a Compound of Shells.
        TopoDS_Shape solidResult;
        
        // Helper to ensure we have a Shell
        TopoDS_Shell aShell;
        if (sewedShape.ShapeType() == TopAbs_SHELL) {
            aShell = TopoDS::Shell(sewedShape);
        } else {
            // If it's a compound, extract the shell
            TopExp_Explorer shellExp(sewedShape, TopAbs_SHELL);
            if (shellExp.More()) {
                aShell = TopoDS::Shell(shellExp.Current());
            }
        }

        if (!aShell.IsNull()) {
            BRepBuilderAPI_MakeSolid mkSolid(aShell);
            if (mkSolid.IsDone()) {
                solidResult = mkSolid.Solid();
            }
        }

        // If MakeSolid failed or didn't run, use the raw sewed shape 
        // (ShapeFix might be able to upgrade it to a solid)
        if (solidResult.IsNull()) {
            solidResult = sewedShape;
        }

        // ---------------------------------------------------------
        // Step 4: Heal and Validate (Orientation fix)
        // ---------------------------------------------------------
        ShapeFix_Shape fixer(solidResult);
        fixer.Perform();
        solidResult = fixer.Shape();

        // ---------------------------------------------------------
        // Step 5: Simplification (UnifySameDomain)
        // ---------------------------------------------------------
        // This removes internal seam edges and merges coplanar faces
        ShapeUpgrade_UnifySameDomain unifier(solidResult, Standard_True, Standard_True, Standard_True);
        unifier.Build();
        
        return toShape(unifier.Shape());

    } catch (...) {
        return nullptr;
    }
}


int mucad_write_step(Shape s, const char* filename)
{
    if (!s || !filename) return -1;
    try {
        STEPControl_Writer writer;
        // STEPControl_AsIs preserves the geometry type best
        IFSelect_ReturnStatus status = writer.Transfer(fromShape(s), STEPControl_AsIs);
        
        if (status != IFSelect_RetDone) return -1;
        
        status = writer.Write(filename);
        return (status == IFSelect_RetDone) ? 0 : -1;
    } catch (...) { return -1; }
}

int mucad_write_stl(Shape s, const char* filename)
{
    if (!s || !filename) return -1;
    try {
        StlAPI_Writer writer;
        // Force triangulation before writing to ensure STL is valid
        // (Linear deflection 0.01 is a reasonable default for visuals)
        BRepMesh_IncrementalMesh mesh(fromShape(s), 0.01);
        
        // Write binary STL (Standard_True for binary)
        Standard_Boolean success = writer.Write(fromShape(s), filename);
        return success ? 0 : -1;
    } catch (...) { return -1; }
}

int mucad_write_iges(Shape s, const char* filename)
{
    if (!s || !filename) return -1;
    try {
        IGESControl_Writer writer;
        // 0 = Metric (Millimeters), 1 = Inches. OCCT default is usually MM.
        const int unit_mm = 0; 
        writer.ChangeGlobalSection().SetUnitFlag(unit_mm);
        
        Standard_Boolean ok = writer.AddShape(fromShape(s));
        if (!ok) return -1;
        
        writer.ComputeModel();
        Standard_Boolean success = writer.Write(filename);
        return success ? 0 : -1;
    } catch (...) { return -1; }
}

int mucad_write_obj(Shape s, const char* filename)
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

        file << "# Generated by MUCAD/OCCT Wrapper\n";
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