// ccad.cpp - implements ccad.h using Open CASCADE

#include "ccad.h"

// Standard C++ headers
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <vector>

// OCCT core geometry & topology
#include <BRepBuilderAPI_FindPlane.hxx>
#include <BRepBuilderAPI_GTransform.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepCheck_Analyzer.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakePrism.hxx>
#include <BRepPrimAPI_MakeRevol.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepTools_WireExplorer.hxx>

#include <TopExp_Explorer.hxx>
#include <TopTools_HSequenceOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>
#include <TopoDS_Wire.hxx>

// OCCT boolean operations
#include <BRepAlgoAPI_Cut.hxx>
#include <BRepAlgoAPI_Fuse.hxx>

// OCCT extrusion, loft & sweep
#include <BRepFill_CompatibleWires.hxx>
#include <BRepOffsetAPI_MakeFilling.hxx>
#include <BRepOffsetAPI_MakePipeShell.hxx>
#include <BRepOffsetAPI_ThruSections.hxx> // needed for loft

// OCCT meshing & file export
#include <BRepMesh_IncrementalMesh.hxx>
#include <IGESControl_Writer.hxx>
#include <STEPControl_Writer.hxx>
#include <StlAPI_Writer.hxx>

// OCCT geometry utilities
#include <Geom_BSplineCurve.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <gp_Ax1.hxx>
#include <gp_Ax2.hxx>
#include <gp_GTrsf.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

// OCCT analysis & repair
#include <BRepAdaptor_CompCurve.hxx>
#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>
#include <Poly_Triangulation.hxx>
#include <Precision.hxx>
#include <ShapeAnalysis_FreeBounds.hxx>
#include <ShapeFix_Shape.hxx>
#include <ShapeFix_Solid.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <Standard_Boolean.hxx>


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
static inline Shape toShape(const TopoDS_Shape &s) {
	return new TopoDS_Shape(s);
}

static inline const TopoDS_Shape &fromShape(Shape s) {
	return *static_cast<const TopoDS_Shape *>(s);
}

// ----------------------------------------------------------------
// Cleanup (all returns)
// ----------------------------------------------------------------
void ccadcpp_free(Shape s) {
	if(s) {
		delete static_cast<TopoDS_Shape *>(s);
	}
}

// ----------------------------------------------------------------
// 2-D Primitives
// ----------------------------------------------------------------
Shape ccadcpp_circle(double cx, double cy, double radius) {
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
		if(!mkEdge.IsDone()) {
			return nullptr;
		}

		TopoDS_Edge edge = mkEdge.Edge();

		// 4. Wrap in a Wire (to match your original intent)
		BRepBuilderAPI_MakeWire mkWire(edge);

		if(!mkWire.IsDone()) {
			return nullptr;
		}

		TopoDS_Wire w = mkWire.Wire();

		return toShape(w);
	} catch(...) { return nullptr; }
}

Shape ccadcpp_rectangle(double x1, double y1, double x2, double y2) {
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
	} catch(...) { return nullptr; }
}

Shape ccadcpp_triangle(double x1, double y1, double x2, double y2, double x3, double y3) {
	try {
		gp_Pnt p1(x1, y1, 0), p2(x2, y2, 0), p3(x3, y3, 0);
		BRepBuilderAPI_MakeWire mkWire;
		mkWire.Add(BRepBuilderAPI_MakeEdge(p1, p2).Edge());
		mkWire.Add(BRepBuilderAPI_MakeEdge(p2, p3).Edge());
		mkWire.Add(BRepBuilderAPI_MakeEdge(p3, p1).Edge());
		return toShape(mkWire.Wire());
	} catch(...) { return nullptr; }
}

Shape ccadcpp_polygon(const double *pts, size_t npts) {
	if(npts < 3) {
		return nullptr;
	}
	try {
		BRepBuilderAPI_MakeWire mkWire;
		for(size_t i = 0; i < npts; ++i) {
			gp_Pnt p1(pts[2 * i], pts[2 * i + 1], 0.0);
			gp_Pnt p2(pts[2 * ((i + 1) % npts)], pts[2 * ((i + 1) % npts) + 1], 0.0);
			mkWire.Add(BRepBuilderAPI_MakeEdge(p1, p2).Edge());
		}
		return toShape(mkWire.Wire());
	} catch(...) { return nullptr; }
}

// ----------------------------------------------------------------
// 3-D Primitives
// ----------------------------------------------------------------
Shape ccadcpp_sphere(double cx, double cy, double cz, double radius) {
	try {
		gp_Pnt center(cx, cy, cz);
		return toShape(BRepPrimAPI_MakeSphere(center, radius).Shape());
	} catch(...) { return nullptr; }
}

Shape ccadcpp_box(double x1, double y1, double z1, double x2, double y2, double z2) {
	try {
		gp_Pnt p1(x1, y1, z1);
		gp_Pnt p2(x2, y2, z2);
		return toShape(BRepPrimAPI_MakeBox(p1, p2).Shape());
	} catch(...) { return nullptr; }
}

// ----------------------------------------------------------------
// Transformations
// ----------------------------------------------------------------
Shape ccadcpp_translate(Shape s, double dx, double dy, double dz) {
	if(!s) {
		return nullptr;
	}
	try {
		gp_Trsf trsf;
		trsf.SetTranslation(gp_Vec(dx, dy, dz));
		BRepBuilderAPI_Transform transf(fromShape(s), trsf);
		return toShape(transf.Shape());
	} catch(...) { return nullptr; }
}

Shape ccadcpp_rotate(Shape s, double ax, double ay, double az,
		      double ux, double uy, double uz, double angle) {
	if(!s) {
		return nullptr;
	}
	try {
		gp_Ax1 axis(gp_Pnt(ax, ay, az), gp_Dir(ux, uy, uz));
		gp_Trsf trsf;
		trsf.SetRotation(axis, angle);
		BRepBuilderAPI_Transform transf(fromShape(s), trsf);
		return toShape(transf.Shape());
	} catch(...) { return nullptr; }
}

Shape ccadcpp_transform(Shape s, const double mat[16]) {
	if(!s || !mat) {
		return nullptr;
	}
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

	} catch(...) {
		return nullptr;
	}
}

Shape ccadcpp_scale(Shape s, double sx, double sy, double sz) {
	if(!s) {
		return nullptr;
	}
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
	} catch(...) {
		return nullptr;
	}
}

Shape ccadcpp_mirror(Shape s, double ax, double ay, double az,
		      double ux, double uy, double uz) {
	if(!s) {
		return nullptr;
	}
	try {
		gp_Pnt p(ax, ay, az);
		gp_Dir d(ux, uy, uz);
		gp_Ax2 plane(p, d); // Defines a coordinate system (Main Axis = Normal)

		gp_Trsf trsf;
		trsf.SetMirror(plane); // Mirror across the XY plane of the given axis

		// Standard Transform is sufficient (Mirror is a rigid body transform)
		BRepBuilderAPI_Transform applier(fromShape(s), trsf);

		return toShape(applier.Shape());
	} catch(...) { return nullptr; }
}

// ----------------------------------------------------------------
// Operations: Extrude, Sweep, Revolve, Loft
// ----------------------------------------------------------------
Shape ccadcpp_extrude(Shape shape2d, double ux, double uy, double uz) {
	if(!shape2d) {
		return nullptr;
	}
	try {
		const TopoDS_Shape &w = fromShape(shape2d);
		TopoDS_Face f;

		if(w.ShapeType() == TopAbs_FACE) {
			f = TopoDS::Face(w);
		} else if(w.ShapeType() == TopAbs_WIRE) {
			// BRepBuilderAPI_MakeFace assumes wire is planar
			f = BRepBuilderAPI_MakeFace(TopoDS::Wire(w)).Face();
		} else {
			return nullptr;
		}

		gp_Vec dir(ux, uy, uz);
		return toShape(BRepPrimAPI_MakePrism(f, dir).Shape());
	} catch(...) { return nullptr; }
}

//// Sweep a profile along an arbitrary path shape
//  profile : the cross-section to sweep (must be a valid shape)
//  path    : the spine of the sweep - must be a wire or a shape that can be
//            interpreted as a wire (e.g. a single edge or a collection of
//            edges).  If the input is not a wire, we try to build one from its
//            edges.  On failure the function returns nullptr.
//
//  Returns the swept solid or nullptr on error.
//
Shape ccadcpp_sweep(Shape profile, Shape path) {
	// --- 1. Basic sanity check ------------------------------------------------
	if(!profile || !path) {
		return nullptr;
	}

	try {
		// --- 2. Convert the user-supplied path into a TopoDS_Wire ----------------
		TopoDS_Shape rawPath = fromShape(path); // generic shape
		TopoDS_Wire pathWire;

		if(rawPath.ShapeType() == TopAbs_WIRE) {
			// Path is already a wire - cast it directly
			pathWire = TopoDS::Wire(rawPath);
		} else {
			// Look for a wire nested inside the shape
			TopExp_Explorer expW(rawPath, TopAbs_WIRE);
			if(expW.More()) {
				pathWire = TopoDS::Wire(expW.Current());
			} else {
				// No wire - try to build one from all edges in the shape
				BRepBuilderAPI_MakeWire mkWire;
				for(TopExp_Explorer expE(rawPath, TopAbs_EDGE); expE.More(); expE.Next()) {
					mkWire.Add(TopoDS::Edge(expE.Current()));
				}

				if(!mkWire.IsDone() || mkWire.Wire().IsNull()) {
					return nullptr; // cannot form a valid spine
				}

				pathWire = mkWire.Wire();
			}
		}

		// --- 3. Build the pipe shell -------------------------------------------
		BRepOffsetAPI_MakePipeShell pipe(pathWire);
		pipe.Add(fromShape(profile)); // cross-section profile

		pipe.Build();
		if(!pipe.IsDone()) {
			return nullptr;
		}

		// --- 4. Return the resulting solid -------------------------------------
		return toShape(pipe.Shape());
	} catch(...) {
		return nullptr;
	}
}

Shape ccadcpp_revolve(Shape s, double ax, double ay, double az,
		       double ux, double uy, double uz, double angle) {
	if(!s) {
		return nullptr;
	}
	try {
		gp_Pnt p(ax, ay, az);
		gp_Dir d(ux, uy, uz);
		gp_Ax1 axis(p, d);
		BRepPrimAPI_MakeRevol revol(fromShape(s), axis, angle);
		return toShape(revol.Shape());
	} catch(...) { return nullptr; }
}


// ----------------------------------------------------------------
// Loft Implementation
// ----------------------------------------------------------------

//// Main loft routine
// Use BRepFill_CompatibleWires. This OCCT tool automatically
// cuts edges (e.g., splits a circle into 4 arcs) to match the topology
// of the target shape (e.g., a square) while preserving sharp vertices.
Shape ccadcpp_loft(const Shape *shapes,
		    size_t n,
		    bool solid,
		    bool ruled) {
	if(!shapes || n < 2) {
		return nullptr;
	}

	try {
		// 1. Collect input Wiree
		TopTools_SequenceOfShape inputWires;

		for(size_t i = 0; i < n; ++i) {
			TopoDS_Shape s = fromShape(shapes[i]);
			if(s.IsNull()) {
				return nullptr;
			}

			// Extract the Wire from the Shape
			if(s.ShapeType() != TopAbs_WIRE) {
				TopExp_Explorer exp(s, TopAbs_WIRE);
				if(exp.More()) {
					s = exp.Current();
				} else {
					// If no wire found (e.g. passing an Edge), try to make one
					TopExp_Explorer expE(s, TopAbs_EDGE);
					if(expE.More()) {
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

		// 2. Make Wires Compatible
		//   This logic aligns the wires and splits edges so that a Circle (1 edge)
		//   can loft to a Square (4 edges) without smoothing the Square's corne//
		BRepFill_CompatibleWires fixer(inputWires);
		fixer.SetPercent(0.1); // Precision for matching vertices
		fixer.Perform();

		TopTools_SequenceOfShape correctedWires;
		bool compatibilityMode = Standard_True;

		if(fixer.IsDone()) {
			// Use the processed wires (Circle is now split into arcs)
			correctedWires = fixer.Shape();
		} else {
			// Fallback: Use original wires if compatibility alg fails.
			// We disable compatibility check in ThruSections so it doesn't crash,
			// though the mesh might be twisted.
			correctedWires = inputWires;
			compatibilityMode = Standard_False;
		}

		// 3. Build the Loft
		// Param 3 (1e-6) is precision
		BRepOffsetAPI_ThruSections generator(solid, ruled, 1e-6);

		// If wires were successfully made compatible, we turn on the check.
		// This creates a cleaner mesh. If not, we turn it off to force creation.
		generator.CheckCompatibility(compatibilityMode);

		for(int i = 1; i <= correctedWires.Length(); i++) {
			generator.AddWire(TopoDS::Wire(correctedWires.Value(i)));
		}

		generator.Build();
		if(!generator.IsDone()) {
			return nullptr;
		}

		TopoDS_Shape result = generator.Shape();

		// 4. Fix Solid topology if requested
		if(solid && result.ShapeType() == TopAbs_SOLID) {
			ShapeFix_Solid fix;
			fix.Init(TopoDS::Solid(result));
			fix.Perform();
			result = fix.Solid();
		}

		return toShape(result);

	} catch(...) {
		return nullptr;
	}
}

// ----------------------------------------------------------------
// Booleans
// ----------------------------------------------------------------
Shape ccadcpp_union(Shape a, Shape b) {
	if(!a || !b) {
		return nullptr;
	}
	try {
		return toShape(BRepAlgoAPI_Fuse(fromShape(a), fromShape(b)).Shape());
	} catch(...) { return nullptr; }
}

Shape ccadcpp_difference(Shape a, Shape b) {
	if(!a || !b) {
		return nullptr;
	}
	try {
		return toShape(BRepAlgoAPI_Cut(fromShape(a), fromShape(b)).Shape());
	} catch(...) { return nullptr; }
}

// ----------------------------------------------------------------
// Splines & Advanced
// ----------------------------------------------------------------


//// Create a wire from 3-node quadratic segments (FEM style).
//
// The array is interpreted as:
// Segment 1: P[0] (Start), P[1] (Mid), P[2] (End)
// Segment 2: P[2] (Start), P[3] (Mid), P[4] (End)
// ...
//
// pts  : Flat array of coords {x0,y0,z0, x1,y1,z1...}
// npts : Total points. Must be odd (2*N_segments + 1).
//        E.g., 3 points for 1 segment, 5 points for 2 segments.
//  Example:
//  double nodes[] = {
//      0.0, 0.0, 0.0,   // Start
//      0.5, 0.2, 0.0,   // Mid 1 (Curve goes through here)
//      1.0, 0.0, 0.0,   // End 1 / Start 2
//      1.5, -0.2, 0.0,  // Mid 2 (Curve goes through here)
//      2.0, 0.0, 0.0    // End 2
//  };
//  // 5 points = 2 segments
//  Shape fem_wire = occt_quadratic_spline_wire(nodes, 5);
// -----------------------------------------------------------------
Shape ccadcpp_quadratic_spline_wire(const double *pts, size_t npts) {
	// We need at least 3 points (1 segment) and an odd number of points
	// to maintain the Start-Mid-End chain continuity.
	if(!pts || npts < 3 || (npts % 2 == 0)) {
		return nullptr;
	}

	try {
		BRepBuilderAPI_MakeWire mkWire;

		// We construct a Quadratic B-Spline:
		// A degree 2 curve through 3 points P0, P1, P2 is defined.
		//
		// To do this explicitly in OCCT without an approximator:
		// We calculate the geometric "Poles" (Control Points) required to make
		// the curve pass through the Mid point.
		//
		// For a rational B-Spline passing through P0(t=0), P1(t=0.5), P2(t=1):
		// The middle Control Pole (C1) is NOT P1.
		// C1 = 2*P1 - 0.5*P0 - 0.5*P2

		for(size_t i = 0; i < npts - 2; i += 2) {
			gp_Pnt P0(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]);
			gp_Pnt Pmid(pts[3 * (i + 1)], pts[3 * (i + 1) + 1], pts[3 * (i + 1) + 2]);
			gp_Pnt P2(pts[3 * (i + 2)], pts[3 * (i + 2) + 1], pts[3 * (i + 2) + 2]);

			gp_Pnt C1;
			C1.SetX(2.0 * Pmid.X() - 0.5 * P0.X() - 0.5 * P2.X());
			C1.SetY(2.0 * Pmid.Y() - 0.5 * P0.Y() - 0.5 * P2.Y());
			C1.SetZ(2.0 * Pmid.Z() - 0.5 * P0.Z() - 0.5 * P2.Z());

			TColgp_Array1OfPnt poles(1, 3);
			poles(1) = P0;
			poles(2) = C1;
			poles(3) = P2;

			TColStd_Array1OfReal knots(1, 2);
			knots(1) = 0.0;
			knots(2) = 1.0;

			TColStd_Array1OfInteger mults(1, 2);
			mults(1) = 3;
			mults(2) = 3;

			Handle(Geom_Curve) curve = new Geom_BSplineCurve(poles, knots, mults, 2);

			TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve).Edge();
			mkWire.Add(edge);
		}

		if(!mkWire.IsDone()) {
			return nullptr;
		}
		return toShape(mkWire.Wire());

	} catch(...) {
		return nullptr;
	}
}

Shape ccadcpp_face_from_wire(Shape wire) {
	if(!wire) {
		return nullptr;
	}
	try {
		return toShape(BRepBuilderAPI_MakeFace(TopoDS::Wire(fromShape(wire))).Face());
	} catch(...) { return nullptr; }
}


// ----------------------------------------------------------------
// Shape Cleanup
// ----------------------------------------------------------------

// HELPER: Validate Shape Validity
// Returns true only if OCCT thinks the geometry is mathematically sound
bool IsShapeValid(const TopoDS_Shape &s) {
	if(s.IsNull()) {
		return false;
	}
	BRepCheck_Analyzer checker(s);
	return checker.IsValid();
}

// HELPER: Check if Wire is Safe to Cap
bool IsWireSafe(const TopoDS_Wire &w) {
	if(w.IsNull()) {
		return false;
	}

	// 1. Check basic topology validity
	if(!IsShapeValid(w)) {
		return false;
	}

	// 2. Check for Self-Intersection (Major cause of crashes)
	// We need a face to run ShapeAnalysis_Wire, so we make a dummy plane or use null
	Handle(ShapeAnalysis_Wire) analyzer = new ShapeAnalysis_Wire();
	analyzer->Load(w);
	analyzer->CheckSelfIntersection();
	if(analyzer->StatusSelfIntersection(ShapeExtend_FAIL)) {
		return false;
	}

	return true;
}

// HELPER: Fill Hole
TopoDS_Shape FillHole(const TopoDS_Wire &wire) {
	// 1. Sanity Check
	if(!IsWireSafe(wire)) {
		return TopoDS_Shape();
	}

	TopoDS_Shape resultFace;

	try {
		// 2. Try Planar Cap (Fastest, Safest)
		BRepBuilderAPI_FindPlane planeFinder(wire, 1.0e-3);
		if(planeFinder.Found()) {
			BRepBuilderAPI_MakeFace mkFace(planeFinder.Plane(), wire);
			if(mkFace.IsDone()) {
				resultFace = mkFace.Face();
			}
		}

		// 3. Try Ruled Cap (Tent) if Planar failed
		if(resultFace.IsNull()) {
			// Calculate Centroid
			GProp_GProps props;
			BRepGProp::LinearProperties(wire, props);
			gp_Pnt centroid = props.CentreOfMass();
			TopoDS_Vertex vCenter = BRepBuilderAPI_MakeVertex(centroid).Vertex();

			// Loft Wire -> Center
			BRepOffsetAPI_ThruSections lofter(Standard_False, Standard_True, 1.0e-3);
			lofter.AddWire(wire);
			lofter.AddVertex(vCenter);
			lofter.Build();

			if(lofter.IsDone()) {
				resultFace = lofter.Shape();
			}
		}

		// 4. Validate the Result
		// Never return a broken face, it will crash the Sewer or Unifier later
		if(!resultFace.IsNull()) {
			if(IsShapeValid(resultFace)) {
				return resultFace;
			}
		}
	} catch(...) {
		return TopoDS_Shape();
	}

	return TopoDS_Shape();
}

//// Simplify to Solid: fill holes and return solid
Shape ccadcpp_simplify_to_solid(const Shape s) {
	if(!s) {
		return nullptr;
	}

	TopoDS_Shape resultShape;

	try {
		TopoDS_Shape input = fromShape(s);

		// 1. PRE-FIX: Aggressive Topology Repair
		ShapeFix_Shape preFixer(input);
		preFixer.SetPrecision(1.0e-4);
		preFixer.SetMaxTolerance(1.0e-2);
		preFixer.Perform();
		input = preFixer.Shape();

		// 2. SEWING
		BRepBuilderAPI_Sewing sewer(1.0e-3);
		sewer.Add(input);
		sewer.Perform();
		TopoDS_Shape sewedShape = sewer.SewedShape();

		// 3. HOLE FILLING
		ShapeAnalysis_FreeBounds freeBounds(sewedShape, 1.0e-3);
		const TopoDS_Compound &closedWires = freeBounds.GetClosedWires();

		bool capsAdded = false;
		TopExp_Explorer wireExp(closedWires, TopAbs_WIRE);

		while(wireExp.More()) {
			TopoDS_Wire wire = TopoDS::Wire(wireExp.Current());

			// Defensive Cap
			TopoDS_Shape cap = FillHole(wire);
			if(!cap.IsNull()) {
				sewer.Add(cap);
				capsAdded = true;
			}
			wireExp.Next();
		}

		// 4. RE-SEW
		if(capsAdded) {
			sewer.Perform();
			sewedShape = sewer.SewedShape();
		}

		// 5. MAKE SOLID
		TopoDS_Compound resultComp;
		BRep_Builder b;
		b.MakeCompound(resultComp);

		TopExp_Explorer shellExp(sewedShape, TopAbs_SHELL);
		bool foundAny = false;

		while(shellExp.More()) {
			TopoDS_Shell sh = TopoDS::Shell(shellExp.Current());
			BRepBuilderAPI_MakeSolid mkSolid(sh);

			if(mkSolid.IsDone()) {
				TopoDS_Solid rawSolid = mkSolid.Solid();

				// Fix orientation (Make sure volume is positive)
				ShapeFix_Solid solidFixer(rawSolid);
				solidFixer.Perform();

				b.Add(resultComp, solidFixer.Solid());
				foundAny = true;
			}
			shellExp.Next();
		}

		if(!foundAny) {
			resultShape = sewedShape; // Fallback to shell
		} else {
			resultShape = resultComp;
		}

		// 6. SIMPLIFY
		try {
			// Pre-check validity before passing to Unifier
			if(IsShapeValid(resultShape)) {
				ShapeUpgrade_UnifySameDomain unifier(resultShape, Standard_True, Standard_True, Standard_False);
				unifier.AllowInternalEdges(Standard_True);
				unifier.Build();
				return toShape(unifier.Shape());
			} else {
				// If shape is invalid, UnifySameDomain WILL crash. Return un-simplified.
				return toShape(resultShape);
			}
		} catch(...) {
			return toShape(resultShape);
		}

	} catch(Standard_Failure const &) {
		return nullptr;
	} catch(...) {
		return nullptr;
	}
}

//// Simplify: basic unify
Shape ccadcpp_simplify(const Shape s) {
	// return s;
	// OSD::SetSignal(Standard_True);
	if(!s) {
		return nullptr;
	}

	try {
		TopoDS_Shape shape = fromShape(s);

		// 1. Fix
		ShapeFix_Shape fixer(shape);
		fixer.Perform();
		shape = fixer.Shape();

		// 2. Validate
		if(!IsShapeValid(shape)) {
			return s;
		}

		// 3. Unify
		ShapeUpgrade_UnifySameDomain unifier(shape, Standard_True, Standard_True, Standard_False);
		unifier.AllowInternalEdges(Standard_True);
		unifier.Build();
		return toShape(unifier.Shape());
		return toShape(shape);

	} catch(...) {
		return s;
	}
}

// ----------------------------------------------------------------
// Analysis
// ----------------------------------------------------------------
int ccadcpp_bounding_sphere(Shape s, double *x, double *y, double *z, double *r) {
	if(!s || !x || !y || !z || !r) {
		return -1;
	}
	try {
		// 1. Compute Axis Aligned Bounding Box (AABB)
		Bnd_Box box;
		// Use triangulation if available, otherwise geometry controls
		// true = use triangulation (more precise if meshed), false = geometry
		BRepBndLib::Add(fromShape(s), box, false);

		if(box.IsVoid()) {
			return -1;
		}

		// 2. Get Min/Max
		double xmin, ymin, zmin, xmax, ymax, zmax;
		box.Get(xmin, ymin, zmin, xmax, ymax, zmax);

		// 3. Compute Center
		*x = (xmin + xmax) * 0.5;
		*y = (ymin + ymax) * 0.5;
		*z = (zmin + zmax) * 0.5;

		// 4. Compute Radius (Half Diagonal)
		// This is a "loose" bounding sphere. Minimal enclosing sphere
		// is harder to compute, but this covers the object safely.
		double dx = xmax - xmin;
		double dy = ymax - ymin;
		double dz = zmax - zmin;

		*r = 0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);

		return 0;
	} catch(...) { return -1; }
}

double ccadcpp_volume(Shape s) {
	if(!s) {
		return 0.0;
	}
	try {
		GProp_GProps props;
		BRepGProp::VolumeProperties(fromShape(s), props);
		return props.Mass(); // Mass = Volume when density is 1 (default)
	} catch(...) { return 0.0; }
}

double ccadcpp_surface_area(Shape s) {
	if(!s) {
		return 0.0;
	}
	try {
		GProp_GProps props;
		BRepGProp::SurfaceProperties(fromShape(s), props);
		return props.Mass(); // Mass = Area for surface props
	} catch(...) { return 0.0; }
}

int ccadcpp_centroid(Shape s, double *cx, double *cy, double *cz) {
	if(!s || !cx || !cy || !cz) {
		return -1;
	}
	try {
		const TopoDS_Shape &shape = fromShape(s);
		GProp_GProps props;

		// Determine best property type based on shape dimension
		// 1. Try Volume (Solid)
		BRepGProp::VolumeProperties(shape, props);

		// If volume is negligible (e.g. it's a Face or Wire), try Surface
		if(props.Mass() < Precision::Confusion()) {
			BRepGProp::SurfaceProperties(shape, props);
		}

		// If Area is negligible (e.g. it's a Wire/Edge), try Linear
		if(props.Mass() < Precision::Confusion()) {
			BRepGProp::LinearProperties(shape, props);
		}

		// If still zero, we can't compute a centroid
		if(props.Mass() < Precision::Confusion()) {
			return -1;
		}

		gp_Pnt center = props.CentreOfMass();
		*cx = center.X();
		*cy = center.Y();
		*cz = center.Z();

		return 0;
	} catch(...) { return -1; }
}

int ccadcpp_inertia(Shape s, double density,
		     double *Ixx, double *Iyy, double *Izz) {
	if(!s || !Ixx || !Iyy || !Izz) {
		return -1;
	}
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
	} catch(...) {
		return -1;
	}
}

// Enumerator
// 0 = Edge, 1 = Wire, 2 = Face, 3 = Solid, 4 = Other
int ccadcpp_shape_type(Shape s) {
    TopoDS_Shape shape = fromShape(s);
    if (shape.IsNull()) {
        return -1;   // invalid or null shape
    }

    switch (shape.ShapeType()) {
        case TopAbs_EDGE:   return 0;
        case TopAbs_WIRE:   return 1;
        case TopAbs_FACE:   return 2;
        case TopAbs_SOLID:  return 3;
        default:            return 4;   // any other TopAbs type
    }
}

// ----------------------------------------------------------------
// I/O and Cleanup
// ----------------------------------------------------------------
int ccadcpp_write_step(Shape s, const char *filename) {
	if(!s || !filename) {
		return -1;
	}
	try {
		STEPControl_Writer writer;
		writer.Transfer(fromShape(s), STEPControl_AsIs);
		return (writer.Write(filename) == IFSelect_RetDone) ? 0 : -1;
	} catch(...) { return -1; }
}

int ccadcpp_write_stl(Shape s, const char *filename, float resolution) {
	if(!s || !filename) {
		return -1;
	}
	try {
		const TopoDS_Shape &shape = fromShape(s);

		// Force triangulation with a reasonable deflection
		BRepMesh_IncrementalMesh mesh(shape, resolution);
		if(!mesh.IsDone()) { // lazy evaluates.
			return -1;   // Meshing failed
		}

		StlAPI_Writer writer;
		// Optionally set to ASCII for debugging (use Standard_False)
		writer.ASCIIMode() = false; // true = ASCII, false = Binary

		Standard_Boolean success = writer.Write(shape, filename);
		return success ? 0 : -1;
	} catch(...) {
		return -1;
	}
}

int ccadcpp_write_iges(Shape s, const char *filename) {
	if(!s || !filename) {
		return -1;
	}
	try {
		IGESControl_Writer writer;
		Standard_Boolean ok = writer.AddShape(fromShape(s));
		if(!ok) {
			return -1;
		}

		writer.ComputeModel();
		Standard_Boolean success = writer.Write(filename);
		return success ? 0 : -1;
	} catch(...) { return -1; }
}

int ccadcpp_write_obj(Shape s, const char *filename, float resolution) {
	if(!s || !filename) {
		return -1;
	}
	try {
		const TopoDS_Shape &shape = fromShape(s);

		// 1. Ensure the shape is triangulated (Meshed)
		//    Deflection controls quality (lower = smoother).
		BRepMesh_IncrementalMesh meshGen(shape, resolution);
		if(!meshGen.IsDone()) {
			return -1;
		}

		std::ofstream file(filename);
		if(!file.is_open()) {
			return -1;
		}

		file << "# Generated by ccadcpp/OCCT Wrapper\n";
		file << std::fixed << std::setprecision(6);

		// OBJ indices are global and 1-based. We must track offsets per face.
		int global_vertex_offset = 1;

		// 2. Iterate over every face to extract its triangulation
		TopExp_Explorer expl(shape, TopAbs_FACE);
		for(; expl.More(); expl.Next()) {
			const TopoDS_Face &face = TopoDS::Face(expl.Current());
			TopLoc_Location loc;
			Handle(Poly_Triangulation) tri = BRep_Tool::Triangulation(face, loc);

			if(tri.IsNull()) {
				continue;
			}

			// Get nodes (vertices)
			// OCCT Poly_Triangulation nodes are 1-based in the array
			Standard_Integer nbNodes = tri->NbNodes();
			for(Standard_Integer i = 1; i <= nbNodes; ++i) {
				gp_Pnt p = tri->Node(i);
				// Apply location transformation if the face is instanced/moved
				p.Transform(loc.Transformation());
				file << "v " << p.X() << " " << p.Y() << " " << p.Z() << "\n";
			}

			// Get triangles (faces)
			Standard_Integer nbTriangles = tri->NbTriangles();
			for(Standard_Integer i = 1; i <= nbTriangles; ++i) {
				const Poly_Triangle &t = tri->Triangle(i);
				Standard_Integer n1, n2, n3;
				t.Get(n1, n2, n3);

				// Check face orientation. If reversed, flip normal by swapping indices
				if(face.Orientation() == TopAbs_REVERSED) {
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
	} catch(...) { return -1; }
}
