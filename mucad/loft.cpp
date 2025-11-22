#include "occt_simple.h"

// OCCT Includes
#include <BRepAdaptor_CompCurve.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeSolid.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepBuilderAPI_MakeWire.hxx>
#include <BRepBuilderAPI_Sewing.hxx>
#include <BRepOffsetAPI_MakeFilling.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <GCPnts_UniformAbscissa.hxx>
#include <Geom_BSplineCurve.hxx>
#include <Precision.hxx>
#include <ShapeFix_Solid.hxx>
#include <ShapeUpgrade_UnifySameDomain.hxx>
#include <TColStd_Array1OfInteger.hxx>
#include <TColStd_Array1OfReal.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopoDS_Wire.hxx>
#include <algorithm>
#include <cmath>
#include <gp_Pnt.hxx>
#include <vector>

// ----------------------------------------------------------------
// HELPER: Quadratic B-Spline Edge Construction
// ----------------------------------------------------------------

TopoDS_Edge internal_make_quadratic_edge(const gp_Pnt &p1, const gp_Pnt &pMid, const gp_Pnt &p2) {
	gp_Pnt C1;
	C1.SetX(2.0 * pMid.X() - 0.5 * p1.X() - 0.5 * p2.X());
	C1.SetY(2.0 * pMid.Y() - 0.5 * p1.Y() - 0.5 * p2.Y());
	C1.SetZ(2.0 * pMid.Z() - 0.5 * p1.Z() - 0.5 * p2.Z());

	TColgp_Array1OfPnt poles(1, 3);
	poles(1) = p1;
	poles(2) = C1;
	poles(3) = p2;

	TColStd_Array1OfReal knots(1, 2);
	knots(1) = 0.0;
	knots(2) = 1.0;

	TColStd_Array1OfInteger mults(1, 2);
	mults(1) = 3;
	mults(2) = 3;

	Handle(Geom_BSplineCurve) curve = new Geom_BSplineCurve(poles, knots, mults, 2);

	return BRepBuilderAPI_MakeEdge(curve).Edge();
}

// ----------------------------------------------------------------
// HELPER: Struct to hold sampled wire data
// ----------------------------------------------------------------
struct SampledWire {
	TopoDS_Wire wire;
	std::vector<gp_Pnt> cloud;
	std::vector<double> params;
	BRepAdaptor_CompCurve adaptor;
	bool isClosed;
	double period;


	SampledWire() : isClosed(false), period(0.0) {}
};

// ----------------------------------------------------------------
// MAIN FUNCTION
// ----------------------------------------------------------------
Shape mucad_loft_robust(const Shape *shapes, size_t n, bool solid, bool ruled) {

	if(!shapes || n < 2) {
		return nullptr;
	}

	try {

		std::vector<SampledWire> wires(n);
		Standard_Integer maxVerts = 0;

		// -------------------------------------------------
		// 1. Pre-process: Validation and Sampling
		// -------------------------------------------------
		for(size_t i = 0; i < n; ++i) {
			const TopoDS_Shape &s = fromShape(shapes[i]);
			if(s.IsNull()) {
				return nullptr;
			}

			if(s.ShapeType() == TopAbs_WIRE) {
				wires[i].wire = TopoDS::Wire(s);
			} else {
				TopExp_Explorer exp(s, TopAbs_WIRE);
				if(exp.More()) {
					wires[i].wire = TopoDS::Wire(exp.Current());
				} else {
					return nullptr; // Input was not a wire
				}
			}

			// Check if wire has ANY edges
			TopExp_Explorer eExp(wires[i].wire, TopAbs_EDGE);
			if(!eExp.More()) {
				return nullptr; // Empty wire
			}

			wires[i].adaptor.Initialize(wires[i].wire);
			wires[i].isClosed = wires[i].wire.Closed();
			if(wires[i].isClosed) {
				wires[i].period = wires[i].adaptor.LastParameter() - wires[i].adaptor.FirstParameter();
			}

			// Count Vertices
			int vCount = 0;
			TopExp_Explorer vExp(wires[i].wire, TopAbs_VERTEX);
			for(; vExp.More(); vExp.Next()) {
				vCount++;
			}
			if(vCount > maxVerts) {
				maxVerts = vCount;
			}
		}


		if(maxVerts == 0) {
			return nullptr;
		}


		// Minimum 50, Max based on vertices, but hard cap at 500 (usually sufficient for alignment)
		int cloudSize = std::max(50, maxVerts * 10);
		if(cloudSize > 500) {
			cloudSize = 500;
		}

		for(size_t i = 0; i < n; ++i) {
			SampledWire &w = wires[i];
			w.cloud.reserve(cloudSize);
			w.params.reserve(cloudSize);

			GCPnts_UniformAbscissa uniform(w.adaptor, cloudSize, Precision::Confusion());

			if(uniform.IsDone()) {

				for(Standard_Integer p = 1; p <= cloudSize; ++p) {
					double u = uniform.Parameter(p);
					w.params.push_back(u);
					w.cloud.push_back(w.adaptor.Value(u));
				}
			} else {
				// Fallback linear sampling
				double first = w.adaptor.FirstParameter();
				double last = w.adaptor.LastParameter();
				for(int p = 0; p < cloudSize; ++p) {
					double u = first + (last - first) * ((double) p / (double) (cloudSize - 1));
					w.params.push_back(u);
					w.cloud.push_back(w.adaptor.Value(u));
				}
			}
		}


		BRepBuilderAPI_Sewing sewer(1.0e-4);

		// -------------------------------------------------
		// 2. Lofting Loop
		// -------------------------------------------------
		for(size_t i = 0; i < n - 1; ++i) {
			SampledWire &w1 = wires[i];
			SampledWire &w2 = wires[i + 1];

			// A. Extract Vertices from W1
			std::vector<gp_Pnt> v1_pnts;
			std::vector<double> v1_params_on_w1;


			BRepTools_WireExplorer exp(w1.wire);
			for(; exp.More(); exp.Next()) {
				TopoDS_Vertex v = exp.CurrentVertex();
				gp_Pnt p = BRep_Tool::Pnt(v);
				v1_pnts.push_back(p);

				// Find approx param
				double u = 0.0;
				double minDist = Precision::Infinite();
				for(size_t k = 0; k < w1.cloud.size(); ++k) {
					double d = p.SquareDistance(w1.cloud[k]);
					if(d < minDist) {
						minDist = d;
						u = w1.params[k];
					}
				}
				v1_params_on_w1.push_back(u);

				// If open wire and this is the last edge, grab the end vertex
				if(!w1.isClosed && exp.Current() == exp.Current()) { // Check for last edge?
										     // BRepTools_WireExplorer doesn't have IsLast(), so we check implicitly via loop.
										     // Better approach for logic:
				}
			}
			// Re-implementing extraction to be 100% safe for Open Wires:
			// Clear and restart
			v1_pnts.clear();
			v1_params_on_w1.clear();

			exp.Init(w1.wire);
			while(exp.More()) {
				TopoDS_Vertex v = exp.CurrentVertex();
				gp_Pnt p = BRep_Tool::Pnt(v);
				v1_pnts.push_back(p);

				// Find param
				double u = 0.0;
				double minDist = Precision::Infinite();
				for(size_t k = 0; k < w1.cloud.size(); ++k) {
					double d = p.SquareDistance(w1.cloud[k]);
					if(d < minDist) {
						minDist = d;
						u = w1.params[k];
					}
				}
				v1_params_on_w1.push_back(u);

				// If open, get the last vertex of the last edge
				TopoDS_Edge e = exp.Current();
				exp.Next();
				if(!exp.More() && !w1.isClosed) {
					TopoDS_Vertex vLast = TopExp::LastVertex(e);
					gp_Pnt pLast = BRep_Tool::Pnt(vLast);
					v1_pnts.push_back(pLast);

					// Find param for last
					double uLast = 0.0;
					double minDLast = Precision::Infinite();
					for(size_t k = 0; k < w1.cloud.size(); ++k) {
						double d = pLast.SquareDistance(w1.cloud[k]);
						if(d < minDLast) {
							minDLast = d;
							uLast = w1.params[k];
						}
					}
					v1_params_on_w1.push_back(uLast);
				}
			}

			// B. Match Vertices of W1 to Cloud of W2
			std::vector<int> w2_indices;
			std::vector<double> w2_params;

			int N1 = (int) v1_pnts.size();
			int N2 = (int) w2.cloud.size();


			if(N1 == 0 || N2 == 0) {
				continue;
			}

			int step = std::max(1, N2 / N1);
			int bestShift = 0;
			double minTotalDist = Precision::Infinite();

			for(int shift = 0; shift < N2; shift += step) {
				double dist = 0;
				for(int k = 0; k < N1; ++k) {
					int targetIdx = (shift + (Standard_Integer) ((double) k * N2 / N1)) % N2;
					dist += v1_pnts[k].SquareDistance(w2.cloud[targetIdx]);
				}
				if(dist < minTotalDist) {
					minTotalDist = dist;
					bestShift = shift;
				}
			}

			int currentSearchIdx = bestShift;
			for(int k = 0; k < N1; ++k) {
				double locMin = Precision::Infinite();
				int bestLocIdx = -1;
				int searchRange = std::max(5, (N2 / N1) * 2); // Ensure minimal search range

				for(int s = 0; s < searchRange; ++s) {
					int idx = (currentSearchIdx + s) % N2;
					double d = v1_pnts[k].SquareDistance(w2.cloud[idx]);
					if(d < locMin) {
						locMin = d;
						bestLocIdx = idx;
					}
				}
				w2_indices.push_back(bestLocIdx);
				w2_params.push_back(w2.params[bestLocIdx]);
				currentSearchIdx = bestLocIdx;
			}

			// C. Build Patches
			// For open wires, we have N-1 segments. For closed, we have N.
			int numSegments = w1.isClosed ? N1 : N1 - 1;

			for(int k = 0; k < numSegments; ++k) {
				int nextK = (k + 1) % N1;

				gp_Pnt p00 = v1_pnts[k];
				gp_Pnt p10 = v1_pnts[nextK];
				gp_Pnt p01 = w2.cloud[w2_indices[k]];
				gp_Pnt p11 = w2.cloud[w2_indices[nextK]];


				double u1_mid;
				double u1_start = v1_params_on_w1[k];
				double u1_end = v1_params_on_w1[nextK];

				if(w1.isClosed) {
					if(u1_end < u1_start) {
						u1_end += w1.period;
					}
					u1_mid = (u1_start + u1_end) * 0.5;
					// Normalize back to range if needed
					// (BRepAdaptor usually handles out-of-range, but safe to normalize)
				} else {
					// Linear average for open segments
					u1_mid = (u1_start + u1_end) * 0.5;
				}
				gp_Pnt pMid_Bottom = w1.adaptor.Value(u1_mid);

				double u2_mid;
				double u2_start = w2_params[k];
				double u2_end = w2_params[nextK];

				if(w2.isClosed) {
					if(u2_end < u2_start) {
						u2_end += w2.period;
					}
					u2_mid = (u2_start + u2_end) * 0.5;
				} else {
					u2_mid = (u2_start + u2_end) * 0.5;
				}
				gp_Pnt pMid_Top = w2.adaptor.Value(u2_mid);

				gp_Pnt pMid_Left = p00.Translated(gp_Vec(p00, p01) * 0.5);
				gp_Pnt pMid_Right = p10.Translated(gp_Vec(p10, p11) * 0.5);


				TopoDS_Edge eBottom = internal_make_quadratic_edge(p00, pMid_Bottom, p10);
				TopoDS_Edge eTop = internal_make_quadratic_edge(p01, pMid_Top, p11);

				// Note: Orientation of eTop needs to be REVERSED relative to eBottom for the filling loop?
				// BRepOffsetAPI_MakeFilling expects edges to form a closed loop.
				// p00 -> p10 (Bottom)
				// p10 -> p11 (Right)
				// p11 -> p01 (Top - Needs to be reversed to flow p11->p01)
				// p01 -> p00 (Left)

				// However, MakeFilling Add() doesn't strictly enforce wire direction if passed as edges,
				// but geometric continuity matters.
				// Let's construct them directionally:
				TopoDS_Edge eRight = ruled ? BRepBuilderAPI_MakeEdge(p10, p11).Edge()
							   : internal_make_quadratic_edge(p10, pMid_Right, p11);

				// Top Edge: Construct P11 -> P01 to maintain loop flow
				TopoDS_Edge eTopRev = internal_make_quadratic_edge(p11, pMid_Top, p01);

				TopoDS_Edge eLeftRev = ruled ? BRepBuilderAPI_MakeEdge(p01, p00).Edge()
							     : internal_make_quadratic_edge(p01, pMid_Left, p00);

				BRepOffsetAPI_MakeFilling fill(2, 10, 2, Standard_False, 1e-5, 1e-4, 0.1, 0.01);
				fill.Add(eBottom, GeomAbs_C0);
				fill.Add(eRight, GeomAbs_C0);
				fill.Add(eTopRev, GeomAbs_C0);
				fill.Add(eLeftRev, GeomAbs_C0);

				fill.Build();


				if(fill.IsDone()) {
					sewer.Add(fill.Shape());
				}
				// Optional: Else log error?
			}
		}

		// -------------------------------------------------
		// 3. Capping (if Solid)
		// -------------------------------------------------
		if(solid) {
			TopoDS_Wire startW = wires[0].wire;
			TopoDS_Wire endW = wires[n - 1].wire;


			if(startW.Closed()) {
				BRepBuilderAPI_MakeFace mkStart(startW, true); // Try planar
				if(mkStart.IsDone()) {
					sewer.Add(mkStart.Face());
				}
				// Fallback logic omitted for brevity, standard MakeFace usually best for caps
			}

			if(endW.Closed()) {
				BRepBuilderAPI_MakeFace mkEnd(endW, true);
				if(mkEnd.IsDone()) {
					sewer.Add(mkEnd.Face());
				}
			}
		}

		// -------------------------------------------------
		// 4. Finalize
		// -------------------------------------------------
		sewer.Perform();
		TopoDS_Shape sewed = sewer.SewedShape();

		TopoDS_Shape result = sewed;


		if(solid) {
			bool isSolid = false;
			TopExp_Explorer shellExp(sewed, TopAbs_SHELL);
			if(shellExp.More()) {
				TopoDS_Shell sh = TopoDS::Shell(shellExp.Current());
				BRepBuilderAPI_MakeSolid mkS(sh);
				if(mkS.IsDone()) {
					result = mkS.Solid();
					isSolid = true;
				}
			}

			// Only run fix if we successfully made a solid
			if(isSolid) {
				ShapeFix_Solid fix;
				fix.Init(TopoDS::Solid(result));
				fix.Perform();
				result = fix.Solid();
			}
		}

		ShapeUpgrade_UnifySameDomain unifier(result, Standard_True, Standard_True, Standard_True);
		unifier.Build();

		return toShape(unifier.Shape());

	} catch(...) {
		return nullptr;
	}
}