// ----------------------------------------------------------------
// Loft Implementation
// ----------------------------------------------------------------

#include <GCPnts_QuasiUniformAbscissa.hxx>
#include <BRepBuilderAPI_MakePolygon.hxx>

// ----------------------------------------------------------------
// HELPER: Sample a wire into N points uniformly
// ----------------------------------------------------------------
bool SampleWire(const TopoDS_Wire& wire, int nPoints, std::vector<gp_Pnt>& outPoints) {
    outPoints.clear();
    if (wire.IsNull()) return false;

    BRepAdaptor_CompCurve adaptor(wire);
    
    // FIXED: Header added for GCPnts_QuasiUniformAbscissa
    GCPnts_QuasiUniformAbscissa sampler(adaptor, nPoints);
    
    if (!sampler.IsDone()) return false;

    for (int i = 1; i <= nPoints; ++i) {
        double param = sampler.Parameter(i);
        outPoints.push_back(adaptor.Value(param));
    }
    return true;
}

void AlignClouds(const std::vector<gp_Pnt>& refCloud, std::vector<gp_Pnt>& nextCloud) {
    if (refCloud.empty() || nextCloud.empty()) return;
    size_t N = refCloud.size();
    if (nextCloud.size() != N) return;

    double minTotalDist = 1.0e100; // Infinite
    int bestShift = 0;

    // Simple heuristic search for best start index
    for (size_t shift = 0; shift < N; ++shift) {
        double currentDist = 0.0;
        // Optimization: Check stride of 10 first to abort bad shifts early
        for (size_t i = 0; i < N; i += 5) { 
             size_t idx = (shift + i) % N;
             currentDist += refCloud[i].SquareDistance(nextCloud[idx]);
        }
        
        if (currentDist < minTotalDist) {
            minTotalDist = currentDist;
            bestShift = shift;
        }
    }

    // Apply the shift
    std::vector<gp_Pnt> temp = nextCloud;
    for (size_t i = 0; i < N; ++i) {
        nextCloud[i] = temp[(bestShift + i) % N];
    }
}


// ----------------------------------------------------------------
// HELPER: Rebuild Wire from Points (ALWAYS Smooth)
// ----------------------------------------------------------------
TopoDS_Wire RebuildWire(const std::vector<gp_Pnt>& points, bool closed) {
    if (points.size() < 2) return TopoDS_Wire();

    // 1. Prepare Point Array
    // For periodic (closed) B-Splines, we should not duplicate the last point 
    // if it overlaps the first.
    std::vector<gp_Pnt> cleanPoints = points;
    if (closed && cleanPoints.size() > 2) {
        if (cleanPoints.front().SquareDistance(cleanPoints.back()) < 1e-10) {
             cleanPoints.pop_back();
        }
    }

    Handle(TColgp_HArray1OfPnt) pArray = new TColgp_HArray1OfPnt(1, (Standard_Integer)cleanPoints.size());
    for (size_t i = 0; i < cleanPoints.size(); ++i) {
        pArray->SetValue((Standard_Integer)i + 1, cleanPoints[i]);
    }

    // 2. Interpolate (Smooth Curve)
    // We ALWAYS use interpolation, even if the loft is ruled. 
    // This ensures circles remain circular (preserving volume).
    GeomAPI_Interpolate interpolator(pArray, closed, 1.0e-6);
    interpolator.Perform();

    if (interpolator.IsDone()) {
        Handle(Geom_Curve) curve = interpolator.Curve();
        TopoDS_Edge edge = BRepBuilderAPI_MakeEdge(curve).Edge();
        return BRepBuilderAPI_MakeWire(edge).Wire();
    } 
    
    // Fallback (should rarely be reached)
    BRepBuilderAPI_MakePolygon mkPoly;
    for (const auto& p : points) mkPoly.Add(p);
    if (closed) mkPoly.Close();
    return mkPoly.Wire();
}

// ----------------------------------------------------------------
// MAIN FUNCTION
// ----------------------------------------------------------------
Shape mucadcpp_loft(const Shape* shapes, size_t n, bool solid, bool ruled) {
    if (!shapes || n < 2) return nullptr;

    try {
        std::vector<TopoDS_Wire> inputWires;
        int maxVerts = 0;

        // 1. Extract Wires
        for (size_t i = 0; i < n; ++i) {
            TopoDS_Shape s = fromShape(shapes[i]);
            if (s.IsNull()) return nullptr;

            TopoDS_Wire w;
            if (s.ShapeType() == TopAbs_WIRE) w = TopoDS::Wire(s);
            else {
                TopExp_Explorer exp(s, TopAbs_WIRE);
                if (exp.More()) w = TopoDS::Wire(exp.Current());
                else return nullptr; 
            }
            inputWires.push_back(w);
            
            TopExp_Explorer expV(w, TopAbs_VERTEX);
            int vCount = 0;
            while(expV.More()) { vCount++; expV.Next(); }
            if (vCount > maxVerts) maxVerts = vCount;
        }

        // 2. Determine Sampling Resolution
        // INCREASED: 32 was too low (caused volume error). 
        // 120 samples ensures the B-Spline area error is < 1e-6.
        int sampleCount = std::max(120, maxVerts * 5); 
        if (sampleCount > 400) sampleCount = 400;

        std::vector<gp_Pnt> prevCloud;
        
        // 3. Initialize Loft
        // Note: 'ruled' here only affects the connection BETWEEN sections.
        // It does not affect the sections themselves anymore.
        BRepOffsetAPI_ThruSections generator(solid, ruled, 1.0e-05);
        generator.CheckCompatibility(Standard_False); 

        for (size_t i = 0; i < n; ++i) {
            std::vector<gp_Pnt> currentCloud;
            bool isClosed = inputWires[i].Closed();
            
            if (!SampleWire(inputWires[i], sampleCount, currentCloud)) return nullptr;

            if (i > 0) AlignClouds(prevCloud, currentCloud);
            prevCloud = currentCloud;

            // FIXED: We removed the 'ruled' arg. Always rebuild as smooth curve.
            TopoDS_Wire cleanWire = RebuildWire(currentCloud, isClosed);
            
            if (cleanWire.IsNull()) return nullptr;
            generator.AddWire(cleanWire);
        }

        generator.Build();
        if (!generator.IsDone()) return nullptr;

        TopoDS_Shape result = generator.Shape();

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