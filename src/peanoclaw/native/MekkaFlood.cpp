#if defined(SWE) || defined(PEANOCLAW_FULLSWOF2D)

#include <algorithm>

#include "peanoclaw/native/MekkaFlood.h"

/*MekkaFlood_SWEKernelScenario::MekkaFlood_SWEKernelScenario(double domainSize) : 
    domainSize(domainSize),
    bathymetryHelper("smith_sandwell.nc", "X2MIN", "Y2MIN","ROSE")
{

}*/

MekkaFlood_SWEKernelScenario::MekkaFlood_SWEKernelScenario(DEM& dem) : dem(dem)
{

}

MekkaFlood_SWEKernelScenario::~MekkaFlood_SWEKernelScenario() {}

void MekkaFlood_SWEKernelScenario::initializePatch(peanoclaw::Patch& patch) {
    // dam coordinates
    //double x0=domainSize*0.5;
    //double y0=domainSize*0.5;
  
    double x_size = dem.upper_right(0) - dem.lower_left(0);
    double y_size = dem.upper_right(1) - dem.lower_left(1);

    double x0=(x_size) * 0.5;
    double y0=(y_size) * 0.8;
    
    // Riemann states of the dam break problem
    double radDam = 0.05*std::min(x_size,y_size);
    double hl = 0.; // 2
    double ul = 0.;
    double vl = 0.;
    double hr = 0.; // 1
    double ur = 0.;
    double vr = 0.;

    double q0 = 0;
    double q1 = 0;
    
    // compute from mesh data
    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();
    const tarch::la::Vector<DIMENSIONS, double> patchPosition = patch.getPosition();
    const tarch::la::Vector<DIMENSIONS, double> meshWidth = patch.getSubcellSize();
 
    int ghostlayerWidth = patch.getGhostlayerWidth();

    // initialize new part only
    tarch::la::Vector<DIMENSIONS, int> subcellIndex;
    for (int yi = 0; yi < patch.getSubdivisionFactor()(1); yi++) {
        for (int xi = 0; xi < patch.getSubdivisionFactor()(0); xi++) {
            subcellIndex(0) = xi;
            subcellIndex(1) = yi;
 
            double X = patchPosition(0) + xi*meshWidth(0);
            double Y = patchPosition(1) + yi*meshWidth(1);
  
            tarch::la::Vector<DIMENSIONS, double> coords = mapMeshToCoordinates(X, Y);
            //double bathymetry = bathymetryHelper.getHeight(coords(0),coords(1));
            double bathymetry = dem(coords(0), coords(1));

            //std::cout << "bathymetry: " << bathymetry << " @ " << coords(0) << " " << coords(1) << std::endl;

            double r = sqrt((X-x0)*(X-x0) + (Y-y0)*(Y-y0));
            double h = hl*(r<=radDam) + hr*(r>radDam);
            double u = hl*ul*(r<=radDam) + hr*ur*(r>radDam);
            double v = hl*vl*(r<=radDam) + hr*vr*(r>radDam);
  
            //bathymetry = 0.0; // BENCHMARK

            patch.setValueUNew(subcellIndex, 0, h);
            patch.setValueUNew(subcellIndex, 1, u);
            patch.setValueUNew(subcellIndex, 2, v);
            patch.setValueUNew(subcellIndex, 3, bathymetry);
            patch.setValueAux(subcellIndex, 0, bathymetry);
 
            patch.setValueUNew(subcellIndex, 4, h * u);
            patch.setValueUNew(subcellIndex, 5, h * v);
        }
    }
}

double MekkaFlood_SWEKernelScenario::computeDemandedMeshWidth(peanoclaw::Patch& patch) {
    double retval = 0.0;

    const tarch::la::Vector<DIMENSIONS, double> patchPosition = patch.getPosition();
    const tarch::la::Vector<DIMENSIONS, double> patchSize = patch.getSize();

    tarch::la::Vector<DIMENSIONS, double> patchCenter;
    patchCenter(0) = patchPosition(0) + patchSize(0)/2.0f;
    patchCenter(1) = patchPosition(1) + patchSize(1)/2.0f;

    double outerRadius = std::max(patchSize(0),patchSize(1))/2.0*sqrt(2);

    //const tarch::la::Vector<DIMENSIONS, double> mekkaPosition = mapCoordinatesToMesh(mekka_lon,mekka_lat);

    // check if mekka is inside or at least near to our current patch
    //double mekka_distance = sqrt((mekkaPosition(0)-patchCenter(0))*(mekkaPosition(0)-patchCenter(0)) + (mekkaPosition(1)-patchCenter(1))*(mekkaPosition(1)-patchCenter(1)));
 
    const tarch::la::Vector<DIMENSIONS, double> meshWidth = patch.getSubcellSize();
   
    // update bathymetry (gets somehow lost)
    double x_size = dem.upper_right(0) - dem.lower_left(0);
    double y_size = dem.upper_right(1) - dem.lower_left(1);
    
    double x0=(x_size) * 0.5;
    double y0=(y_size) * 0.5;
    double radDam = 0.05*std::min(x_size,y_size);
    bool isInsideCircle = false;
    tarch::la::Vector<DIMENSIONS, int> subcellIndex;
    for (int yi = -1; yi < patch.getSubdivisionFactor()(1)+1; yi++) {
        for (int xi = -1; xi < patch.getSubdivisionFactor()(0)+1; xi++) {
            subcellIndex(0) = xi;
            subcellIndex(1) = yi;
 
            double X = patchPosition(0) + xi*meshWidth(0);
            double Y = patchPosition(1) + yi*meshWidth(1);
  
            tarch::la::Vector<DIMENSIONS, double> coords = mapMeshToCoordinates(X, Y);
            //double bathymetry = bathymetryHelper.getHeight(coords(0),coords(1));
            double bathymetry = dem(0.0, 0.0);

            double r = sqrt((X-x0)*(X-x0) + (Y-y0)*(Y-y0));

            //bathymetry = 0.0; // BENCHMARK
            
            if (yi >= 0 && xi >= 0 
                && yi < patch.getSubdivisionFactor()(1) && xi < patch.getSubdivisionFactor()(0)
               ) {
                //patch.setValueAux(subcellIndex, 0,  bathymetry);
            }
            //patch.setValueUOld(subcellIndex, 3, bathymetry);

            isInsideCircle |=(r <radDam);
        }
    }
    
    // try to adapt to bathymetry
    // loop starts at one, due to central finite differences
    double max_curvature = 0.0; // second derivative
    for (int yi = 1; yi < patch.getSubdivisionFactor()(1)-1; yi++) {
        for (int xi = 1; xi < patch.getSubdivisionFactor()(0)-1; xi++) {
            subcellIndex(0) = xi;
            subcellIndex(1) = yi;
            double bathemetry_11= patch.getValueAux(subcellIndex, 0);

            subcellIndex(0) = xi-1;
            subcellIndex(1) = yi;
            double bathemetry_01= patch.getValueAux(subcellIndex, 0);
 
            subcellIndex(0) = xi+1;
            subcellIndex(1) = yi;
            double bathemetry_21= patch.getValueAux(subcellIndex, 0);
 
            double curvature_x = std::abs((bathemetry_21 + 2*bathemetry_11 + bathemetry_01)/meshWidth(0));

            subcellIndex(0) = xi;
            subcellIndex(1) = yi-1;
            double bathemetry_10= patch.getValueAux(subcellIndex, 0);
 
            subcellIndex(0) = xi;
            subcellIndex(1) = yi+1;
            double bathemetry_12= patch.getValueAux(subcellIndex, 0);

            double curvature_y = std::abs((bathemetry_12 + 2*bathemetry_11 + bathemetry_10)/meshWidth(1));

            max_curvature = std::max(max_curvature, std::max(curvature_x,curvature_y));
        }
    }

    /*if (mekka_distance > outerRadius) {
        retval = domainSize/6/243;
    } else {
        //retval = 1e-3;
        retval = domainSize/6/243;
    }*/

    /*if (meshWidth(0) < 1e-4 || max_curvature > 1e1) {
        retval = meshWidth(0);
    } else {
        //std::cout << "max_curvature is: " << max_curvature << std::endl;
        retval = meshWidth(0)/3.0;
    }*/

    /*if (isInsideCircle && meshWidth(0) > 3e-3) {
        retval = meshWidth(0)/3.0;
    } else {
        retval = meshWidth(0);
    }*/

    retval = std::min(x_size,y_size) / (27.0 * 6.0);
    return retval;
}

// box sizes:
// - complete saudi arabia:
//  longitude: 32 - 60E
//  latitude: 13 - 40N
//
// - mekka area


tarch::la::Vector<DIMENSIONS, double> MekkaFlood_SWEKernelScenario::mapCoordinatesToMesh(double longitude, double latitude) {
    tarch::la::Vector<DIMENSIONS, double> mesh;
 
    /*// put mekkah in the center - adjust bei 0*5 * scale / domainSize

    double scale = 0.1;

    // peano coordinates
    //float mesh_y = (latitude-13.0f)/27.0f * domainSize;
    //float mesh_x = (longitude-32.0f)/28.0f * domainSize;
 
    float mesh_y = (latitude-39.8167f)/1.0f * domainSize / scale + 0.5;
    float mesh_x = (longitude-21.4167f)/1.0f * domainSize / scale + 0.5;

    mesh(0) = mesh_x;
    mesh(1) = mesh_y;*/
	
    double lower_left_0 = dem.lower_left(0);
    double lower_left_1 = dem.lower_left(1);

    double ws_0;
    double ws_1;
    dem.pixelspace_to_worldspace(ws_0, ws_1, longitude, latitude);

    // TODO_: the offset change is just due to a problem with plotting, meh ...
    mesh(0) = ws_0 - lower_left_0;
    mesh(1) = ws_1 - lower_left_1;
    return mesh;
}

tarch::la::Vector<DIMENSIONS, double> MekkaFlood_SWEKernelScenario::mapMeshToCoordinates(double x, double y) {
    tarch::la::Vector<DIMENSIONS, double> coords;

    // put mekkah in the center - adjust bei 0*5 * scale / domainSize

    /*double scale = 0.5;

    //double latitude = x / domainSize * 27.0 + 13.0;
    //double longitude = y / domainSize * 28.0 + 32.0;
 
    double latitude = (y-0.5)* (scale / domainSize) * 1.0 + 21.4167;
    double longitude = (x-0.5)* (scale / domainSize) * 1.0 + 39.8167;
    coords(0) = longitude;
    coords(1) = latitude;*/
 
    double lower_left_0 = dem.lower_left(0);
    double lower_left_1 = dem.lower_left(1);

    double ps_0 = 0.0;
    double ps_1 = 0.0;
    dem.worldspace_to_pixelspace(ps_0, ps_1, x+lower_left_0, y+lower_left_1);

    coords(0) = ps_0;
    coords(1) = ps_1;
    return coords;
}

#endif
