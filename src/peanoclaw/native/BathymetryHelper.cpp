// configure accordingly for dataset

// values for smith_sandwell
#define LON_NAME "X2MIN"
#define LAT_NAME "Y2MIN"
#define DATA_NAME "ROSE"

// values for ferret_listing
//#define LON_NAME "ETOPO05_X1_4321"
//#define LAT_NAME "ETOPO05_Y1_2160"
//#define DATA_NAME "ROSE"

#include "tarch/la/Vector.h"

#include "BathymetryHelper.h"

BathymetryHelper::BathymetryHelper(const std::string& nc_filename,
                 const std::string& lon_name,
                 const std::string& lat_name,
                 const std::string& dataset_name
) {
    int retval = 0;

    /* Open the file. NC_NOWRITE tells netCDF we want read-only access
     * to the file.*/
    if ((retval = nc_open(nc_filename.c_str(), NC_NOWRITE, &ncid))) {
        std::cout << "could not open file: " << nc_filename << std::endl;
    }

    /* There are a number of inquiry functions in netCDF which can be
     * used to learn about an unknown netCDF file. NC_INQ tells how
     * many netCDF variables, dimensions, and global attributes are in
     * the file; also the dimension id of the unlimited dimension, if
     * there is one. */
    int ndims_in, nvars_in, ngatts_in, unlimdimid_in;
    if ((retval = nc_inq(ncid, &ndims_in, &nvars_in, &ngatts_in,
        &unlimdimid_in))) {
        printf("could not get information about file\n");
    } else {
       printf("file has %i dimensions, %i variables, %i attributes, %i unlimited dimension_id\n", ndims_in, nvars_in, ngatts_in, unlimdimid_in);
    }

    /* Get the varids of the latitude and longitude coordinate
     * variables. */
    int lat_varid, lon_varid;
    if ((retval = nc_inq_varid(ncid, lat_name.c_str(), &lat_varid))) {
        printf("could not get variable id of latitude variable\n");
    }
  
    if ((retval = nc_inq_varid(ncid, lon_name.c_str(), &lon_varid))) {
        printf("could not get variable id of longitude variable\n");
    }
  
    if ((retval = nc_inq_varid(ncid, dataset_name.c_str(), &data_varid))) {
        printf("could not get variable id of data variable\n");
    }
  
    nc_type lat_type;
    int lat_ndims, lat_dimids[NC_MAX_VAR_DIMS], lat_natts;
    if ((retval = nc_inq_var(ncid, lat_varid, 0, &lat_type, &lat_ndims, lat_dimids, &lat_natts))) {
        printf("could not get information about latitude variable\n");
    } else {
       printf("latitude variable has %i dimensions, %i attributes\n", lat_ndims, lat_natts);
       int dim = 0;
       /*for (dim=0; dim <lat_ndims; dim++) {
           size_t dimension_length = 0;
           nc_inq_dimlen(ncid, lat_dimids[dim], &dimension_length);
           printf("latitude variable: dimension %i has size %zi\n", dim, dimension_length);
       }*/
       nc_inq_dimlen(ncid, lat_dimids[0], &lat_length);
    }
  
    nc_type lon_type;
    int lon_ndims, lon_dimids[NC_MAX_VAR_DIMS], lon_natts;
    if ((retval = nc_inq_var(ncid, lon_varid, 0, &lon_type, &lon_ndims, lon_dimids, &lon_natts))) {
        printf("could not get information about longitude variable\n");
    } else {
       printf("latitude variable has %i dimensions, %i attributes\n", lon_ndims, lon_natts);
       int dim = 0;
       nc_inq_dimlen(ncid, lon_dimids[0], &lon_length);
    }

    printf("latitude size %zi, longitude size %zi\n", lat_length, lon_length);

    latitude_values = new float[lat_length];
    longitude_values = new float[lon_length];

    nc_get_var_float(ncid, lat_varid, latitude_values);
    nc_get_var_float(ncid, lon_varid, longitude_values);

    size_t data_index[2];
    data_index[0] = 0;
    data_index[1] = 0;

    float data_value = 0.0;
    nc_get_var1_float(ncid, data_varid, data_index, &data_value);
 
    printf("first longitude value %f first latitude value %f => data value %f\n", longitude_values[0], latitude_values[0], data_value);

    /*ptrdiff_t lon_closest_index = -1.0;
    ptrdiff_t lat_closest_index = -1.0;
    lon_closest_index = getClosestIndex1D(370.0, longitude_values, lon_length);
    if (lon_closest_index != -1) {
        printf("closest index for %f: %zi = %f\n", 370.0, lon_closest_index, longitude_values[lon_closest_index]);
    }
 
    lon_closest_index = getClosestIndex1D(0.0, longitude_values, lon_length);
    if (lon_closest_index != -1) {
        printf("closest index for %f: %zi = %f\n", 0.0, lon_closest_index, longitude_values[lon_closest_index]);
    }
 
    lon_closest_index = getClosestIndex1D(10.0, longitude_values, lon_length);
    if (lon_closest_index != -1) {
        printf("closest index for %f: %zi = %f\n", 10.0, lon_closest_index, longitude_values[lon_closest_index]);
    }*/
}
 
BathymetryHelper::~BathymetryHelper() {
    delete[] latitude_values;
    delete[] longitude_values;

    /* Close the file, freeing all resources. */
    int retval = 0;
    if ((retval = nc_close(ncid))) {
        printf("could not close file\n");
    }
}

// binary search
tarch::la::Vector<DIMENSIONS,size_t> BathymetryHelper::getClosestInterval1D(float probe, float *data, const size_t size) {
    size_t search_left = 0;
    size_t search_size = size;
 
    float difference, data_value;
    size_t search_index = -1;
 
    //std::cout << "start search_size " << search_size << std::endl;

    // difference != 0.0
    while (search_size > 1) {
        search_index = search_left + (search_size / 2);
        data_value = data[search_index];
        difference = probe - data_value;

        if (difference > 0) {
            // use right interval
            search_left = search_index + 1;
            search_size = search_size / 2;
        } else {
            // use left interval
            search_size = search_size / 2;
        }
    }

    // in case we have only one item left: lets see if we have to choose our left or right neighbor
    tarch::la::Vector<DIMENSIONS, size_t> interval;
    if (search_size == 1) {
        data_value = data[search_left];
        difference = probe - data_value;
        if (difference > 0.0) {
            interval(0) = search_left;
            interval(1) = search_left+1;
            //std::cout << "correction: right neighbor: probe=" << probe << " difference=" << difference << " data_value=" << data_value << std::endl;
        } else {
            interval(0) = search_left-1;
            interval(1) = search_left;
            //std::cout << "correction: left neighbor: probe=" << probe << " difference=" << difference << " data_value=" << data_value << std::endl;
        }
    } else {
        interval(0) = search_left;
        interval(1) = search_left + search_size - 1;
    }
    return interval;
}

float BathymetryHelper::getHeight(double longitude, double latitude) {
    tarch::la::Vector<DIMENSIONS, size_t> longitude_interval = getClosestInterval1D(longitude, longitude_values, lon_length);
    tarch::la::Vector<DIMENSIONS, size_t> latitude_interval = getClosestInterval1D(latitude, latitude_values, lat_length);
 
    //std::cout << "longitude index interval " << longitude_interval(0) << " " << longitude_interval(1) << std::endl;
    //std::cout << "latitude index interval " << latitude_interval(0) << " " << latitude_interval(1) << std::endl;
  
    tarch::la::Vector<DIMENSIONS, double> longitude_points;
    longitude_points(0) = longitude_values[longitude_interval(0)];
    longitude_points(1) = longitude_values[longitude_interval(1)];
 
    //std::cout << "longitude interval " << longitude_points(0) << " " << longitude_points(1) << std::endl;
 
    tarch::la::Vector<DIMENSIONS, double> latitude_points;
    latitude_points(0) = latitude_values[latitude_interval(0)];
    latitude_points(1) = latitude_values[latitude_interval(1)];
 
    //std::cout << "latitude interval " << latitude_points(0) << " " << latitude_points(1) << std::endl;

    float longitude_distance = longitude_points(1) - longitude_points[0];
    float latitude_distance = latitude_points(1) - latitude_points(0);

    //std::cout << "longitude distance " << longitude_distance << " latitude distance " << latitude_distance << std::endl;
 
    size_t data_index[2];

    // NOTE: dataset is organized as latitude / longitude

    float left_bottom_value = 0.0f;
    data_index[1] = longitude_interval(0);
    data_index[0] = latitude_interval(0);
    nc_get_var1_float(ncid, data_varid, data_index, &left_bottom_value);
  
    float right_bottom_value = 0.0f;
    data_index[1] = longitude_interval(1);
    data_index[0] = latitude_interval(0);
    nc_get_var1_float(ncid, data_varid, data_index, &right_bottom_value);
  
    // 1D interpolation along bottom
    float bottom_interpolated = ((longitude_points(1) - longitude) * left_bottom_value) / longitude_distance + 
                                ((longitude - longitude_points(0)) * right_bottom_value) / longitude_distance;

    //std::cout << "bottom interpolated: " << bottom_interpolated << " with values " << left_bottom_value << " " << right_bottom_value << std::endl;

    float left_top_value = 0.0f;
    data_index[1] = longitude_interval(0);
    data_index[0] = latitude_interval(1);
    nc_get_var1_float(ncid, data_varid, data_index, &left_top_value);
  
    float right_top_value = 0.0f;
    data_index[1] = longitude_interval(1);
    data_index[0] = latitude_interval(1);
    nc_get_var1_float(ncid, data_varid, data_index, &right_top_value);
 
    // 1D interpolation along top
    float top_interpolated = ((longitude_points(1) - longitude) * left_top_value) / longitude_distance + 
                             ((longitude - longitude_points(0)) * right_top_value) / longitude_distance;

    //std::cout << "top interpolated: " << top_interpolated << std::endl;

    // complete the bilinear interpolation along latitude
    float data_value = ((latitude_points(1) - latitude) * bottom_interpolated) / latitude_distance + 
                       ((latitude - latitude_points(0)) * top_interpolated) / latitude_distance;
 
    //std::cout << "bilinear interpolated: " << data_value << std::endl;

    //printf("value(%f,%f) = data(%zi,%zi) = %f\n", latitude, longitude, data_index[0], data_index[1], data_value);

    return data_value;
}

// useful for testing
#if 0
// based on simple_xy and sfc_pres_temp_rd examples
int main(int argc, char **argv) {
    int retval = 0;

    if (argc != 2) {
        printf("wrong arguments\n");
        return -1;
    }
 
    BathymetryHelper bathymetryHelper(argv[1], LON_NAME, LAT_NAME, DATA_NAME);

    // Munich
    float lat = 48.0;
    float lon = 11.0;
 
    // Thuwal, KSA
    //float lat = 22.2833;
    //float lon = 39.1;
 
    // Mekkah, KSA
    //float lat = 21.4167;
    //float lon = 39.8167;

    float height = bathymetryHelper.getHeight(lon,lat);
    std::cout << "height = " << height << std::endl;
    return 0;
}
#endif
