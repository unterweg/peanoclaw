// BASED ON FULLSWOF2D SOLVER

#ifndef _MEKKAFLOOD_SOLVER_H_
#define _MEKKAFLOOD_SOLVER_H_

#include <limits> 

class MekkaFlood_solver {
    protected:

    private:
        // this structs are not allocated and serve just as links to temp data
        struct SchemeArrays {
            float *h;
            float *u;
            float *v;
            float *q1;
            float *q2;
            float *z;
        };

    public:
        struct InputArrays {
             // 0 = stride between elements, 1 stride between rows, 2 stride between patches

            float *h;
            float *u;
            float *v;
            float *z;
        };

        struct TempArrays {
            // set by MUSCL reconstruction
            float *h1r;
            float *h1l;

            float *u1r;
            float *u1l;

            float *v1r;
            float *v1l;

            float *z1r;
            float *z1l;
            float *delta_z1; // maybe we can rules this one away its just: z[i+1][j]-z[i][j];
            float *delzc1;
            float *delz1;

            float *h2r;
            float *h2l;
 
            float *u2r;
            float *u2l;

            float *v2r;
            float *v2l;

            float *z2r;
            float *z2l;
            float *delta_z2; // maybe we can rules this one away its just: z[i+1][j]-z[i][j];
            float *delzc2;
            float *delz2;


            // flux computation
            float *f1;
            float *f2;
            float *f3;

            float *g1;
            float *g2;
            float *g3;

            // scheme
            float *q1;
            float *q2;

            float *hs;
            float *us;
            float *vs;
            float *qs1;
            float *qs2;

            float *hsa;
            float *usa;
            float *vsa;
            float *qsa1;
            float *qsa2;

            float *Vin1;
            float *Vin2;
            float *Vin_tot;

            // friction
            float *Fric_tab;

            // set by hydrostatic reconstruction
            float *h1right;
            float *h1left;
            float *h2right;
            float *h2left;

            // rain:
            float *Tab_rain;
        };


        struct Constants {
            Constants(int nx, int ny, float dx, float dy) : 
                NXCELL(nx), 
                NYCELL(ny), 
                DX(dx),
                DY(dy)
            {
                GRAV = 9.81;
                GRAV_DEM = 4.905;
                CONST_CFL_X = 0.5;
                CONST_CFL_Y = 0.5;
                HE_CA = 1e-8;
                VE_CA = 1e-8;
                MAX_CFL_X = 0.;
                MAX_CFL_Y = 0.;
                NB_CHAR = 256;
                ZERO = 0.;
                IE_CA = 1.e-8;
                EPSILON = 1.e-9;
                Ratio_Close_cell = 1e-3;
                MAX_SCAL = std::numeric_limits<float>::max();
                RainIntensity = 0.001;
                
                FRICCOEF = 0.0; // TODO get real value for this
                CFL_FIX = 0.5;
            }

            float GRAV; // 9.81

            float GRAV_DEM; // 4.905
            float CONST_CFL_X; // 0.5
            float CONST_CFL_Y; // 0.5
            float HE_CA; //1e-12;
            float VE_CA; //1e-12;

            float MAX_CFL_X; // 0.
            float MAX_CFL_Y; // 0.

            int NB_CHAR; // 256 // TODO: what is its use?
            float ZERO; // 0
            float IE_CA; // 1.e-8;
            float EPSILON; // 1.e-13

            float Ratio_Close_cell; // 1e-3
            float MAX_SCAL; // DBL_MAX

            float RainIntensity; // 0.001;

            float FRICCOEF;
            float CFL_FIX;

            const int NXCELL;
            const int NYCELL;
            const float DX;
            const float DY;
        };

        // index helper
        static inline unsigned int linearizeIndex(int dim, unsigned int* index, unsigned int* strides) {
            unsigned int result = 0;
            for (int i=0; i < dim; i++) {
                result += index[i] * strides[i];
            }
            return result;
        }

        static void initializeStrideinfo(const Constants& constants, int dim, unsigned int* strideinfo);
        static void allocateInput(int nr_patches, int dim, unsigned int* strideinfo, InputArrays& input);
        static void allocateTemp(int nr_patches, int dim, unsigned int* strideinfo, TempArrays& temp);
        static void freeInput(InputArrays& input);
        static void freeTemp(TempArrays& temp);

        MekkaFlood_solver();
        virtual ~MekkaFlood_solver();

        // returns used timestep
        static float calcul(const int patchid, int dim, unsigned int* strideinfo, InputArrays& input, TempArrays& temp, const Constants& constants, float dt_max);

        static void boundary(const int patchid, int dim, unsigned int* strideinfo, SchemeArrays& input, TempArrays& temp, const Constants& constants,
                             float time_tmp);
        
        // minmod slope limiter
        static float lim_minmod(float a, float b);

        // muscl reconstruction
        static void rec_muscl_init(const int patchid, int dim, unsigned int* strideinfo, InputArrays& input, TempArrays& temp, const Constants& constants);
        static void rec_muscl(const int patchid, int dim, unsigned int* strideinfo, SchemeArrays& input, TempArrays& temp, const Constants& constants);

        // hydrostatic reconstruction
        static void rec_hydro(float hg, float hd,float dz, float& hg_rec, float& hd_rec);

        // HLL flux
        static void flux_hll(const Constants& constants, float h_L,float u_L,float v_L,float h_R,float u_R,float v_R, float& f1, float& f2, float& f3, float& cfl);

        // rain
        static void rain(const int patchid, int dim, unsigned int* strideinfo, SchemeArrays& input, TempArrays& temp, const Constants& constants,
                         float time);

        // friction: (actually no friction)
        // CAREFUL: input.q1 == output.q1 and input.q2 == output.q2
        static void friction(float uold, float vold, float hnew, float q1new, float q2new, float dt, float cf, float& q1mod, float& q2mod);

        // infiltration: (actually no infiltration)
        static void infiltration(const int patchid, int dim, unsigned int* strideinfo, SchemeArrays& input, TempArrays& temp, const Constants& constants,
                                 float dt);

        // general scheme:
        static void maincalcflux(const int patchid, int dim, unsigned int* strideinfo, TempArrays& temp, const Constants& constants,
                                 float cflfix, float dt_max, float& dt); // dt is both input AND output
        static void maincalcscheme(const int patchid, int dim, unsigned int* strideinfo, SchemeArrays& input, SchemeArrays& output, TempArrays& temp, const Constants& constants,
                                   float tps, float dt, int verif); // there was originally a "n" input parameter which was not used here.
};

#endif // _MEKKAFLOOD_SOLVER_H_
