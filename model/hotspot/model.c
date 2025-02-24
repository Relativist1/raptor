/*
 * model file for HARM3D data
 *
 * Please note that most of the code for the harm3d model was adapted
 * from GRMONTY (Dolence et al., 2009).
 *
 * GRMONTY is released under the GNU GENERAL PUBLIC LICENSE.
 * Modifications were made in August 2016 by T. Bronzwaer and
 * J. Davelaar.
 * Modifications were made in August 2022 by J. Davelaar
 */

#include "definitions.h"
#include "functions.h"
#include "global_vars.h"
#include "model_definitions.h"
#include "model_functions.h"
#include "model_global_vars.h"

// GLOBAL VARS
//////////////
double ****p;

int N1, N2, N3;

double gam;

double R0, Rin, Rout, a, hslope;
double startx[NDIM], stopx[NDIM], dx[NDIM];
double Rh;
double L_unit, T_unit;
double RHO_unit, U_unit, B_unit;
double Ne_unit, Thetae_unit;
double r_isco;
double keplerian_factor = 0.5;
double infall_factor = 0.5;

double Ne_unit = 1.e5;
double Te_unit = 1.e11;

// FUNCTIONS
////////////
#define MULOOP for(int mu=0;mu<NDIM;mu++)
#define MUNULOOP for(int mu=0;mu<NDIM;mu++) for(int nu=0;nu<NDIM;nu++)
void init_model() {

    set_units(M_UNIT);

    a = 0.;
    double R0 = 0.;
    

    Rh = 1 + sqrt(1. - a * a) ;
    Rin = Rh + 1.e-5; 
    Rout = 500;

    double z1 = 1. + pow(1. - a * a, 1. / 3.) * (pow(1. + a, 1. / 3.) + pow(1. - a, 1. / 3.));
    double z2 = sqrt(3. * a * a + z1 * z1);
    r_isco = 3. + z2 - copysign(sqrt((3. - z1) * (3. + z1 + 2. * z2)), a);
    hslope = 1;

    startx[0] = 0.0;
    startx[1] = log(Rin);
    startx[2] = 0.0;
    startx[3] = 0.0;
    stopx[0] = 0.0;
    stopx[1] = log(Rout);
    stopx[2] = M_PI;
    stopx[3] = 2*M_PI;

}


void get_model_ucon(double X[NDIM], double bl_Ucon[NDIM])
{
    double r, th;
    double bl_gcov[NDIM][NDIM], bl_gcon[NDIM][NDIM];

    bl_coord(X, &r, &th);

    gcon_bl(r, th, bl_gcon);
    gcov_bl(r, th, bl_gcov);

    double omegaK, omegaFF, omega, uphi;
    double K, ur, ut;
    double a2 = a*a;    

  if (r < Rh) {
    double bl_Ucov[NDIM];
    bl_Ucov[0] = -1;
    bl_Ucov[1] = 0.;
    bl_Ucov[2] = 0.;
    bl_Ucov[3] = 0.;
   bl_raise(bl_gcon, bl_Ucon, bl_Ucov);

  } else if (r < r_isco) {
    double omegaK_isco = 1. / (pow(r_isco, 3./2) + a);
  
    // Get conserved quantities at the ISCO...
    double bl_Ucon_isco[NDIM], bl_Ucov_isco[NDIM];
    bl_Ucon_isco[0] = 1.0;
    bl_Ucon_isco[1] = 0.0;
    bl_Ucon_isco[2] = 0.0;
    bl_Ucon_isco[3] = omegaK_isco;

    double bl_gcov_isco[NDIM][NDIM];

    gcov_bl(r_isco, th, bl_gcov_isco);

    normalize(bl_Ucon_isco, bl_gcov_isco);
    bl_lower(bl_gcov_isco, bl_Ucov_isco, bl_Ucon_isco);

    double e = bl_Ucov_isco[0];
    double l = bl_Ucov_isco[3];

    // ...then set the infall velocity and find omega
    double bl_Ucon_tmp[NDIM], bl_Ucov_tmp[NDIM];
    double K_con = bl_gcon[0][0] * e * e + 2.0 * bl_gcon[0][3] * e * l + bl_gcon[3][3] * l * l;
    double urk_precut = -(1.0 + K_con) / bl_gcon[1][1];
    double urk = -sqrt(fmax(0.0, urk_precut));
    bl_Ucov_tmp[0] = e;
    bl_Ucov_tmp[1] = urk;
    bl_Ucov_tmp[2] = 0.0;
    bl_Ucov_tmp[3] = l;

    bl_raise(bl_gcon, bl_Ucon_tmp, bl_Ucov_tmp);

    omegaK = bl_Ucon_tmp[3] / bl_Ucon_tmp[0];

    omegaFF = bl_gcon[0][3] / bl_gcon[0][0];
    // Compromise
    omega = omegaK + (1 - keplerian_factor)*(omegaFF - omegaK);

    // Then set the infall rate
    double urFF = -sqrt(fmax(0.0, -(1.0 + bl_gcon[0][0]) * bl_gcon[1][1]));
    ur = bl_Ucon_tmp[1] + infall_factor * (urFF - bl_Ucon_tmp[1]);


    // Finally, get Ucon in BL coordinates
    K = bl_gcov[0][0] + 2*omega*bl_gcov[0][3] + omega*omega*bl_gcov[3][3];
    ut = sqrt(fmax(0.0, -(1. + ur*ur*bl_gcov[1][1]) / K));
    bl_Ucon[0] = ut;
    bl_Ucon[1] = ur;
    bl_Ucon[2] = 0.;
    bl_Ucon[3] = omega * ut;
  } else {
    // Outside r_isco, Keplerian
    omegaK = 1. / (pow(r, 3./2) + a);
    omegaFF = bl_gcon[0][3] / bl_gcon[0][0];

    // Compromise
    omega = omegaK + (1 - keplerian_factor)*(omegaFF - omegaK);
    // Set infall rate
    ur = infall_factor * -sqrt(fmax(0.0, -(1.0 + bl_gcon[0][0]) * bl_gcon[1][1]));

    // Get the normal observer velocity for Ucon/Ucov, in BL coordinates
    K = bl_gcov[0][0] + 2*omega*bl_gcov[0][3] + omega*omega*bl_gcov[3][3];
    ut = sqrt(fmax(0.0, -(1. + ur*ur*bl_gcov[1][1]) / K));
    bl_Ucon[0] = ut;
    bl_Ucon[1] = ur;
    bl_Ucon[2] = 0.;
    bl_Ucon[3] = omega * ut;
  }
}

void get_model_B_spatial(double r, double th, double B_vector[4])
{

    double gcov_c[NDIM][NDIM];
    gcov_bl(r, th, gcov_c);
  
  double gdet3 = gcov_c[1][1]*gcov_c[2][2]*(gcov_c[3][3]*gcov_c[0][0] - gcov_c[0][3]*gcov_c[0][3]);
  double g = (sqrt(fabs(gdet3)));


  int bfield = 2;  
  double theta = th;
  double B_1,B_2,B_3;

  if(bfield==1){
  //TOROIDAL FIELD

      B_1 = 0.0;
      B_2 = 0.0;
      B_3 = 1.;
  }

  else  if(bfield==2){
  //Vertical

      B_1 = r*cos(th);
      B_2 = -sin(th);
      B_3 = 0.0 ;
  }

  else if (bfield==3){
  //Radial

      B_1 = -(sin(th) * sign(cos(th)));
      B_2 = 0.;
//       fprintf(stderr, "%g", g);
      B_3 = 0.;
  }

  B_vector[0] = 0;
  B_vector[1] = B_1/g;
  B_vector[2] = B_2/g;
  B_vector[3] = B_3/g;
  // fprintf(stderr, "%f, %f, %f\n", B_1, B_2, B_3);
}


void get_model_bcon(double X[NDIM], double Bcon[NDIM], double Ucon[NDIM], double Ucov[NDIM])
{

  double B_spatial[4], B_for_norm[4];
  double B_1, B_2, B_3;
  double r, th;

  bl_coord(X, &r, &th);
  get_model_B_spatial(r, th, B_spatial);

  B_1 = B_spatial[1];
  B_2 = B_spatial[2];
  B_3 = B_spatial[3];

  // Normalize to the value measured at the photon sphere (i.e sufficiently far away from horizon)
  get_model_B_spatial(3., M_PI/2., B_for_norm);

  double norm = sqrt(B_for_norm[1]*B_for_norm[1] + B_for_norm[2]*B_for_norm[2]*r*r + B_for_norm[3]*B_for_norm[3]*r*r);

  double b00 = 30; // 5G for M87*

  B_1 = b00*B_1/norm; 
  B_2 = b00*B_2/norm;
  B_3 = b00*B_3/norm;

  Bcon[0] = B_1  * Ucov[1] + B_2 * Ucov[2] + B_3 * Ucov[3] ;
  Bcon[1] = (B_1 + Bcon[0] * Ucon[1]) / Ucon[0];
  Bcon[2] = (B_2 + Bcon[0] * Ucon[2]) / Ucon[0];
  Bcon[3] = (B_3 + Bcon[0] * Ucon[3]) / Ucon[0];
  
  }

double get_spot(double X[NDIM]){
    double r, th;
    double Xbl[4];
    bl_coord(X, &r, &th);

    double tau = 0.;

    double r_current2 = r;

    double rplus = 1. + sqrt(1. - a * a);
    double rmin  = 1. - sqrt(1. - a * a);

    double ph = X[3] - (a / (2. * sqrt(1. - a * a)) * log((r_current2 - rplus) / (r_current2 - rmin))); //probably dont need this. Come back later on to this.
    double ne = 0;
    double theta_hs = M_PI/2;
    double r_hs = 6.;
    double omega_hs = 1./(pow(r_hs, 1.5) + a);
    double P = 2*M_PI/omega_hs;
    double phi_hs = omega_hs*tau;
    
    double Rspot = 1.;
   
    double Xspot[NDIM];
    Xspot[0] = 0.;
    Xspot[1] = sqrt(r_hs*r_hs + a*a)*cos(phi_hs);
    Xspot[2] = sqrt(r_hs*r_hs + a*a)*sin(phi_hs) - sqrt(r_hs*r_hs + a*a)*sin(M_PI/2);
    Xspot[3] = 0.;
    
    // fprintf(stderr, "spot %f %f %f\n", Xspot[1], Xspot[2], Xspot[3]);
    double Xn[NDIM];
    Xn[0] = 0;
    Xn[1] = sqrt(r*r + a*a)*cos(ph)*sin(th);
    Xn[2] = sqrt(r*r + a*a)*sin(ph)*sin(th);
    Xn[3] = sqrt(r*r + a*a)*cos(th);

    // fprintf(stderr, "BL: %f %f %f %f\n",Xspot[1] - Xn[1], Xn[1], Xspot[2], Xn[2]);

    double pos2 = fabs((Xspot[1] - Xn[1])*(Xspot[1] - Xn[1]) + 
                       (Xspot[2] - Xn[2])*(Xspot[2] - Xn[2]) + 
                       (Xspot[3] - Xn[3])*(Xspot[3] - Xn[3]));

    double Ne = 0;

    // Spot within 4*Rspot and z>0
    if (pos2<16*Rspot*Rspot &&  r*cos(th)>0){
        Ne = exp(-pos2/2./Rspot/Rspot);
         // fprintf(stderr, "BL: %f %f %f %f\n",Ne, pos2, Xspot[1], Xspot[2]);
    }
    return Ne;
}




int get_fluid_params(double X[NDIM], struct GRMHD *modvar) {

    int i, j, k;   
    double del[NDIM];
    double rho, uu;
    double Bp[4], V_u[NDIM], Vfac, VdotV, UdotBp;
    double g_uu[NDIM][NDIM], g_dd[NDIM][NDIM], coeff[4];
    double bsq, beta, beta_trans, b2, trat, Th_unit, two_temp_gam;
    double Ucov[NDIM], Ucon[NDIM];

    if (X[1] < startx[1] || X[1] > stopx[1] || X[2] < startx[2] ||
       X[2] > stopx[2]) {
       // fprintf(stderr, "%f %f %f %f %f %f\n", X[1], X[2], startx[1], stopx[1], startx[2], stopx[2]);
       (*modvar).n_e = 0.;
      return 0;
    }
    // fprintf(stderr, "%f %f %f %f %f %f\n", X[1], X[2], startx[1], stopx[1], startx[2], stopx[2]);
       
    double r, th;
    bl_coord(X, &r, &th);

    //////if RIAF///////

    double zc=r*cos(th);
    double rc=r*sin(th);
    double HH = 0.5;
    double nth0 = 4.;
    double ne_spot=nth0*exp(-zc*zc/2./rc/rc/HH/HH)*pow(r,-1.5)*Ne_unit;
    ////////////////////
    // double nth0 = 8.;
    // double ne_spot = get_spot(X)*nth0*Ne_unit;
  
    double Te_bg = 6;
    double Te = Te_bg * get_spot(X) * pow(r, -0.84)* Te_unit;
    // double Te = Te_bg*pow(r, -0.84) * Te_unit;
    double theta_E = BOLTZMANN_CONSTANT*Te/ELECTRON_MASS/SPEED_OF_LIGHT/SPEED_OF_LIGHT;

    // double thetaE_0 = 200;
    // double theta_E =  thetaE_0 * get_spot(X) * pow(r, -0.84);
    
    // if (ne_spot>0) {
    //      fprintf(stderr, "BL: %f %f\n", get_spot(X),  theta_E);
    // }

    (*modvar).n_e = ne_spot ;
    (*modvar).theta_e = theta_E;

    double bl_gcov[NDIM][NDIM], bl_gcon[NDIM][NDIM];

    gcon_bl(r, th, bl_gcon);
    gcov_bl(r, th, bl_gcov);

    double bl_Ucon[NDIM], bl_Ucov[NDIM] ;
    get_model_ucon(X, bl_Ucon);
    lower_index(X, bl_Ucon, bl_Ucov);

    double KS_Utemp[4];

    BL_to_KS(X, bl_Ucon, Ucon);

    // vec_from_ks(X, KS_Utemp, Ucon); //Probably need this from KS to native coords.

    for (i = 0; i < NDIM; i++){
        (*modvar).U_u[i] = Ucon[i];
    }

    lower_index(X, (*modvar).U_u, (*modvar).U_d);
 
    double bl_Bcon[4],bl_Bcov[4], dummy_B_ks[4], ucov_dummy[4];
    double Bcov[NDIM], Bcon[NDIM];

    get_model_bcon(X, bl_Bcon, bl_Ucon, bl_Ucov);

    lower_index(X, bl_Bcon, bl_Bcov);

    BL_to_KS(X, bl_Bcon, Bcon);

    // vec_from_ks(X, dummy_B_ks, Bcon);   //Probably need this from KS to native coords.

    for (i = 0; i < NDIM; i++){
        (*modvar).B_u[i] = Bcon[i];
    }

    lower_index(X, (*modvar).B_u, (*modvar).B_d);

    bsq =     bl_Bcov[0] * bl_Bcon[0] +
              bl_Bcov[1] * bl_Bcon[1] +
              bl_Bcov[2] * bl_Bcon[2] +
              bl_Bcov[3] * bl_Bcon[3];

    (*modvar).B = sqrt(bsq) * B_unit + 1e-40;

    // (*modvar).B = 5 * get_spot(X) * B_unit;

    // fprintf(stderr, "Bcon: %f, %f, %f, %f\n", bl_Bcon[0], bl_Bcon[1], bl_Bcon[2], bl_Bcon[3]);
    // fprintf(stderr, "Bcov: %f, %f, %f, %f\n", Bcov[0], Bcov[1], Bcov[2], Bcov[3]);

// #if (DEBUG)

    if (uu < 0)
        fprintf(stderr, "U %e %e\n", uu, p[UU][i][j][k]);
    ;

    if ((*modvar).theta_e < 0)
        fprintf(stderr, "Te %e\n", (*modvar).theta_e);
    if ((*modvar).B < 0)
        fprintf(stderr, "B %e\n", (*modvar).B);
    if ((*modvar).n_e < 0)
        fprintf(stderr, "Ne %e %e\n", (*modvar).n_e, p[KRHO][i][j][k]);
// #endif


    return 1;
}



void set_units(double M_unit_) {

    L_unit = GGRAV * MBH / (SPEED_OF_LIGHT * SPEED_OF_LIGHT);
    T_unit = L_unit / SPEED_OF_LIGHT;
  
    RHO_unit = Ne_unit * (PROTON_MASS + ELECTRON_MASS);

    M_unit_ = RHO_unit * pow(L_unit, 3);

    U_unit = RHO_unit * SPEED_OF_LIGHT * SPEED_OF_LIGHT;
    B_unit = SPEED_OF_LIGHT * sqrt(4. * M_PI * RHO_unit);
}



void Xtoijk(double X[NDIM], int *i, int *j, int *k, double del[NDIM]) {
    double phi;

    phi = fmod(X[3], stopx[3]);
    if (phi < 0.)
        phi = stopx[3] + phi;

    *i = (int)((X[1] - startx[1]) / dx[1] - 0.5 + 1000) - 1000;
    *j = (int)((X[2] - startx[2]) / dx[2] - 0.5 + 1000) - 1000;
    *k = (int)(phi / dx[3] + 1000 - 0.5) - 1000;

    if (*i < 0) {
        *i = 0;
        del[1] = 0.;
    } else if (*i > N1 - 2) {
        *i = N1 - 2;
        del[1] = 1.;
    } else {
        del[1] = (X[1] - ((*i + 0.5) * dx[1] + startx[1])) / dx[1];
    }

    if (*j < 0) {
        *j = 0;
        del[2] = 0.;
    } else if (*j > N2 - 2) {
        *j = N2 - 2;
        del[2] = 1.;
    } else {
        del[2] = (X[2] - ((*j + 0.5) * dx[2] + startx[2])) / dx[2];
    }

    if (*k < 0) {
        *k = 0;
        del[3] = 0.;
    } else if (*k > N3 - 2) {
        *k = N3 - 2;
        del[3] = 1.;
    } else {
        del[3] = (phi - ((*k + 0.5) * dx[3])) / dx[3];
    }
    if (del[3] < 0)
        fprintf(stderr, "%e %e %e %d\n", del[3], phi, (*k + 0.5) * dx[3], *k);

    return;
}

static void *malloc_rank1(int n1, int alloc_size) {
    void *A;

    if ((A = (void *)malloc(n1 * alloc_size)) == NULL) {
        fprintf(stderr, "malloc failure in malloc_rank1\n");
        exit(123);
    }

    return A;
}

double ***malloc_rank3(int n1, int n2, int n3) {
    double ***A;
    double *space;
    int i, j;

    space = malloc_rank1(n1 * n2 * n3, sizeof(double));

    A = malloc_rank1(n1, sizeof(double *));

    for (i = 0; i < n1; i++) {
        A[i] = malloc_rank1(n2, sizeof(double *));
        for (j = 0; j < n2; j++) {
            A[i][j] = &(space[n3 * (j + n2 * i)]);
        }
    }

    return A;
}

void init_storage(void) {
    int i, j, k;

    p = (double ****)malloc(NPRIM * sizeof(double ***));
    for (i = 0; i < NPRIM; i++) {
        p[i] = (double ***)malloc(N1 * sizeof(double **));
        for (j = 0; j < N1; j++) {
            p[i][j] = (double **)malloc(N2 * sizeof(double *));
            for (k = 0; k < N2; k++) {
                p[i][j][k] = (double *)malloc(N3 * sizeof(double));
            }
        }
    }

    return;
}

void compute_spec_user(struct Camera *intensityfield,
                       double energy_spectrum[num_frequencies][nspec]) {

    return;
}
