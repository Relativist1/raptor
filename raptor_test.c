//
//  raptor_test.c
//  Simple test version of RAPTOR for Linux compilation
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

// Constants matching the C implementation
#define DIM 4
#define MAX_STEPS 10000
#define CUTOFF_INNER 1.0
#define MSUN 1.989e33
#define KPCTOCM 3.086e21
#define ELECTRON_MASS 9.1093837015e-28
#define SPEED_OF_LIGHT 2.99792458e10
#define BOLTZMANN_CONSTANT 1.380649e-16
#define PLANCK_CONSTANT 6.62607015e-27
#define ELECTRON_CHARGE 4.80320425e-10
#define M_PI 3.14159265358979323846

// Model parameters structure
typedef struct {
    double MBH;
    double M_UNIT;
    double TIME_INIT;
    double INCLINATION;
    double R_HIGH;
    double R_LOW;
    double FREQS_PER_DEC;
    double FREQ_MIN;
    double FREQ_MAX;
    double SOURCE_DIST;
    int IMG_WIDTH;
    int IMG_HEIGHT;
    double CAM_SIZE_X;
    double CAM_SIZE_Y;
    double STEPSIZE;
    int max_level;
} ModelParameters;

// Camera structure for storing results
typedef struct {
    double alpha;
    double beta;
    double IQUV[4];  // Stokes parameters
    double tau;
    double tauF;
} CameraPixel;

// Geodesic state structure
typedef struct {
    double position[4];    // t, r, theta, phi
    double momentum[4];    // k_t, k_r, k_theta, k_phi
} GeodesicState;

// Kerr metric components
void compute_metric_components(double pos[4], double a, double M, double g[4]) {
    double r = pos[1];
    double theta = pos[2];
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r*r + a*a*cos_theta*cos_theta;
    double Delta = r*r - 2*M*r + a*a;
    
    g[0] = -(1.0 - 2*M*r/rho2);  // g_tt
    g[1] = rho2/Delta;            // g_rr
    g[2] = rho2;                  // g_theta_theta
    g[3] = (r*r + a*a + 2*M*r*a*a*sin_theta*sin_theta/rho2)*sin_theta*sin_theta; // g_phi_phi
}

// Christoffel symbols for Kerr metric (simplified)
void compute_christoffel_symbols(double pos[4], double a, double M, double Gamma[4][4][4]) {
    double r = pos[1];
    double theta = pos[2];
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r*r + a*a*cos_theta*cos_theta;
    double Delta = r*r - 2*M*r + a*a;
    
    // Initialize to zero
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                Gamma[i][j][k] = 0.0;
            }
        }
    }
    
    // Non-zero Christoffel symbols for Kerr metric (simplified)
    double drho2_dr = 2*r;
    double drho2_dtheta = -2*a*a*cos_theta*sin_theta;
    double dDelta_dr = 2*r - 2*M;
    
    // Gamma^r_tt
    Gamma[1][0][0] = -0.5 * (1.0/Delta) * (2*M/rho2 - 4*M*r*drho2_dr/(rho2*rho2));
    
    // Gamma^r_rr
    Gamma[1][1][1] = 0.5 * (1.0/Delta) * dDelta_dr;
    
    // Gamma^r_theta_theta
    Gamma[1][2][2] = -0.5 * (1.0/Delta) * drho2_dr;
    
    // Gamma^r_phi_phi
    Gamma[1][3][3] = -0.5 * (1.0/Delta) * (2*M*a*a*sin_theta*sin_theta/rho2 - 
                                            4*M*r*a*a*sin_theta*sin_theta*drho2_dr/(rho2*rho2));
    
    // Gamma^theta_rr
    Gamma[2][1][1] = 0.5 * (1.0/rho2) * drho2_dr;
    
    // Gamma^theta_theta_theta
    Gamma[2][2][2] = 0.5 * (1.0/rho2) * drho2_dtheta;
    
    // Gamma^theta_phi_phi
    Gamma[2][3][3] = -0.5 * (1.0/rho2) * (2*M*r*a*a*2*sin_theta*cos_theta/rho2 - 
                                           2*M*r*a*a*sin_theta*sin_theta*drho2_dtheta/(rho2*rho2));
}

// Geodesic equation: d^2x^mu/dlambda^2 + Gamma^mu_nu_sigma * (dx^nu/dlambda) * (dx^sigma/dlambda) = 0
void geodesic_equation(double pos[4], double mom[4], double a, double M, double acceleration[4]) {
    double Gamma[4][4][4];
    compute_christoffel_symbols(pos, a, M, Gamma);
    
    // Initialize acceleration to zero
    for (int mu = 0; mu < 4; mu++) {
        acceleration[mu] = 0.0;
    }
    
    // Compute acceleration using Christoffel symbols
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int sigma = 0; sigma < 4; sigma++) {
                acceleration[mu] -= Gamma[mu][nu][sigma] * mom[nu] * mom[sigma];
            }
        }
    }
}

// RK4 integrator for geodesic equations
void rk4_step(GeodesicState* state, double dt, double a, double M) {
    double k1_pos[4], k1_mom[4], k2_pos[4], k2_mom[4], k3_pos[4], k3_mom[4], k4_pos[4], k4_mom[4];
    double k1_acc[4], k2_acc[4], k3_acc[4], k4_acc[4];
    
    // Step 1
    for (int i = 0; i < 4; i++) {
        k1_pos[i] = state->momentum[i] * dt;
    }
    geodesic_equation(state->position, state->momentum, a, M, k1_acc);
    for (int i = 0; i < 4; i++) {
        k1_mom[i] = k1_acc[i] * dt;
    }
    
    // Step 2
    double pos2[4], mom2[4];
    for (int i = 0; i < 4; i++) {
        pos2[i] = state->position[i] + 0.5 * k1_pos[i];
        mom2[i] = state->momentum[i] + 0.5 * k1_mom[i];
    }
    for (int i = 0; i < 4; i++) {
        k2_pos[i] = mom2[i] * dt;
    }
    geodesic_equation(pos2, mom2, a, M, k2_acc);
    for (int i = 0; i < 4; i++) {
        k2_mom[i] = k2_acc[i] * dt;
    }
    
    // Step 3
    double pos3[4], mom3[4];
    for (int i = 0; i < 4; i++) {
        pos3[i] = state->position[i] + 0.5 * k2_pos[i];
        mom3[i] = state->momentum[i] + 0.5 * k2_mom[i];
    }
    for (int i = 0; i < 4; i++) {
        k3_pos[i] = mom3[i] * dt;
    }
    geodesic_equation(pos3, mom3, a, M, k3_acc);
    for (int i = 0; i < 4; i++) {
        k3_mom[i] = k3_acc[i] * dt;
    }
    
    // Step 4
    double pos4[4], mom4[4];
    for (int i = 0; i < 4; i++) {
        pos4[i] = state->position[i] + k3_pos[i];
        mom4[i] = state->momentum[i] + k3_mom[i];
    }
    for (int i = 0; i < 4; i++) {
        k4_pos[i] = mom4[i] * dt;
    }
    geodesic_equation(pos4, mom4, a, M, k4_acc);
    for (int i = 0; i < 4; i++) {
        k4_mom[i] = k4_acc[i] * dt;
    }
    
    // Update state
    for (int i = 0; i < 4; i++) {
        state->position[i] += (k1_pos[i] + 2.0*k2_pos[i] + 2.0*k3_pos[i] + k4_pos[i]) / 6.0;
        state->momentum[i] += (k1_mom[i] + 2.0*k2_mom[i] + 2.0*k3_mom[i] + k4_mom[i]) / 6.0;
    }
}

// Integrate geodesic for a single pixel
int integrate_geodesic(double alpha, double beta, double* lightpath, double dt, double a, double M, double cutoff_inner) {
    GeodesicState state;
    
    // Set initial position (observer at infinity)
    state.position[0] = 0.0;  // t
    state.position[1] = 1e6;  // r
    state.position[2] = M_PI/2.0;  // theta
    state.position[3] = 0.0;  // phi
    
    // Set initial momentum based on alpha, beta (impact parameters)
    double r_obs = 1e6;
    double sin_theta_obs = sin(M_PI/2.0);
    double cos_theta_obs = cos(M_PI/2.0);
    
    // Convert alpha, beta to initial momentum components (simplified)
    state.momentum[0] = 1.0;
    state.momentum[1] = -alpha/r_obs;
    state.momentum[2] = -beta/(r_obs*sin_theta_obs);
    state.momentum[3] = 0.0;
    
    // Normalize momentum to null condition
    double g[4];
    compute_metric_components(state.position, a, M, g);
    double norm = 0.0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            norm += g[i] * state.momentum[i] * state.momentum[j];
        }
    }
    double norm_factor = sqrt(-norm);
    for (int i = 0; i < 4; i++) {
        state.momentum[i] /= norm_factor;
    }
    
    int steps = 0;
    
    // Store initial state
    lightpath[steps * 9 + 0] = state.position[0];  // t
    lightpath[steps * 9 + 1] = state.position[1];  // r
    lightpath[steps * 9 + 2] = state.position[2];  // theta
    lightpath[steps * 9 + 3] = state.position[3];  // phi
    lightpath[steps * 9 + 4] = state.momentum[0];  // k_t
    lightpath[steps * 9 + 5] = state.momentum[1];  // k_r
    lightpath[steps * 9 + 6] = state.momentum[2];  // k_theta
    lightpath[steps * 9 + 7] = state.momentum[3];  // k_phi
    lightpath[steps * 9 + 8] = 0.0;  // affine parameter
    steps++;
    
    // Integrate geodesic
    while (steps < MAX_STEPS && state.position[1] > cutoff_inner) {
        rk4_step(&state, dt, a, M);
        
        // Store state
        lightpath[steps * 9 + 0] = state.position[0];
        lightpath[steps * 9 + 1] = state.position[1];
        lightpath[steps * 9 + 2] = state.position[2];
        lightpath[steps * 9 + 3] = state.position[3];
        lightpath[steps * 9 + 4] = state.momentum[0];
        lightpath[steps * 9 + 5] = state.momentum[1];
        lightpath[steps * 9 + 6] = state.momentum[2];
        lightpath[steps * 9 + 7] = state.momentum[3];
        lightpath[steps * 9 + 8] = steps * dt;
        steps++;
    }
    
    return steps;
}

// Simple radiative transfer (simplified)
void compute_radiative_transfer(double* lightpath, int steps, double frequency, 
                               double* IQUV, double* tau, double* tauF) {
    // Initialize Stokes parameters
    IQUV[0] = 0.0;  // I
    IQUV[1] = 0.0;  // Q
    IQUV[2] = 0.0;  // U
    IQUV[3] = 0.0;  // V
    *tau = 0.0;
    *tauF = 0.0;
    
    // Simple emission model (mock)
    for (int step = 0; step < steps; step++) {
        double r = lightpath[step * 9 + 1];
        if (r > 10.0) {  // Only emit from outer regions
            double intensity = 1.0 / (r * r);  // Simple 1/r^2 scaling
            IQUV[0] += intensity;
            *tau += 0.1;  // Simple optical depth
        }
    }
    
    *tauF = *tau * pow(frequency, -2.0);  // Faraday rotation
}

// Main function for testing
int main(int argc, const char * argv[]) {
    printf("RAPTOR Test Implementation\n");
    printf("==========================\n");
    
    // Set up model parameters
    ModelParameters params = {
        .MBH = 6.2e9 * MSUN,           // M87 black hole mass
        .M_UNIT = 1e9 * MSUN,          // Mass unit
        .TIME_INIT = 0.0,
        .INCLINATION = 17.0,           // degrees
        .R_HIGH = 1.0,
        .R_LOW = 1.0,
        .FREQS_PER_DEC = 4,
        .FREQ_MIN = 1e10,              // 10 GHz
        .FREQ_MAX = 1e15,              // 1000 GHz
        .SOURCE_DIST = 16.4 * KPCTOCM, // M87 distance
        .IMG_WIDTH = 64,               // Reduced for testing
        .IMG_HEIGHT = 64,              // Reduced for testing
        .CAM_SIZE_X = 40.0,            // GM/c^2
        .CAM_SIZE_Y = 40.0,            // GM/c^2
        .STEPSIZE = 0.1,
        .max_level = 4
    };
    
    printf("Model parameters:\n");
    printf("  Black hole mass: %.2e Msun\n", params.MBH / MSUN);
    printf("  Image size: %dx%d pixels\n", params.IMG_WIDTH, params.IMG_HEIGHT);
    printf("  Field of view: %.1f x %.1f GM/c^2\n", params.CAM_SIZE_X, params.CAM_SIZE_Y);
    printf("\n");
    
    int totalPixels = params.IMG_WIDTH * params.IMG_HEIGHT;
    double* lightpaths = malloc(totalPixels * MAX_STEPS * 9 * sizeof(double));
    int* steps_array = malloc(totalPixels * sizeof(int));
    CameraPixel* results = malloc(totalPixels * sizeof(CameraPixel));
    
    if (!lightpaths || !steps_array || !results) {
        printf("Memory allocation failed!\n");
        return 1;
    }
    
    printf("Starting geodesic integration...\n");
    
    // Process each pixel
    for (int i = 0; i < params.IMG_WIDTH; i++) {
        for (int j = 0; j < params.IMG_HEIGHT; j++) {
            int pixel = i * params.IMG_HEIGHT + j;
            
            // Calculate alpha, beta (impact parameters)
            double alpha = (i - params.IMG_WIDTH/2) * params.CAM_SIZE_X / params.IMG_WIDTH;
            double beta = (j - params.IMG_HEIGHT/2) * params.CAM_SIZE_Y / params.IMG_HEIGHT;
            
            // Integrate geodesic for this pixel
            int lightpath_offset = pixel * MAX_STEPS * 9;
            int steps = integrate_geodesic(alpha, beta, 
                                         &lightpaths[lightpath_offset], 
                                         params.STEPSIZE, 0.9, params.MBH, CUTOFF_INNER);
            
            steps_array[pixel] = steps;
            
            // Simple radiative transfer
            double IQUV[4];
            double tau, tauF;
            compute_radiative_transfer(&lightpaths[lightpath_offset], steps, 
                                     params.FREQ_MIN, IQUV, &tau, &tauF);
            
            // Store results
            results[pixel].alpha = alpha;
            results[pixel].beta = beta;
            for (int k = 0; k < 4; k++) {
                results[pixel].IQUV[k] = IQUV[k];
            }
            results[pixel].tau = tau;
            results[pixel].tauF = tauF;
            
            if (pixel % 100 == 0) {
                printf("Processed pixel %d/%d (alpha=%.2f, beta=%.2f, steps=%d)\n", 
                       pixel, totalPixels, alpha, beta, steps);
            }
        }
    }
    
    printf("\nGeodesic integration completed!\n");
    
    // Print some sample results
    printf("\nSample results:\n");
    printf("Central pixel (0,0):\n");
    CameraPixel center = results[totalPixels/2];
    printf("  I = %.6e, Q = %.6e, U = %.6e, V = %.6e\n", 
           center.IQUV[0], center.IQUV[1], center.IQUV[2], center.IQUV[3]);
    printf("  tau = %.6e, tauF = %.6e\n", center.tau, center.tauF);
    
    // Calculate total flux
    double total_flux = 0.0;
    for (int pixel = 0; pixel < totalPixels; pixel++) {
        total_flux += results[pixel].IQUV[0];
    }
    printf("\nTotal integrated flux: %.6e\n", total_flux);
    
    // Clean up
    free(lightpaths);
    free(steps_array);
    free(results);
    
    printf("\nTest completed successfully!\n");
    return 0;
}