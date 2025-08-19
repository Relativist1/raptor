//
//  raptor_metal.m
//  RAPTOR Metal GPU Implementation
//
//  GPU-accelerated version of the RAPTOR ray-tracing code
//  for MacBook Pro M2 using Metal Performance Shaders
//

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <MetalPerformanceShaders/MetalPerformanceShaders.h>

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

// GPU kernel for geodesic integration
static NSString* geodesicKernelSource = @"
#include <metal_stdlib>
using namespace metal;

struct GeodesicState {
    double4 position;    // t, r, theta, phi
    double4 momentum;    // k_t, k_r, k_theta, k_phi
};

struct ModelParams {
    double MBH;
    double M_UNIT;
    double a;  // black hole spin
    double R_ISCO;
};

// Kerr metric components
double4 compute_metric_components(double4 pos, double a, double M) {
    double r = pos.y;
    double theta = pos.z;
    double sin_theta = sin(theta);
    double cos_theta = cos(theta);
    double rho2 = r*r + a*a*cos_theta*cos_theta;
    double Delta = r*r - 2*M*r + a*a;
    
    double4 g = double4(0.0);
    g.x = -(1.0 - 2*M*r/rho2);  // g_tt
    g.y = rho2/Delta;            // g_rr
    g.z = rho2;                  // g_theta_theta
    g.w = (r*r + a*a + 2*M*r*a*a*sin_theta*sin_theta/rho2)*sin_theta*sin_theta; // g_phi_phi
    
    return g;
}

// Christoffel symbols for Kerr metric
void compute_christoffel_symbols(double4 pos, double a, double M, 
                                thread double4x4x4& Gamma) {
    double r = pos.y;
    double theta = pos.z;
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
    
    // Non-zero Christoffel symbols for Kerr metric
    // This is a simplified version - full implementation would include all components
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
double4 geodesic_equation(double4 pos, double4 mom, double a, double M) {
    double4x4x4 Gamma;
    compute_christoffel_symbols(pos, a, M, Gamma);
    
    double4 acceleration = double4(0.0);
    
    // Compute acceleration using Christoffel symbols
    for (int mu = 0; mu < 4; mu++) {
        for (int nu = 0; nu < 4; nu++) {
            for (int sigma = 0; sigma < 4; sigma++) {
                acceleration[mu] -= Gamma[mu][nu][sigma] * mom[nu] * mom[sigma];
            }
        }
    }
    
    return acceleration;
}

// RK4 integrator for geodesic equations
void rk4_step(thread GeodesicState& state, double dt, double a, double M) {
    double4 k1_pos, k1_mom, k2_pos, k2_mom, k3_pos, k3_mom, k4_pos, k4_mom;
    double4 k1_acc, k2_acc, k3_acc, k4_acc;
    
    // Step 1
    k1_pos = state.momentum * dt;
    k1_acc = geodesic_equation(state.position, state.momentum, a, M);
    k1_mom = k1_acc * dt;
    
    // Step 2
    double4 pos2 = state.position + 0.5 * k1_pos;
    double4 mom2 = state.momentum + 0.5 * k1_mom;
    k2_pos = mom2 * dt;
    k2_acc = geodesic_equation(pos2, mom2, a, M);
    k2_mom = k2_acc * dt;
    
    // Step 3
    double4 pos3 = state.position + 0.5 * k2_pos;
    double4 mom3 = state.momentum + 0.5 * k2_mom;
    k3_pos = mom3 * dt;
    k3_acc = geodesic_equation(pos3, mom3, a, M);
    k3_mom = k3_acc * dt;
    
    // Step 4
    double4 pos4 = state.position + k3_pos;
    double4 mom4 = state.momentum + k3_mom;
    k4_pos = mom4 * dt;
    k4_acc = geodesic_equation(pos4, mom4, a, M);
    k4_mom = k4_acc * dt;
    
    // Update state
    state.position += (k1_pos + 2.0*k2_pos + 2.0*k3_pos + k4_pos) / 6.0;
    state.momentum += (k1_mom + 2.0*k2_mom + 2.0*k3_mom + k4_mom) / 6.0;
}

kernel void integrate_geodesics(device const double* alpha_beta,
                               device double* lightpaths,
                               device int* steps_array,
                               constant ModelParams& params,
                               constant double& dt,
                               constant int& max_steps,
                               constant double& cutoff_inner,
                               uint pixel_id [[thread_position_in_grid]]) {
    
    double alpha = alpha_beta[pixel_id * 2];
    double beta = alpha_beta[pixel_id * 2 + 1];
    
    // Initialize geodesic state
    GeodesicState state;
    
    // Set initial position (observer at infinity)
    state.position = double4(0.0, 1e6, M_PI/2.0, 0.0);  // t, r, theta, phi
    
    // Set initial momentum based on alpha, beta (impact parameters)
    double r_obs = 1e6;
    double sin_theta_obs = sin(M_PI/2.0);
    double cos_theta_obs = cos(M_PI/2.0);
    
    // Convert alpha, beta to initial momentum components
    // This is a simplified version - full implementation would use proper coordinate transformations
    state.momentum = double4(1.0, -alpha/r_obs, -beta/(r_obs*sin_theta_obs), 0.0);
    
    // Normalize momentum to null condition
    double4 g = compute_metric_components(state.position, params.a, params.MBH);
    double norm = 0.0;
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            norm += g[i] * state.momentum[i] * state.momentum[j];
        }
    }
    state.momentum = state.momentum / sqrt(-norm);
    
    int steps = 0;
    int lightpath_offset = pixel_id * max_steps * 9;
    
    // Store initial state
    lightpaths[lightpath_offset + steps * 9 + 0] = state.position.x;  // t
    lightpaths[lightpath_offset + steps * 9 + 1] = state.position.y;  // r
    lightpaths[lightpath_offset + steps * 9 + 2] = state.position.z;  // theta
    lightpaths[lightpath_offset + steps * 9 + 3] = state.position.w;  // phi
    lightpaths[lightpath_offset + steps * 9 + 4] = state.momentum.x;  // k_t
    lightpaths[lightpath_offset + steps * 9 + 5] = state.momentum.y;  // k_r
    lightpaths[lightpath_offset + steps * 9 + 6] = state.momentum.z;  // k_theta
    lightpaths[lightpath_offset + steps * 9 + 7] = state.momentum.w;  // k_phi
    lightpaths[lightpath_offset + steps * 9 + 8] = 0.0;  // affine parameter
    steps++;
    
    // Integrate geodesic
    while (steps < max_steps && state.position.y > cutoff_inner) {
        rk4_step(state, dt, params.a, params.MBH);
        
        // Store state
        lightpaths[lightpath_offset + steps * 9 + 0] = state.position.x;
        lightpaths[lightpath_offset + steps * 9 + 1] = state.position.y;
        lightpaths[lightpath_offset + steps * 9 + 2] = state.position.z;
        lightpaths[lightpath_offset + steps * 9 + 3] = state.position.w;
        lightpaths[lightpath_offset + steps * 9 + 4] = state.momentum.x;
        lightpaths[lightpath_offset + steps * 9 + 5] = state.momentum.y;
        lightpaths[lightpath_offset + steps * 9 + 6] = state.momentum.z;
        lightpaths[lightpath_offset + steps * 9 + 7] = state.momentum.w;
        lightpaths[lightpath_offset + steps * 9 + 8] = steps * dt;
        steps++;
    }
    
    steps_array[pixel_id] = steps;
}
";

// GPU kernel for radiative transfer
static NSString* radiativeTransferKernelSource = @"
#include <metal_stdlib>
using namespace metal;

struct StokesParameters {
    double I, Q, U, V;
};

// Synchrotron emission coefficients for kappa distribution
double j_synch_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
    double w = (kappa-3)*theta_e/kappa;
    double nuc = ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    
    // Simplified emission coefficient
    double factor = pow(X, 1.0/3.0) * exp(-pow(X, 0.84) / 30.0);
    return factor * n_e * pow(theta_e, 2.0) * pow(B, 1.5);
}

// Absorption coefficient
double alpha_synch_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {
    double w = (kappa-3)*theta_e/kappa;
    double nuc = ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);
    double X = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));
    
    // Simplified absorption coefficient
    double factor = pow(X, -2.0/3.0) * exp(-pow(X, 0.84) / 30.0);
    return factor * n_e * pow(theta_e, 1.5) * pow(B, 2.5);
}

kernel void compute_radiative_transfer(device const double* lightpaths,
                                     device const int* steps_array,
                                     device const double* frequencies,
                                     device double* IQUV_results,
                                     device double* tau_results,
                                     device double* tauF_results,
                                     constant double& B_field,
                                     constant double& n_e,
                                     constant double& theta_e,
                                     uint pixel_id [[thread_position_in_grid]]) {
    
    int steps = steps_array[pixel_id];
    if (steps == 0) return;
    
    int lightpath_offset = pixel_id * MAX_STEPS * 9;
    
    // Initialize Stokes parameters
    StokesParameters stokes = {0.0, 0.0, 0.0, 0.0};
    double tau = 0.0;
    double tauF = 0.0;
    
    // For each frequency
    for (int f = 0; f < num_frequencies; f++) {
        double nu = frequencies[f];
        
        // Reset for this frequency
        stokes.I = 0.0;
        stokes.Q = 0.0;
        stokes.U = 0.0;
        stokes.V = 0.0;
        tau = 0.0;
        tauF = 0.0;
        
        // Integrate along light path
        for (int step = 0; step < steps; step++) {
            int offset = lightpath_offset + step * 9;
            
            double r = lightpaths[offset + 1];
            double theta = lightpaths[offset + 2];
            double B_local = B_field * pow(r, -2.0);  // Simplified B field scaling
            double n_e_local = n_e * pow(r, -2.0);    // Simplified density scaling
            double theta_e_local = theta_e * pow(r, -1.0); // Simplified temperature scaling
            
            // Compute emission and absorption coefficients
            double j_nu = j_synch_kappa(theta_e_local, n_e_local, nu, B_local, M_PI/4.0);
            double alpha_nu = alpha_synch_kappa(theta_e_local, n_e_local, nu, B_local, M_PI/4.0);
            
            // Simple radiative transfer integration
            double ds = 1.0;  // Simplified step size
            double dtau = alpha_nu * ds;
            
            // Update optical depth
            tau += dtau;
            
            // Update Stokes parameters (simplified)
            double intensity = j_nu / alpha_nu * (1.0 - exp(-dtau));
            stokes.I += intensity * exp(-tau);
            
            // Update cumulative optical depth for Faraday rotation
            tauF += dtau * pow(nu, -2.0);
        }
        
        // Store results for this frequency
        int result_offset = (pixel_id * num_frequencies + f) * 4;
        IQUV_results[result_offset + 0] = stokes.I;
        IQUV_results[result_offset + 1] = stokes.Q;
        IQUV_results[result_offset + 2] = stokes.U;
        IQUV_results[result_offset + 3] = stokes.V;
        
        tau_results[pixel_id * num_frequencies + f] = tau;
        tauF_results[pixel_id * num_frequencies + f] = tauF;
    }
}
";

@interface RAPTORMetal : NSObject

@property (nonatomic, strong) id<MTLDevice> device;
@property (nonatomic, strong) id<MTLCommandQueue> commandQueue;
@property (nonatomic, strong) id<MTLLibrary> library;
@property (nonatomic, strong) id<MTLComputePipelineState> geodesicPipeline;
@property (nonatomic, strong) id<MTLComputePipelineState> radiativeTransferPipeline;

- (instancetype)init;
- (void)setupPipelines;
- (void)runRayTracing:(ModelParameters)params 
            alphaBeta:(double*)alphaBeta 
            frequencies:(double*)frequencies
            numFrequencies:(int)numFrequencies
            results:(CameraPixel*)results;

@end

@implementation RAPTORMetal

- (instancetype)init {
    self = [super init];
    if (self) {
        self.device = MTLCreateSystemDefaultDevice();
        if (!self.device) {
            NSLog(@"Metal is not supported on this device");
            return nil;
        }
        
        self.commandQueue = [self.device newCommandQueue];
        [self setupPipelines];
    }
    return self;
}

- (void)setupPipelines {
    NSError *error = nil;
    
    // Create library from source
    self.library = [self.device newLibraryWithSource:geodesicKernelSource options:nil error:&error];
    if (!self.library) {
        NSLog(@"Failed to create library: %@", error);
        return;
    }
    
    // Create geodesic pipeline
    id<MTLFunction> geodesicFunction = [self.library newFunctionWithName:@"integrate_geodesics"];
    self.geodesicPipeline = [self.device newComputePipelineStateWithFunction:geodesicFunction error:&error];
    if (!self.geodesicPipeline) {
        NSLog(@"Failed to create geodesic pipeline: %@", error);
        return;
    }
    
    // Create radiative transfer pipeline
    id<MTLLibrary> rtLibrary = [self.device newLibraryWithSource:radiativeTransferKernelSource options:nil error:&error];
    if (!rtLibrary) {
        NSLog(@"Failed to create RT library: %@", error);
        return;
    }
    
    id<MTLFunction> rtFunction = [rtLibrary newFunctionWithName:@"compute_radiative_transfer"];
    self.radiativeTransferPipeline = [self.device newComputePipelineStateWithFunction:rtFunction error:&error];
    if (!self.radiativeTransferPipeline) {
        NSLog(@"Failed to create RT pipeline: %@", error);
        return;
    }
}

- (void)runRayTracing:(ModelParameters)params 
            alphaBeta:(double*)alphaBeta 
            frequencies:(double*)frequencies
            numFrequencies:(int)numFrequencies
            results:(CameraPixel*)results {
    
    int totalPixels = params.IMG_WIDTH * params.IMG_HEIGHT;
    
    // Create buffers
    id<MTLBuffer> alphaBetaBuffer = [self.device newBufferWithBytes:alphaBeta 
                                                             length:totalPixels * 2 * sizeof(double) 
                                                            options:MTLResourceStorageModeShared];
    
    id<MTLBuffer> frequenciesBuffer = [self.device newBufferWithBytes:frequencies 
                                                               length:numFrequencies * sizeof(double) 
                                                              options:MTLResourceStorageModeShared];
    
    id<MTLBuffer> lightpathsBuffer = [self.device newBufferWithLength:totalPixels * MAX_STEPS * 9 * sizeof(double) 
                                                               options:MTLResourceStorageModeShared];
    
    id<MTLBuffer> stepsBuffer = [self.device newBufferWithLength:totalPixels * sizeof(int) 
                                                          options:MTLResourceStorageModeShared];
    
    id<MTLBuffer> IQUVBuffer = [self.device newBufferWithLength:totalPixels * numFrequencies * 4 * sizeof(double) 
                                                          options:MTLResourceStorageModeShared];
    
    id<MTLBuffer> tauBuffer = [self.device newBufferWithLength:totalPixels * numFrequencies * sizeof(double) 
                                                        options:MTLResourceStorageModeShared];
    
    id<MTLBuffer> tauFBuffer = [self.device newBufferWithLength:totalPixels * numFrequencies * sizeof(double) 
                                                         options:MTLResourceStorageModeShared];
    
    // Create command buffer
    id<MTLCommandBuffer> commandBuffer = [self.commandQueue commandBuffer];
    
    // Step 1: Integrate geodesics
    id<MTLComputeCommandEncoder> geodesicEncoder = [commandBuffer computeCommandEncoder];
    [geodesicEncoder setComputePipelineState:self.geodesicPipeline];
    [geodesicEncoder setBuffer:alphaBetaBuffer offset:0 atIndex:0];
    [geodesicEncoder setBuffer:lightpathsBuffer offset:0 atIndex:1];
    [geodesicEncoder setBuffer:stepsBuffer offset:0 atIndex:2];
    
    // Set constants
    ModelParams modelParams = {params.MBH, params.M_UNIT, 0.9, 6.0}; // a=0.9, R_ISCO=6
    [geodesicEncoder setBytes:&modelParams length:sizeof(ModelParams) atIndex:3];
    
    double dt = params.STEPSIZE;
    [geodesicEncoder setBytes:&dt length:sizeof(double) atIndex:4];
    
    int maxSteps = MAX_STEPS;
    [geodesicEncoder setBytes:&maxSteps length:sizeof(int) atIndex:5];
    
    double cutoff = CUTOFF_INNER;
    [geodesicEncoder setBytes:&cutoff length:sizeof(double) atIndex:6];
    
    // Dispatch geodesic computation
    MTLSize gridSize = MTLSizeMake(totalPixels, 1, 1);
    MTLSize threadgroupSize = MTLSizeMake(64, 1, 1);
    [geodesicEncoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
    [geodesicEncoder endEncoding];
    
    // Step 2: Compute radiative transfer
    id<MTLComputeCommandEncoder> rtEncoder = [commandBuffer computeCommandEncoder];
    [rtEncoder setComputePipelineState:self.radiativeTransferPipeline];
    [rtEncoder setBuffer:lightpathsBuffer offset:0 atIndex:0];
    [rtEncoder setBuffer:stepsBuffer offset:0 atIndex:1];
    [rtEncoder setBuffer:frequenciesBuffer offset:0 atIndex:2];
    [rtEncoder setBuffer:IQUVBuffer offset:0 atIndex:3];
    [rtEncoder setBuffer:tauBuffer offset:0 atIndex:4];
    [rtEncoder setBuffer:tauFBuffer offset:0 atIndex:5];
    
    // Set plasma parameters (simplified)
    double B_field = 100.0;  // Gauss
    double n_e = 1e6;        // cm^-3
    double theta_e = 10.0;   // dimensionless temperature
    [rtEncoder setBytes:&B_field length:sizeof(double) atIndex:6];
    [rtEncoder setBytes:&n_e length:sizeof(double) atIndex:7];
    [rtEncoder setBytes:&theta_e length:sizeof(double) atIndex:8];
    
    // Dispatch radiative transfer computation
    [rtEncoder dispatchThreads:gridSize threadsPerThreadgroup:threadgroupSize];
    [rtEncoder endEncoding];
    
    // Commit and wait for completion
    [commandBuffer commit];
    [commandBuffer waitUntilCompleted];
    
    // Copy results back
    double* IQUVData = [IQUVBuffer contents];
    double* tauData = [tauBuffer contents];
    double* tauFData = [tauFBuffer contents];
    
    for (int pixel = 0; pixel < totalPixels; pixel++) {
        results[pixel].alpha = alphaBeta[pixel * 2];
        results[pixel].beta = alphaBeta[pixel * 2 + 1];
        
        // For simplicity, use first frequency results
        int offset = pixel * numFrequencies * 4;
        results[pixel].IQUV[0] = IQUVData[offset + 0];  // I
        results[pixel].IQUV[1] = IQUVData[offset + 1];  // Q
        results[pixel].IQUV[2] = IQUVData[offset + 2];  // U
        results[pixel].IQUV[3] = IQUVData[offset + 3];  // V
        
        results[pixel].tau = tauData[pixel * numFrequencies];
        results[pixel].tauF = tauFData[pixel * numFrequencies];
    }
}

@end

// Main function for testing
int main(int argc, const char * argv[]) {
    @autoreleasepool {
        NSLog(@"RAPTOR Metal GPU Implementation");
        NSLog(@"Initializing Metal device...");
        
        RAPTORMetal *raptor = [[RAPTORMetal alloc] init];
        if (!raptor) {
            NSLog(@"Failed to initialize RAPTOR Metal");
            return 1;
        }
        
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
            .IMG_WIDTH = 256,
            .IMG_HEIGHT = 256,
            .CAM_SIZE_X = 40.0,            // GM/c^2
            .CAM_SIZE_Y = 40.0,            // GM/c^2
            .STEPSIZE = 0.1,
            .max_level = 4
        };
        
        int numFrequencies = 20;
        int totalPixels = params.IMG_WIDTH * params.IMG_HEIGHT;
        
        // Allocate arrays
        double *alphaBeta = malloc(totalPixels * 2 * sizeof(double));
        double *frequencies = malloc(numFrequencies * sizeof(double));
        CameraPixel *results = malloc(totalPixels * sizeof(CameraPixel));
        
        // Initialize alpha, beta (impact parameters)
        for (int i = 0; i < params.IMG_WIDTH; i++) {
            for (int j = 0; j < params.IMG_HEIGHT; j++) {
                int pixel = i * params.IMG_HEIGHT + j;
                alphaBeta[pixel * 2] = (i - params.IMG_WIDTH/2) * params.CAM_SIZE_X / params.IMG_WIDTH;
                alphaBeta[pixel * 2 + 1] = (j - params.IMG_HEIGHT/2) * params.CAM_SIZE_Y / params.IMG_HEIGHT;
            }
        }
        
        // Initialize frequencies
        for (int f = 0; f < numFrequencies; f++) {
            frequencies[f] = params.FREQ_MIN * pow(10.0, (double)f / params.FREQS_PER_DEC);
        }
        
        NSLog(@"Starting GPU ray tracing...");
        NSLog(@"Image size: %dx%d pixels", params.IMG_WIDTH, params.IMG_HEIGHT);
        NSLog(@"Number of frequencies: %d", numFrequencies);
        
        // Run ray tracing
        [raptor runRayTracing:params 
                   alphaBeta:alphaBeta 
                  frequencies:frequencies 
                numFrequencies:numFrequencies 
                      results:results];
        
        NSLog(@"Ray tracing completed!");
        
        // Print some results
        NSLog(@"Sample results for central pixel:");
        CameraPixel center = results[totalPixels/2];
        NSLog(@"I = %.6e, Q = %.6e, U = %.6e, V = %.6e", 
              center.IQUV[0], center.IQUV[1], center.IQUV[2], center.IQUV[3]);
        NSLog(@"tau = %.6e, tauF = %.6e", center.tau, center.tauF);
        
        // Clean up
        free(alphaBeta);
        free(frequencies);
        free(results);
        
        NSLog(@"That's all folks!");
    }
    return 0;
}