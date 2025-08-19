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
static NSString* geodesicKernelSource = @"#include <metal_stdlib>\n"
"using namespace metal;\n"
"\n"
"struct GeodesicState {\n"
"    double4 position;    // t, r, theta, phi\n"
"    double4 momentum;    // k_t, k_r, k_theta, k_phi\n"
"};\n"
"\n"
"struct ModelParams {\n"
"    double MBH;\n"
"    double M_UNIT;\n"
"    double a;  // black hole spin\n"
"    double R_ISCO;\n"
"};\n"
"\n"
"// Kerr metric components\n"
"double4 compute_metric_components(double4 pos, double a, double M) {\n"
"    double r = pos.y;\n"
"    double theta = pos.z;\n"
"    double sin_theta = sin(theta);\n"
"    double cos_theta = cos(theta);\n"
"    double rho2 = r*r + a*a*cos_theta*cos_theta;\n"
"    double Delta = r*r - 2*M*r + a*a;\n"
"    \n"
"    double4 g = double4(0.0);\n"
"    g.x = -(1.0 - 2*M*r/rho2);  // g_tt\n"
"    g.y = rho2/Delta;            // g_rr\n"
"    g.z = rho2;                  // g_theta_theta\n"
"    g.w = (r*r + a*a + 2*M*r*a*a*sin_theta*sin_theta/rho2)*sin_theta*sin_theta; // g_phi_phi\n"
"    \n"
"    return g;\n"
"}\n"
"\n"
"// Christoffel symbols for Kerr metric\n"
"void compute_christoffel_symbols(double4 pos, double a, double M, \n"
"                                thread double4x4x4& Gamma) {\n"
"    double r = pos.y;\n"
"    double theta = pos.z;\n"
"    double sin_theta = sin(theta);\n"
"    double cos_theta = cos(theta);\n"
"    double rho2 = r*r + a*a*cos_theta*cos_theta;\n"
"    double Delta = r*r - 2*M*r + a*a;\n"
"    \n"
"    // Initialize to zero\n"
"    for (int i = 0; i < 4; i++) {\n"
"        for (int j = 0; j < 4; j++) {\n"
"            for (int k = 0; k < 4; k++) {\n"
"                Gamma[i][j][k] = 0.0;\n"
"            }\n"
"        }\n"
"    }\n"
"    \n"
"    // Non-zero Christoffel symbols for Kerr metric\n"
"    // This is a simplified version - full implementation would include all components\n"
"    double drho2_dr = 2*r;\n"
"    double drho2_dtheta = -2*a*a*cos_theta*sin_theta;\n"
"    double dDelta_dr = 2*r - 2*M;\n"
"    \n"
"    // Gamma^r_tt\n"
"    Gamma[1][0][0] = -0.5 * (1.0/Delta) * (2*M/rho2 - 4*M*r*drho2_dr/(rho2*rho2));\n"
"    \n"
"    // Gamma^r_rr\n"
"    Gamma[1][1][1] = 0.5 * (1.0/Delta) * dDelta_dr;\n"
"    \n"
"    // Gamma^r_theta_theta\n"
"    Gamma[1][2][2] = -0.5 * (1.0/Delta) * drho2_dr;\n"
"    \n"
"    // Gamma^r_phi_phi\n"
"    Gamma[1][3][3] = -0.5 * (1.0/Delta) * (2*M*a*a*sin_theta*sin_theta/rho2 - \n"
"                                            4*M*r*a*a*sin_theta*sin_theta*drho2_dr/(rho2*rho2));\n"
"    \n"
"    // Gamma^theta_rr\n"
"    Gamma[2][1][1] = 0.5 * (1.0/rho2) * drho2_dr;\n"
"    \n"
"    // Gamma^theta_theta_theta\n"
"    Gamma[2][2][2] = 0.5 * (1.0/rho2) * drho2_dtheta;\n"
"    \n"
"    // Gamma^theta_phi_phi\n"
"    Gamma[2][3][3] = -0.5 * (1.0/rho2) * (2*M*r*a*a*2*sin_theta*cos_theta/rho2 - \n"
"                                           2*M*r*a*a*sin_theta*sin_theta*drho2_dtheta/(rho2*rho2));\n"
"}\n"
"\n"
"// Geodesic equation: d^2x^mu/dlambda^2 + Gamma^mu_nu_sigma * (dx^nu/dlambda) * (dx^sigma/dlambda) = 0\n"
"double4 geodesic_equation(double4 pos, double4 mom, double a, double M) {\n"
"    double4x4x4 Gamma;\n"
"    compute_christoffel_symbols(pos, a, M, Gamma);\n"
"    \n"
"    double4 acceleration = double4(0.0);\n"
"    \n"
"    // Compute acceleration using Christoffel symbols\n"
"    for (int mu = 0; mu < 4; mu++) {\n"
"        for (int nu = 0; nu < 4; nu++) {\n"
"            for (int sigma = 0; sigma < 4; sigma++) {\n"
"                acceleration[mu] -= Gamma[mu][nu][sigma] * mom[nu] * mom[sigma];\n"
"            }\n"
"        }\n"
"    }\n"
"    \n"
"    return acceleration;\n"
"}\n"
"\n"
"// RK4 integrator for geodesic equations\n"
"void rk4_step(thread GeodesicState& state, double dt, double a, double M) {\n"
"    double4 k1_pos, k1_mom, k2_pos, k2_mom, k3_pos, k3_mom, k4_pos, k4_mom;\n"
"    double4 k1_acc, k2_acc, k3_acc, k4_acc;\n"
"    \n"
"    // Step 1\n"
"    k1_pos = state.momentum * dt;\n"
"    k1_acc = geodesic_equation(state.position, state.momentum, a, M);\n"
"    k1_mom = k1_acc * dt;\n"
"    \n"
"    // Step 2\n"
"    double4 pos2 = state.position + 0.5 * k1_pos;\n"
"    double4 mom2 = state.momentum + 0.5 * k1_mom;\n"
"    k2_pos = mom2 * dt;\n"
"    k2_acc = geodesic_equation(pos2, mom2, a, M);\n"
"    k2_mom = k2_acc * dt;\n"
"    \n"
"    // Step 3\n"
"    double4 pos3 = state.position + 0.5 * k2_pos;\n"
"    double4 mom3 = state.momentum + 0.5 * k2_mom;\n"
"    k3_pos = mom3 * dt;\n"
"    k3_acc = geodesic_equation(pos3, mom3, a, M);\n"
"    k3_mom = k3_acc * dt;\n"
"    \n"
"    // Step 4\n"
"    double4 pos4 = state.position + k3_pos;\n"
"    double4 mom4 = state.momentum + k3_mom;\n"
"    k4_pos = mom4 * dt;\n"
"    k4_acc = geodesic_equation(pos4, mom4, a, M);\n"
"    k4_mom = k4_acc * dt;\n"
"    \n"
"    // Update state\n"
"    state.position += (k1_pos + 2.0*k2_pos + 2.0*k3_pos + k4_pos) / 6.0;\n"
"    state.momentum += (k1_mom + 2.0*k2_mom + 2.0*k3_mom + k4_mom) / 6.0;\n"
"}\n"
"\n"
"kernel void integrate_geodesics(device const double* alpha_beta,\n"
"                               device double* lightpaths,\n"
"                               device int* steps_array,\n"
"                               constant ModelParams& params,\n"
"                               constant double& dt,\n"
"                               constant int& max_steps,\n"
"                               constant double& cutoff_inner,\n"
"                               uint pixel_id [[thread_position_in_grid]]) {\n"
"    \n"
"    double alpha = alpha_beta[pixel_id * 2];\n"
"    double beta = alpha_beta[pixel_id * 2 + 1];\n"
"    \n"
"    // Initialize geodesic state\n"
"    GeodesicState state;\n"
"    \n"
"    // Set initial position (observer at infinity)\n"
"    state.position = double4(0.0, 1e6, M_PI/2.0, 0.0);  // t, r, theta, phi\n"
"    \n"
"    // Set initial momentum based on alpha, beta (impact parameters)\n"
"    double r_obs = 1e6;\n"
"    double sin_theta_obs = sin(M_PI/2.0);\n"
"    double cos_theta_obs = cos(M_PI/2.0);\n"
"    \n"
"    // Convert alpha, beta to initial momentum components\n"
"    // This is a simplified version - full implementation would use proper coordinate transformations\n"
"    state.momentum = double4(1.0, -alpha/r_obs, -beta/(r_obs*sin_theta_obs), 0.0);\n"
"    \n"
"    // Normalize momentum to null condition\n"
"    double4 g = compute_metric_components(state.position, params.a, params.MBH);\n"
"    double norm = 0.0;\n"
"    for (int i = 0; i < 4; i++) {\n"
"        for (int j = 0; j < 4; j++) {\n"
"            norm += g[i] * state.momentum[i] * state.momentum[j];\n"
"        }\n"
"    }\n"
"    state.momentum = state.momentum / sqrt(-norm);\n"
"    \n"
"    int steps = 0;\n"
"    int lightpath_offset = pixel_id * max_steps * 9;\n"
"    \n"
"    // Store initial state\n"
"    lightpaths[lightpath_offset + steps * 9 + 0] = state.position.x;  // t\n"
"    lightpaths[lightpath_offset + steps * 9 + 1] = state.position.y;  // r\n"
"    lightpaths[lightpath_offset + steps * 9 + 2] = state.position.z;  // theta\n"
"    lightpaths[lightpath_offset + steps * 9 + 3] = state.position.w;  // phi\n"
"    lightpaths[lightpath_offset + steps * 9 + 4] = state.momentum.x;  // k_t\n"
"    lightpaths[lightpath_offset + steps * 9 + 5] = state.momentum.y;  // k_r\n"
"    lightpaths[lightpath_offset + steps * 9 + 6] = state.momentum.z;  // k_theta\n"
"    lightpaths[lightpath_offset + steps * 9 + 7] = state.momentum.w;  // k_phi\n"
"    lightpaths[lightpath_offset + steps * 9 + 8] = 0.0;  // affine parameter\n"
"    steps++;\n"
"    \n"
"    // Integrate geodesic\n"
"    while (steps < max_steps && state.position.y > cutoff_inner) {\n"
"        rk4_step(state, dt, params.a, params.MBH);\n"
"        \n"
"        // Store state\n"
"        lightpaths[lightpath_offset + steps * 9 + 0] = state.position.x;\n"
"        lightpaths[lightpath_offset + steps * 9 + 1] = state.position.y;\n"
"        lightpaths[lightpath_offset + steps * 9 + 2] = state.position.z;\n"
"        lightpaths[lightpath_offset + steps * 9 + 3] = state.position.w;\n"
"        lightpaths[lightpath_offset + steps * 9 + 4] = state.momentum.x;\n"
"        lightpaths[lightpath_offset + steps * 9 + 5] = state.momentum.y;\n"
"        lightpaths[lightpath_offset + steps * 9 + 6] = state.momentum.z;\n"
"        lightpaths[lightpath_offset + steps * 9 + 7] = state.momentum.w;\n"
"        lightpaths[lightpath_offset + steps * 9 + 8] = steps * dt;\n"
"        steps++;\n"
"    }\n"
"    \n"
"    steps_array[pixel_id] = steps;\n"
"}\n";

// GPU kernel for radiative transfer
static NSString* radiativeTransferKernelSource = @"#include <metal_stdlib>\n"
"using namespace metal;\n"
"\n"
"struct StokesParameters {\n"
"    double I, Q, U, V;\n"
"};\n"
"\n"
"// Synchrotron emission coefficients for kappa distribution\n"
"double j_synch_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {\n"
"    double w = (kappa-3)*theta_e/kappa;\n"
"    double nuc = ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);\n"
"    double X = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));\n"
"    \n"
"    // Simplified emission coefficient\n"
"    double factor = pow(X, 1.0/3.0) * exp(-pow(X, 0.84) / 30.0);\n"
"    return factor * n_e * pow(theta_e, 2.0) * pow(B, 1.5);\n"
"}\n"
"\n"
"// Absorption coefficient\n"
"double alpha_synch_kappa(double theta_e, double n_e, double nu, double B, double theta_B) {\n"
"    double w = (kappa-3)*theta_e/kappa;\n"
"    double nuc = ELECTRON_CHARGE * B / (2 * M_PI * ELECTRON_MASS * SPEED_OF_LIGHT);\n"
"    double X = nu / (nuc * pow(w * kappa, 2) * sin(theta_B));\n"
"    \n"
"    // Simplified absorption coefficient\n"
"    double factor = pow(X, -2.0/3.0) * exp(-pow(X, 0.84) / 30.0);\n"
"    return factor * n_e * pow(theta_e, 1.5) * pow(B, 2.5);\n"
"}\n"
"\n"
"kernel void compute_radiative_transfer(device const double* lightpaths,\n"
"                                     device const int* steps_array,\n"
"                                     device const double* frequencies,\n"
"                                     device double* IQUV_results,\n"
"                                     device double* tau_results,\n"
"                                     device double* tauF_results,\n"
"                                     constant double& B_field,\n"
"                                     constant double& n_e,\n"
"                                     constant double& theta_e,\n"
"                                     uint pixel_id [[thread_position_in_grid]]) {\n"
"    \n"
"    int steps = steps_array[pixel_id];\n"
"    if (steps == 0) return;\n"
"    \n"
"    int lightpath_offset = pixel_id * MAX_STEPS * 9;\n"
"    \n"
"    // Initialize Stokes parameters\n"
"    StokesParameters stokes = {0.0, 0.0, 0.0, 0.0};\n"
"    double tau = 0.0;\n"
"    double tauF = 0.0;\n"
"    \n"
"    // For each frequency\n"
"    for (int f = 0; f < num_frequencies; f++) {\n"
"        double nu = frequencies[f];\n"
"        \n"
"        // Reset for this frequency\n"
"        stokes.I = 0.0;\n"
"        stokes.Q = 0.0;\n"
"        stokes.U = 0.0;\n"
"        stokes.V = 0.0;\n"
"        tau = 0.0;\n"
"        tauF = 0.0;\n"
"        \n"
"        // Integrate along light path\n"
"        for (int step = 0; step < steps; step++) {\n"
"            int offset = lightpath_offset + step * 9;\n"
"            \n"
"            double r = lightpaths[offset + 1];\n"
"            double theta = lightpaths[offset + 2];\n"
"            double B_local = B_field * pow(r, -2.0);  // Simplified B field scaling\n"
"            double n_e_local = n_e * pow(r, -2.0);    // Simplified density scaling\n"
"            double theta_e_local = theta_e * pow(r, -1.0); // Simplified temperature scaling\n"
"            \n"
"            // Compute emission and absorption coefficients\n"
"            double j_nu = j_synch_kappa(theta_e_local, n_e_local, nu, B_local, M_PI/4.0);\n"
"            double alpha_nu = alpha_synch_kappa(theta_e_local, n_e_local, nu, B_local, M_PI/4.0);\n"
"            \n"
"            // Simple radiative transfer integration\n"
"            double ds = 1.0;  // Simplified step size\n"
"            double dtau = alpha_nu * ds;\n"
"            \n"
"            // Update optical depth\n"
"            tau += dtau;\n"
"            \n"
"            // Update Stokes parameters (simplified)\n"
"            double intensity = j_nu / alpha_nu * (1.0 - exp(-dtau));\n"
"            stokes.I += intensity * exp(-tau);\n"
"            \n"
"            // Update cumulative optical depth for Faraday rotation\n"
"            tauF += dtau * pow(nu, -2.0);\n"
"        }\n"
"        \n"
"        // Store results for this frequency\n"
"        int result_offset = (pixel_id * num_frequencies + f) * 4;\n"
"        IQUV_results[result_offset + 0] = stokes.I;\n"
"        IQUV_results[result_offset + 1] = stokes.Q;\n"
"        IQUV_results[result_offset + 2] = stokes.U;\n"
"        IQUV_results[result_offset + 3] = stokes.V;\n"
"        \n"
"        tau_results[pixel_id * num_frequencies + f] = tau;\n"
"        tauF_results[pixel_id * num_frequencies + f] = tauF;\n"
"    }\n"
"}\n";

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