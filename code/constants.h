#ifndef CONSTANTS_H
#define CONSTANTS_H


// --- SETTINGS ---
// Solve the Navier-Stokes equations, making the fluid move under forces and pressure differences. 
#define FLOW
// Use the Shan-Chen method to simulate two (partially) immiscible fluids.
#undef DUALCOMPONENT
// Solve the Advection-Diffusion Equation for the temperature field. 
#undef TEMPERATURE
// Allow phase change using an Enthalpy method.
#undef PHASECHANGE


// --- BOUNDARY CONDITIONS ---
#define XPERIODIC
#undef WEST_NOSLIP
#undef EAST_NOSLIP
#undef WEST_VELOCITY_PRESCRIBED
#undef WEST_PRESSURE_PRESCRIBED
#undef EAST_VELOCITY_PRESCRIBED
#undef EAST_PRESSURE_PRESCRIBED

#undef YPERIODIC
#define SOUTH_NOSLIP
#define NORTH_NOSLIP
#undef SOUTH_VELOCITY_PRESCRIBED
#undef SOUTH_PRESSURE_PRESCRIBED
#undef NORTH_VELOCITY_PRESCRIBED
#undef NORTH_PRESSURE_PRESCRIBED


// --- INITIAL CONDITIONS ---
#define DENSITY_UNITY
#undef DENSITY_CONSTANT_HORIZONTAL_GRADIENT
#define VELOCITY_ZERO

#ifdef DENSITY_CONSTANT_HORIZONTAL_GRADIENT
    static double rho_ini_left = 1.0;
    static double rho_ini_right = 0.99997;
#endif


// --- PARAMETERS ---
// Number of timesteps
#define NTIME   10
// Store macroscopic quantities after NSTORE timesteps        
#define NSTORE  1
// Print progress percentage after NLOG timesteps
#define NLOG    1
// Output macroscopic quantities to h5 files
#define OUTPUT

// Number of cells in the x-direction
#define NX      32
// Number of cells in the y-direction
#define NY      64

#ifdef FLOW
    #ifdef DUALCOMPONENT
        static double tau_lava = 1.0;
        static double tau_air = 0.625;
    #else
        static double tau = 1.0;
    #endif
    static double Fx_body = 1e-5;
    static double Fy_body = 0.0;
#endif

#ifdef DUALCOMPONENT
    static double G = 4.0;
#endif

#ifdef TEMPERATURE
    #ifdef PHASECHANGE
        static double tau_g_liquid = 0.53;
        static double tau_g_solid = 0.53;
        static double c_liquid = 0.95;
        static double c_solid = 0.95;
        static double L_f = 1.0;
    #else
        static double tau_g = 0.53;
    #endif
#endif


// --- D2Q9 STENCIL ---
#define NP 9
#define INDEX_2D(i,j)       NY*(i) + (j)
#define INDEX_3D(i,j,p)     NY*NP*(i) + NP*(j) + (p)
static double cs2 = 1.0/3.0;
static int cx_i[NP] = {0, 1, -1, 0, 0, 1, -1, 1, -1};
static int cy_i[NP] = {0, 0, 0, 1, -1, 1, 1, -1, -1};
static double cx[NP] = {0.0, 1.0, -1.0, 0.0, 0.0, 1.0, -1.0, 1.0, -1.0};
static double cy[NP] = {0.0, 0.0, 0.0, 1.0, -1.0, 1.0, 1.0, -1.0, -1.0};
static double w[NP] = {4.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/9.0, 1.0/36.0, 1.0/36.0, 1.0/36.0, 1.0/36.0};
static int p_bounceback[NP] = {0, 2, 1, 4, 3, 8, 7, 6, 5};

#endif