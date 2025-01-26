#ifndef CONSTANTS_H
#define CONSTANTS_H


// --- TOOLS ---
static double cs2 = 1.0/3.0;
#define NP 9
#define INDEX_2D(i,j)       NY*(i) + (j)
#define INDEX_3D(i,j,p)     NY*NP*(i) + NP*(j) + (p)
#define M_PI 3.14159265358979323846


// --- SETTINGS ---
// Solve the Navier-Stokes equations, making the fluid move under forces and pressure differences. 
#undef FLOW
// Use the Shan-Chen method to simulate two (partially) immiscible fluids.
#undef DUALCOMPONENT
// Solve the Advection-Diffusion Equation for the temperature field. 
#define TEMPERATURE
// Allow phase change using an Enthalpy method.
#define PHASECHANGE


// --- BOUNDARY CONDITIONS ---
#define XPERIODIC
#undef WEST_NOSLIP
#undef EAST_NOSLIP
#undef WEST_PRESSURE_WETNODE_NEBB
#undef EAST_PRESSURE_WETNODE_NEBB

#ifdef WEST_PRESSURE_WETNODE_NEBB
    static double rho_left = 1.0;
#endif

#ifdef EAST_PRESSURE_WETNODE_NEBB
    static double rho_right = 1.0 - 32e-5*3.0;
#endif

#undef YPERIODIC
#define SOUTH_NOSLIP_HALFWAY_BOUNCEBACK
#define NORTH_NOSLIP_HALFWAY_BOUNCEBACK

#ifdef SOUTH_NOSLIP_HALFWAY_BOUNCEBACK
    static double T_bottom = -1.0/0.95;
#endif

#ifdef NORTH_NOSLIP_HALFWAY_BOUNCEBACK
    static double T_top = 0.0;
#endif


// --- INITIAL CONDITIONS ---
#define DENSITY_UNITY
#undef DENSITY_CONSTANT_HORIZONTAL_GRADIENT
#define VELOCITY_ZERO
#undef VELOCITY_SINE
#define TEMPERATURE_CONSTANT_VALUE
#undef TEMPERATURE_CONSTANT_VERTICAL_GRADIENT
#define PHI_CONSTANT_VALUE

#ifdef DENSITY_CONSTANT_HORIZONTAL_GRADIENT
    static double rho_ini_left = 1.0;
    static double rho_ini_right = 1.0 - 32e-5*3.0;
#endif

#ifdef VELOCITY_SINE
    static double u_ini_amplitude = 0.0001;
    static double v_ini_amplitude = 0.0001;
    static double u_ini_frequency = 1.0;
    static double v_ini_frequency = 1.0;
#endif

#ifdef TEMPERATURE_CONSTANT_VALUE
    static double T_ini = 0.0;
#endif

#ifdef PHI_CONSTANT_VALUE
    static double phi_ini = 1.0;
#endif


// --- PARAMETERS ---
// Number of timesteps
#define NTIME   1000000
// Store macroscopic quantities after NSTORE timesteps        
#define NSTORE  10000
// Print progress percentage after NLOG timesteps
#define NLOG    100
// Output macroscopic quantities to h5 files
#define OUTPUT
// Animate a macroscopic quantity
#undef ANIMATE

// Animate the x-component of the velocity
#undef ANIMATE_U
// Animate the temperature
#define ANIMATE_T
// Number of pixels per cell, must be at least 2
#define CELL_SIZE 8
// Render a frame of the macroscopic quantity after NFRAME steps
#define NFRAME  100
// Title of animation window
#define ANIMATION_TITLE "Rayleigh-Benard"

// Number of cells in the x-direction
#define NX 4
// Number of cells in the y-direction
#define NY 2048

// FLOW
static double tau = 1.0;
static double Fx_body = 0.0;
static double Fy_body = 0.0;

// DUALCOMPONENT
static double tau_lava = 1.0;
static double tau_air = 0.625;
static double G = 4.0;

// TEMPERATURE
static double tau_g = 0.5017801724872908;
static double alpha = 207e-6;
static double g = 0.009661835748792274;

// PHASECHANGE
static double tau_g_liquid = 0.50498;
static double tau_g_solid = 0.50498;
static double c_liquid = 0.95;
static double c_solid = 0.95;
static double L_f = 1.0;
static double T_m = 0.0;

#endif