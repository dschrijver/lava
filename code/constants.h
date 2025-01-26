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
#define FLOW
// Use the Shan-Chen method to simulate two (partially) immiscible fluids.
#undef DUALCOMPONENT
// Solve the Advection-Diffusion Equation for the temperature field. 
#define TEMPERATURE
// Allow phase change using an Enthalpy method.
#define PHASECHANGE


// --- BOUNDARY CONDITIONS ---
#undef XPERIODIC
#undef WEST_NOSLIP
#undef EAST_NOSLIP
#define WEST_PRESSURE_WETNODE_NEBB
#define EAST_PRESSURE_WETNODE_NEBB
#define WEST_TEMPERATURE_WETNODE_NEBB_CONSTANT_VALUE
#define EAST_TEMPERATURE_WETNODE_NEBB_CONSTANT_GRADIENT

#ifdef WEST_PRESSURE_WETNODE_NEBB
    static double rho_left = 1.0;
#endif

#ifdef EAST_PRESSURE_WETNODE_NEBB
    static double rho_right = 1.0 - 200e-5*3.0;
#endif

#ifdef WEST_TEMPERATURE_WETNODE_NEBB_CONSTANT_VALUE
    static double T_left = 0.0526;
#endif

#undef YPERIODIC
#define SOUTH_NOSLIP_HALFWAY_BOUNCEBACK
#define NORTH_NOSLIP_HALFWAY_BOUNCEBACK

#ifdef SOUTH_NOSLIP_HALFWAY_BOUNCEBACK
    static double T_bottom = -1.0;
#endif

#ifdef NORTH_NOSLIP_HALFWAY_BOUNCEBACK
    static double T_top = -1.0;
#endif


// --- INITIAL CONDITIONS ---
#undef DENSITY_UNITY
#define DENSITY_CONSTANT_HORIZONTAL_GRADIENT
#define VELOCITY_ZERO
#undef VELOCITY_SINE
#undef TEMPERATURE_CONSTANT_VALUE
#undef TEMPERATURE_CONSTANT_VERTICAL_GRADIENT
#define TEMPERATURE_CHANNEL
#undef PHI_CONSTANT_VALUE
#define PHI_CHANNEL

#ifdef DENSITY_CONSTANT_HORIZONTAL_GRADIENT
    static double rho_ini_left = 1.0;
    static double rho_ini_right = 1.0 - 200e-5*3.0;
#endif

#ifdef VELOCITY_SINE
    static double u_ini_amplitude = 0.0001;
    static double v_ini_amplitude = 0.0001;
    static double u_ini_frequency = 1.0;
    static double v_ini_frequency = 1.0;
#endif

#ifdef TEMPERATURE_CONSTANT_VALUE
    static double T_ini = 0.0526;
#endif

#ifdef TEMPERATURE_CHANNEL
    static double T_solid = -1.0;
    static double T_liquid = 0.0526;
#endif

#ifdef PHI_CONSTANT_VALUE
    static double phi_ini = 1.0;
#endif

#ifdef PHI_CHANNEL
    static double channel_width = 100.0;
#endif


// --- PARAMETERS ---
// Number of timesteps
#define NTIME   100000
// Store macroscopic quantities after NSTORE timesteps        
#define NSTORE  1000
// Print progress percentage after NLOG timesteps
#define NLOG    1
// Output macroscopic quantities to h5 files
#undef OUTPUT
// Animate a macroscopic quantity
#define ANIMATE

// Animate the x-component of the velocity
#undef ANIMATE_U
// Animate the temperature
#undef ANIMATE_T
// Animate the liquid fraction
#define ANIMATE_PHI
// Number of pixels per cell, must be at least 2
#define CELL_SIZE 3
// Render a frame of the macroscopic quantity after NFRAME steps
#define NFRAME  6
// Title of animation window
#define ANIMATION_TITLE "Channel"

// Number of cells in the x-direction
#define NX 200
// Number of cells in the y-direction
#define NY 300

// FLOW
static double tau = 1.0;
static double Fx_body = 0.0;
static double Fy_body = 0.0;

// DUALCOMPONENT
static double tau_lava = 1.0;
static double tau_air = 0.625;
static double G = 4.0;

// TEMPERATURE
static double tau_g = 1.0;
static double alpha = 0.0;
static double g = 0.009661835748792274;

// PHASECHANGE
static double tau_g_liquid = 1.0;
static double tau_g_solid = 1.0;
static double c_liquid = 0.95;
static double c_solid = 0.95;
static double L_f = 0.1;
static double T_m = 0.0;

#endif