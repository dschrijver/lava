#ifndef CONSTANTS_H
#define CONSTANTS_H




/*────────────────────────────────────────────────────────────────────────────*/
/*                                 PARAMETERS                                 */
/*────────────────────────────────────────────────────────────────────────────*/

/*───────────────────────────── Time Control ────────────────────────────────*/
// Number of timesteps.
#define NTIME       100000
// Print progress percentage after NLOG timesteps
#define NLOG        1

/*───────────────────────────── Domain Size ─────────────────────────────────*/
// Number of cells in the x-direction.
#define NX          200
// Number of cells in the y-direction.
#define NY          300

/*───────────────────────────── Output ──────────────────────────────────────*/
// Output macroscopic quantities to h5 files.
#undef OUTPUT
// Store macroscopic quantities in an h5 file after NSTORE timesteps.    
#define NSTORE      1000

/*───────────────────────────── Animation ───────────────────────────────────*/
// Animate a macroscopic quantity.
#define ANIMATION
// Animate the x-component of the velocity.
#undef ANIMATE_U
// Animate the temperature.
#undef ANIMATE_T
// Animate the liquid fraction.
#define ANIMATE_PHI
// Number of pixels per cell, must be at least 2.
#define CELL_SIZE   3
// Render a frame of the macroscopic quantity after NFRAME steps.
#define NFRAME      3
// Title of animation window.
#define ANIMATION_TITLE "Channel"



/*────────────────────────────────────────────────────────────────────────────*/
/*                                SOLVERS                                     */
/*────────────────────────────────────────────────────────────────────────────*/

// Solve the Navier-Stokes equations, making the fluid move under forces and pressure differences. 
#define FLOW
#ifdef FLOW
    static double Fx_body = 0.0;
    static double Fy_body = 0.0;
#endif
// Use the Shan-Chen method to simulate two (partially) immiscible fluids.
#undef DUALCOMPONENT
#ifdef DUALCOMPONENT
    static double G = 4.0;
#endif
// Solve the Advection-Diffusion Equation for the temperature field. 
#define TEMPERATURE
#ifdef TEMPERATURE
    static double alpha = 0.0;
    static double g = 0.009661835748792274;
#endif
// Allow phase change using an Enthalpy method. The constant TEMPERATURE must be defined too.
#define PHASECHANGE
#ifdef PHASECHANGE
    static double c_liquid = 0.95;
    static double c_solid = 0.95;
    static double L_f = 0.1;
    static double T_melt = 0.0;
#endif

/*───────────────────────────── Between Solvers ─────────────────────────────*/
#if defined(FLOW) && defined(DUALCOMPONENT)
    static double tau_lava = 1.0;
    static double tau_air = 0.625;
#elif defined(FLOW) && !defined(DUALCOMPONENT)
    static double tau = 1.0;
#endif

#if defined(TEMPERATURE) && defined(PHASECHANGE)
    static double tau_g_liquid = 1.0;
    static double tau_g_solid = 1.0;
#elif defined(TEMPERATURE) && !defined(PHASECHANGE)
    static double tau_g = 1.0;
#endif




/*────────────────────────────────────────────────────────────────────────────*/
/*                          BOUNDARY CONDITIONS                               */
/*────────────────────────────────────────────────────────────────────────────*/

/*─────────────────────── General Boundary Conditions ───────────────────────*/
// Fluid is periodic in the X-direction
#undef XPERIODIC
// Fluid is periodic in the Y-direction
#undef YPERIODIC
// Implement no-slip SOUTH walls using the Halfway Bounce-Back Method.
#define SOUTH_NOSLIP_HALFWAY_BOUNCEBACK
#if defined (SOUTH_NOSLIP_HALFWAY_BOUNCEBACK) && defined(TEMPERATURE)
    static double T_bottom = -1.0;
#endif
// Implement no-slip NORTH walls using the Halfway Bounce-Back Method.
#define NORTH_NOSLIP_HALFWAY_BOUNCEBACK
#if defined (NORTH_NOSLIP_HALFWAY_BOUNCEBACK) && defined(TEMPERATURE)
    static double T_top = -1.0;
#endif

/*───────────────────── Hydrodynamic Boundary Conditions ────────────────────*/
// Apply a constant pressure on the WEST boundary using the Non-equilibrium Bounce-Back Method.
#define WEST_PRESSURE_NEBB
#ifdef WEST_PRESSURE_NEBB
    static double rho_left = 1.0;
#endif
// Apply a constant pressure on the EAST boundary using the Non-equilibrium Bounce-Back Method.
#define EAST_PRESSURE_NEBB
#ifdef EAST_PRESSURE_NEBB
    static double rho_right = 1.0 - 200e-5*3.0;
#endif

/*─────────────────────── Thermal Boundary Conditions ───────────────────────*/
// Apply a constant temperature on the WEST wall using the Non-equilibrium Bounce-Back Method.
#define WEST_TEMPERATURE_NEBB_CONSTANT_VALUE 
#ifdef WEST_TEMPERATURE_NEBB_CONSTANT_VALUE
    static double T_left = 0.0526;
#endif
// Apply a constant temperature gradient on the EAST wall using the Non-equilibrium Bounce-Back Method.
#define EAST_TEMPERATURE_NEBB_CONSTANT_GRADIENT




/*────────────────────────────────────────────────────────────────────────────*/
/*                             INITIAL CONDITIONS                             */
/*────────────────────────────────────────────────────────────────────────────*/

/*───────────────────── Density Initial Conditions ──────────────────────────*/
// Initialize the density as 1.
#undef DENSITY_UNITY
// Initialize the density as a constant gradient over the x-direction.
#define DENSITY_CONSTANT_HORIZONTAL_GRADIENT
#ifdef DENSITY_CONSTANT_HORIZONTAL_GRADIENT
    static double rho_ini_left = 1.0;
    static double rho_ini_right = 1.0 - 200e-5*3.0;
#endif

/*───────────────────── Velocity Initial Conditions ─────────────────────────*/
// Initialize the velocity as zero.
#undef VELOCITY_ZERO
// Initialize both velocity components as a sine wave. 
#undef VELOCITY_SINE
#ifdef VELOCITY_SINE
    static double u_ini_amplitude = 0.0001;
    static double v_ini_amplitude = 0.0001;
    static double u_ini_frequency = 1.0;
    static double v_ini_frequency = 1.0;
#endif
#define VELOCITY_CHANNEL
#ifdef VELOCITY_CHANNEL
    static double dp_dx = 1e-5;
#endif

/*───────────────────── Temperature Initial Conditions ──────────────────────*/
// Initialize the temperature as a constant value.
#undef TEMPERATURE_CONSTANT_VALUE
#ifdef TEMPERATURE_CONSTANT_VALUE
    static double T_ini = 0.0526;
#endif
// Initialize the temperature as a constant gradient over the x-direction.
#undef TEMPERATURE_CONSTANT_VERTICAL_GRADIENT
#ifdef TEMPERATURE_CONSTANT_VERTICAL_GRADIENT
    static double T_ini_bottom = 294.0;
    static double T_ini_top = 293.0;
#endif
// Initialize the temperature to below the melting temperature near the channel walls, 
// and above the melting temperature in the middle of the channel.
#define TEMPERATURE_CHANNEL
#ifdef TEMPERATURE_CHANNEL
    static double T_ini_solid = -1.0;
    static double T_ini_liquid = 0.0526;
#endif

/*───────────────────── Liquid Fraction Initial Conditions ───────────────────*/
// Initialize the liquid fraction as a constant value. 
#undef PHI_CONSTANT_VALUE
#ifdef PHI_CONSTANT_VALUE
    static double phi_ini = 1.0;
#endif
// Initialize the liquid fraction to be 0 near the channel walls,
// and 1 in the middle of the channel.
#define PHI_CHANNEL
#ifdef PHI_CHANNEL
    static double channel_width = 100.0;
#endif




/*────────────────────────────────────────────────────────────────────────────*/
/*                                   TOOLS                                    */
/*────────────────────────────────────────────────────────────────────────────*/
static double cs2 = 1.0/3.0;
#define NP 9
#define INDEX_2D(i,j)       NY*(i) + (j)
#define INDEX_3D(i,j,p)     NY*NP*(i) + NP*(j) + (p)
#define M_PI 3.14159265358979323846

#endif