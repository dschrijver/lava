#ifndef MACROSCOPIC_H
#define MACROSCOPIC_H


/**
 * @brief Compute the density and velocity fields of the fluid.
 */
void calculate_hydrodynamic_macroscopic_quantities();

/**
 * @brief Compute the temperature of the fluid.
 */
void calculate_thermal_macroscopic_quantities();

/**
 * @brief Compute liquid fraction using the enthalpy.
 */
void calculate_enthalpy_and_liquid_fraction();

#endif