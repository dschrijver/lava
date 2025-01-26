#ifndef COLLIDE_H
#define COLLIDE_H


#include "../constants.h"

#ifdef FLOW

/**
 * @brief Perform the collision step on the hydrodynamic populations f.
 * 
 * @note This function also computes the source term from the forces on the fluid.
 */
void collide_hydrodynamic_populations();

#endif

#ifdef TEMPERATURE

/**
 * @brief Perform the collision step on the thermal populations g.
 */
void collide_thermal_populations();

#endif

#endif