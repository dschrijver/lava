#ifndef STREAM_H
#define STREAM_H


#include "../constants.h"

#ifdef FLOW

/**
 * @brief Perform the streaming step on the hydrodynamic populations.
 * 
 * @note This function applies the boundary conditions during the streaming step.
 */
void stream_hydrodynamic_populations();

#endif

#ifdef TEMPERATURE

/**
 * @brief Perform the streaming step on the thermal populations.
 * 
 * @note This function applies the boundary conditions during the streaming step.
 */
void stream_thermal_populations();

#endif

#endif