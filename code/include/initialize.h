#ifndef INITIALIZE_H
#define INITIALIZE_H


#include "../constants.h"

/**
 * @brief Initialize macroscopic quantities and populations.
 */
void initialize();

#ifdef FLOW

/**
 * @brief Shift initial velocity by local force.
 */
void shift_velocity();

#endif

#endif