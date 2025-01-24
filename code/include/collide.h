#ifndef COLLIDE_H
#define COLLIDE_H


/**
 * @brief Perform the collision step on the hydrodynamic populations f.
 * 
 * @note This function also computes the source term from the forces on the fluid.
 */
void collide_hydrodynamic_populations();

#endif