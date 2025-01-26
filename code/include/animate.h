#ifndef ANIMATE_H
#define ANIMATE_H


/**
 * @brief Initialize the Raylib animation.
 */
void initialize_animation();

/**
 * @brief Close the animation window.
 */
void close_animation();

/**
 * @brief Render a frame with the current macroscopic quantity.
 */
int render_frame();

#endif