#include <raylib.h>
#include <stdio.h>  // sprintf()

#include "../constants.h"
#include "../include/globals.h"
#include "../include/tools.h"
#include "../include/animate.h"


static double time_last_frame;


void initialize_animation() {

    int screen_width = NX*CELL_SIZE;
    int screen_height = NY*CELL_SIZE;
    InitWindow(screen_width, screen_height, ANIMATION_TITLE);

    time_last_frame = GetTime();

}


void close_animation() {
    
    CloseWindow();

}


int render_frame() {

    double *macroscopic_quantity;
    double fraction;
    char value;
    Color cell_color = {0,0,0,255};
    double time_this_frame;
    char frame_string[32];

#if defined ANIMATE_T
    macroscopic_quantity = T;
#elif defined ANIMATE_U
    macroscopic_quantity = u;
#elif defined ANIMATE_PHI
    macroscopic_quantity = phi;
#endif

    double min_value = calculate_minimum_macroscopic_quantity(macroscopic_quantity);
    double max_value = calculate_maximum_macroscopic_quantity(macroscopic_quantity);

    PollInputEvents(); // Poll input events (SUPPORT_CUSTOM_FRAME_CONTROL)
    if (WindowShouldClose()) {
        return 0;
    }

    BeginDrawing();

        for (int i = 0; i < NX; i++) {
            for (int j = 0; j < NY; j++) {

                fraction = (macroscopic_quantity[INDEX_2D(i,j)]-min_value)/(max_value-min_value);
                value = (char) (fraction*255.0);
                cell_color.r = value;
                cell_color.b = 255 - value;
                DrawRectangle(i*CELL_SIZE, (NY-1-j)*CELL_SIZE, CELL_SIZE, CELL_SIZE, cell_color);

            }
        }

        time_this_frame = GetTime();
        sprintf(frame_string, "FPS = %.2f", 1.0/(time_this_frame-time_last_frame));
        time_last_frame = time_this_frame;
        DrawText(frame_string, 20, 20, 20, BLACK);

    EndDrawing();
    SwapScreenBuffer();

    return 1;

}

