#include <stdio.h>  // printf(), fflush(), stdout()

#include "constants.h"
#include "include/memory.h"
#include "include/initialize.h"
#include "include/output.h"
#include "include/collide.h"
#include "include/stream.h"
#include "include/macroscopic.h"
#include "include/main.h"


int main(void) {

    allocate_memory();

    initialize();

#ifdef OUTPUT
    output_data(0);
#endif

    shift_velocity();

    // --- MAIN LOOP ---
    for (int t = 0; t < NTIME; t++) {
        
#ifdef OUTPUT
        if (((t % NSTORE) == 0) && (t > 0)) {
            output_data(t);
        }
#endif

#ifdef FLOW
        collide_hydrodynamic_populations();
#endif

#ifdef FLOW
        stream_hydrodynamic_populations();
#endif

#ifdef FLOW
        calculate_hydrodynamic_macroscopic_quantities();
#endif

        // --- LOG PROGRESS ---
        if (((t+1) % NLOG) == 0) {
            printf("\rProgress: %.1f%%", ((double)(t+1))/((double)NTIME)*100.0);
            fflush(stdout);
        }
    }

    output_data(NTIME);

    free_memory();

    printf("\nDone!\n");

    return 0;
}