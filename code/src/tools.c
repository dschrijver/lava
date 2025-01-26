#include "../constants.h"
#include "../include/tools.h"


double calculate_minimum_macroscopic_quantity(double *macroscopic_quantity) {

    double value;

    double result = macroscopic_quantity[INDEX_2D(0,0)];

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {

            value = macroscopic_quantity[INDEX_2D(i,j)];
            if (value < result) result = value;

        }
    }

    return result;

}


double calculate_maximum_macroscopic_quantity(double *macroscopic_quantity) {

    double value;

    double result = macroscopic_quantity[INDEX_2D(0,0)];

    for (int i = 0; i < NX; i++) {
        for (int j = 0; j < NY; j++) {
            
            value = macroscopic_quantity[INDEX_2D(i,j)];
            if (value > result) result = value;

        }
    }

    return result;

}