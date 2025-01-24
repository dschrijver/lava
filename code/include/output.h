#ifndef OUTPUT_H
#define OUTPUT_H


/**
 * @brief Output macroscopic quantities to an h5 file.
 * 
 * @param t Time at which the macroscopic quantities are defined.
 * @param n_output The name of the output file is given by data_{n_output}.h5.
 */
void output_data(int t);

#endif