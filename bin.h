#ifndef __BIN_H__
#define __BIN_H__
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"

typedef struct
{
    int* indeces;
    int bin_size;
    int capacity;
}bin_t;

// This function will allocate initial memory for each bin
void init_grid(int num_bins, bin_t *bin_list);

// This function calculates the size of the bins, usually rectangular
void set_grid_size(int &x, int &y, int num_bins);

// This function bins particles to bins, should be called only before the simulation
void bin_particles(int n, particle_t *particles, int num_bins, bin_t *bin_list, double bin_x, double bin_y, int num_rows);

// Deallocate the memory of each bin
void clear_grid(int num_bins, bin_t *bin_list);

// Add particle i to bin j, dynamically adjust the bin size if necessary
void add_particle(bin_t *bin_list, int i, int j);

// Remove particle i from bin j
void remove_particle(bin_t *bin_list, int i, int j);

// This clears and reallocates memory for one col of the grid from start to end
void clear_bin_col(bin_t *bin_list, int start, int end, int col, int lda);
#endif