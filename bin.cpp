#include "bin.h"

void init_grid(int num_bins, bin_t *bin_list)
{
    for(int i = 0; i < num_bins; i++)
    {
        bin_list[i].capacity = 20;
        bin_list[i].bin_size = 0;
        bin_list[i].indeces = (int*)malloc(bin_list[i].capacity*sizeof(int));
    }
       
}

void clear_grid(int num_bins, bin_t *bin_list)
{
    for(int i = 0; i < num_bins; i++)
    {   
        free(bin_list[i].indeces);
        bin_list[i].indeces = NULL;
    }
        
}

void clear_bin_col(bin_t *bin_list, int start, int end, int col, int lda)
{   
    for(int r = start; r < end; r++)
    {
        int index = r + col * lda;
        free(bin_list[index].indeces);
        bin_list[index].indeces = (int*)malloc(sizeof(int)*bin_list[index].capacity);
        bin_list[index].bin_size = 0;
    }
}

// The list of bins will be separated into y*x
void set_grid_size(int &x, int &y, int num_bins)
{
    int m = static_cast<int>(sqrt(num_bins));
    for(int i = m; i > 0; i--)
    {
        if(num_bins % i == 0) 
        {
            y = i;
            break;
        }
    }

    x = num_bins / y;
}

void bin_particles(int n, particle_t *particles, int num_bins, bin_t *bin_list, double bin_x, double bin_y, int num_rows)
{
    int i;
    for(i = 0; i < num_bins; i++)
    {
        bin_list[i].bin_size = 0;
    }

    // Binning particles to bins
    for(i = 0; i < n; i++)
    {   
        int x = particles[i].x / bin_x, y = particles[i].y / bin_y;
        int index = y + x*num_rows;
        add_particle(bin_list, i, index);
    }
}

// This function adds a particle to the end of a bin
void add_particle(bin_t *bin_list, int i, int j)
{ 
    if (bin_list[j].bin_size == bin_list[j].capacity)
    {
        // Need to allocate more memory here
        bin_list[j].indeces = (int*) realloc(bin_list[j].indeces, 2*bin_list[j].capacity*sizeof(int));
        bin_list[j].capacity  *= 2;
    }

    bin_list[j].indeces[bin_list[j].bin_size] = i;
    bin_list[j].bin_size ++;
}

void remove_particle(bin_t *bin_list, int i, int j)
{
    for (int k = 0; k < bin_list[j].bin_size; k++)
    {
        if (bin_list[j].indeces[k] == i) // Need to remove this particle
        {
            for (int l = k; l < bin_list[j].bin_size; l++) 
            {
                bin_list[j].indeces[l] = bin_list[j].indeces[l+1];
            }
            break;
        }
        if (k == bin_list[j].bin_size - 1) 
        {
            printf("Failed to remove particle %d from bin %d\n", i, j);
            printf("bin %d has partiles: \n", j);
            for(int p = 0; p < bin_list[j].bin_size; p ++)
                printf("%d \n", bin_list[j].indeces[p]);
            printf("\n");

            exit(1);
        }
    }

    bin_list[j].bin_size --;
}