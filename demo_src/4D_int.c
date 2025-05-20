#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vofi.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Test functions for 4D implicit surfaces
double hyperplane_func(const double x[], void *par) {
    return x[3] - 0.5;  // w = 0.5 divides the hypercube
}

double hypersphere_func(const double x[], void *par) {
    double radius = 0.4;
    double center[4] = {0.5, 0.5, 0.5, 0.5}; // Center of domain
    return (x[0]-center[0])*(x[0]-center[0]) + 
           (x[1]-center[1])*(x[1]-center[1]) + 
           (x[2]-center[2])*(x[2]-center[2]) + 
           (x[3]-center[3])*(x[3]-center[3]) - radius*radius;
}

// 4D torus (with major radius R and minor radius r)
double torus4d_func(const double x[], void *par) {
    double R = 0.35; // Major radius
    double r = 0.15; // Minor radius
    double center[4] = {0.5, 0.5, 0.5, 0.5};
    
    // Distance from center to point in the xy plane
    double xy_dist = sqrt((x[0]-center[0])*(x[0]-center[0]) + 
                         (x[1]-center[1])*(x[1]-center[1]));
    
    // Distance from center to point in the zw plane
    double zw_dist = sqrt((x[2]-center[2])*(x[2]-center[2]) + 
                         (x[3]-center[3])*(x[3]-center[3]));
    
    // Distance from the circular "spine" of the torus
    double spine_dist = sqrt((xy_dist-R)*(xy_dist-R) + (zw_dist-R)*(zw_dist-R));
    
    return spine_dist - r;
}


int main() {
    // Mesh parameters
    int nx = 20;
    int ny = 20;
    int nz = 20;
    int nw = 20;
    double domain_size = 1.0;
    double dx = domain_size / nx;
    
    // Cell counts
    int total_cells = nx * ny * nz * nw;
    int empty_cells = 0;
    int full_cells = 0;
    int interface_cells = 0;
    
    printf("=== 4D Cell Type Test on %dx%dx%dx%d Mesh ===\n\n", nx, ny, nz, nw);
    printf("Domain: [0,%g] x [0,%g] x [0,%g] x [0,%g]\n", 
           domain_size, domain_size, domain_size, domain_size);
    printf("Cell size: %g x %g x %g x %g\n\n", dx, dx, dx, dx);
    
    // Select test function
    printf("Select test function:\n");
    printf("1. Hyperplane (w = 0.5)\n");
    printf("2. Hypersphere (radius = 0.4)\n"); 
    printf("3. 4D Torus\n");
    
    int choice = 2; // Default to hypersphere
    printf("Using function %d\n\n", choice);
    
    integrand func;
    switch(choice) {
        case 1: func = hyperplane_func; 
                printf("Using hyperplane function\n"); 
                break;
        case 2: func = hypersphere_func;
                printf("Using hypersphere function\n"); 
                break;
        case 3: func = torus4d_func;
                printf("Using 4D torus function\n"); 
                break;
        default: func = hypersphere_func;
                printf("Using default hypersphere function\n");
    }
    
    // Process the mesh
    printf("Processing %d cells...\n", total_cells);
    
    double x0[4], h0[4] = {dx, dx, dx, dx};
    int cell_type;
    
    // We'll take a 3D slice at w=0.5 for visualization
    int w_slice = nw/2;
    char slice_visual[nx][ny];
    
    // Initialize slice visualization
    for (int i = 0; i < nx; i++)
        for (int j = 0; j < ny; j++)
            slice_visual[i][j] = ' ';
    
    // Process all cells
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                for (int l = 0; l < nw; l++) {
                    // Set cell origin
                    x0[0] = i * dx;
                    x0[1] = j * dx;
                    x0[2] = k * dx;
                    x0[3] = l * dx;
                    
                    // Get cell type
                    cell_type = vofi_get_cell_type(func, NULL, x0, h0, 4);
                    
                    // Count cells by type
                    if (cell_type == 0) empty_cells++;
                    else if (cell_type == 1) full_cells++;
                    else interface_cells++;
                    
                    // Mark cells in the visualization slice
                    if (l == w_slice) {
                        if (cell_type == 0) slice_visual[i][j] = ' ';
                        else if (cell_type == 1) slice_visual[i][j] = '#';
                        else slice_visual[i][j] = 'I';
                    }
                }
            }
        }
    }
    
    // Print results
    printf("\nResults:\n");
    printf("Total cells: %d\n", total_cells);
    printf("Empty cells: %d (%.2f%%)\n", empty_cells, 100.0*empty_cells/total_cells);
    printf("Full cells: %d (%.2f%%)\n", full_cells, 100.0*full_cells/total_cells);
    printf("Interface cells: %d (%.2f%%)\n", interface_cells, 100.0*interface_cells/total_cells);
    
    // Print visualization of the slice at w=0.5
    printf("\n2D Slice at w=0.5 (z-middle):\n");
    printf("(' '=empty, '#'=full, 'I'=interface)\n\n");
    
    int z_slice = nz/2;  // Middle slice in z
    for (int j = ny-1; j >= 0; j--) {
        printf("%2d | ", j);
        for (int i = 0; i < nx; i++) {
            printf("%c ", slice_visual[i][j]);
        }
        printf("|\n");
    }
    printf("   ");
    for (int i = 0; i < nx; i++) {
        printf("--");
    }
    printf("\n    ");
    for (int i = 0; i < nx; i++) {
        printf("%d ", i % 10);
    }
    printf("\n");
    
        // Add this near the end of your main() function, before returning
    // filepath: /home/libat/github/vofi/demo_src/4D_int.c
    
    // Save data to file for Python visualization
    FILE *outfile = fopen("4d_cell_data.txt", "w");
    if (!outfile) {
        fprintf(stderr, "Error: Could not open output file\n");
        return 1;
    }
    
    // Write mesh parameters
    fprintf(outfile, "%d %d %d %d\n", nx, ny, nz, nw);
    fprintf(outfile, "%lf\n", domain_size);
    fprintf(outfile, "%d\n", choice); // Which function is used
    
    // Write all cell data
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            for (int k = 0; k < nz; k++) {
                for (int l = 0; l < nw; l++) {
                    // Set cell origin
                    x0[0] = i * dx;
                    x0[1] = j * dx;
                    x0[2] = k * dx;
                    x0[3] = l * dx;
                    
                    // Get cell type
                    cell_type = vofi_get_cell_type(func, NULL, x0, h0, 4);
                    
                    // Write: i j k l cell_type
                    fprintf(outfile, "%d %d %d %d %d\n", i, j, k, l, cell_type);
                }
            }
        }
    }
    
    fclose(outfile);
    printf("\nData saved to '4d_cell_data.txt' for Python visualization\n");

    return 0;
}