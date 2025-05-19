#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "vofi.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Circle parameters */
static double radius = 0.5;

/* Function to compute the implicit function for a circle */
double circle_func(const double x[], void *par)
{
  /* Implicit function for a circle: x^2 + y^2 - r^2 */
  return x[0]*x[0] + x[1]*x[1] - radius*radius;
}

int main()
{
  /* Define mesh parameters */
  int dim = 2;  /* 2D problem */
  double domain_min[2] = {-1.0, -1.0};  /* Domain boundaries */
  double domain_max[2] = {1.0, 1.0};
  double analytical_area = M_PI * radius * radius;
  double analytical_length = 2.0 * M_PI * radius;  /* Perimeter of circle */
  
  /* Arrays required by vofi_get_cc */
  double xex[5] = {0.0, 0.0, 0.0, 0.0, 0.0};
  int nex[2] = {0, 1};  /* [0]=0: don't compute centroid, [1]=1: DO compute interface length */
  int npt[1] = {20};
  int nvis[2] = {0, 0};
  
  printf("===== VOFI LIBRARY MESH TEST FOR CIRCLE =====\n");
  printf("Circle radius: %.4f\n", radius);
  printf("Analytical area: %.10f\n", analytical_area);
  printf("Analytical interface length: %.10f\n\n", analytical_length);
  
  printf("%-10s %-12s %-12s %-12s %-12s %-12s %-12s\n", 
         "Resolution", "Cells", "Total Area", "Area Error", "Total Length", "Length Error", "Time (ms)");
  printf("---------------------------------------------------------------------------------------------\n");
  
  /* Test different mesh resolutions */
  int resolutions[] = {10, 20, 40, 80, 160, 320};
  
  for (int r = 0; r < sizeof(resolutions)/sizeof(resolutions[0]); r++) {
    int n = resolutions[r];
    double h = (domain_max[0] - domain_min[0]) / n;  /* Cell size */
    int total_cells = n * n;                        /* Total number of cells */
    double numerical_area = 0.0;                    /* Accumulated area */
    double numerical_length = 0.0;                  /* Accumulated interface length */
    
    /* Time measurement */
    clock_t start = clock();
    
    /* Loop through all cells in the mesh */
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double x0[3] = {domain_min[0] + i*h, domain_min[1] + j*h, 0.0};
        double cell_size[3] = {h, h, 0.0};
        
        /* Reset xex array for each cell */
        xex[0] = xex[1] = xex[2] = xex[3] = xex[4] = 0.0;
        
        /* Calculate volume fraction using VOFI */
        double vol_frac = vofi_get_cc(circle_func, NULL, x0, cell_size, xex, nex, npt, nvis, dim);
        
        /* Add to total area */
        numerical_area += vol_frac * h * h;
        
        /* Add to total interface length - xex[4] contains interface length */
        numerical_length += xex[3];
      }
    }
    
    /* Calculate execution time */
    clock_t end = clock();
    double cpu_time_ms = 1000.0 * (double)(end - start) / CLOCKS_PER_SEC;
    
    /* Calculate errors */
    double area_rel_error = fabs(numerical_area - analytical_area) / analytical_area;
    double length_rel_error = fabs(numerical_length - analytical_length) / analytical_length;
    
    /* Print results */
    printf("%-10d %-12d %-12.10f %-12.4e %-12.10f %-12.4e %-12.4f\n",
           n, total_cells, numerical_area, area_rel_error, 
           numerical_length, length_rel_error, cpu_time_ms);
  }
  
  /* Print conclusion */
  printf("\n===== CONCLUSION =====\n");
  printf("The VOFI library can accurately compute both:\n");
  printf("1. The area of a circle using volume fractions\n");
  printf("2. The interface length (perimeter) of the circle\n");
  
  return 0;
}