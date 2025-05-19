#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "vofi.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Sphere parameters */
static double radius = 0.5;
static double center_x = 0.0;
static double center_y = 0.0;
static double center_z = 0.0;

/* Function to compute the implicit function for a sphere */
double sphere_func(const double x[], void *par)
{
  /* Implicit function for a sphere: (x-x0)^2 + (y-y0)^2 + (z-z0)^2 - r^2 */
  return (x[0]-center_x)*(x[0]-center_x) + 
         (x[1]-center_y)*(x[1]-center_y) + 
         (x[2]-center_z)*(x[2]-center_z) - radius*radius;
}

/* Function to calculate distance from center */
double calculate_distance(double x, double y, double z) {
  return sqrt((x - center_x)*(x - center_x) + 
              (y - center_y)*(y - center_y) + 
              (z - center_z)*(z - center_z));
}

int main()
{
  /* Define mesh parameters */
  int dim = 3;  /* 3D problem */
  double domain_min[3] = {-1.0, -1.0, -1.0};  /* Domain boundaries */
  double domain_max[3] = {1.0, 1.0, 1.0};
  double analytical_area = 4.0 * M_PI * radius * radius;  /* Surface area of sphere */
  double analytical_centroid[3] = {center_x, center_y, center_z};  /* Expected centroid */
  
  /* Arrays required by vofi_get_cc */
  double xex[8] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Now needs 8 elements for 3D
  int nex[3] = {0, 1, 1};  /* [0]=0: no volume centroid, [1]=1: DO compute interface area, [2]=1: DO compute interface centroid */
  int npt[4] = {5, 5, 10, 10};  // More refinement for 3D
  int nvis[2] = {0, 0};
  
  printf("===== VOFI LIBRARY TEST FOR 3D INTERFACE CENTROID =====\n");
  printf("Sphere radius: %.4f\n", radius);
  printf("Sphere center: (%.4f, %.4f, %.4f)\n", center_x, center_y, center_z);
  printf("Analytical centroid: (%.10f, %.10f, %.10f)\n\n", 
         analytical_centroid[0], analytical_centroid[1], analytical_centroid[2]);
  
  printf("%-10s %-12s %-15s %-15s %-15s %-15s %-15s %-15s\n", 
         "Resolution", "Cells", "Centroid X", "Centroid Y", "Centroid Z", 
         "Error X", "Error Y", "Error Z");
  printf("-------------------------------------------------------------------------------------------\n");
  
  /* Test different mesh resolutions */
  int resolutions[] = {10, 20, 30, 40};
  int selected_resolution = 20;  // Resolution to use for detailed analysis
  FILE *fp = NULL;
  
  for (int r = 0; r < sizeof(resolutions)/sizeof(resolutions[0]); r++) {
    int n = resolutions[r];
    double h = (domain_max[0] - domain_min[0]) / n;  /* Cell size */
    int total_cells = n * n * n;                     /* Total number of cells */
    double weighted_centroid[3] = {0.0, 0.0, 0.0};   /* For accumulating area-weighted centroid */
    double total_area = 0.0;                         /* Total interface area */
    
    /* Statistical variables for analysis */
    double min_radius = 1e10, max_radius = 0.0, sum_radius = 0.0;
    int interface_cell_count = 0;
    
    /* Open CSV file for output if this is the selected resolution */
    if (n == selected_resolution) {
      fp = fopen("interface_data_3d.csv", "w");
      fprintf(fp, "cell_x,cell_y,cell_z,volume_fraction,interface_area,centroid_x,centroid_y,centroid_z,distance,expected_radius\n");
    }
    
    /* Time measurement */
    clock_t start = clock();
    
    /* Loop through all cells in the mesh */
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < n; k++) {
          double x0[3] = {domain_min[0] + i*h, 
                          domain_min[1] + j*h, 
                          domain_min[2] + k*h};
          double cell_size[3] = {h, h, h};
          
          /* Reset xex array for each cell */
          for (int m = 0; m < 8; m++) xex[m] = 0.0;
          
          /* Calculate volume fraction using VOFI with interface centroid */
          double vol_frac = vofi_get_cc(sphere_func, NULL, x0, cell_size, xex, nex, npt, nvis, dim);
          
          /* If this cell has an interface segment, accumulate weighted centroid */
          if (xex[3] > 0.0) {  /* If interface area > 0 */
              /* Add weighted contribution of this cell's interface centroid */
              weighted_centroid[0] += xex[5] * xex[3];  /* x-coord * segment area */
              weighted_centroid[1] += xex[6] * xex[3];  /* y-coord * segment area */
              weighted_centroid[2] += xex[7] * xex[3];  /* z-coord * segment area */
              total_area += xex[3];
              
              /* Perform radial analysis */
              double distance = calculate_distance(xex[5], xex[6], xex[7]);
              min_radius = (distance < min_radius) ? distance : min_radius;
              max_radius = (distance > max_radius) ? distance : max_radius;
              sum_radius += distance;
              
              interface_cell_count++;
              
              /* Write data to CSV file if this is the selected resolution */
              if (n == selected_resolution && fp != NULL) {
                  fprintf(fp, "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f\n",
                          x0[0], x0[1], x0[2], vol_frac, xex[3], 
                          xex[5], xex[6], xex[7], distance, radius);
              }
          }
        }
      }
    }
    
    /* Close CSV file if open */
    if (n == selected_resolution && fp != NULL) {
        fclose(fp);
        fp = NULL;
    }
    
    /* Calculate final centroid */
    double calculated_centroid[3] = {0.0, 0.0, 0.0};
    if (total_area > 0.0) {
        calculated_centroid[0] = weighted_centroid[0] / total_area;
        calculated_centroid[1] = weighted_centroid[1] / total_area;
        calculated_centroid[2] = weighted_centroid[2] / total_area;
    }
    
    /* Calculate radial statistics */
    double avg_radius = (interface_cell_count > 0) ? sum_radius / interface_cell_count : 0.0;
    
    /* Calculate errors */
    double centroid_error_x = fabs(calculated_centroid[0] - analytical_centroid[0]);
    double centroid_error_y = fabs(calculated_centroid[1] - analytical_centroid[1]);
    double centroid_error_z = fabs(calculated_centroid[2] - analytical_centroid[2]);
    
    /* Timing */
    clock_t end = clock();
    double cpu_time = ((double) (end - start)) / CLOCKS_PER_SEC;
    
    /* Print results */
    printf("%-10d %-12d %-15.10f %-15.10f %-15.10f %-15.10e %-15.10e %-15.10e\n",
           n, total_cells, 
           calculated_centroid[0], calculated_centroid[1], calculated_centroid[2],
           centroid_error_x, centroid_error_y, centroid_error_z);
           
    /* Print detailed analysis for selected resolution */
    if (n == selected_resolution) {
        fp = fopen("interface_analysis_3d.txt", "w");
        if (fp != NULL) {
            /* Interface statistics */
            fprintf(fp, "===== SPHERE INTERFACE ANALYSIS =====\n");
            fprintf(fp, "Resolution: %d x %d x %d\n", n, n, n);
            fprintf(fp, "Total cells: %d\n", total_cells);
            fprintf(fp, "Interface cells: %d (%.2f%%)\n", interface_cell_count, 
                    100.0 * interface_cell_count / total_cells);
            fprintf(fp, "Computational time: %.2f seconds\n\n", cpu_time);
            
            /* Centroid information */
            fprintf(fp, "===== CENTROID ANALYSIS =====\n");
            fprintf(fp, "Calculated centroid: (%.10f, %.10f, %.10f)\n", 
                    calculated_centroid[0], calculated_centroid[1], calculated_centroid[2]);
            fprintf(fp, "Expected centroid: (%.10f, %.10f, %.10f)\n", 
                    analytical_centroid[0], analytical_centroid[1], analytical_centroid[2]);
            fprintf(fp, "Absolute errors: (%.10e, %.10e, %.10e)\n\n", 
                    centroid_error_x, centroid_error_y, centroid_error_z);
            
            /* Radius information */
            fprintf(fp, "===== RADIUS ANALYSIS =====\n");
            fprintf(fp, "Expected radius: %.6f\n", radius);
            fprintf(fp, "Min radius: %.6f (error: %.2f%%)\n", 
                    min_radius, 100.0 * fabs(min_radius - radius) / radius);
            fprintf(fp, "Max radius: %.6f (error: %.2f%%)\n", 
                    max_radius, 100.0 * fabs(max_radius - radius) / radius);
            fprintf(fp, "Average radius: %.6f (error: %.2f%%)\n\n", 
                    avg_radius, 100.0 * fabs(avg_radius - radius) / radius);
            
            /* Surface area information */
            fprintf(fp, "===== AREA ANALYSIS =====\n");
            fprintf(fp, "Calculated surface area: %.6f\n", total_area);
            fprintf(fp, "Expected surface area: %.6f\n", analytical_area);
            fprintf(fp, "Relative error: %.2f%%\n", 
                    100.0 * fabs(total_area - analytical_area) / analytical_area);
            
            fclose(fp);
            printf("\nDetailed analysis written to interface_analysis_3d.txt\n");
            printf("Interface data written to interface_data_3d.csv\n");
        }
    }
  }
  
  /* Test with off-center sphere */
  printf("\n===== OFF-CENTER SPHERE TEST =====\n");
  center_x = 0.2; center_y = -0.1; center_z = 0.3;
  double x0[3] = {0.0, 0.0, 0.0};
  double cell_size[3] = {2.0, 2.0, 2.0};  // Large cell containing sphere
  
  /* Reset xex array */
  for (int k = 0; k < 8; k++) xex[k] = 0.0;
  
  /* Calculate with interface centroid */
  double vol_frac = vofi_get_cc(sphere_func, NULL, x0, cell_size, xex, nex, npt, nvis, dim);
  
  printf("Volume fraction: %.10f\n", vol_frac);
  printf("Interface area: %.10f\n", xex[3]);
  printf("Interface centroid: (%.10f, %.10f, %.10f)\n", xex[5], xex[6], xex[7]);
  printf("Expected centroid: (%.10f, %.10f, %.10f)\n", center_x, center_y, center_z);
  double dist_error = sqrt(pow(xex[5]-center_x, 2) + pow(xex[6]-center_y, 2) + pow(xex[7]-center_z, 2));
  printf("Distance error: %.10e\n", dist_error);
  
  /* Print conclusion */
  printf("\n===== CONCLUSION =====\n");
  printf("The VOFI library now computes:\n");
  printf("1. The volume fraction of a cell (%.10f)\n", vol_frac);
  printf("2. The interface area in a cell (%.10f)\n", xex[3]);
  printf("3. The interface centroid in a cell (%.10f, %.10f, %.10f)\n", xex[5], xex[6], xex[7]);
  printf("\n");
  printf("Analysis files generated:\n");
  printf("1. interface_data_3d.csv - Raw interface data for visualization\n");
  printf("2. interface_analysis_3d.txt - Statistical analysis of interface properties\n");
  
  return 0;
}