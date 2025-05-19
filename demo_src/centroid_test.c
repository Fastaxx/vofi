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
static double center_x = 0.0;  // Can be changed to test off-center circles
static double center_y = 0.0;

/* Function to compute the implicit function for a circle */
double circle_func(const double x[], void *par)
{
  /* Implicit function for a circle: (x-x0)^2 + (y-y0)^2 - r^2 */
  return (x[0]-center_x)*(x[0]-center_x) + (x[1]-center_y)*(x[1]-center_y) - radius*radius;
}

/* Function to calculate angle in degrees (0-360) from x,y coordinates */
double calculate_angle(double x, double y) {
  double angle = atan2(y - center_y, x - center_x) * 180.0 / M_PI;
  if (angle < 0) angle += 360.0;
  return angle;
}

/* Function to calculate radial distance from center */
double calculate_radius(double x, double y) {
  return sqrt((x - center_x)*(x - center_x) + (y - center_y)*(y - center_y));
}

int main()
{
  /* Define mesh parameters */
  int dim = 2;  /* 2D problem */
  double domain_min[2] = {-1.0, -1.0};  /* Domain boundaries */
  double domain_max[2] = {1.0, 1.0};
  double analytical_length = 2.0 * M_PI * radius;  /* Perimeter of circle */
  double analytical_centroid[2] = {center_x, center_y};  /* Expected centroid */
  
  /* Arrays required by vofi_get_cc - now extended for centroid */
  double xex[7] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // Now needs 7 elements
  int nex[3] = {0, 1, 1};  /* [0]=0: no volume centroid, [1]=1: DO compute interface length, [2]=1: DO compute interface centroid */
  int npt[1] = {20};
  int nvis[2] = {0, 0};
  
  printf("===== VOFI LIBRARY TEST FOR INTERFACE CENTROID =====\n");
  printf("Circle radius: %.4f\n", radius);
  printf("Circle center: (%.4f, %.4f)\n", center_x, center_y);
  printf("Analytical centroid: (%.10f, %.10f)\n\n", analytical_centroid[0], analytical_centroid[1]);
  
  printf("%-10s %-12s %-15s %-15s %-15s %-15s\n", 
         "Resolution", "Cells", "Centroid X", "Centroid Y", "Error X", "Error Y");
  printf("----------------------------------------------------------------------------------\n");
  
  /* Test different mesh resolutions */
  int resolutions[] = {10, 20, 40, 80, 160};
  int selected_resolution = 20;  // Resolution to use for detailed analysis
  FILE *fp = NULL;
  
  for (int r = 0; r < sizeof(resolutions)/sizeof(resolutions[0]); r++) {
    int n = resolutions[r];
    double h = (domain_max[0] - domain_min[0]) / n;  /* Cell size */
    int total_cells = n * n;                        /* Total number of cells */
    double weighted_centroid[2] = {0.0, 0.0};       /* For accumulating length-weighted centroid */
    double total_length = 0.0;                      /* Total interface length */
    
    /* Statistical variables for analysis */
    double min_radius = 1e10, max_radius = 0.0, sum_radius = 0.0, sum_radius_squared = 0.0;
    int interface_cell_count = 0;
    
    /* Create array for angular distribution analysis (36 bins of 10 degrees) */
    double angle_bins[36] = {0};
    double length_bins[36] = {0};
    
    /* Open CSV file for output if this is the selected resolution */
    if (n == selected_resolution) {
      fp = fopen("interface_data.csv", "w");
      fprintf(fp, "cell_x,cell_y,volume_fraction,interface_length,centroid_x,centroid_y,radial_distance,angle,expected_radius\n");
    }
    
    /* Time measurement */
    clock_t start = clock();
    
    /* Loop through all cells in the mesh */
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < n; j++) {
        double x0[3] = {domain_min[0] + i*h, domain_min[1] + j*h, 0.0};
        double cell_size[3] = {h, h, 0.0};
        
        /* Reset xex array for each cell */
        for (int k = 0; k < 7; k++) xex[k] = 0.0;
        
        /* Calculate volume fraction using VOFI with interface centroid */
        double vol_frac = vofi_get_cc(circle_func, NULL, x0, cell_size, xex, nex, npt, nvis, dim);
        
        /* If this cell has an interface segment, accumulate weighted centroid */
        if (xex[3] > 0.0) {  /* If interface length > 0 */
            /* Add weighted contribution of this cell's interface centroid */
            weighted_centroid[0] += xex[5] * xex[3];  /* x-coord * segment length */
            weighted_centroid[1] += xex[6] * xex[3];  /* y-coord * segment length */
            total_length += xex[3];
            
            /* Perform radial analysis */
            double radial_dist = calculate_radius(xex[5], xex[6]);
            min_radius = (radial_dist < min_radius) ? radial_dist : min_radius;
            max_radius = (radial_dist > max_radius) ? radial_dist : max_radius;
            sum_radius += radial_dist;
            sum_radius_squared += radial_dist * radial_dist;
            
            /* Perform angular analysis */
            double angle = calculate_angle(xex[5], xex[6]);
            int bin_index = (int)(angle / 10.0);
            if (bin_index >= 36) bin_index = 0;  /* Handle edge case */
            angle_bins[bin_index] += 1;
            length_bins[bin_index] += xex[3];
            
            interface_cell_count++;
            
            /* Write data to CSV file if this is the selected resolution */
            if (n == selected_resolution && fp != NULL) {
                fprintf(fp, "%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.2f,%.6f\n",
                        x0[0], x0[1], vol_frac, xex[3], xex[5], xex[6], 
                        radial_dist, angle, radius);
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
    double calculated_centroid[2] = {0.0, 0.0};
    if (total_length > 0.0) {
        calculated_centroid[0] = weighted_centroid[0] / total_length;
        calculated_centroid[1] = weighted_centroid[1] / total_length;
    }
    
    /* Calculate radial statistics */
    double avg_radius = sum_radius / interface_cell_count;
    double variance = (sum_radius_squared / interface_cell_count) - (avg_radius * avg_radius);
    double std_dev = sqrt(variance);
    
    /* Calculate errors */
    double centroid_error_x = fabs(calculated_centroid[0] - analytical_centroid[0]);
    double centroid_error_y = fabs(calculated_centroid[1] - analytical_centroid[1]);
    
    /* Print results */
    printf("%-10d %-12d %-15.10f %-15.10f %-15.10e %-15.10e\n",
           n, total_cells, calculated_centroid[0], calculated_centroid[1], 
           centroid_error_x, centroid_error_y);
           
    /* Print detailed analysis for selected resolution */
    if (n == selected_resolution) {
        /* Open analysis file */
        fp = fopen("interface_analysis.txt", "w");
        if (fp == NULL) {
            printf("Error opening interface_analysis.txt file!\n");
        } else {
            /* Radial distance analysis */
            fprintf(fp, "===== RADIAL DISTANCE ANALYSIS =====\n");
            fprintf(fp, "Expected radius: %.6f\n", radius);
            fprintf(fp, "Min radius: %.6f\n", min_radius);
            fprintf(fp, "Max radius: %.6f\n", max_radius);
            fprintf(fp, "Average radius: %.6f\n", avg_radius);
            fprintf(fp, "Standard deviation: %.6f\n", std_dev);
            fprintf(fp, "Relative error: %.6f%%\n\n", 100.0 * fabs(avg_radius - radius) / radius);
            
            /* Angular distribution analysis */
            fprintf(fp, "===== ANGULAR DISTRIBUTION ANALYSIS =====\n");
            fprintf(fp, "Angle (degrees), Cell Count, Interface Length, Expected Length\n");
            double expected_length_per_bin = analytical_length / 36.0;
            for (int i = 0; i < 36; i++) {
                fprintf(fp, "%d-%d, %.0f, %.6f, %.6f\n", 
                        i*10, (i+1)*10, angle_bins[i], length_bins[i], expected_length_per_bin);
            }
            
            /* Length distribution analysis */
            fprintf(fp, "\n===== LENGTH DISTRIBUTION ANALYSIS =====\n");
            fprintf(fp, "Total interface length: %.6f\n", total_length);
            fprintf(fp, "Expected interface length: %.6f\n", analytical_length);
            fprintf(fp, "Relative error: %.6f%%\n", 100.0 * fabs(total_length - analytical_length) / analytical_length);
            fprintf(fp, "Average length per cell: %.6f\n", total_length / interface_cell_count);
            
            fclose(fp);
            printf("\nDetailed analysis written to interface_analysis.txt\n");
            printf("Interface data written to interface_data.csv\n");
        }
    }
  }
  
  /* Test with a single cell that has the interface */
  printf("\n===== SINGLE CELL TEST =====\n");
  double x0[3] = {center_x - radius/2, center_y - radius/2, 0.0};  /* Cell containing part of the circle */
  double cell_size[3] = {radius, radius, 0.0};
  
  /* Reset xex array */
  for (int k = 0; k < 7; k++) xex[k] = 0.0;
  
  /* Calculate with interface centroid */
  double vol_frac = vofi_get_cc(circle_func, NULL, x0, cell_size, xex, nex, npt, nvis, dim);
  
  printf("Volume fraction: %.10f\n", vol_frac);
  printf("Interface length: %.10f\n", xex[3]);
  printf("Interface centroid: (%.10f, %.10f)\n", xex[5], xex[6]);
  
  /* Additional analysis for single cell */
  double radial_dist = calculate_radius(xex[5], xex[6]);
  double angle = calculate_angle(xex[5], xex[6]);
  printf("Radial distance: %.10f (expected: %.10f)\n", radial_dist, radius);
  printf("Angular position: %.2f degrees\n", angle);
  
  /* Print conclusion */
  printf("\n===== CONCLUSION =====\n");
  printf("The VOFI library now computes:\n");
  printf("1. The volume fraction of a cell (%.10f)\n", vol_frac);
  printf("2. The interface length in a cell (%.10f)\n", xex[3]);
  printf("3. The interface centroid in a cell (%.10f, %.10f)\n", xex[5], xex[6]);
  printf("\n");
  printf("Analysis files generated:\n");
  printf("1. interface_data.csv - Raw interface data for visualization\n");
  printf("2. interface_analysis.txt - Statistical analysis of interface properties\n");
  printf("\nTo visualize the results, use the provided Python script.\n");
  
  return 0;
}