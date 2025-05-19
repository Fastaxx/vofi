import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import matplotlib.patches as patches

# Load the data
data = pd.read_csv("interface_data.csv")

# Define the center and radius of the circle
center_x, center_y = 0, 0
radius = a=0.5

# Create figure with 3 subplots
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(18, 6))

# Calculate cell size from data
unique_cell_x = sorted(list(set(data['cell_x'])))
if len(unique_cell_x) >= 2:
    cell_size = unique_cell_x[1] - unique_cell_x[0]
else:
    cell_size = 0.1  # Default if can't determine

# Plot 1: Interface centroids with cells explicitly shown
ax1.set_title('Interface Centroids and Cell Geometry')
theta = np.linspace(0, 2*np.pi, 100)
x = center_x + radius * np.cos(theta)
y = center_y + radius * np.sin(theta)
ax1.plot(x, y, 'k--', label='Analytical circle')

# Plot the square cells as rectangles
for _, row in data.iterrows():
    rect = patches.Rectangle((row['cell_x'], row['cell_y']), 
                            cell_size, cell_size, 
                            linewidth=1, edgecolor='gray', 
                            facecolor='lightblue', alpha=0.2)
    ax1.add_patch(rect)

# Plot the cell centers
ax1.scatter(data['cell_x'] + cell_size/2, data['cell_y'] + cell_size/2, 
           color='blue', s=10, marker='+', label='Cell centers')

# Plot the interface centroids
ax1.scatter(data['centroid_x'], data['centroid_y'], 
           s=data['interface_length']*30, c='red', alpha=0.7, 
           label='Interface centroids')

ax1.set_aspect('equal')
ax1.grid(True, alpha=0.3)
ax1.legend()

# Plot 2: Radial distance histogram
ax2.set_title('Radial Distance Distribution')
ax2.hist(data['radial_distance'], bins=20, alpha=0.7)
ax2.axvline(radius, color='r', linestyle='--', label=f'Expected radius ({radius})')
ax2.set_xlabel('Distance from center')
ax2.set_ylabel('Frequency')
ax2.grid(True, alpha=0.3)
ax2.legend()

# Plot 3: Angular distribution
ax3.set_title('Angular Distribution of Interface')
angles = data['angle'].values
length = data['interface_length'].values
bins = np.linspace(0, 360, 37)  # 36 bins of 10 degrees each
hist, _ = np.histogram(angles, bins=bins, weights=length)
center = (bins[:-1] + bins[1:]) / 2
ax3.bar(center, hist, width=10, alpha=0.7)
ax3.set_xlabel('Angle (degrees)')
ax3.set_ylabel('Total interface length')
ax3.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('interface_analysis.png')
plt.show()

# Create a focused view of cells and centroids
plt.figure(figsize=(10, 10))
plt.title('Detailed View of Cells and Interface Centroids')

# Plot the analytical circle
plt.plot(x, y, 'k--', label='Analytical circle')

# Plot all cells with interface
for _, row in data.iterrows():
    rect = patches.Rectangle((row['cell_x'], row['cell_y']), 
                            cell_size, cell_size, 
                            linewidth=1, edgecolor='black', 
                            facecolor='lightblue', alpha=0.2)
    plt.gca().add_patch(rect)

# Plot cell centers
plt.scatter(data['cell_x'] + cell_size/2, data['cell_y'] + cell_size/2, 
           color='blue', s=15, marker='+', label='Cell centers')

# Plot interface centroids with size proportional to interface length
plt.scatter(data['centroid_x'], data['centroid_y'], 
           s=data['interface_length']*50, c='red', alpha=0.7, 
           label='Interface centroids')

# Connect interface centroids with lines to show the approximated interface
sorted_data = data.sort_values(by='angle')
plt.plot(sorted_data['centroid_x'], sorted_data['centroid_y'], 
         'r-', alpha=0.3, label='Approximated interface')

plt.grid(True, alpha=0.3)
plt.axis('equal')
plt.legend()
plt.savefig('detailed_interface_view.png')
plt.show()