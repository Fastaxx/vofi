import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Load the data
data = pd.read_csv("interface_data_3d.csv")

# Create a 3D figure
fig = plt.figure(figsize=(15, 12))
ax = fig.add_subplot(111, projection='3d')

# Calculate cell size from data
unique_cell_x = sorted(list(set(data['cell_x'])))
if len(unique_cell_x) >= 2:
    cell_size = unique_cell_x[1] - unique_cell_x[0]
else:
    cell_size = 0.1  # Default if can't determine

# Get the sphere center and radius
center_x = np.median(data['centroid_x'])
center_y = np.median(data['centroid_y'])
center_z = np.median(data['centroid_z'])
radius = data['expected_radius'].iloc[0]

# ===== PLOT CELLS WITH INTERFACE =====
# Sample cells to avoid overcrowding the plot (adjust the number based on your data)
sample_size = min(50, len(data))
sampled_data = data.sample(sample_size)

# Create cell wireframes for sampled cells
for _, row in sampled_data.iterrows():
    # Define the 8 vertices of the cube
    x = row['cell_x']
    y = row['cell_y']
    z = row['cell_z']
    
    # Define the vertices of the cell
    vertices = [
        [x, y, z],
        [x+cell_size, y, z],
        [x+cell_size, y+cell_size, z],
        [x, y+cell_size, z],
        [x, y, z+cell_size],
        [x+cell_size, y, z+cell_size],
        [x+cell_size, y+cell_size, z+cell_size],
        [x, y+cell_size, z+cell_size]
    ]
    
    # Define the faces of the cube using indices of vertices
    faces = [
        [0, 1, 2, 3],  # bottom
        [4, 5, 6, 7],  # top
        [0, 1, 5, 4],  # front
        [2, 3, 7, 6],  # back
        [0, 3, 7, 4],  # left
        [1, 2, 6, 5]   # right
    ]
    
    # Create a Poly3DCollection for the cell
    poly3d = [[vertices[face[j]] for j in range(4)] for face in faces]
    collection = Poly3DCollection(poly3d, linewidths=1, alpha=0.1, facecolor='lightblue', edgecolor='gray')
    ax.add_collection3d(collection)

# ===== PLOT CELL CENTROIDS =====
# Plot cell centers
cell_centers = ax.scatter(
    sampled_data['cell_x'] + cell_size/2, 
    sampled_data['cell_y'] + cell_size/2, 
    sampled_data['cell_z'] + cell_size/2,
    color='blue', 
    marker='^', 
    s=30, 
    alpha=0.6,
    label='Cell centers'
)

# ===== PLOT INTERFACE CENTROIDS =====
# Plot interface centroids - now larger and more distinct
interface_centroids = ax.scatter(
    data['centroid_x'], 
    data['centroid_y'], 
    data['centroid_z'],
    c='red',
    marker='o',
    s=data['interface_area']*1000,  # Make them larger
    alpha=0.8,
    edgecolor='black',  # Add black outline for better visibility
    label='Interface centroids'
)

# ===== PLOT ANALYTICAL SPHERE =====
# Plot the analytical sphere
u = np.linspace(0, 2 * np.pi, 30)
v = np.linspace(0, np.pi, 30)
x = radius * np.outer(np.cos(u), np.sin(v))
y = radius * np.outer(np.sin(u), np.sin(v))
z = radius * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot analytical sphere surface
ax.plot_surface(
    x + center_x, 
    y + center_y, 
    z + center_z,
    color='green',
    alpha=0.1,
    label='Analytical sphere'
)

# ===== PLOT CENTER POINT =====
# Add the sphere center point
ax.scatter([center_x], [center_y], [center_z], 
           color='black', marker='*', s=200, label='Sphere center')

# ===== VISUALIZATION SETTINGS =====
# Set equal aspect ratio
ax.set_box_aspect([1, 1, 1])
ax.set_xlabel('X', fontsize=12)
ax.set_ylabel('Y', fontsize=12)
ax.set_zlabel('Z', fontsize=12)
ax.set_title('3D Interface Centroids and Cell Geometry', fontsize=16)

# Add custom legend with larger markers
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='*', color='w', markerfacecolor='black', 
           markersize=15, label='Sphere center'),
    Line2D([0], [0], marker='^', color='w', markerfacecolor='blue', 
           markersize=10, label='Cell centers'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='red', 
           markersize=15, label='Interface centroids'),
    Line2D([0], [0], linestyle='none', marker='s', color='w', 
           markerfacecolor='lightblue', markersize=10, alpha=0.3, label='Mesh cells'),
    Line2D([0], [0], color='green', alpha=0.5, label='Analytical sphere')
]
ax.legend(handles=legend_elements, loc='best', fontsize=12)

# Add grid
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('sphere_interface_visualization.png', dpi=300)
plt.show()

# ===== ADDITIONAL ANALYSIS PLOTS =====
# Create a figure with 2 subplots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

# Plot 1: Distance distribution histogram
ax1.hist(data['distance'], bins=30, alpha=0.7, color='skyblue', edgecolor='black')
ax1.axvline(radius, color='r', linestyle='--', linewidth=2, 
           label=f"Expected radius: {radius:.3f}")
ax1.set_xlabel('Distance from center', fontsize=12)
ax1.set_ylabel('Frequency', fontsize=12)
ax1.set_title('Distribution of Interface Centroid Distances', fontsize=14)
ax1.legend(fontsize=12)
ax1.grid(True, alpha=0.3)

# Plot 2: 3D error scatter plot
# Calculate error vectors
error_x = data['centroid_x'] - center_x
error_y = data['centroid_y'] - center_y
error_z = data['centroid_z'] - center_z

# Normalize the errors by the radius to get a better view
error_magnitude = np.sqrt(error_x**2 + error_y**2 + error_z**2)
max_error = error_magnitude.max()

# Create a scatter plot of errors
scatter = ax2.scatter(error_x, error_y, c=error_magnitude, 
                      cmap='viridis', s=50, alpha=0.7)
ax2.set_xlabel('Error in X', fontsize=12)
ax2.set_ylabel('Error in Y', fontsize=12)
ax2.set_title('Distribution of Centroid Position Errors', fontsize=14)
ax2.grid(True, alpha=0.3)
ax2.set_aspect('equal')

# Add a reference circle showing r=radius
theta = np.linspace(0, 2*np.pi, 100)
ax2.plot(0.05*np.cos(theta), 0.05*np.sin(theta), 'r--', 
         label='5% radius error', alpha=0.8)

# Add a colorbar
cbar = plt.colorbar(scatter, ax=ax2)
cbar.set_label('Error magnitude', fontsize=12)

ax2.legend(fontsize=12)

plt.tight_layout()
plt.savefig('sphere_error_analysis.png', dpi=300)
plt.show()

# Create an additional visualization showing cell slices
# This will show a 2D slice of the 3D data for better clarity
fig, ax = plt.subplots(figsize=(10, 10))

# Choose a slice near the center in Z direction
z_values = sorted(list(set(data['cell_z'])))
z_median_idx = len(z_values) // 2
selected_z = z_values[z_median_idx]
slice_tolerance = cell_size / 2  # Allow for some tolerance

# Filter data for the selected Z slice
slice_data = data[abs(data['cell_z'] - selected_z) < slice_tolerance]

# Plot cells in the slice
for _, row in slice_data.iterrows():
    rect = plt.Rectangle(
        (row['cell_x'], row['cell_y']), 
        cell_size, cell_size,
        linewidth=1, 
        edgecolor='gray', 
        facecolor='lightblue', 
        alpha=0.2
    )
    ax.add_patch(rect)

# Plot cell centers
ax.scatter(
    slice_data['cell_x'] + cell_size/2, 
    slice_data['cell_y'] + cell_size/2,
    color='blue', 
    marker='+', 
    s=30, 
    label='Cell centers'
)

# Plot interface centroids
ax.scatter(
    slice_data['centroid_x'], 
    slice_data['centroid_y'],
    c='red', 
    s=slice_data['interface_area']*1000, 
    alpha=0.7,
    label='Interface centroids'
)

# Draw circle representing the intersection of the sphere with this plane
circle_radius = np.sqrt(max(0, radius**2 - (selected_z - center_z)**2))
theta = np.linspace(0, 2*np.pi, 100)
ax.plot(
    center_x + circle_radius * np.cos(theta),
    center_y + circle_radius * np.sin(theta),
    'g--',
    label=f'Sphere slice at z={selected_z:.3f}'
)

ax.set_aspect('equal')
ax.grid(True, alpha=0.3)
ax.set_xlabel('X', fontsize=12)
ax.set_ylabel('Y', fontsize=12)
ax.set_title(f'2D Slice at Z = {selected_z:.3f}', fontsize=14)
ax.legend(fontsize=12)

plt.tight_layout()
plt.savefig('sphere_slice_visualization.png', dpi=300)
plt.show()