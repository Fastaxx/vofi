import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from mpl_toolkits.mplot3d import Axes3D

# Read the data file
with open('4d_cell_data.txt', 'r') as f:
    # Read mesh parameters
    nx, ny, nz, nw = map(int, f.readline().split())
    domain_size = float(f.readline())
    func_choice = int(f.readline())
    
    # Calculate cell size
    dx = domain_size / nx
    
    # Initialize data array
    cell_types = np.zeros((nx, ny, nz, nw), dtype=int)
    
    # Read cell data
    for line in f:
        i, j, k, l, cell_type = map(int, line.split())
        cell_types[i, j, k, l] = cell_type

# Get function name
function_names = {
    1: "Hyperplane (w = 0.5)",
    2: "Hypersphere (radius = 0.4)",
    3: "4D Torus"
}
func_name = function_names.get(func_choice, "Unknown Function")

# Count cell types
empty_cells = np.sum(cell_types == 0)
full_cells = np.sum(cell_types == 1)
interface_cells = np.sum(cell_types == -1)
total_cells = nx * ny * nz * nw

print(f"\nResults for {func_name}:")
print(f"Total cells: {total_cells}")
print(f"Empty cells: {empty_cells} ({100.0*empty_cells/total_cells:.2f}%)")
print(f"Full cells: {full_cells} ({100.0*full_cells/total_cells:.2f}%)")
print(f"Interface cells: {interface_cells} ({100.0*interface_cells/total_cells:.2f}%)")

# Create mesh grid for plotting
x = np.linspace(dx/2, domain_size-dx/2, nx)
y = np.linspace(dx/2, domain_size-dx/2, ny)
z = np.linspace(dx/2, domain_size-dx/2, nz)
w = np.linspace(dx/2, domain_size-dx/2, nw)

# Create visualizations
fig = plt.figure(figsize=(15, 10))
plt.suptitle(f"4D Cell Classification - {func_name}", fontsize=16)

# Create custom colormap: empty=white, full=blue, interface=red
cmap = ListedColormap(['white', 'blue', 'red'])
norm = plt.Normalize(-1, 1)

# Choose middle slices for visualization
w_slice = nw//2
z_slice = nz//2

# 2D XY slice
ax1 = fig.add_subplot(221)
slice_xy = cell_types[:,:,z_slice,w_slice]
im = ax1.imshow(slice_xy.T, origin='lower', cmap=cmap, norm=norm, 
               extent=[0, domain_size, 0, domain_size])
ax1.set_title(f"XY Slice at z=0.5, w=0.5")
ax1.set_xlabel("x")
ax1.set_ylabel("y")

# 3D visualization at w=0.5
ax2 = fig.add_subplot(222, projection='3d')
x_grid, y_grid, z_grid = np.meshgrid(x, y, z, indexing='ij')
types_3d = cell_types[:,:,:,w_slice]

# Plot full cells (blue)
mask_full = types_3d == 1
if np.any(mask_full):
    ax2.scatter(x_grid[mask_full], y_grid[mask_full], z_grid[mask_full], 
               c='blue', marker='o', alpha=0.7, label='Full')

# Plot interface cells (red)
mask_interface = types_3d == -1
if np.any(mask_interface):
    ax2.scatter(x_grid[mask_interface], y_grid[mask_interface], z_grid[mask_interface], 
               c='red', marker='s', alpha=0.9, label='Interface')

ax2.set_title(f"3D View at w=0.5")
ax2.set_xlabel("x")
ax2.set_ylabel("y")
ax2.set_zlabel("z")
ax2.legend()

# XZ slice
ax3 = fig.add_subplot(223)
y_slice = ny//2
slice_xz = cell_types[:,y_slice,:,w_slice]
ax3.imshow(slice_xz.T, origin='lower', cmap=cmap, norm=norm,
          extent=[0, domain_size, 0, domain_size])
ax3.set_title(f"XZ Slice at y=0.5, w=0.5")
ax3.set_xlabel("x")
ax3.set_ylabel("z")

# YZ slice
ax4 = fig.add_subplot(224)
x_slice = nx//2
slice_yz = cell_types[x_slice,:,:,w_slice]
ax4.imshow(slice_yz.T, origin='lower', cmap=cmap, norm=norm,
          extent=[0, domain_size, 0, domain_size])
ax4.set_title(f"YZ Slice at x=0.5, w=0.5")
ax4.set_xlabel("y")
ax4.set_ylabel("z")

plt.colorbar(im, ax=[ax1, ax3, ax4], label='Cell Type: Empty (0), Full (1), Interface (-1)')
plt.tight_layout(rect=[0, 0, 1, 0.95])

# Save the visualization
plt.savefig(f'4d_visualization_{func_name.replace(" ", "_")}.png', dpi=300)
plt.show()