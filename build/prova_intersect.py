import trimesh
import numpy as np

def intersect_with_cube(input_stl, output_stl, cube_size=10.0, scale_factor=np.sqrt(2.0)):
    # Load the input STL file
    input_mesh = trimesh.load(input_stl)
    theta = np.pi/6
    # Define the cube geometry
    # cube = trimesh.creation.box(extents=(cube_size, cube_size*scale_factor*0.99, cube_size*scale_factor*0.99))
    # cube = trimesh.creation.box(extents=(cube_size, cube_size/np.cos(theta), cube_size*np.sqrt(2)/np.cos(np.pi/4-theta) ))
    cube = trimesh.creation.box(extents=(cube_size, cube_size*5, cube_size*5))

    # Perform the intersection
    intersected_mesh = input_mesh.intersection(cube)
    
    # Save the resulting mesh to a new STL file
    intersected_mesh.export(output_stl)

# Example usage
input_stl_file = 'G50_5x5.stl'
output_stl_file = 'output2.stl'
theta = np.pi/6
LL = 10.
scale_f = 1/np.cos(theta)
intersect_with_cube(input_stl_file, output_stl_file)
