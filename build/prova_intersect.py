import trimesh
import numpy as np

def intersect_with_cube(input_stl, output_stl, cube_size=14.0, scale_factor=np.sqrt(2.0)):
    # Load the input STL file
    input_mesh = trimesh.load(input_stl)
    
    # Scale the input mesh by the specified scale factor
    # input_mesh.apply_scale(scale_factor)
    
    # Define the cube geometry
    cube = trimesh.creation.box(extents=(cube_size, cube_size*scale_factor, cube_size*scale_factor))
    
    # Position the cube at the center of the input mesh's bounding box
    # cube.apply_translation(input_mesh.bounds.mean(axis=0))
    
    # Perform the intersection
    intersected_mesh = input_mesh.intersection(cube)
    
    # Save the resulting mesh to a new STL file
    intersected_mesh.export(output_stl)

# Example usage
input_stl_file = 'myTPMS.stl'
output_stl_file = 'output2.stl'
intersect_with_cube(input_stl_file, output_stl_file)
