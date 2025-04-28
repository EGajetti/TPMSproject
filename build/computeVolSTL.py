import numpy as np
import trimesh
import sys


def intersect_with_cube(input_stl, cube_size=10.0):
    # Load the input STL file
    input_mesh = trimesh.load(input_stl)

    # Define the cube geometry
    cube = trimesh.creation.box(extents=(cube_size, cube_size, cube_size))
      
    # Perform the intersection
    intersected_mesh = input_mesh.intersection(cube)
    
    # Save the resulting mesh to a new STL file
    # intersected_mesh.export(output_stl)
    return intersected_mesh

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python3 computeVolSTL.py <input_STL>")
        sys.exit(1)

    input_STL = str(sys.argv[1])
    input_STL = str(sys.argv[1])
    LL = 10.0
    int_mesh = intersect_with_cube(input_STL, cube_size=10.0)

    por = int_mesh.volume/LL**3
    print(f'The porosity of {input_STL} is {por:.4f}')