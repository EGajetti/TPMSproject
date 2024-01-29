#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkMassProperties.h>
#include <vtkQuadricDecimation.h>
#include <vtkAppendPolyData.h>
#include <vtkStaticCleanPolyData.h>

#include "Definition.h"
#include "Utils.h"
#include "Tpms.h"

#ifdef USE_FLYING_EDGES
#include <vtkFlyingEdges3D.h>
#else
#include <vtkMarchingCubes.h>
#endif

using namespace std;



int main(int argc, char* argv[])
{
	string* var_value;
	char* out_file;

	if (argc == 3) {
		var_value = readConfiguration(argv[1]);
		out_file = argv[2];
	}
	else {
		var_value = readConfiguration("configuration.txt");
		out_file = (char*) "myTPMS.stl";
	}

	clock_t t0 = clock();


	// Input

	int nFinal = stoi(var_value[0]);

	char TPMSname = var_value[1][0];
	string type = var_value[2];

	// int numCellX = stoi(var_value[3]);
	// int numCellY = stoi(var_value[4]);
	// int numCellZ = stoi(var_value[5]);

	int numCellX = 2;
	int numCellY = 2;
	int numCellZ = 2;

	float tarSize = stof(var_value[3]);

	double* origin = convertOrigin(var_value[4]);

	float rvalue = stof(var_value[5]);

	// Thickness of the walls
	// Upper wall (+1.0 to have a thicker wall, then it will cut via blockMesh)
	float thick1 = stof(var_value[6]) + 1.0;
	// Lower wall
	float thick2 = stof(var_value[7]) + 1.0;

	// Saving TPMS to stl file and display the TPMS

	bool saveSTL = stoi(var_value[8]);
	bool graph = stoi(var_value[9]);


	// Vtk objects

	vtkNew<vtkImageData> volume;
#ifdef USE_FLYING_EDGES
	vtkNew<vtkFlyingEdges3D> surface;
#else
	vtkNew<vtkMarchingCubes> surface;
#endif

	// Generation of the final TPMS lattice object

	Tpms tpms_final(nFinal, tarSize, numCellX, numCellY, numCellZ, TPMSname, origin, rvalue);
	tpms_final.SetVtkObjects(volume, surface);
	tpms_final.TpmsSet(type);

	
	vtkNew<vtkAppendPolyData> appendTPMS = tpms_final.TpmsAppend(tarSize, origin, thick1, thick2);

	// Saving to .stl file

	if (saveSTL) {
		tpms_final.TpmsWriteToSTL(out_file,appendTPMS);
	}


	// Printing execution time

	clock_t t1 = clock();
	printTime(t0, t1);


	// Graphic

#ifdef GRAPHICAL
	if (graph)
		renderSurface(surface, appendTPMS);
#endif // GRAPHICAL


	return EXIT_SUCCESS;

}
