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
	string out_file;

	if (argc == 3) {
		var_value = readConfiguration(argv[1]);
		out_file = (string) argv[2];
	}
	else {
		var_value = readConfiguration("configuration.txt");
		out_file = "myTPMS.stl";
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

	// Saving TPMS to stl file and display the TPMS

	bool saveSTL = stoi(var_value[6]);
	bool graph = stoi(var_value[7]);


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

	
	vtkNew<vtkTransformPolyDataFilter> finalTPMS = tpms_final.TpmsTransform();

	vtkNew<vtkLinearSubdivisionFilter> boxRefined = tpms_final.TpmsBox(tarSize, origin);

	vtkNew<vtkStaticCleanPolyData> solidTPMS = tpms_final.TpmsSolid(finalTPMS, boxRefined, tarSize);
	vtkNew<vtkStaticCleanPolyData> fluidTPMS = tpms_final.TpmsFluid(finalTPMS, boxRefined, tarSize);

	// Saving to .stl file

	if (saveSTL) {
		string solidName = out_file + "_solid.stl";
		string fluidName = out_file + "_fluid.stl";

		tpms_final.TpmsWriteToSTL(solidName,solidTPMS);
		tpms_final.TpmsWriteToSTL(fluidName,fluidTPMS);
	}


	// Printing execution time

	clock_t t1 = clock();
	printTime(t0, t1);


	// Graphic

#ifdef GRAPHICAL
	if (graph)
		renderSurface(solidTPMS);
#endif // GRAPHICAL


	return EXIT_SUCCESS;

}
