#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkMassProperties.h>
#include <vtkQuadricDecimation.h>
#include <vtkAppendPolyData.h>
#include <vtkStaticCleanPolyData.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkSTLWriter.h>


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

	// int nFinal = stoi(var_value[0]);
	int nFinal = 100;

	char TPMSname = var_value[0][0];
	string type = var_value[1];

	// int numCellX = stoi(var_value[3]);
	// int numCellY = stoi(var_value[4]);
	// int numCellZ = stoi(var_value[5]);

	int numCellX = 1;
	int numCellY = 1;
	int numCellZ = 1;

	float tarSize = stof(var_value[2]);

	// double* origin = convertOrigin(var_value[3]);
	// double trasla = (tarSize + 1.0)/nFinal;
	double trasla = 50./nFinal*tarSize/2.;
	// double origin[3] = {-numCellX*tarSize/2.0 - trasla, -numCellY*tarSize/2.0 - trasla, -numCellZ*tarSize/2.0 - trasla};
	double origin[3] = {-numCellX*tarSize/2.0 - trasla, -numCellY*tarSize/2.0 - trasla, -numCellZ*tarSize/2.0 - trasla};

	float rvalue = stof(var_value[3]);

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

	
	// vtkNew<vtkTransformPolyDataFilter> finalTPMS = tpms_final.TpmsTransform();
	vtkNew<vtkQuadricDecimation> finalTPMS = tpms_final.TpmsQuadricDecimation();

	// vtkNew<vtkLinearSubdivisionFilter> boxRefined = tpms_final.TpmsBox(tarSize, origin);
	// vtkNew<vtkStaticCleanPolyData> fluidTPMS = tpms_final.TpmsFluid(finalTPMS, boxRefined, tarSize);


	// Saving to .stl file
	tpms_final.TpmsWriteToSTL(out_file,finalTPMS);



	// Printing execution time

	clock_t t1 = clock();
	printTime(t0, t1);


	return EXIT_SUCCESS;

}
