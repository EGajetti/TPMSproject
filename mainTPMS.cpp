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
	// TPMS topology
	char TPMSname = var_value[0][0];
	// sheet or solid type
	string type = var_value[1];
	// Target size
	float tarSize = stof(var_value[2]);
	// f(x,y,z) = rvalue
	float rvalue = stof(var_value[3]);
	// Num cells along x,y,z directions
	int numCellX = stoi(var_value[4]);
	int numCellY = stoi(var_value[5]);
	int numCellZ = stoi(var_value[6]);

	// double* origin = convertOrigin(var_value[3]);
	// double trasla = (tarSize + 1.0)/nFinal;
	// To have the TPMS centered in (0,0,0)
	double trasla = 50./nFinal*tarSize/2.;
	// double origin[3] = {-numCellX*tarSize/2.0 - trasla, -numCellY*tarSize/2.0 - trasla, -numCellZ*tarSize/2.0 - trasla};
	double origin[3] = {-numCellX*tarSize/2.0 - trasla, -numCellY*tarSize/2.0 - trasla, -numCellZ*tarSize/2.0 - trasla};

	// Boolean to rotate the TPMS
	bool ruota = stoi(var_value[7]);	

	// VTK objects
	vtkNew<vtkImageData> volume;
#ifdef USE_FLYING_EDGES
	vtkNew<vtkFlyingEdges3D> surface;
#else
	vtkNew<vtkMarchingCubes> surface;
#endif

	// Generation of the  TPMS lattice object
	Tpms tpms_final(nFinal, tarSize, numCellX, numCellY, numCellZ, TPMSname, origin, rvalue);
	tpms_final.SetVtkObjects(volume, surface);
	tpms_final.TpmsSet(type);

	vtkNew<vtkQuadricDecimation> finalTPMS = tpms_final.TpmsQuadricDecimation();

	if (ruota) {
		double* angles = convertOrigin(var_value[8]);
		vtkNew<vtkTransformPolyDataFilter> rotateTPMS = tpms_final.TpmsTransform(finalTPMS, angles);
		// Saving to .stl file
		tpms_final.TpmsWriteToSTL(out_file,rotateTPMS);

	}
	else {
		// Saving to .stl file
		tpms_final.TpmsWriteToSTL(out_file,finalTPMS);
	}

	// vtkNew<vtkLinearSubdivisionFilter> boxRefined = tpms_final.TpmsBox(tarSize, origin);
	// vtkNew<vtkStaticCleanPolyData> fluidTPMS = tpms_final.TpmsFluid(finalTPMS, boxRefined, tarSize);


	// Printing execution time

	clock_t t1 = clock();
	printTime(t0, t1);


	return EXIT_SUCCESS;

}
