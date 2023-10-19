#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkMassProperties.h>

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
		out_file = "myTPMSCpp.stl";
	}

	clock_t t0 = clock();


	// Input

	int nFinal = stoi(var_value[0]);
	int n = stoi(var_value[1]);

	char TPMSname = var_value[2][0];
	string type = var_value[3];
	float volFrac = stof(var_value[4]);

	int numCellX = stoi(var_value[5]);
	int numCellY = stoi(var_value[6]);
	int numCellZ = stoi(var_value[7]);

	float tarSize = stof(var_value[8]);

	double* origin = convertOrigin(var_value[9]);


	// Calibration of volume fraction

	float rstart = stof(var_value[10]);
	float rStep = stof(var_value[11]);
	float volTol = stof(var_value[12]);


	// Saving TPMS to stl file and display the TPMS

	bool saveSTL = stoi(var_value[13]);
	bool graph = stoi(var_value[14]);


	// Vtk objects

	vtkNew<vtkImageData> volume;
	vtkNew<vtkMassProperties> massProperties;
#ifdef USE_FLYING_EDGES
	vtkNew<vtkFlyingEdges3D> surface;
#else
	vtkNew<vtkMarchingCubes> surface;
#endif


	// Generation of the final TPMS lattice object

	Tpms tpms_final(nFinal, tarSize, numCellX, numCellY, numCellZ, TPMSname, origin, rstart);
	tpms_final.SetVtkObjects(volume, surface, massProperties);
	tpms_final.TpmsSet();

	double stlVol = tpms_final.TpmsVolume();
	double stlArea = tpms_final.TpmsArea();

	double volFracFinal = stlVol / (tarSize * tarSize * tarSize);

	cout << "Final volume fraction TPMS: " << volFracFinal << endl;


	// Saving to .stl file

	if (saveSTL) {
		tpms_final.TpmsWriteToSTL(out_file);
	}


	// Printing execution time

	clock_t t1 = clock();
	printTime(t0, t1);


	// Graphic

#ifdef GRAPHICAL
	if (graph)
		renderSurface(surface);
#endif // GRAPHICAL


	return EXIT_SUCCESS;

}

//For problem closing borders
//#include <vtkFillHolesFilter.h> //For problem closing borders
/*
vtkNew<vtkFillHolesFilter> fillHolesFilter;
fillHolesFilter->SetInputData(surface->GetOutput());
fillHolesFilter->SetHoleSize(1.0);
fillHolesFilter->Update();
*/