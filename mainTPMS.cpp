#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkMassProperties.h>
#include <vtkQuadricDecimation.h>


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
		out_file = "myTPMSFinal.stl";
	}

	clock_t t0 = clock();


	// Input

	int nFinal = stoi(var_value[0]);

	char TPMSname = var_value[1][0];
	string type = var_value[2];

	int numCellX = stoi(var_value[3]);
	int numCellY = stoi(var_value[4]);
	int numCellZ = stoi(var_value[5]);

	float tarSize = stof(var_value[6]);

	double* origin = convertOrigin(var_value[7]);


	// Calibration of volume fraction

	float rstart = stof(var_value[8]);


	// Saving TPMS to stl file and display the TPMS

	bool saveSTL = stoi(var_value[9]);
	bool graph = stoi(var_value[10]);


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
	tpms_final.TpmsSet(type);

	double stlVol = tpms_final.TpmsVolume();
	double stlArea = tpms_final.TpmsArea();

	// Reducing mesh size
	vtkNew<vtkQuadricDecimation> decimate = tpms_final.TpmsQuadricDecimation(surface);

	// Intersecting to cut the coarse edges
	vtkNew<vtkPolyDataBooleanFilter> intersezione = tpms_final.TpmsIntersect(decimate);


	double volFracFinal = stlVol / (tarSize * tarSize * tarSize);

	cout << "Volume TPMS: " << volFracFinal << endl;



	vtkNew<vtkSTLWriter> writer;
	writer->SetInputData(decimate->GetOutput());
	// writer->SetInputData(intersezione->GetOutput());
	writer->SetFileName("Gyroid.stl");
	writer->SetFileTypeToBinary();
	writer->Update();


	// Saving to .stl file

	if (saveSTL) {
		tpms_final.TpmsWriteToSTL(out_file,intersezione);
	}


	// Printing execution time

	clock_t t1 = clock();
	printTime(t0, t1);


	// Graphic

#ifdef GRAPHICAL
	if (graph)
		renderSurface(surface, decimate);
#endif // GRAPHICAL


	return EXIT_SUCCESS;

}
