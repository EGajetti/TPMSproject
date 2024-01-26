#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkMassProperties.h>
#include <vtkQuadricDecimation.h>
#include <vtkCubeSource.h>
#include <vtkAppendPolyData.h>
#include <vtkStaticCleanPolyData.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkSphereSource.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkTriangleFilter.h>
#include <vtkLinearSubdivisionFilter.h>

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

	int numCellX = stoi(var_value[3]);
	int numCellY = stoi(var_value[4]);
	int numCellZ = stoi(var_value[5]);

	float tarSize = stof(var_value[6]);

	double* origin = convertOrigin(var_value[7]);


	// Calibration of volume fraction

	float rstart = stof(var_value[8]);

	// Thickness of the walls
	// Upper wall (+1.0 to have a thicker wall, then it will cut via blockMesh)
	float thick1 = stof(var_value[9]) + 1.0;
	// Lower wall
	float thick2 = stof(var_value[10]) + 1.0;

	// Saving TPMS to stl file and display the TPMS

	bool saveSTL = stoi(var_value[11]);
	bool graph = stoi(var_value[12]);


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

	// Reduce the meshing to improve the grid quality
	vtkNew<vtkQuadricDecimation> decimate = tpms_final.TpmsQuadricDecimation();

	// Create two boxes (upper and lower walls)
	double originCubo1[3];
	for (int i = 0; i < 2; i++)
		originCubo1[i] = origin[i] + tarSize/2.0;
	originCubo1[2] = origin[2] + tarSize + thick1/2.0;

	// vtkNew<vtkCubeSource> cubo1;
	// cubo1->SetXLength(tarSize + 1.0);
	// cubo1->SetYLength(tarSize + 1.0);
	// cubo1->SetZLength(thick1);
	// cubo1->SetCenter(originCubo1);
	// cubo1->Update();

	// double originCubo2[3];
	// for (int i = 0; i < 2; i++)
	// 	originCubo2[i] = origin[i] + tarSize/2.0;
	// originCubo2[2] = origin[2] - thick2/2.0;

	// vtkNew<vtkCubeSource> cubo2;
	// cubo2->SetXLength(tarSize + 1.0);
	// cubo2->SetYLength(tarSize + 1.0);
	// cubo2->SetZLength(thick2);
	// cubo2->SetCenter(originCubo2);
	// cubo2->Update();

	// vtkNew<vtkAppendPolyData> appendi;
	// appendi->AddInputData(cubo1->GetOutput());
	// appendi->AddInputData(cubo2->GetOutput());
	// appendi->AddInputData(decimate->GetOutput());
	// appendi->Update();

	// vtkNew<vtkStaticCleanPolyData> cleanAppend;
	// cleanAppend->SetInputConnection(appendi->GetOutputPort());
	// cleanAppend->SetTolerance(1e-3);
	// cleanAppend->Update();

	// vtkNew<vtkTransformPolyDataFilter> TPMSsalva = tpms_final.TpmsTransform(cleanAppend);

	vtkNew<vtkTransformPolyDataFilter> TPMStrasforma = tpms_final.TpmsTransform(decimate);

	vtkNew<vtkCubeSource> box;
	box->SetXLength(tarSize);
	box->SetYLength(tarSize);
	box->SetZLength(tarSize);
	box->SetCenter(0.0, 0.0, 0.0);
	box->Update();


	auto triang2 = vtkSmartPointer<vtkTriangleFilter>::New();
	triang2->SetInputData(box->GetOutput());
	triang2->Update();

	vtkNew<vtkLinearSubdivisionFilter> boxRefined;
	boxRefined->SetInputData(triang2->GetOutput());
	boxRefined->SetNumberOfSubdivisions(5);
	boxRefined->Update();

   vtkNew<vtkBooleanOperationPolyDataFilter> booleanOperation;
   booleanOperation->SetOperationToIntersection();
   booleanOperation->SetInputData(0, TPMStrasforma->GetOutput());
   booleanOperation->SetInputData(1, boxRefined->GetOutput());
   booleanOperation->Update();

	// vtkNew<vtkStaticCleanPolyData> cleanedTPMS;
	// cleanedTPMS->SetInputConnection(booleanOperation->GetOutputPort());
	// cleanedTPMS->SetTolerance(1e-3);
	// cleanedTPMS->Update();

	// vtkNew<vtkSTLWriter> writerTPMS;
	// writerTPMS->SetInputData(cleanedTPMS->GetOutput());
	// writerTPMS->SetFileName("myTPMS_cleaned.stl");
	// writerTPMS->SetFileTypeToBinary();
	// writerTPMS->Update();


	// double volFracFinal = stlVol / (tarSize * tarSize * tarSize);

	// cout << "Volume TPMS: " << volFracFinal << endl;


	// Saving to .stl file

	if (saveSTL) {
		tpms_final.TpmsWriteToSTL(out_file,booleanOperation);
	}


	// Printing execution time

	clock_t t1 = clock();
	printTime(t0, t1);


	// Graphic

#ifdef GRAPHICAL
	if (graph)
		renderSurface(surface, booleanOperation);
#endif // GRAPHICAL


	return EXIT_SUCCESS;

}
