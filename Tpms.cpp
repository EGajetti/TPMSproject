/*
 * This file is part of {{ tpmsGenerator }}.
 *
 * .Package for creating a Tpms object
 *  Autors: E. Gajetti, U. Follo
 */

 // cpp file of Tpms class

#include "Tpms.h"

using namespace std;

// -------- Costruttori e distruttore ---------//

Tpms::Tpms() {
	nPoints = 100;
	scaleVtk = 1;
	numCell = 1;
	typeTpms = 'G';
	rStart = -0.07;
	isoValue = 0.;
}

Tpms::Tpms(int npoints, float scalevtk, int numcell, char typetpms, double origin[3], float rstart) {
	nPoints = npoints;
	scaleVtk = scalevtk;
	numCell = numcell;
	typeTpms = typetpms;
	rStart = rstart;
	isoValue = 0.;
	for (int i = 0; i < 3; i++)
		Origin[i] = origin[i];
}

Tpms::~Tpms() {}




//--------- Metodi della classe --------------//

#ifdef USE_FLYING_EDGES
void Tpms::SetVtkObjects(vtkImageData* volume, vtkFlyingEdges3D* surface, vtkMassProperties* massproperties) {
	Volume = volume;
	Surface = surface;
	massProperties = massproperties;
}
#else


void Tpms::SetVtkObjects(vtkImageData* volume, vtkMarchingCubes* surface, vtkMassProperties* massproperties) {
	Volume = volume;
	Surface = surface;
	massProperties = massproperties;
}
#endif


void Tpms::TpmsSet() {

	Volume->ReleaseData();

	int extent[6] = { -1, nPoints * numCell + 1, -1, nPoints * numCell + 1, -1, nPoints * numCell + 1 };
	Volume->SetExtent(extent);
	Volume->SetOrigin(Origin);

	double spacing[3] = { 1. / nPoints / numCell * scaleVtk, 1. / nPoints / numCell * scaleVtk, 1. / nPoints / numCell * scaleVtk };
	Volume->SetSpacing(spacing);

	Volume->AllocateScalars(VTK_FLOAT, 1);

	TpmsGenerator(nPoints, numCell, typeTpms, rStart, Volume);
}


void Tpms::TpmsUpdate(float rstep) {
	for (int z = 0; z < nPoints * numCell; z++)
		for (int y = 0; y < nPoints * numCell; y++)
			for (int x = 0; x < nPoints * numCell; x++) {
				float* a = static_cast<float*> (Volume->GetScalarPointer(x, y, z));
				*a = *a + rstep;
			}
}


double Tpms::TpmsVolume() {
	Surface->RemoveAllInputs();
	Surface->SetInputData(Volume);
	Surface->ComputeNormalsOn();
	Surface->SetValue(0, isoValue);
	Surface->Update();
	massProperties->SetInputConnection(Surface->GetOutputPort());
	massProperties->Update();
	double stlVol = massProperties->GetVolume();
	cout << "Volume fraction evaluated: " << stlVol << endl;
	return stlVol;
}


double Tpms::TpmsArea() {
	Surface->RemoveAllInputs();
	Surface->SetInputData(Volume);
	Surface->ComputeNormalsOn();
	Surface->SetValue(0, isoValue);
	Surface->Update();
	massProperties->RemoveAllInputs();
	massProperties->SetInputConnection(Surface->GetOutputPort());
	massProperties->Update();
	double stlArea = massProperties->GetSurfaceArea();
	return stlArea;
}


void Tpms::TpmsWriteToSTL(const char* filename) {
	vtkNew<vtkSTLWriter> writer;
	writer->SetInputData(Surface->GetOutput());
	writer->SetFileName(filename);
	writer->SetFileTypeToBinary();
	writer->Update();
}