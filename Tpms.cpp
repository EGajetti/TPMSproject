/*
 * This file is part of {{ tpmsGenerator }}.
 *
 * .Package for creating a Tpms object
 *  Authors: E. Gajetti, U. Follo
 */

 // cpp file of Tpms class

#include "Tpms.h"

using namespace std;

// -------- Costruttori e distruttore ---------//

Tpms::Tpms() {
	nPoints = 100;
	scaleVtk = 1;
	numCellX = 1;
	numCellY = 1;
	numCellZ = 1;
	typeTpms = 'G';
	rStart = -0.07;
	isoValue = 0.;
}

Tpms::Tpms(int npoints, float scalevtk, int numcellx, int numcelly, int numcellz, char typetpms, double origin[3], float rstart) {
	nPoints = npoints;
	scaleVtk = scalevtk;
	numCellX = numcellx;
	numCellY = numcelly;
	numCellZ = numcellz;
	typeTpms = typetpms;
	rStart = rstart;
	isoValue = 0.;
	for (int i = 0; i < 3; i++)
		Origin[i] = origin[i];
}

Tpms::~Tpms() {}



//--------- Class methods --------------//

#ifdef USE_FLYING_EDGES
void Tpms::SetVtkObjects(vtkImageData* volume, vtkFlyingEdges3D* surface, vtkMassProperties* massproperties, vtkCleanPolyData* cleanpoly) {
	Volume = volume;
	Surface = surface;
	massProperties = massproperties;
	cleanPoly = cleanpoly;
}
#else


void Tpms::SetVtkObjects(vtkImageData* volume, vtkMarchingCubes* surface, vtkMassProperties* massproperties, vtkCleanPolyData* cleanpoly) {
	Volume = volume;
	Surface = surface;
	massProperties = massproperties;
	cleanPoly = cleanpoly;
}
#endif


void Tpms::TpmsSet() {

	Volume->ReleaseData();

	// int extent[6] = { -2, nPoints * numCellX + 2, -2, nPoints * numCellY + 2, -2, nPoints * numCellZ + 2 };
	int extent[6] = { -1, nPoints * numCellX , -1, nPoints * numCellY , -1, nPoints * numCellZ };

	Volume->SetExtent(extent);
	Volume->SetOrigin(Origin);

	//double spacing[3] = { 1. / nPoints / numCellX * scaleVtk, 1. / nPoints / numCellY * scaleVtk, 1. / nPoints / numCellZ * scaleVtk };
	double spacing[3] = { 1. / nPoints * scaleVtk, 1. / nPoints * scaleVtk, 1. / nPoints * scaleVtk };

	Volume->SetSpacing(spacing);

	Volume->AllocateScalars(VTK_FLOAT, 1);

	TpmsGenerator(nPoints, numCellX, numCellY, numCellZ, typeTpms, rStart, Volume);
}


void Tpms::TpmsUpdate(float rstep) {
	for (int z = 0; z < nPoints * numCellX; z++)
		for (int y = 0; y < nPoints * numCellY; y++)
			for (int x = 0; x < nPoints * numCellZ; x++) {
				float* a = static_cast<float*> (Volume->GetScalarPointer(x, y, z));
				*a = *a + rstep;
			}
}

void Tpms::TpmsClean() {
	Surface->RemoveAllInputs();
	Surface->SetInputData(Volume);
	Surface->ComputeNormalsOn();
	Surface->SetValue(0, isoValue);
	Surface->Update();
	cleanPoly->SetInputConnection(Surface->GetOutputPort());
	cleanPoly->ConvertLinesToPointsOn();
	cleanPoly->ConvertLinesToPointsOn();
	cleanPoly->ConvertPolysToLinesOn();
	cleanPoly->ConvertStripsToPolysOn();
	cleanPoly->Update();
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

vtkNew<vtkQuadricDecimation> Tpms::TpmsQuadricDecimation() {
	vtkNew<vtkQuadricDecimation> decimate;
	float reduction = 0.9;
  	decimate->SetInputData(Surface->GetOutput());
  	decimate->SetTargetReduction(reduction);
  	decimate->VolumePreservationOn();
  	decimate->Update();
	return decimate;
}

void Tpms::TpmsWriteToSTL(const char* filename, vtkQuadricDecimation* decimate) {
	vtkNew<vtkSTLWriter> writer;
	// writer->SetInputData(Surface->GetOutput());
	//writer->SetInputData(cleanPoly->GetOutput());
	writer->SetInputData(decimate->GetOutput());
	writer->SetFileName(filename);
	writer->SetFileTypeToBinary();
	writer->Update();
}