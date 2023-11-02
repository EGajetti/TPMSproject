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


void Tpms::TpmsSet(string type) {

	Volume->ReleaseData();
	// Better to have a sligthly larger extent, so to cut coarse edges a posteriori
	int extent[6] = { -2, nPoints * numCellX + 2, -2, nPoints * numCellY + 2, -2, nPoints * numCellZ + 2 };
	// int extent[6] = { -1, nPoints * numCellX , -1, nPoints * numCellY , -1, nPoints * numCellZ };

	Volume->SetExtent(extent);
	Volume->SetOrigin(Origin);

	// Using the first one will result in deformed cells
	//double spacing[3] = { 1. / nPoints / numCellX * scaleVtk, 1. / nPoints / numCellY * scaleVtk, 1. / nPoints / numCellZ * scaleVtk };
	double spacing[3] = { 1. / nPoints * scaleVtk, 1. / nPoints * scaleVtk, 1. / nPoints * scaleVtk };

	Volume->SetSpacing(spacing);

	Volume->AllocateScalars(VTK_FLOAT, 1);

	if (type == "solid"){
		TpmsSolidGenerator(nPoints, numCellX, numCellY, numCellZ, typeTpms, rStart, Volume);

	}
	else if (type == "sheet"){
		TpmsSheetGenerator(nPoints, numCellX, numCellY, numCellZ, typeTpms, rStart, Volume);
	}
	else{
		cout << "Invalid TPMS type" << endl;
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
	cout << "Evaluated volume: " << stlVol << endl;
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

vtkNew <vtkPolyDataNormals> Tpms::TpmsNormals() {
	vtkNew<vtkPolyDataNormals> normals;
	normals->SetInputConnection(Surface->GetOutputPort());
	normals->FlipNormalsOn();
	normals->Update();
	return normals;
}

vtkNew <vtkStaticCleanPolyData> Tpms::TpmsClean() {
	vtkNew<vtkPolyDataNormals> normals = TpmsNormals();
	vtkNew<vtkStaticCleanPolyData> cleaned;
	cleaned->SetInputConnection(normals->GetOutputPort());
	cleaned->SetTolerance(1e-3);
	cleaned->Update();
	return cleaned;
}

vtkNew<vtkQuadricDecimation> Tpms::TpmsQuadricDecimation(){
	vtkNew <vtkStaticCleanPolyData> cleaned = TpmsClean();
	vtkNew<vtkQuadricDecimation> decimate;
	float reduction = 0.8;
  	decimate->SetInputData(cleaned->GetOutput());
  	decimate->SetTargetReduction(reduction);
  	decimate->VolumePreservationOn();
  	decimate->Update();
	return decimate;
}



// vtkNew<vtkQuadricDecimation> Tpms::TpmsQuadricDecimation() {
// 	vtkNew<vtkQuadricDecimation> decimate;
// 	float reduction = 0.8;
//   	decimate->SetInputData(Surface->GetOutput());
//   	decimate->SetTargetReduction(reduction);
//   	decimate->VolumePreservationOn();
//   	decimate->Update();
// 	return decimate;
// }

// vtkNew <vtkPolyDataNormals> Tpms::TpmsNormals(){
// 	vtkNew<vtkQuadricDecimation> decimate = TpmsQuadricDecimation();
// 	vtkNew<vtkPolyDataNormals> normals;
// 	normals->SetInputConnection(decimate->GetOutputPort());
// 	normals->FlipNormalsOn();
// 	normals->Update();
// 	return normals;
// }

// vtkNew <vtkStaticCleanPolyData> Tpms::TpmsClean(){
// 	vtkNew<vtkPolyDataNormals> normals = TpmsNormals();
// 	vtkNew<vtkStaticCleanPolyData> cleaned;
// 	cleaned->SetInputConnection(normals->GetOutputPort());
// 	cleaned->SetTolerance(1e-2);
// 	cleaned->Update();
// 	return cleaned;
// }

void Tpms::TpmsWriteToSTL(const char* filename, vtkQuadricDecimation* decimate) {
// void Tpms::TpmsWriteToSTL(const char* filename, vtkStaticCleanPolyData* cleaned) {
	vtkNew<vtkSTLWriter> writer;
	writer->SetInputData(decimate->GetOutput());
	// writer->SetInputData(cleaned->GetOutput());
	writer->SetFileName(filename);
	writer->SetFileTypeToBinary();
	writer->Update();
}