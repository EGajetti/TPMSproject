/*
 * This file is part of {{ tpmsProject }}.
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
	numCellX = 2;
	numCellY = 2;
	numCellZ = 2;
	typeTpms = 'G';
	rValue = -0.07;
	isoValue = 0.;
}

Tpms::Tpms(int npoints, float scalevtk, int numcellx, int numcelly, int numcellz, char typetpms, double origin[3], float rvalue) {
	nPoints = npoints;
	scaleVtk = scalevtk;
	numCellX = numcellx;
	numCellY = numcelly;
	numCellZ = numcellz;
	typeTpms = typetpms;
	rValue = rvalue;
	isoValue = 0.;
	for (int i = 0; i < 3; i++)
		Origin[i] = origin[i];
}

Tpms::~Tpms() {}



//--------- Class methods --------------//

#ifdef USE_FLYING_EDGES
void Tpms::SetVtkObjects(vtkImageData* volume, vtkFlyingEdges3D* surface) {
	Volume = volume;
	Surface = surface;
}
#else


void Tpms::SetVtkObjects(vtkImageData* volume, vtkMarchingCubes* surface) {
	Volume = volume;
	Surface = surface;
}
#endif


void Tpms::TpmsSet(string type) {

	Volume->ReleaseData();
	// Better to have a sligthly larger extent, so to cut coarse edges a posteriori
	// int extent[6] = { -2*(int)scaleVtk - 1, nPoints * numCellX + 2*(int)scaleVtk + 1, 
	// 				  -2*(int)scaleVtk - 1, nPoints * numCellY + 2*(int)scaleVtk + 1,
	// 				  -2*(int)scaleVtk - 1, nPoints * numCellZ + 2*(int)scaleVtk + 1};
	int extent[6] = { -50, nPoints * numCellX + 50, -50, nPoints * numCellY + 50, -50, nPoints * numCellZ + 50};
	// int extent[6] = { 0, nPoints * numCellX -1, 0, nPoints * numCellY -1, 0, nPoints * numCellZ - 1};

	Volume->SetExtent(extent);
	Volume->SetOrigin(Origin);

	// Using the first one will result in deformed cells
	double spacing[3] = { 1. / nPoints * scaleVtk, 1. / nPoints * scaleVtk, 1. / nPoints * scaleVtk };

	Volume->SetSpacing(spacing);

	Volume->AllocateScalars(VTK_FLOAT, 1);

	if (type == "solid"){
		TpmsSolidGenerator(nPoints, numCellX, numCellY, numCellZ, typeTpms, rValue, Volume, scaleVtk);

	}
	else if (type == "sheet"){
		TpmsSheetGenerator(nPoints, numCellX, numCellY, numCellZ, typeTpms, rValue, Volume, scaleVtk);
	}
	else{
		cout << "Invalid TPMS type" << endl;
	}

	Surface->RemoveAllInputs();
	Surface->SetInputData(Volume);
	Surface->ComputeNormalsOn();
	Surface->SetValue(0, isoValue);
	Surface->Update();
	
}



double Tpms::TpmsVolume(vtkBooleanOperationPolyDataFilter* intersectTPMS, float tarSize) {

	vtkNew<vtkMassProperties> massProperties;

	double cubeVolume = (tarSize*1.01)*(tarSize*1.01)*(tarSize*1.01);

	massProperties->SetInputConnection(intersectTPMS->GetOutputPort());
	massProperties->Update();
	double stlVol = massProperties->GetVolume();

	double porosity = 1.0 - stlVol/cubeVolume;

	cout << "Evaluated volume: " << stlVol << endl;
	cout << "Fluid porosity: " << porosity << endl;
	return porosity;
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
	float reduction = 0.6;
  	decimate->SetInputData(cleaned->GetOutput());
  	decimate->SetTargetReduction(reduction);
  	decimate->VolumePreservationOn();
  	decimate->Update();
	return decimate;
}


// vtkNew<vtkTransformPolyDataFilter> Tpms::TpmsTransform() {
// 	vtkNew<vtkQuadricDecimation> decimateCubo = TpmsQuadricDecimation();
// 	vtkNew<vtkTransform> trasformazione;
// 	trasformazione->Translate(-numCellX*scaleVtk/2.0, -numCellY*scaleVtk/2.0, -numCellZ*scaleVtk/2.0);
// 	vtkNew<vtkTransformPolyDataFilter> trasformaCubo;
// 	trasformaCubo->SetTransform(trasformazione);
// 	trasformaCubo->SetInputData(decimateCubo->GetOutput());
// 	trasformaCubo->Update();
// 	return trasformaCubo;
// }

vtkNew<vtkStaticCleanPolyData> Tpms::TpmsCleanBoolOper(vtkBooleanOperationPolyDataFilter* boolTPMS){
	vtkNew<vtkStaticCleanPolyData> finalTPMS;
	finalTPMS->SetInputData(boolTPMS->GetOutput());
	finalTPMS->SetTolerance(5e-4);
	finalTPMS->ConvertLinesToPointsOn();
	finalTPMS->ConvertLinesToPointsOn();
	finalTPMS->Update();
	return finalTPMS;
}


// vtkNew<vtkQuadricDecimation> Tpms::TpmsFluid(vtkTransformPolyDataFilter* translateTPMS, vtkLinearSubdivisionFilter* boxRefined, float tarSize){
//     vtkNew<vtkBooleanOperationPolyDataFilter> subtractTPMS;
//     subtractTPMS->SetOperationToDifference();
// 	subtractTPMS->SetInputData(0, boxRefined->GetOutput());
//     subtractTPMS->SetInputData(1, translateTPMS->GetOutput());
//     subtractTPMS->Update();
// 	vtkNew<vtkStaticCleanPolyData> fluidTPMS = TpmsCleanBoolOper(subtractTPMS);
// 	vtkNew <vtkQuadricDecimation> fluidTPMS_red;
// 	float reduction = 0.2;
//   	fluidTPMS_red->SetInputData(fluidTPMS->GetOutput());
//   	fluidTPMS_red->SetTargetReduction(reduction);
//   	fluidTPMS_red->VolumePreservationOn();
//   	fluidTPMS_red->Update();
// 	return fluidTPMS_red;
// 	}	


void Tpms::TpmsWriteToSTL(const string filename, vtkTransformPolyDataFilter* finalTPMS) {	
	vtkNew<vtkSTLWriter> writer;
	writer->SetInputData(finalTPMS->GetOutput());
	writer->SetFileName(filename.c_str());
	writer->SetFileTypeToBinary();
	writer->Update();
}

void Tpms::TpmsWriteToSTL(const string filename, vtkQuadricDecimation* finalTPMS) {	
	vtkNew<vtkSTLWriter> writer;
	writer->SetInputData(finalTPMS->GetOutput());
	writer->SetFileName(filename.c_str());
	writer->SetFileTypeToBinary();
	writer->Update();
}



