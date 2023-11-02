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
	int extent[6] = { -10, nPoints * numCellX + 10, -10, nPoints * numCellY + 10, -10, nPoints * numCellZ + 10 };
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

vtkNew<vtkQuadricDecimation> Tpms::TpmsQuadricDecimation(vtkFlyingEdges3D* intersectTPMS) {
	vtkNew<vtkQuadricDecimation> decimate;
	float reduction = 0.9;
  	decimate->SetInputData(Surface->GetOutput());
	// decimate->SetInputData(intersectTPMS->GetOutput());
  	decimate->SetTargetReduction(reduction);
  	decimate->VolumePreservationOn();
  	decimate->Update();
	return decimate;
}

vtkNew<vtkPolyDataBooleanFilter> Tpms::TpmsIntersect(vtkQuadricDecimation* decimate){

	// Creation of the box to intersect
	vtkNew<vtkCubeSource> cubo;
	// Placing the center at the center or the TPMS (which goes from 0 to numCell*scaleVtk)
	cubo->SetCenter(numCellX*scaleVtk/2.0, numCellY*scaleVtk/2.0, numCellZ*scaleVtk/2.0);
	cout << numCellX*scaleVtk/2.0 << endl;
	cout << numCellX*scaleVtk << endl;
	cubo->SetXLength(static_cast<float> (numCellX*scaleVtk));
	cubo->SetYLength(static_cast<float>  (numCellY*scaleVtk));
	cubo->SetZLength(static_cast<float>  (numCellZ*scaleVtk));
	cubo->Update();

	// Creating the boolean filter
	vtkNew<vtkPolyDataBooleanFilter> intersezione;
	intersezione->SetInputConnection(0, cubo->GetOutputPort());
	intersezione->SetInputData(1, decimate->GetOutput());
	// For some reason, difference works as intersection and viceversa
	intersezione->SetOperModeToDifference();
	intersezione->Update();
	return intersezione;
}



// void Tpms::TpmsWriteToSTL(const char* filename, vtkQuadricDecimation* decimate) {
	void Tpms::TpmsWriteToSTL(const char* filename, vtkPolyDataBooleanFilter* intersezione) {
	vtkNew<vtkSTLWriter> writer;
	// writer->SetInputData(decimate->GetOutput());
	writer->SetInputData(intersezione->GetOutput());
	writer->SetFileName(filename);
	writer->SetFileTypeToBinary();
	writer->Update();
}
