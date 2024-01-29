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
	numCellX = 1;
	numCellY = 1;
	numCellZ = 1;
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
	int extent[6] = { -5, nPoints * numCellX + 5, -5, nPoints * numCellY + 5, -5, nPoints * numCellZ + 5 };
	// int extent[6] = { -1, nPoints * numCellX , -1, nPoints * numCellY , -1, nPoints * numCellZ };

	Volume->SetExtent(extent);
	Volume->SetOrigin(Origin);

	// Using the first one will result in deformed cells
	//double spacing[3] = { 1. / nPoints / numCellX * scaleVtk, 1. / nPoints / numCellY * scaleVtk, 1. / nPoints / numCellZ * scaleVtk };
	double spacing[3] = { 1. / nPoints * scaleVtk, 1. / nPoints * scaleVtk, 1. / nPoints * scaleVtk };

	Volume->SetSpacing(spacing);

	Volume->AllocateScalars(VTK_FLOAT, 1);

	if (type == "solid"){
		TpmsSolidGenerator(nPoints, numCellX, numCellY, numCellZ, typeTpms, rValue, Volume);

	}
	else if (type == "sheet"){
		TpmsSheetGenerator(nPoints, numCellX, numCellY, numCellZ, typeTpms, rValue, Volume);
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
	float reduction = 0.5;
  	decimate->SetInputData(cleaned->GetOutput());
  	decimate->SetTargetReduction(reduction);
  	decimate->VolumePreservationOn();
  	decimate->Update();
	return decimate;
}

	vtkNew<vtkTransformPolyDataFilter> Tpms::TpmsTransform() {
	vtkNew<vtkQuadricDecimation> decimateCubo = TpmsQuadricDecimation();
	vtkNew<vtkTransform> trasformazione;
	trasformazione->Translate(-numCellX*scaleVtk/2.0, -numCellY*scaleVtk/2.0, -numCellZ*scaleVtk/2.0);
	vtkNew<vtkTransformPolyDataFilter> trasformaCubo;
	trasformaCubo->SetTransform(trasformazione);
	trasformaCubo->SetInputData(decimateCubo->GetOutput());
	trasformaCubo->Update();
	return trasformaCubo;
}

vtkNew<vtkBooleanOperationPolyDataFilter> Tpms::TpmsIntersecting(float tarSize, double* origin){
	vtkNew<vtkCubeSource> box;
	box->SetXLength(tarSize + 0.1*tarSize);
	box->SetYLength(tarSize + 0.1*tarSize);
	box->SetZLength(tarSize + 0.1*tarSize);
	box->SetCenter(origin);
	box->Update();

	// convert the box in a triangular grid
	vtkNew<vtkTriangleFilter> boxTriang;
	boxTriang->SetInputData(box->GetOutput());
	boxTriang->Update();

	vtkNew<vtkLinearSubdivisionFilter> boxRefined;
	boxRefined->SetInputData(boxTriang->GetOutput());
	boxRefined->SetNumberOfSubdivisions(5);
	boxRefined->Update();

	vtkNew<vtkTransformPolyDataFilter> translateTPMS = TpmsTransform();
    vtkNew<vtkBooleanOperationPolyDataFilter> intersectTPMS;
    intersectTPMS->SetOperationToIntersection();
    intersectTPMS->SetInputData(0, translateTPMS->GetOutput());
    intersectTPMS->SetInputData(1, boxRefined->GetOutput());
    intersectTPMS->Update();

	return intersectTPMS;
	}

vtkNew<vtkAppendPolyData> Tpms::TpmsAppend(float tarSize, double* origin, float thick1, float thick2){
	// Create two boxes (upper and lower walls)
	double originCubo1[3];
	for (int i = 0; i < 2; i++)
		originCubo1[i] = origin[i];
	originCubo1[2] = origin[2] + tarSize/2.0 + thick1/2.0;

	vtkNew<vtkCubeSource> cubo1;
	cubo1->SetXLength(tarSize + 0.1*tarSize);
	cubo1->SetYLength(tarSize + 0.1*tarSize);
	cubo1->SetZLength(thick1);
	cubo1->SetCenter(originCubo1);
	cubo1->Update();

	double originCubo2[3];
	for (int i = 0; i < 2; i++)
		originCubo2[i] = origin[i];
	originCubo2[2] = origin[2] - tarSize/2.0 - thick2/2.0;

	vtkNew<vtkCubeSource> cubo2;
	cubo2->SetXLength(tarSize + 0.1*tarSize);
	cubo2->SetYLength(tarSize + 0.1*tarSize);
	cubo2->SetZLength(thick2);
	cubo2->SetCenter(originCubo2);
	cubo2->Update();

	// Create the TPMS
	vtkNew<vtkBooleanOperationPolyDataFilter> intersectTPMS = TpmsIntersecting(tarSize, origin);

	vtkNew<vtkAppendPolyData> appendTPMS;
	appendTPMS->AddInputData(cubo1->GetOutput());
	appendTPMS->AddInputData(cubo2->GetOutput());
	appendTPMS->AddInputData(intersectTPMS->GetOutput());
	appendTPMS->Update();

	return appendTPMS;
}


void Tpms::TpmsWriteToSTL(const char* filename, vtkAppendPolyData* appendTPMS) {	

	vtkNew<vtkSTLWriter> writer;
	writer->SetInputData(appendTPMS->GetOutput());
	writer->SetFileName(filename);
	writer->SetFileTypeToBinary();
	writer->Update();
}


