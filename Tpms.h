#ifndef TPMS_H
#define TPMS_H

#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkMassProperties.h>
#include <vtkSTLWriter.h>
#include <vtkQuadricDecimation.h>
#include <vtkStaticCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkTransform.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkCubeSource.h>
#include <vtkTriangleFilter.h>
#include <vtkLinearSubdivisionFilter.h>

#include "Definition.h"
#include "Utils.h"


#ifdef USE_FLYING_EDGES
#include <vtkFlyingEdges3D.h>
#else
#include <vtkMarchingCubes.h>
#endif




/**
 * \class TPMS
 *
 * \ingroup tpmsGenerator
 * \brief Concrete class for generating a TPMS
 *
 */
class Tpms {
public:

	/// Standard constructor
	Tpms();
	/// Constructor with initialized itemList
	Tpms(int npoints, float scalevtk, int numcellx, int numcelly, int numcellz, char typetpms, double origin[3], float rvalue);
	/// Destructor
	~Tpms();

	/**
	 *  \brief Set the properties of the vtk obbjects
	*/
	void TpmsSet(string type);


	/**
	*  \brief Get the TPMS volume
	*/
	double TpmsVolume(vtkBooleanOperationPolyDataFilter* intersectTPMS, float tarSize);
	
	/**
	 * \brief Flip normals orientiation
	*/
	vtkNew<vtkPolyDataNormals> TpmsNormals();

	/**
	 * \brief Cleaning the mesh after isosurface
	*/
	vtkNew <vtkStaticCleanPolyData> TpmsClean();

	/**
	 * \brief Reducing mesh size
	*/
	vtkNew <vtkQuadricDecimation> TpmsQuadricDecimation();

	/**
	 * \brief Create a box L x L x L
	 * @param tarSize Target size of the sube length
	 * @param origin The origin of the TPMS (centre)
	*/
	vtkNew<vtkLinearSubdivisionFilter> TpmsBox(float tarSize, double* origin);

	/**
	 * \brief Transform the geometry (i.e. rotating and translating)
	*/
	vtkNew<vtkTransformPolyDataFilter> TpmsTransform();

	/**
	 * \brief Clean the geometry after intersection and difference
	 * @param boolTPMS The boolean operation performed, i.e. intersection (solid domain) or difference (fluid domain) 
	*/
	vtkNew<vtkStaticCleanPolyData> TpmsCleanBoolOper(vtkBooleanOperationPolyDataFilter* boolTPMS);

	/**
	 * \brief Intersecting the TPMS with the cube L x L x L
	 * @param translateTPMS TPMS with the center in origin
	 * @param boxRefined The cube (refined and triangulized)
	 * @param tarSize the desired size of the cube
	*/
	vtkNew<vtkStaticCleanPolyData> TpmsSolid(vtkTransformPolyDataFilter* translateTPMS, vtkLinearSubdivisionFilter* boxRefined, float tarSize);
	
	/**
	 * \brief Subtracting the TPMS from a cube L x L x L
	 * @param translateTPMS TPMS with the center in origin
	 * @param boxRefined The cube (refined and triangulized)
	 * @param tarSize the desired size of the cube
	*/
	vtkNew<vtkQuadricDecimation> TpmsFluid(vtkTransformPolyDataFilter* translateTPMS, vtkLinearSubdivisionFilter* boxRefined, float tarSize);

	/**
	*  \brief Write the Tpms to the stl file
	*  @param filename Output filename
	*/
	void TpmsWriteToSTL(const string filename, vtkTransformPolyDataFilter* finalTPMS);
	void TpmsWriteToSTL(const string filename, vtkStaticCleanPolyData* finalTPMS);
	void TpmsWriteToSTL(const string filename, vtkQuadricDecimation* finalTPMS);



	/**
	*  \brief Set the pointers to the vtk object of the class
	*  @param volume vtkImageData object used for the TPMS stucture
	*  @param surface vtkMarchingCubes object used for generating the isosurface of the TPMS
	*/
#ifdef USE_FLYING_EDGES
	void SetVtkObjects(vtkImageData* volume, vtkFlyingEdges3D* surface);
#else
	void SetVtkObjects(vtkImageData* volume, vtkMarchingCubes* surface);
#endif


protected:

	int nPoints;
	float scaleVtk;
	int numCellX;
	int numCellY;
	int numCellZ;
	char typeTpms;
	double isoValue;
	float rValue;
	double Origin[3] = { 0,0,0 };
	vtkImageData* Volume;
	vtkMassProperties* massProperties;
#ifdef USE_FLYING_EDGES
	vtkFlyingEdges3D* Surface;
#else
	vtkMarchingCubes* Surface;
#endif

};

#endif