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
#include <vtkIntersectionPolyDataFilter.h>
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
	double TpmsVolume();

	/**
	*  \brief Get the TPMS area
	*/
	double TpmsArea();

	/**
	 * \brief Cleaning the mesh after isosurface
	*/
	vtkNew <vtkStaticCleanPolyData> TpmsClean();

	/**
	 * \brief Reducing mesh size
	*/
	vtkNew <vtkQuadricDecimation> TpmsQuadricDecimation();

	/**
	 * \brief Flip normals orientiation
	*/
	vtkNew<vtkPolyDataNormals> TpmsNormals();


	/**
	*  \brief Write the Tpms to the stl file
	*  @param filename Output filename
	*/
	void TpmsWriteToSTL(const char* filename, vtkAppendPolyData* appendTPMS);

	/**
	 * \brief Transform the geometry (i.e. rotating and translating)
	*/
	vtkNew<vtkTransformPolyDataFilter> TpmsTransform();

	/**
	 * \brief Intersecting the TPMS with a cube LxLxL
	*/
	vtkNew<vtkBooleanOperationPolyDataFilter> TpmsIntersecting(float tarSize, double* origin);

	/**
	 * \brief Union with walls
	*/
	vtkNew<vtkAppendPolyData> TpmsAppend(float tarSize, double* origin, float thick1, float thick2);

	/**
	*  \brief Set the pointers to the vtk object of the class
	*  @param volume vtkImageData object used for the TPMS stucture
	*  @param surface vtkMarchingCubes object used for generating the isosurface of the TPMS
	*  @param massproperties vtkMassProperties object used to calculate area and volume of the TPMS
	*/
#ifdef USE_FLYING_EDGES
	void SetVtkObjects(vtkImageData* volume, vtkFlyingEdges3D* surface, vtkMassProperties* massproperties);
#else
	void SetVtkObjects(vtkImageData* volume, vtkMarchingCubes* surface, vtkMassProperties* massproperties);
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