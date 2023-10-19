#ifndef TPMS_H
#define TPMS_H

#include <vtkImageData.h>
#include <vtkNew.h>
#include <vtkMassProperties.h>
#include <vtkSTLWriter.h>
#include <vtkCleanPolyData.h>

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
	Tpms(int npoints, float scalevtk, int numcellx, int numcelly, int numcellz, char typetpms, double origin[3], float rstart);
	/// Destructor
	~Tpms();

	/**
	 *  \brief Set the properties of the vtk obbjects
	*/
	void TpmsSet();

	/**
	 *  \brief Update the TPMS adding the rstep to the scalar matrix
	 *  @param rstep Float to add to the scalar of the matrix
	*/
	void TpmsUpdate(float rstep);

	/**
	 *  \brief Clean the output mesh
	*/
	void TpmsClean();

	/**
	*  \brief Get the TPMS volume
	*/
	double TpmsVolume();

	/**
	*  \brief Get the TPMS area
	*/
	double TpmsArea();

	/**
	*  \brief Write the Tpms to the stl file
	*  @param filename Output filename
	*/
	void TpmsWriteToSTL(const char* filename);

	/**
	*  \brief Set the pointers to the vtk object of the class
	*  @param volume vtkImageData object used for the TPMS stucture
	*  @param surface vtkMarchingCubes object used for generating the isosurface of the TPMS
	*  @param massproperties vtkMassProperties object used to calculate area and volume of the TPMS
	*/
#ifdef USE_FLYING_EDGES
	void SetVtkObjects(vtkImageData* volume, vtkFlyingEdges3D* surface, vtkMassProperties* massproperties, vtkCleanPolyData* cleanpoly);
#else
	void SetVtkObjects(vtkImageData* volume, vtkMarchingCubes* surface, vtkMassProperties* massproperties, vtkCleanPolyData* cleanpoly);
#endif


protected:

	int nPoints;
	float scaleVtk;
	int numCellX;
	int numCellY;
	int numCellZ;
	char typeTpms;
	double isoValue;
	float rStart;
	double Origin[3] = { 0,0,0 };
	vtkImageData* Volume;
	vtkMassProperties* massProperties;
	vtkCleanPolyData* cleanPoly;
#ifdef USE_FLYING_EDGES
	vtkFlyingEdges3D* Surface;
#else
	vtkMarchingCubes* Surface;
#endif

};

#endif