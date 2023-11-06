/*
 * This file is part of {{ tpmsGenerator }}.
 *
 * .Package for creating a Tpms object
 *  Autors: E. Gajetti, U. Follo
 */

#pragma once
#ifndef UTILS_H
#define UTILS_H

#include <iostream>
#include <vtkImageData.h>
#include "Definition.h"

 //Libraries needed for showing the isosurface
#ifdef GRAPHICAL
#include <vtkNew.h>
#include <vtkFlyingEdges3D.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkProperty.h>
#include <vtkQuadricDecimation.h>
#include <vtkPolyDataBooleanFilter.h>

#endif // GRAPHICAL

using namespace std;



/**
* \brief Functions for setting scalar values in a VTKImageDATA given a function:
*
* @param npoints	Dimension of the VTKImagedata
* @param numcellx	Number of cells along x
* @param numcelly	Number of cells along y
* @param numcellz	Number of cells along z
* @param type       Char selecting the TPMS type
* @param rstart     Start parameter
* @param volume	    VTKImageData
*/
void TpmsSolidGenerator(const int npoints, const int numcellx, const int numcelly, const int numcellz, char type, const float rstart, vtkImageData* volume);

/**
* \brief Functions for setting scalar values in a VTKImageDATA given a function:
*
* @param npoints	Dimension of the VTKImagedata
* @param numcellx	Number of cells along x
* @param numcelly	Number of cells along y
* @param numcellz	Number of cells along z
* @param type       Char selecting the TPMS type
* @param rstart     Start parameter
* @param volume	    VTKImageData
*/
void TpmsSheetGenerator(const int npoints, const int numcellx, const int numcelly, const int numcellz, char type, const float rstart, vtkImageData* volume);



/**
 * \brief Function for reading the configuration file
 *
 * @param filename	Name of the configuration file
*/
string* readConfiguration(const char* filename);


/**
 * \brief Function for converting the origin from the configuration file
 *
 * @param origin	String of the origin
*/
double* convertOrigin(string origin);


/**
 * \brief Function for printing elapsed time in seconds
 *
 * @param start	Start time
 * @param end	Stop time
*/
void printTime(clock_t start, clock_t end);


/**
* \brief Function wich prints the isosurface.
*
* @param surface	vtkFlyingEdges3D object which contains the isosurface
*/
#ifdef GRAPHICAL
// void renderSurface(vtkFlyingEdges3D* surface, vtkQuadricDecimation* decimate);
void renderSurface(vtkFlyingEdges3D* surface, vtkPolyDataBooleanFilter* intersezione);
#endif // GRAPHICAL


//-------------- Old functions no more used --------------//

/**
* \brief Function for reading a txt file.The line structur needs to be:
* 			int int int float
*
* @param filename	Name of the input file
* @param extent[6]	Boundaries of the isosurface colume
* @param origin[3]	Origin of the isosurface volume
* @param spacing[3]	Thickness of the grid of the isosurface volume
*/
void readVolumeFromFile(const char* filename,
	int extent[6], const double origin[3], const double spacing[3],
	vtkImageData* volume);


/**
* \brief Function for writin on a txt file an vtkImageData object.
*
* @param filename	Name of the output file
* @param volume		vtkImageData object containing the isosurface
*/
void writeVolumeToFile(const char* filename, vtkImageData* volume);


#endif // !utils
