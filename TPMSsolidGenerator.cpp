#include <cmath>
#include "Utils.h"

using namespace std;


// Implementation of the function defined in Utils.h

void TpmsSolidGenerator(const int npoints, const int numcellx, const int numcelly, const int numcellz, char type, const float rvalue, vtkImageData* volume, const float scaleVtk)
{
	// int dimx = npoints * numcellx + 2*(int)scaleVtk + 1;
	// int dimy = npoints * numcelly + 2*(int)scaleVtk + 1;
	// int dimz = npoints * numcellz + 2*(int)scaleVtk + 1;
	int dimx = npoints * numcellx + 50;
	int dimy = npoints * numcelly + 50;
	int dimz = npoints * numcellz + 50;
	int dimension = max({dimx, dimy, dimz});

	const float pi2 = 2 * 3.14159265358979323846 / npoints;
	float scal = 0.;

	float* cosv = new float[dimension];
	for (int i = 0; i < dimension; i++)
		cosv[i] = cos(pi2 * i);

	float* senv = new float[dimension];
	for (int i = 0; i < dimension; i++)
		senv[i] = sin(pi2 * i);

	float* cos2 = new float[dimension];
	float* sen2 = new float[dimension];
	
	if ((type == 'I') || (type == 'S') || (type == 'F') || (type == 'K') || (type == 'L'))
		for (int i = 0; i < dimension; i++)
		{
			cos2[i] = cos(2 * pi2 * i);
			sen2[i] = sin(2 * pi2 * i);
		}

	int temp = 0;

	switch (type)
	{

	// SCHWARZ_PRIMITIVE
	case 'P': {
		for (int z = 0; z < dimz; z++)
			for (int y = 0; y < dimy; y++)
				for (int x = 0; x < dimx; x++) {
					scal = (cosv[x] + cosv[y] + cosv[z]) - rvalue;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


	// GYROID
	case 'G': {
		for (int z = 0; z < dimz; z++)
			for (int y = 0; y < dimy; y++)
				for (int x = 0; x < dimx; x++) {
					scal = cosv[x] * senv[y] + cosv[y] * senv[z] + cosv[z] * senv[x] - rvalue;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


	// DIAMOND
	case 'D': {
		for (int z = 0; z < dimz; z++)
			for (int y = 0; y < dimy; y++)
				for (int x = 0; x < dimx; x++) {
					scal = senv[x] * senv[y] * senv[z] + senv[x] * cosv[y] * cosv[z] + cosv[x] * senv[y] * cosv[z] + cosv[x] * cosv[y] * senv[z] - rvalue;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


	// IWP
	case 'I': {
		for (int z = 0; z < dimz; z++)
			for (int y = 0; y < dimy; y++)
				for (int x = 0; x < dimx; x++) {
					scal = 2. * (cosv[x] * cosv[y] + cosv[y] * cosv[z] + cosv[z] * cosv[x]) - (cos2[x] + cos2[y] + cos2[z]) - rvalue;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


	// FISCHER_KOCH_S
	case 'K': {
		for (int z = 0; z < dimz; z++)
			for (int y = 0; y < dimy; y++)
				for (int x = 0; x < dimx; x++) {
					scal = cos2[x] * senv[y] * cosv[z] + cosv[x] * cos2[y] * senv[z] + senv[x] * cosv[y] * cos2[z] - rvalue;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


	// F_RD
	case 'F': {
		for (int z = 0; z < dimz; z++)
			for (int y = 0; y < dimy; y++)
				for (int x = 0; x < dimx; x++) {
					scal = -(4. * (cosv[x] * cosv[y] * cosv[z]) - (cos2[x] * cos2[y] + cos2[y] * cos2[z] + cos2[z] * cos2[x])) - rvalue;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;

	// SplitP
	case 'S': {
		for (int z = 0; z < dimz; z++)
			for (int y = 0; y < dimy; y++)
				for (int x = 0; x < dimx; x++) {
					scal = -(1.1*(sen2[x]*senv[z]*cosv[y] + sen2[y]*senv[x]*cosv[z] + sen2[z]*senv[y]*cosv[x])
							- 0.2*(cos2[x]*cos2[y] + cos2[y]*cos2[z] + cos2[z]*cos2[x]) 
							- 0.4*(cos2[x] + cos2[y] + cos2[z]) - rvalue);
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;

	// Lidinoid
	case 'L': {
		for (int z = 0; z < dimz; z++)
			for (int y = 0; y < dimy; y++)
				for (int x = 0; x < dimx; x++) {
					scal = sen2[x]*cosv[y]*senv[z] + sen2[y]*cosv[z]*senv[x] + sen2[z]*cosv[x]*senv[y] 
						- cos2[x]*cos2[y] - cos2[y]*cos2[z] - cos2[z]*cos2[x] + 0.3 - rvalue;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;
	
	}

	delete[] cosv;
	delete[] senv;
	delete[] cos2;
	delete[] sen2;

}
