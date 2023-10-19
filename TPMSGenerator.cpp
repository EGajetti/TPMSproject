#include <cmath>
#include "Utils.h"

using namespace std;


// Implementation of the function defined in Utils.h

void TpmsGenerator(const int npoints, const int numcell, char type, const float rstart, vtkImageData* volume)
{

	int dimension = npoints * numcell;

	const float pi = 2 * 3.14159265358979323846 / npoints;
	float scal = 0.;

	float* cosv = new float[dimension];
	for (int i = 0; i < dimension; i++)
		cosv[i] = cos(pi * i);

	float* senv = new float[dimension];
	for (int i = 0; i < dimension; i++)
		senv[i] = sin(pi * i);

	float* cosv_t = new float[dimension];
	if ((type == 'I') || (type == 'S') || (type == 'F'))
		for (int i = 0; i < dimension; i++)
			cosv_t[i] = cos(2 * pi * i);

	int temp = 0;

	switch (type)
	{

		// SCHWARZ_PRIMITIVE
	case 'P': {
		for (int z = 0; z < dimension; z++)
			for (int y = 0; y < dimension; y++)
				for (int x = 0; x < dimension; x++) {
					scal = -(cosv[x] + cosv[y] + cosv[z]) - rstart;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


		// SHOEN_GYROID
	case 'G': {
		for (int z = 0; z < dimension; z++)
			for (int y = 0; y < dimension; y++)
				for (int x = 0; x < dimension; x++) {
					scal = cosv[x] * senv[y] + cosv[y] * senv[z] + cosv[z] * senv[x] - rstart;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


		// SCHWARZ_DIAMOND
	case 'D': {
		for (int z = 0; z < dimension; z++)
			for (int y = 0; y < dimension; y++)
				for (int x = 0; x < dimension; x++) {
					scal = senv[x] * senv[y] * senv[z] + senv[x] * cosv[y] * cosv[z] + cosv[x] * senv[y] * cosv[z] + cosv[x] * cosv[y] * senv[z] - rstart;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


		// SHOEN_IWP
	case 'I': {
		for (int z = 0; z < dimension; z++)
			for (int y = 0; y < dimension; y++)
				for (int x = 0; x < dimension; x++) {
					scal = 2. * (cosv[x] * cosv[y] + cosv[y] * cosv[z] + cosv[z] * cosv[x]) - (cosv_t[x] + cosv_t[y] + cosv_t[z]) - rstart;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


		// FISCHER_KOCH_S
	case 'S': {
		for (int z = 0; z < dimension; z++)
			for (int y = 0; y < dimension; y++)
				for (int x = 0; x < dimension; x++) {
					scal = cosv_t[x] * senv[y] * cosv[z] + cosv[x] * cosv_t[y] * senv[z] + senv[x] * cosv[y] * cosv_t[z] - rstart;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


		// F_RD
	case 'F': {
		for (int z = 0; z < dimension; z++)
			for (int y = 0; y < dimension; y++)
				for (int x = 0; x < dimension; x++) {
					scal = -(4. * (cosv[x] * cosv[y] * cosv[z]) - (cosv_t[x] * cosv_t[y] + cosv_t[y] * cosv_t[z] + cosv_t[z] * cosv_t[x])) - rstart;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;


		// F_RD
	case 'N': {
		for (int z = 0; z < dimension; z++)
			for (int y = 0; y < dimension; y++)
				for (int x = 0; x < dimension; x++) {
					scal = 3. * (cosv[x] + cosv[y] + cosv[z]) + 4. * cosv[x] * cosv[y] * cosv[z] - rstart;
					float* a = static_cast<float*> (volume->GetScalarPointer(x, y, z));
					*a = scal;
					temp++;
				}
	} break;

	}

	delete[] cosv;
	delete[] senv;
	delete[] cosv_t;

}
