
// Phase2UDM.h : Defines the exported functions for the DLL application.
//

#include "pch.h"
#include "Phase2UDM.h"
#include <iostream>
#include <fstream>
#include <windows.h>

// This is the constructor of a class that has been exported.
// see RS3DLL.h for the class definition

//Constitutive functions
Elastic::Elastic() : MatUserDefined()
{
	isTangential = false;
	nStateVars = 1;

	_dE = 0.0;
	_dNu = 0.0;
	_dKappa = 0.0;
	_dLambda = 0.0;
	_dM = 0.0;
}


Elastic::~Elastic(void)
{

}

Matrix Elastic::ComputeElasticMatrix(Vector stress, Vector strain, Vector& stateVariables, int& failureType)
{
	Matrix Emat(6, 6);
	Emat = 0.;

	assert(abs((1 + _dNu) * (1 - 2 * _dNu)) > 10e-15);
	assert(abs(1 - _dNu) > 10e-15);
	assert(abs(2 * (1 - _dNu)) > 10e-15);

	double dE = _dE, dNu = _dNu;
	double c44 = dE / (2.0 + 2.0 * dNu);
	double c12 = dE * dNu / ((1.0 + dNu) * (1.0 - 2.0 * dNu));
	double c11 = c12 + 2.0 * c44;

	Emat(0, 0) = c11;
	Emat(0, 1) = c12;
	Emat(0, 2) = c12;
	Emat(1, 0) = c12;
	Emat(1, 1) = c11;
	Emat(1, 2) = c12;
	Emat(2, 0) = c12;
	Emat(2, 1) = c12;
	Emat(2, 2) = c11;
	Emat(3, 3) = c44;
	Emat(4, 4) = c44;
	Emat(5, 5) = c44;

	return Emat;
}

Matrix Elastic::ComputePlasticMatrix(Vector stress, Vector strain, Vector& stateVariables, int& failureType)
{
	Matrix Emat = ComputeElasticMatrix(stress, strain, stateVariables, failureType);

	return Emat;
}

void Elastic::UpdateStress(Vector& stress, Vector& strain, Vector& stateVariables, Vector& dStrain, int& failureType, Vector GaussPointInfo)

{
	if (stateVariables.size() == 0) stateVariables.resize(1);
	Matrix D = ComputeElasticMatrix(stress, strain, stateVariables, failureType);
	//size_t sz = dStrain.size();
	Vector ds = D * dStrain;
	stress += ds;
}



////functions needed to integrate with Phase2


void Elastic::SetConstitutiveParameters(Vector& parameters)
{
	_dE = parameters[0];
	_dNu = parameters[1];
	_dKappa = parameters[2];
	_dLambda = parameters[3];
	_dM = parameters[4];
}

MatUserDefined* Elastic::MakeCopy()
{
	MatUserDefined* newMat = new (Elastic)(*this);
	return newMat;

}


///Properties needed for different types of simulation (consolidation, dynamic, etc..)

double Elastic::GetElasticModulus() { return _dE; }

double Elastic::GetBulkModulus() { return _dE / (3.0 * (1.0 - 2.0 * _dNu)); }

double Elastic::GetShearModulus() { return _dE / (2.0 * (1.0 + _dNu)); }

double Elastic::GetPoissonRatio() { return _dNu; }

