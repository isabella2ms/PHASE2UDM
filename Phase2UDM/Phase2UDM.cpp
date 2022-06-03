
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

	_dCK = 0.0;
	_dn = 0.0;
	_dNu = 0.0;
	_dRf = 0.0;
	_dc = 0.0;
	_dfi = 0.0;
	_dE = 0.0;
	_de = 0.0;
	_dG = 0.0;
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

	_dE = _dCK * pow(stress[2], _dn);
	_dG = _dE / 2*(1 - _dNu);
	double num, den;

	num = _dRf * (stress[0] - stress[2])*(1 - sin(_dfi));
	den = 2 * _dc*cos(_dfi) + 2 * stress[2] * sin(_dfi);
	_de = _dE * pow((1 - num / den),2);

	double dE = _dE, dNu = _dNu, dG = _dG, de =_de;
	double c1 = 1/dG;
	double c2 = -dNu/de;
	double c3 = 1/de;

	Emat(0, 0) = c3;
	Emat(0, 1) = c2;
	Emat(0, 2) = c2;
	Emat(1, 0) = c2;
	Emat(1, 1) = c3;
	Emat(1, 2) = c2;
	Emat(2, 0) = c2;
	Emat(2, 1) = c2;
	Emat(2, 2) = c3;
	Emat(3, 3) = c1;
	Emat(4, 4) = c1;
	Emat(5, 5) = c1;

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
	_dCK = parameters[0];
	_dn = parameters[1];
	_dNu = parameters[2];
	_dRf = parameters[3];
	_dc = parameters[4];
	_dfi = parameters[5];

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

