#pragma once
#include "MatUserDefined.h"
#include "Matrix.h"

#include <iostream>
#include <string>
#include <vector>


class Elastic : public MatUserDefined {
public:

	//! Default Constructor. 
	Elastic(void);

	//! Default Destructor. 
	~Elastic(void);

	//////Constitutive model

	virtual void SetConstitutiveParameters(Vector& parameters);

	virtual Matrix ComputeElasticMatrix(Vector stress, Vector strain, Vector& stateVariables, int& failureType);

	virtual Matrix ComputePlasticMatrix(Vector stress, Vector strain, Vector& stateVariables, int& failureType);

	virtual void UpdateStress(Vector& stress, Vector& strain, Vector& stateVariables, Vector& dStrain, int& failureType, Vector GaussPointInfo);

	virtual MatUserDefined* MakeCopy();

	virtual double GetElasticModulus();
	virtual double GetBulkModulus();
	virtual double GetShearModulus();
	virtual double GetPoissonRatio();

private:
	//! Young's Modulus of Elasticity of the material.
	double _dE;

	//! Poisson's ratio of the material.
	double _dNu;

	//! Coeficiente Kappa
	double _dKappa;

	//! Coeficiente Lambda
	double _dLambda;

	//! Inclinação M
	double _dM;

};

extern "C" __declspec(dllexport) MatUserDefined * GetMaterial(char* names)
{
	MatUserDefined* mat = NULL;
	if (strcmp(names, "Elastic2") == 0) {
		mat = new Elastic();
	}
	return mat;
}


extern "C" __declspec(dllexport) int GetNumberOfMaterial()
{
	return 1;
}

extern "C" __declspec(dllexport) char* GetMaterialName(int id)
{
	char* cmatName;
	std::string matName;

	if (id == 0) { matName = "Elastic2"; }
	else {
		matName = "No material was found";
	}
	cmatName = const_cast<char*>(matName.c_str());
	return cmatName;
}


extern "C" __declspec(dllexport) int GetNumberOfConstitutiveParameter(char* matName)
{
	int nVars = 0;
	if (strcmp(matName, "Elastic2") == 0) { nVars = 5; }
	return nVars;
}

extern "C" __declspec(dllexport) char* GetNameOfConstitutiveParameter(char* matName, int paramId)
{
	char* cparamName;
	std::string paramName;
	
	if (strcmp(matName, "Elastic2") == 0) {
		switch (paramId) {
		case 0:
			paramName = "Young";
			break;
		case 1:
			paramName = "Poisson";
			break;
		case 2:
			paramName = "Kappa";
			break;
		case 3:
			paramName = "Lambda";
			break;
		case 4:
			paramName = "M";
			break;
		default:
			paramName = "No param found";
			break;

		}
	}
	else {
		paramName = "No material was found";
	}

	//cparamName = const_cast<char*>(paramName.c_str());

	cparamName = new char[paramName.size() + 1];
	std::copy(paramName.begin(), paramName.end(), cparamName);
	cparamName[paramName.size()] = '\0';

	return cparamName;
}

