//----------------------------------------------------------------------
//
// Rocscience.  All rights reserved.
//
//----------------------------------------------------------------------

#pragma once
#include "matrix.h"

class MatUserDefined
{
public: 
	 //! Default Constructor. 
	MatUserDefined(void);

	//! Default Destructor. 
	~MatUserDefined(void); 

//////Constitutive model

	virtual void SetConstitutiveParameters(Vector& parameters){};

	virtual Matrix ComputeElasticMatrix(Vector stress, Vector strain, Vector& stateVariables, int& failureType){Matrix dummy(0,0); return dummy;};

	virtual Matrix ComputePlasticMatrix(Vector stress, Vector strain, Vector& stateVariables, int& failureType){Matrix dummy(0,0); return dummy;};

	virtual void UpdateStress(Vector& stress, Vector& strain, Vector& stateVariables, Vector& dStrain, int& failureType, Vector GaussPointInfo) {};

///Properties needed for different types of simulation (consolidation, dynamic, etc..)
	
	virtual double GetElasticModulus() { return -1.; }
	virtual double GetBulkModulus() { return -1; }
	virtual double GetShearModulus() { return -1; }
	virtual double GetPoissonRatio() { return -1; }
	virtual double GetElasticModulus(Vector& stress, Vector& strain, Vector& stateVariables, Vector GaussPointInfo) { return GetElasticModulus(); }
	virtual double GetBulkModulus(Vector& stress, Vector& strain, Vector& stateVariables, Vector GaussPointInfo) { return GetBulkModulus(); }
	virtual double GetShearModulus(Vector& stress, Vector& strain, Vector& stateVariables, Vector GaussPointInfo) { return GetShearModulus(); }
	virtual double GetPoissonRatio(Vector& stress, Vector& strain, Vector& stateVariables, Vector GaussPointInfo) { return GetPoissonRatio(); }

/////////////////////////////////////////////////////////////////////////////////////////////////////////
	virtual int GetNumberOfStateVariables();
	
	virtual bool GetIsTangential();

	virtual void SetNumberOfStateVariables(int n);

	virtual void SetIsTangential(bool b);
	
	virtual MatUserDefined* MakeCopy();
	virtual void GetJointProp_Strength_UDM(double& c1, double& p1, double& t1, double& sigmaNN, Vector& stress, Vector& strain, Vector& stateVariables, Vector GaussPointInfo) { return; }

protected:

	bool isTangential;

	int nStateVars;
}; 

