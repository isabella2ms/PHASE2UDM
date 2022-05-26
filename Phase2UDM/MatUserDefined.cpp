#include "MatUserDefined.h"


MatUserDefined::MatUserDefined()
{
}


//! Default Destructor. 
MatUserDefined::~MatUserDefined()
{
}



MatUserDefined* MatUserDefined::MakeCopy() 
{
	MatUserDefined* newMat = new (MatUserDefined)(*this);
	return newMat;

}

int MatUserDefined::GetNumberOfStateVariables()
{return nStateVars;};

bool MatUserDefined::GetIsTangential()
{return isTangential;};

void MatUserDefined::SetNumberOfStateVariables(int n)
{nStateVars = n;}

void MatUserDefined::SetIsTangential(bool b)
{isTangential = b;}
