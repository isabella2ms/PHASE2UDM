//----------------------------------------------------------------------
//
// Rocscience.  All rights reserved.
//
//----------------------------------------------------------------------
#ifndef _REALTYPE_H
#define _REALTYPE_H

#include <float.h>

#ifdef _OPENMP
#include <omp.h>
#endif

typedef double Real;  
#define REAL_MAX DBL_MAX
#define REAL_MIN DBL_MIN
#define REAL_EPSILON DBL_EPSILON

typedef ptrdiff_t isize_t;  // assumed signed equivalent of size_t? 
                            // LONG_PTR, INT_PTR, windef.h::SSIZE_T?

////////////////////////////////////////////////////////////////////////////
#pragma warning (disable:4100) //  unreferenced formal parameter
#pragma warning (disable:4514) //  unreferenced inline function has been removed
#pragma warning (disable:4710) //  function 'ABC' not inlined
#pragma warning (disable:4702) //  unreachable code
////////////////////////////////////////////////////////////////////////////


#endif
