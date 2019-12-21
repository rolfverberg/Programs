#if !defined(CONVERTFVTOTECPLOT_H__INCLUDED_)
#define CONVERTFVTOTECPLOT_H__INCLUDED_

#ifndef FALSE
#define FALSE	0
#endif
#ifndef TRUE
#define TRUE	1
#endif
#ifndef NULL
#define NULL	0
#endif

#ifdef _DEBUG
#ifndef DEBUG
#define DEBUG	1
#endif
#ifndef ALLOCDEBUG
#define ALLOCDEBUG	1
#endif
#endif

#if !defined(__EXTERN_VARIABLES__)
#define __EXTERN_VARIABLES__

#endif // !defined(__EXTERN_VARIABLES__)

#include "../../Source/clFieldview.h"
#include "../../Source/clTecplot.h"

#endif // CONVERTFVTOTECPLOT_H__INCLUDED_
