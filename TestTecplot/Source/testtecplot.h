#ifndef TESTTECPLOT_H__INCLUDED_
#define TESTTECPLOT_H__INCLUDED_

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

#ifndef __EXTERN_VARIABLES__
#define __EXTERN_VARIABLES__

#endif // !defined(__EXTERN_VARIABLES__)

#include <stdio.h>
#include "../../Source/clTecplot.h"
#include "../../Source/clString.h"

enum constants {NAME_LENGTH=128,LINE_LENGTH=256};

#endif // TESTTECPLOT_H__INCLUDED_
