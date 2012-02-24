#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstGamma1_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstGamma1.c
* \author       Bill Hill
* \date         November 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2012],
* The University Court of the University of Edinburgh,
* Old College, Edinburgh, UK.
* 
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be
* useful but WITHOUT ANY WARRANTY; without even the implied
* warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
* PURPOSE.  See the GNU General Public License for more
* details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
* Boston, MA  02110-1301, USA.
* \brief	Simple test for AlgGammaP() and AlgGammaLog() functions.
* \ingroup	binAlgTst
*/
#include <float.h>
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

int             main(int argc, char **argv)
{
  int           ok = 0;
  double        a = 0.0,
                g,
                gln,
                x;
  AlgError      algErr = ALG_ERR_NONE;

  switch(argc)
  {
    case 2:
      if((sscanf(*(argv + 1), "%lg", &x) == 1) && (x >= 0.0))
      {
        ok = 1;
      }
      break;
    case 3:
      if((sscanf(*(argv + 1), "%lg", &x) == 1) && (x >= 0.0) &&
         (sscanf(*(argv + 2), "%lg", &a) == 1) && (a >= DBL_EPSILON))
      {
        ok = 1;
      }
      break;
    default:
      ok = 0;
      break;
  }
  if(ok == 0)
  {
    (void )fprintf(stderr,
    "Usage: %s <x> [<a>]\n%s\n",
    *argv,
    "A test for AlgGammaP() and AlgGammaLog() which computes either the\n"
    "gamma function P(x) from a single value, or the incomplete gamma\n"
    "function P(a,x) from a pair of values.\n");
  }
  else
  {
    if(a >= DBL_EPSILON)
    {
      g = AlgGammaP(a, x, &algErr);
    }
    else
    {
      gln = AlgGammaLog(x, &algErr);
      g = exp(gln);
    }
    if(algErr == ALG_ERR_NONE)
    {
      printf("%g\n", g);
    }
    else
    {
      (void )fprintf(stderr, "%s: Failed to compute gamma function.\n", *argv);
    }
  }
  return(!ok);
}
