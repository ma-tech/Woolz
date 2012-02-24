#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgTstConvolve1_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binAlgTst/AlgTstConvolve1.c
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
* \brief	Simple test for AlgConvolve().
* \ingroup	binAlgTst
*/
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

int             main(int argc, char *argv[])
{
  int           idx;
  AlgPadType    pad  = ALG_PAD_NONE;
  AlgError      errCode = ALG_ERR_NONE;
  static double dat[25],
                krn[5] = {0.50, 0.87, 1.00, 0.87, 0.50},
                cnv[25];
  const int     datSz = 25,
                krnSz = 5;

  for(idx = 0; idx < krnSz; ++idx)
  {
    krn[idx] /= 3.74;
  }
  for(idx = 0; idx < datSz; ++idx)
  {
    dat[idx] = 1.0;
  } 
  dat[0] = 2.0;
  dat[13] = 2.0;
  dat[24] = 2.0;
  errCode = AlgConvolve(25, cnv, 5, krn, 25, dat, pad);
  if(errCode == ALG_ERR_NONE)
  {
    for(idx = 0; idx < datSz; ++idx)
    {
      (void )printf("% 8g    % 8g\n", dat[idx], cnv[idx]);
    }   
  }     
  else  
  {     
    (void )printf("AlgConvolve() returned error code %d\n", (int )errCode);
  }     
  return((int )errCode);
}

