#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _AlgTstMixtureMLG1_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         AlgTstMixtureMLG1.c
* \author       Bill Hill
* \date         November 2010
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2010 Medical research Council, UK.
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
* \brief	Simple test for AlgMixtureMLG().
* \ingroup	binAlgTst
*/
#include <stdio.h>
#include <unistd.h>
#include <Alc.h>
#include <Alg.h>

int             main(int argc, char **argv)
{
  int           idx,
                ok,
                done,
                ioCnt,
                iMax = 1000,
                nCls,
                nDbn,
                nItn;
  double        dataOrg,
                dataItv,
                tol = 0.001,
                logL;
  FILE          *fP = NULL;
  AlgError      algErr = ALG_ERR_NONE;
  static char   ioBuf[256];
  static double freq[256],
                alpha[4],
                mu[4],
                sigma[4];

  ok = 1;
  /* Parse command line. */
  if(argc != 2)
  {
    (void )fprintf(stderr,
                   "Usage: %s <file>\n",
                   *argv);
    ok = 0;
  }
  if(ok)
  {
    if((fP = fopen(*(argv + 1), "r")) == NULL)
    {
      (void )fprintf(stderr,
                     "%s: Failed to open input file %s\n",
                     *argv, *(argv + 1));
      ok = 0;
    }
  }
  if(ok)
  {
    /* Read initial estimates of alpha, mu and sigma until blank line. */
    idx = 0;
    done = 0;
    while(!done && (idx < 256) && (fgets(ioBuf, 256, fP) != NULL))
    {
      ioBuf[255] = '\0';
      if(sscanf(ioBuf, "%lg %lg %lg\n",
                alpha + idx, mu + idx, sigma + idx) != 3)
      {
        done = 1;
        nDbn = idx;
      }
      else
      {
        ++idx;
      }
    }
    if(done == 0)
    {
      (void )fprintf(stderr,
                     "%s: Failed to read input distribution parameters.\n",
                     *argv);
      ok = 0;
    }
  }
  if(ok)
  {
    /* Read data origin and step size followed by blank line. */
    if((fgets(ioBuf, 256, fP) == NULL) ||
       (sscanf(ioBuf, "%lg %lg", &dataOrg, &dataItv) != 2))
    {
      (void )fprintf(stderr,
                     "%s: Failed to read data origin and interval.\n",
                     *argv);
      ok = 0;
    }
    else
    {
      /* Consume blank line. */
      (void )fgets(ioBuf, 256, fP);
    }
  }
  if(ok)
  {
    /* Read data f(x) until end of file */
    idx = 0;
    while((fgets(ioBuf, 256, fP) != NULL) &&
          (sscanf(ioBuf, "%lg", freq + idx) == 1))
    {
      ++idx;
    }
    nCls = idx;
  }
  if(ok)
  {
    /* Compute ML Gaussian mixture */
    algErr = AlgMixtureMLG(nDbn, nCls, dataOrg, dataItv, freq,
                           alpha, mu, sigma,
                           tol, iMax, &logL, &nItn);
    if(algErr != ALG_ERR_NONE)
    {
      (void )fprintf(stderr,
                     "%s: Failed to compute mixture of Gaussians (%d).\n",
                     *argv, (int )algErr);
      ok = 0;
    }
  }
  if(ok)
  {
    /* Output new estimate of Gaussian mu, sigma and alpha */
    for(idx = 0; idx < nDbn; ++idx)
    {
      (void )printf("%8g %8g %8g %d %8g\n",
                    alpha[idx], mu[idx], sigma[idx], nItn, logL);

    }
  }
  if(fP)
  {
    (void )fclose(fP);
  }
  exit(!ok);
}

