#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _AlgDPSearch_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         libAlg/AlgDPSearch.c
* \author       Richard Baldock
* \date         August 2005
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
* \brief        A 1D dynamic programming search procedure
*		assuming a rectangular search region and a given non-local
*		cost function.
* \ingroup      AlgDPSearch
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <Alc.h>

/*! 
* \ingroup      AlgDPSearch
* \brief        Use dynamic programming to establish an optima path given local
 and non-local costs.
*
* \return       zero
* \param    imax	number of points on the path
* \param    jmax	number of locations per path point
* \param    local_cost	local cost for each path point
* \param    optimal_cost	return for optimal path cost through each point
* \param    optimal_path	return for optimal path through each point.
* \param    non_local_cost	non-local cost function calculated in terms
 of current and previous position and current optimal paths.
* \par      Source:
*                AlgDPSearch.c
*/
int AlgDPSearch(
  int 		imax,
  int		jmax,
  double	**local_cost,
  double	**optimal_cost,
  int 		**optimal_path,
  double 	(*non_local_cost)(int, int, int, int **))
{
  /* local variables */
  int 		i, j, jp;
  double	cost, min_cost;

  /* initialise start cost and select initial
     optimal path to minimise the non-local cost */
  for(j=0; j < jmax; j++)
  {
    min_cost = cost = (*non_local_cost)(0,j,0,optimal_path);
    optimal_path[0][j] = 0;
    for(jp=1; jp < jmax; jp++)
    {
      cost = (*non_local_cost)(0,j,jp,optimal_path);
      if( cost < min_cost )
      {
	min_cost = cost;
	optimal_path[0][j] = jp;
      }
    }
    optimal_cost[0][j] = min_cost + local_cost[0][j];
  }

  /* determine the optimal paths */
  for(i=1; i < imax; i++)
  {
    for(j=0; j < jmax; j++)
    {
      min_cost = cost = optimal_cost[i-1][0]
	+ (*non_local_cost)(i,j,0,optimal_path);
      optimal_path[i][j] = 0;
      for(jp=1; jp < jmax; jp++)
      {
	cost = optimal_cost[i-1][jp]
	  + (*non_local_cost)(i,j,jp,optimal_path);
	if( cost < min_cost )
	{
	  min_cost = cost;
	  optimal_path[i][j] = jp;
	}
      }
      optimal_cost[i][j] = min_cost + local_cost[i][j];
    }
  }

  return( 0 );
}

/*!
* \return	zero
* \ingroup	AlgDPSearch
* \brief
* \param    imax	number of points on the path
* \param    jmax	number of locations per path point
* \param    optimal_cost	return for optimal path cost through each point
* \param    optimal_path	return for optimal path through each point.
* \param    non_local_cost	non-local cost function calculated in terms
*/
int AlgDPTotalCosts(
  int 		imax,
  int		jmax,
  double	**optimal_cost,
  int 		**optimal_path,
  double 	(*non_local_cost)(int, int, int, int **))
{
  int		i, j, jp;
  double	cost, min_cost, *tmp;

  /* now determine the total optimal-costs for each point */
  tmp = (double *) AlcMalloc(sizeof(double) * jmax);
  for(i=imax-1; i > 0; i--)
  {
    for(j=0; j < jmax; j++)
    {
      cost = optimal_cost[i][0]
	- optimal_cost[i-1][optimal_path[i][0]]
	+ (*non_local_cost)(i,0,j,optimal_path)
	- (*non_local_cost)(i,0,optimal_path[i][0],
			    optimal_path);
      min_cost = cost;
      for(jp=1; jp < jmax; jp++)
      {
	cost = optimal_cost[i][jp]
	  - optimal_cost[i-1][optimal_path[i][jp]]
	  + (*non_local_cost)(i,jp,j,optimal_path)
	  - (*non_local_cost)(i,jp,
			      optimal_path[i][jp],
			      optimal_path);
	if( cost < min_cost )
	{
	  min_cost = cost;
	}
      }
      tmp[j] = min_cost;
    }
    for(j=0; j < jmax; j++)
    {
      optimal_cost[i-1][j] += tmp[j];
    }
  }

  AlcFree( tmp );
  return( 0 );
}
