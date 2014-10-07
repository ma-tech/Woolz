#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzNearbyDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzNearbyDomain.c
* \author       Bill Hill
* \date         October 2014
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2014],
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
* \brief	Computes a domain which is newarby some some location
* 		in a reference domain.
* \ingroup	wlznearbydomain "WlzNearbyDomain"
*/

/*!
\ingroup BinWlz
\defgroup wlznearbydomain WlzNearbyDomain
\par Name
WlzNearbyDomain - computes the portion of a domain object which is nearby
                  given positions
\par Synopsis
\verbatim
WlzNearbyDomain [-d#] [-m#] [-p<pos>] [-o<out file>] [-h]
                [<Reference object file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-d</b></td>
    <td>Maximum distance from the given location(s).</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Distance function:
      <table width="500" border="0">
      <tr> <td>1</td> <td>octagonal (2D and 3D) - default</td></tr>
      <tr> <td>4</td> <td>4-connected (2D)</td></tr>
      <tr> <td>6</td> <td>6-connected (3D)</td></tr>
      <tr> <td>8</td> <td>8-connected (2D)</td></tr>
      <tr> <td>18</td> <td>18-connected (3D)</td></tr>
      <tr> <td>26</td> <td>26-connected (3D)</td></tr>
    </td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Locations given as x,y pairs or x,y,z triples for two and
        three dimensions respectively.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
</table>
\par Description
Computes a nearby domain in which all pixels/voxels are near to
the given locations within the reference domain.
If the given location(s) are not in the reference domain,
then an empty object will be output.
By default all files are read from the standard input and written to
the standard output.
\par Examples
\verbatim
WlzNearbyDomain -o out.wlz -d 10 -p 10,20,30;40,50,60 ref.wlz
\endverbatim
Creates a new domain object wrtponding to all voxels in the reference
spatial domain object (read from the file ref.wlz) which are less than
or equal to 10 voxels from the locations 10,20,30 and 40,50,60. This
object is written to the file out.wlz.
\par File
\ref WlzNearbyDomain.c "WlzNearbyDomain.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzDomainNearby "WlzDomainNearby(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok = 1,
  		option,
		con = 0,
		nPos = 0,
		usage = 0;
  char		*refFileStr,
		*outFileStr;
  double	dMax = 0.0;
  FILE		*fP = NULL;
  WlzObject	*outObj = NULL,
		*refObj = NULL;
  WlzVertexP	pos;
  WlzDistanceType dFn;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "hd:m:p:o:";
  const char    fileStrDef[] = "-";
  const WlzConnectType defDFn = WLZ_OCTAGONAL_DISTANCE;

  /* Parse the argument list and check for input files. */
  opterr = 0;
  dFn = defDFn;
  pos.v = NULL;
  outFileStr = (char *)fileStrDef;
  refFileStr = (char *)fileStrDef;
  while(ok && !usage && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'd':
	if((sscanf(optarg, "%lg", &dMax) != 1) || (dMax < 0.0))
	{
	  usage = 1;
	}
        break;
      case 'm':
	if(sscanf(optarg, "%d", &con) != 1)
	{
	  usage = 1;
	}
	else
	{
	  switch(con)
	  {
	    case 1:
	      dFn = WLZ_OCTAGONAL_DISTANCE;
	      break;
	    case 4:
	      dFn = WLZ_4_DISTANCE;
	      break;
	    case 6:
	      dFn = WLZ_6_DISTANCE;
	      break;
	    case 8:
	      dFn = WLZ_8_DISTANCE;
	      break;
	    case 18:
	      dFn = WLZ_18_DISTANCE;
	      break;
	    case 26:
	      dFn = WLZ_26_DISTANCE;
	      break;
	    default:
	      usage = 1;
	      break;
	  }
	}
	break;
      case 'p':
	if(optarg == NULL)
	{
	  usage = 1;
	}
	else
	{
	  /* Count the number of locations which will be number of location
	   * seperators + 1. */
	  {
	    char *s;

	    nPos = 1;
	    s = optarg;
	    while(*s)
	    {
	      if(*s++ == ';')
	      {
		++nPos;
	      }
	    }
	  }
	  if((pos.d3 = (WlzDVertex3 *)
	               AlcMalloc(sizeof(WlzDVertex3) * nPos)) == NULL)
	  {
	    (void )fprintf(stderr,
	        "%s: Failed to allocate location buffer.\n",
		argv[0]);
	    ok = 0;
	  }
	  if(ok)
	  {
	    /* Parse the locations using zero for any missing coordinate
	     * values. */
	    int		i;
	    char	*s;

	    i = 0;
	    s = strtok(optarg, ";");
	    while(s && !usage && (i < nPos))
	    {
	      WlzDVertex3 *p;

	      p = pos.d3 + i;
	      WLZ_VTX_3_ZERO(*p);
	      if(sscanf(s, "%lg,%lg,%lg",
	                &(p->vtX), &(p->vtY), &(p->vtZ)) < 2)
	      {
	        usage = 1;
	      }
	      else
	      {
	        ++i;
	        s = strtok(NULL, ";");
	      }
	    }
	    nPos = i;
	  }
	}
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  if(ok)
  {
    ok = !usage;
  }
  if(ok && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
      ok = 0;
    }
    else
    {
      refFileStr = argv[optind];
    }
  }
  /* Read the reference object. */
  if(ok)
  {
    if((refFileStr == NULL) ||
       (*refFileStr == '\0') ||
       ((fP = (strcmp(refFileStr, "-")?
              fopen(refFileStr, "r"): stdin)) == NULL) ||
       ((refObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read reference object from file %s.\n",
                     argv[0], refFileStr);
    }
    if(fP && strcmp(refFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Check object type and parse the promote the default distance function
   * if required. */
  if(ok)
  {
    switch(refObj->type)
    {
      case WLZ_2D_DOMAINOBJ:
	{
	  int i;
	  WlzDVertex3 p;

	  for(i = 0; i < nPos; ++i)
	  {
	    p = pos.d3[i];
	    WLZ_VTX_2_SET(pos.d2[i], p.vtX, p.vtY);
	  }
	}
	break;
      case WLZ_3D_DOMAINOBJ:
        break;
      default:
	ok = 0;
	(void )fprintf(stderr,
	    "%s: Reference object must be either a 2 or 3D spatial\n"
	    "domain object.\n",
	    argv[0]);
	break;
    }
  }
  /* Compute the distance transform object. */
  if(ok)
  {
    outObj = WlzAssignObject(
    	     WlzDomainNearby(refObj, nPos, pos, dFn, dMax, &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to compute nearby domain object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  /* Output the nearby domain object. */
  if(ok)
  {
    if((fP = (strcmp(outFileStr, "-")?
             fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, outObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
                     argv[0], errMsgStr);
    }
  }
  if(fP && strcmp(outFileStr, "-"))
  {
    (void )fclose(fP); fP = NULL;
  }
  (void )WlzFreeObj(refObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-d#] [-m#] [-p<pos>] [-o<output file>] [-h]\n"
    "\t\t[<Reference object file>]\n"
    "Computes a nearby domain in which all pixels/voxels are near to the\n"
    "given locations within the reference domain. If the given location(s)\n"
    "are not in the reference domain, then an empty object will be output.\n"
    "Version: %s\n"
    "Options:\n"
    "  -d  Maximum distance from the given location(s).\n"
    "  -m  Distance function:\n"
    "              1:  octagonal   (2D and 3D) - default\n"
    "              4:  4-connected (2D)\n"
    "              6:  6-connected (3D)\n"
    "              8:  8-connected (2D)\n"
    "             18: 18-connected (3D)\n"
    "             26: 26-connected (3D)\n"
    "  -p  Locations given as x,y pairs or x,y,z triples for two and\n"
    "      three dimensions respectively.\n"
    "  -o  Output object file.\n"
    "  -h  Help - prints this usage message\n"
    "By default all files are read from the standard input and written to\n"
    "the standard output.\n"
    "Example:\n"
    "  %s -o out.wlz -d 10 -p 10,20,30;40,50,60 ref.wlz\n"
    "which creates a domain corresponding to all voxels in the reference\n"
    "spatial domain object (read from the file ref.wlz) which are less than\n"
    "or equal to 10 voxels from the locations 10,20,30 and 40,50,60. This\n"
    "object is written to the file out.wlz\n",
    argv[0],
    WlzVersion(),
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
