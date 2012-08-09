#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzPolarSample_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzPolarSample.c
* \author       Bill Hill
* \date         March 1999
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
* \brief	Polar resampling of domain objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzpolarsample "WlzPolarSample"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzpolarsample WlzPolarSample
\par Name
WlzPolarSample - Applies a polar resampling filter.
\par Synopsis
\verbatim
WlzPolarSample [-o<output object file>] [-i] [-h]
           [<input object file>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-i</b></td>
    <td>inner bounding circle of objects domain used. </td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>
By default the input object is read from the standard input, re-sampled
using the outer bounding circle of the object's domain and  written  to
the standard output.
\par Description
 Polar resamples the input 2D domain object.

\par Examples
\verbatim
       # An example of using WlzPolarSample to resample an object.

       WlzPolarSample -o outfile.wlz infile.wlz

\endverbatim

\par File
\ref WlzPolarSample.c "WlzPolarSample.c"
\par See Also
\ref WlzPolarSample "WlzPolarSample(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             main(int argc, char **argv)
{
  int		option,
		angleSamples,
		outFlag = 1,
		ok = 1,
		usage = 0;
  double	angleInc;
  WlzIVertex2	center,
  		objSz;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  WlzObject	*inObj = NULL,
		*outObj = NULL;
  char 		*outObjFileStr,
  		*inObjFileStr;
  const char	*errMsg;
  static char	optList[] = "o:ih",
		outObjFileStrDef[] = "-",
  		inObjFileStrDef[] = "-";

  opterr = 0;
  outObjFileStr = outObjFileStrDef;
  inObjFileStr = inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'i':
        outFlag = 0;
	break;
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
     (outObjFileStr == NULL) || (*outObjFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
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
      inObjFileStr = *(argv + optind);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inObjFileStr == NULL) ||
       (*inObjFileStr == '\0') ||
       ((fP = (strcmp(inObjFileStr, "-")?
	      fopen(inObjFileStr, "r"): stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s (%s).\n",
		     *argv, inObjFileStr, errMsg);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(ok)
  {
    if(inObj->type == WLZ_EMPTY_OBJ)
    {
      center.vtX = 0;
      center.vtY = 0;
      angleInc = 1.0;
      angleSamples = 1;
    }
    else if(inObj->type == WLZ_2D_DOMAINOBJ)
    {
      if(inObj->domain.core == NULL)
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: input object read from %s has null domain\n",
		       *argv, inObjFileStr);
      }
      else if((inObj->domain.core->type != WLZ_INTERVALDOMAIN_INTVL) &&
      	      (inObj->domain.core->type != WLZ_INTERVALDOMAIN_RECT))
      {
        ok = 0;
	(void )fprintf(stderr,
		       "%s: input object read from %s has invalid domain\n",
		       *argv, inObjFileStr);
      }
      else
      {
        center.vtX = (inObj->domain.i->kol1 + inObj->domain.i->lastkl) / 2;
        center.vtY = (inObj->domain.i->line1 + inObj->domain.i->lastln) / 2;
	objSz.vtX = inObj->domain.i->lastkl - inObj->domain.i->kol1 + 1;
	objSz.vtY = inObj->domain.i->lastln - inObj->domain.i->line1 + 1;
	if(outFlag)
	{
	  angleSamples = WLZ_MAX(objSz.vtX, objSz.vtY) / 2;
	}
	else
	{
	  angleSamples = WLZ_MIN(objSz.vtX, objSz.vtY) / 2;
	}
	if(angleSamples > 0)
	{
	  angleInc = (2.0 * WLZ_M_PI) / angleSamples;
	}
	else
	{
	  angleInc = 1.0;
	}
      }
    }
    else
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: input object read from %s has invalid type\n",
		     *argv, inObjFileStr);
    }
  }
  if(ok)
  {
    if((outObj = WlzPolarSample(inObj, center, angleInc, 1,
    				angleSamples, outFlag, &errNum)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to polar resample object (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outObjFileStr, "-")?
              fopen(outObjFileStr, "w"):
	      stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      		     "%s: failed to write output object (%s)\n",
		     *argv, errMsg);
    }
    if(fP && strcmp(outObjFileStr, "-"))
    {
      fclose(fP);
    }
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-o<out object>] [-i] [-h] [<in object>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -o  Output object file name.\n"
    "  -i  Resampling circle inside object, default is to use the outer\n"
    "      bounding circle.\n"
    "  -h  Help, prints this usage message.\n"
    "Polar resamples the input Woolz object.\n"
    "Objects are read from stdin and written to stdout unless the filenames\n"
    "are given.\n",
    *argv,
    " -o polar.wlz rect.wlz\n"
    "The input Woolz object is read from rect.wlz, resampled and written\n"
    "to polar.wlz.\n");
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
