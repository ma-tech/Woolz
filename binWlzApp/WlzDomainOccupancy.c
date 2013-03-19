#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzDomainOccupancy_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzDomainOccupancy.c
* \author       Bill Hill
* \date         November 2012
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
* \brief	Computes the occupancy of a domain with respect to it's
* 		planes.
* \ingroup	BinWlzApp
*/

/*!
\ingroup BinWlzApp
\defgroup wlzdomainoccupancy WlzDomainOccupancy
\par Name
WlzDomainOccupancy - computes the occupancy of a domain with respect to it's
   	             planes.
\par Synopsis
\verbatim
WlzDomainOccupancy [-h] [-a] [-o<out>] [-b<sec bib>] [-c] [-p<sep>]
                   [-s<sec trm>] [-x] [-y] [-z] [<in obj>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr>
    <td><b>-a</b></td>
    <td>Output plane occupancy area or number of domains rather than the
        default binary values (ie 0 or 1).</td>
  </tr>
  <tr>
    <td><b>-b</b></td>
    <td>Planes defined by an MAPaint style bib file.</td>
  </tr>
  <tr>
    <td><b>-c</b></td>
    <td>Output plane occupancy numbers using comma rather than white space
        separation.</td>
  </tr>
  <tr>
    <td><b>-p</b></td>
    <td>Plane separation distance (default 1.0).</td>
  </tr>
  </tr>
  <tr>
    <td><b>-r</b></td>
    <td>Output plane occupany as a single row of plane coordinates.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Planes defined by a Woolz section transform object.</td>
  </tr>
  <tr>
    <td><b>-x</b></td>
    <td>Planes are columns in the domain.</td>
  </tr>
  <tr>
    <td><b>-y</b></td>
    <td>Planes are lines in the domain.</td>
  </tr>
  <tr>
    <td><b>-z</b></td>
    <td>Planes are object planes in the domain.</td>
  </tr>
</table>
\par Description
Computes the occupancy of a domain or compound array object with
respect to the planes through it.
By default these planes are simply the x-y planes orthogonal to
the z-axis but a different set of planes may be specified by: an
MAPaint style bib file, a Woolz section transform object or by
selecting one of the orthogonal axis plane propagation directions
(x, y or z). The occupancy is output as two columns, the first
containing the plane coordinate (with respect to a section transform)
and the second the occupancy. The occupancy may bean area, the number
of domains or a binary value (ie 0 or 1).
By default the domain is read from the standard input and the occupancy
is written to the standard output.
\par Examples
\verbatim
WlzDomainOccupancy -c heart.wlz
\endverbatim
Outputs a space seperated list of numbers with two columns.
The first column containing the (z) plane coordinate and
the second the occupancy of the domain (read from heart.wlz) for that plane.
\par File
\ref WlzDomainOccupancy.c "WlzDomainOccupancy.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref Wlz3DSectionOcc "Wlz3DSectionOcc(3)"
\ref Wlz3DViewTransformObj "Wlz3DViewTransformObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <unistd.h>
#include <float.h>
#include <Wlz.h>
#include <bibFile.h>

static WlzThreeDViewStruct 	*WlzDomainOccReadBib(
				  char *fStr,
				  WlzErrorNum *dstErr);

int             main(int argc, char **argv)
{
  int		option,
		nOcc = 0,
  		ok = 1,
  		usage = 0,
		useArea = 0,
		useBib = 0,
		useCommas = 0,
		useRow = 0,
		useSec = 0,
		useXYZ = 0;
  char		*bibFileStr,
		*secFileStr,
		*inDomStr,
		*outFileStr;
  double	first = 0.0,
  		last = 0.0,
		sep = 1.0;
  int		*occ = NULL;
  WlzObject	*inObj = NULL;
  WlzThreeDViewStruct *vs = NULL;
  const char    *errMsg = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char	inFileStrDef[] = "-",
  		outFileStrDef[] = "-",
		optList[] = "achrxyzb:o:p:s:";

  opterr = 0;
  bibFileStr = outFileStr = secFileStr = inDomStr = inFileStrDef;
  outFileStr = outFileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
	useArea = 1;
	useRow = 0;
        break;
      case 'b':
	useBib = 1;
	useSec = 0;
	useXYZ = 0;
	bibFileStr = optarg;
        break;
      case 'c':
	useCommas = 1;
        break;
      case 'h':
	usage = 1;
        break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'p':
        if(sscanf(optarg, "%lg", &sep) != 1)
	{
	  usage = 1;
	}
	break;
      case 'r':
	useArea = 0;
	useRow = 1;
        break;
      case 's':
	useBib = 0;
	useXYZ = 0;
	useSec = 1;
	secFileStr = optarg;
        break;
      case 'x': /* FALLTHROUGH */
      case 'y': /* FALLTHROUGH */
      case 'z':
	useBib = 0;
	useSec = 0;
	useXYZ = option;
        break;
      default:
        usage = 1;
	break;
    }
  }
  if(usage == 0)
  {
    if(optind + 1 == argc )
    {
      inDomStr = *(argv + optind);
    }
    else if(optind == argc)
    {
      inDomStr = inFileStrDef;
    }
    else
    {
      usage = 1;
    }
  }
  ok = usage == 0;
  if(ok)
  {
    FILE  	*fP = NULL;

    if((fP = (strcmp(inDomStr, "-")? fopen(inDomStr, "r"): stdin)) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
    else
    {
      inObj = WlzAssignObject(
	      WlzReadObj(fP, &errNum), NULL);
      if(strcmp(inDomStr, "-"))
      {
	(void)fclose(fP);
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      "%s: Failed to read input domain object from file %s (%s).\n",
      *argv, inDomStr, errMsg);
    }
  }
  if(ok)
  {
    if(useBib)
    {
      vs = WlzAssign3DViewStruct(
           WlzDomainOccReadBib(bibFileStr, &errNum), NULL);
    }
    else if(useSec)
    {
      WlzObject *vsObj = NULL;
    
      if(strcmp(secFileStr, "-"))
      {
	vsObj = WlzAssignObject(
		WlzReadObj(stdin, &errNum), NULL);
      }
      else
      {
	FILE    *fP = NULL;

	if((fP = fopen(secFileStr, "r")) == NULL)
	{
	  errNum = WLZ_ERR_READ_EOF;
	}
	else
	{
	  vsObj = WlzAssignObject(
		  WlzReadObj(fP, &errNum), NULL);
	  (void )fclose(fP);
	}
      }
      if(vsObj == NULL)
      {
	errNum = WLZ_ERR_READ_EOF;
      }
      else if(vsObj->type != WLZ_3D_VIEW_STRUCT)
      {
	 errNum = WLZ_ERR_OBJECT_TYPE;
      }
      else if(vsObj->domain.core == NULL)
      {
	 errNum = WLZ_ERR_DOMAIN_NULL;
      }
      else
      {
	vs = WlzAssign3DViewStruct(vsObj->domain.vs3d, NULL);
      }
      WlzFreeObj(vsObj);
    }
    else
    {
      vs = WlzAssign3DViewStruct(
           WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum), NULL);
      if(vs)
      {
	vs->view_mode = WLZ_UP_IS_UP_MODE;
	vs->up.vtX = 0.0;
	vs->up.vtY = 0.0;
	vs->up.vtZ = -1.0;
	switch(useXYZ)
	{
	  case 'x':
	    vs->phi   = ALG_M_PI_2;
	    vs->theta = 0.0;
	    vs->zeta  = ALG_M_PI_2;
	    break;
	  case 'y':
	    vs->phi   = ALG_M_PI_2;
	    vs->theta = ALG_M_PI_2;
	    vs->zeta  = ALG_M_PI_2;
	    break;
	  case 'z': /* FALLTHROUGH */
	  default:
	    vs->phi   = 0.0;
	    vs->theta = 0.0;
	    vs->zeta  = 0.0;
	    break;
	}
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      "%s: Failed to create plane definition view structure (%s).\n",
      *argv, errMsg);
    }
  }
  if(ok)
  {
    errNum = Wlz3DSectionOcc(inObj, vs, sep, &first, &last,
    			     &nOcc, &occ);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      "%s: Failed to compute occupancy of planes (%s).\n",
      *argv, errMsg);
    }
  }
  if(ok)
  {
    FILE  	*fP = NULL;

    if((fP = (strcmp(outFileStr, "-")? fopen(outFileStr, "w"):
                                       stdout)) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
    else
    {
      int 	idN;
      char	fs;

      fs = (useCommas)? ',': ' ';
      if(useRow)
      {
	int	k = 0,
		nK = 0;

        for(idN = 0; idN < nOcc; ++idN)
	{
	  if(occ[idN])
	  {
	    ++nK;
	  }
	}
	if(nK > 0)
	{
	  for(idN = 0; idN < nOcc; ++idN)
	  {
	    if(occ[idN])
	    {
	      if(++k < nK)
	      {
	        (void )fprintf(fP, "%g%c", first + (idN * sep), fs);
	      }
	      else
	      {
	        (void )fprintf(fP, "%g", first + (idN * sep));
	      }
	    }
	  }
	}
	(void )fprintf(fP, "\n");
      }
      else
      {
	for(idN = 0; idN < nOcc; ++idN)
	{
	  (void )fprintf(fP, "%g%c%d\n",
		 first + (idN * sep), fs,
		 (useArea)? occ[idN]: occ[idN] != 0);
	}
      }
    }
    if(strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      "%s: Failed to write occupancies to file %s (%s).\n",
      *argv, outFileStr, errMsg);
    }
  }
  AlcFree(occ);
  (void )WlzFreeObj(inObj);
  (void )WlzFree3DViewStruct(vs);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-h] [-a] [-o<out>] [-b<sec bib>] [-c] [-p<sep>]\n"
    "\t\t[-s<sec trm>] [-x] [-y] [-z] [<in obj>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "-h Prints usage information.\n"
    "-a Output plane occupancy area or number of domains rather than the\n"
    "   default binary\n"
    "   values (ie 0 or 1).\n"
    "-b Planes defined by an MAPaint style bib file.\n"
    "-c Output plane occupancy numbers using comma rather than white space\n"
    "   separation.\n"
    "-p Plane separation distance (default 1.0).\n"
    "-r Output plane occupany as a single row of plane coordinates.\n"
    "-s Planes defined by a Woolz section transform object.\n"
    "-x Planes are columns in the domain.\n"
    "-y Planes are lines in the domain.\n"
    "-z Planes are object planes in the domain.\n"
    "Computes the occupancy of a domain or compound array object with\n"
    "respect to the planes through it.\n"
    "By default these planes are simply the x-y planes orthogonal to\n"
    "the z-axis but a different set of planes may be specified by: an\n"
    "MAPaint style bib file, a Woolz section transform object or by\n"
    "selecting one of the orthogonal axis plane propagation directions\n"
    "(x, y or z). The occupancy is output as two columns, the first\n"
    "containing the plane coordinate (with respect to a section transform)\n"
    "and the second the occupancy. The occupancy may bean area, the number\n"
    "of domains or a binary value (ie 0 or 1).\n"
    "By default the domain is read from the standard input and the occupancy\n"
    "is written to the standard output.\n",
    *argv,
    "-c heart.wlz\n"
    "Outputs a space seperated list of numbers with two columns. The first\n"
    "column containing the (z) plane coordinate and the second the\n"
    "occupancy of the domain (read from heart.wlz) for that plane.\n");
  }
  return(errNum);
}

/*!
* \return	New 3D view data structure read from an MAPaint style
* 		bib file.
* \ingroup	BinWlzApp
* \brief	Reads a view data structure read from an MAPaint style bib
* 		file.
* \param	fStr			File string.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzThreeDViewStruct *WlzDomainOccReadBib(char *fStr, WlzErrorNum *dstErr)
{
  int		tI;
  char		*errMsg = NULL;
  char		vsBuf[64];
  FILE		*fP = NULL;
  BibFileRecord *bibRec;
  BibFileField  *bibFld;
  WlzThreeDViewStruct *vs = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((fP = (strcmp(fStr, "-")?
           fopen(fStr, "rb"): stdin)) == NULL)
  {
    errNum = WLZ_ERR_READ_EOF;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    vs = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    BibFileError  bibErr;

    /* Read the bib file until the first section view entry. */
    bibErr = BibFileRecordRead(&bibRec, &errMsg, fP);
    while((bibErr == BIBFILE_ER_NONE) &&
          (strncmp(bibRec->name, "Wlz3DSectionViewParams", 22)))
    {
      BibFileRecordFree(&bibRec);
      bibErr = BibFileRecordRead(&bibRec, &errMsg, fP);
    }
    if(bibErr == BIBFILE_ER_NONE)
    {
      /* Having found the record, read and parse it. */
      (void )BibFileFieldParseFmt(bibRec->field,
               (void *)&(vs->fixed.vtX), "%lg ,%*lg ,%*lg", "FixedPoint",
               (void *)&(vs->fixed.vtY), "%*lg ,%lg ,%*lg", "FixedPoint",
               (void *)&(vs->fixed.vtZ), "%*lg ,%*lg ,%lg", "FixedPoint",
               (void *)&(vs->dist), "%lg", "Distance",
               (void *)&(vs->phi), "%lg", "Pitch",
               (void *)&(vs->theta), "%lg", "Yaw",
               (void *)&(vs->zeta), "%lg", "Roll",
               (void *)&(vs->scale), "%lg", "Scale",
               (void *)&(vs->up.vtX), "%lg ,%*lg ,%*lg", "UpVector",
               (void *)&(vs->up.vtY), "%*lg ,%lg ,%*lg", "UpVector",
               (void *)&(vs->up.vtZ), "%*lg ,%*lg ,%lg", "UpVector",
               (void *)vsBuf, "%s", "ViewMode",
               NULL);
      /* Doesn't read the view mode correctly. */
      bibFld = bibRec->field;
      while(bibFld)
      {
        if(strncmp(bibFld->name, "ViewMode", 8) == 0)
        {
          strcpy(vsBuf, bibFld->value);
          break;
        }
      }
      bibFld = bibFld->next;
      /* Convert angles to radians and get mode */
      vs->phi *= WLZ_M_PI /180.0;
      vs->theta *= WLZ_M_PI /180.0;
      vs->zeta *= WLZ_M_PI /180.0;
      if(WlzStringMatchValue(&tI, optarg,
                             "up-is-up", WLZ_UP_IS_UP_MODE,
                             "statue", WLZ_STATUE_MODE,
                             "absolute", WLZ_ZETA_MODE,
                             NULL))
      {
        vs->view_mode = (WlzThreeDViewMode )tI;
      }
    }
    BibFileRecordFree(&bibRec);
  }
  if(fP && strcmp(fStr, "-"))
  {
    (void )fclose(fP);
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFree3DViewStruct(vs);
    vs = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(vs);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

