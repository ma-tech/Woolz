#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzDrawDomainObj_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzDrawDomainObj.c
* \author       Bill Hill
* \date         September 2008
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2008 Medical research Council, UK.
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
* \brief	Reads a string of drawing commands and from them
* 		draws a 3D domain object. All drawing commands are
* 		on a section plane which is defined by the view
* 		transform.
* \ingroup	BinWlzApp
* \todo         -
* \bug          None known.
*/

/*!
\ingroup BinWlzApp
\defgroup wlzdrawdomainobj WlzDrawDomainObj
\par Name
WlzDrawDomainObj - reads drawing commands from either a string given on
the command line, a file given on the command line or the standard input
(in this order of precidence).
The drawing commands are used to create a drawn section which is placed
into 3D using the view transform.
\par Synopsis
\verbatim
WlzDrawDomainObj [-2] [-a<pitch,yaw,roll>] [-f <fx,fy,fz>]
                 [-d <dist> [-b <view bib file>]
                 [-m <mode>] [-u<ux,uy,uz>]
                 [-g <ox,oy,oz>] [-h] [-o<out>]
                 [-s <cmd str>] [<cmd str file>>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-2</b></td>
    <td>Ignore the view struct and keep as a 2D object.</td>
  </tr>
  <tr>
    <td><b>-a</b></td>
    <td>Viewing angles: pitch (phi), yaw (theta) and roll (zeta),
        default 0.0,0.0,0.0.</td>
  </tr>
  <tr>
    <td><b>-f</b></td>
    <td>Fixed point position, default 0.0,0.0,0.0.</td>
  </tr>
  <tr>
    <td><b>-d</b></td>
    <td>Distance parameter, default 0.0.</td>
  </tr>
  <tr>
    <td><b>-b</b></td>
    <td>Bib file with view parameters, e.g. saved from MAPaint.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Viewing mode, one of: up-is-up, statue or absolute, default
        is up-is-up.</td>
  </tr>
  <tr>
    <td><b>-u</b></td>
    <td>Up vector, default 0.0,0.0,1.0.</td>
  </tr>
  <tr>
    <td><b>-g</b></td>
    <td>Origin of the drawing with respect to the 2D Woolz object
        cut using the view transform, default 0.0,0.0.</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Drawing command string.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
</table>
\par Description
Reads drawing commands from either a string given on
the command line, a file given on the command line or the standard input
(in this order of precidence).
The drawing commands are used to create a drawn section which is placed
into 3D using the view transform.

The command string must have the following syntax:
\verbatim
  <command string> = <init command>[<command>]+<end command>
  <init command> = <ident>:<version>;
  <end command> = END:;
  <ident> = WLZ_DRAW_DOMAIN
  <version> = 1
  <command> = <command name>:[<parameter>[,<parameter>]*];
  <command name> = PEN | LINE | CIRCLE
\endverbatim
In addition to the init and end commands, the following
drawing commands are recognised:
\verbatim
  CIRCLE:<action>:<radius>,<x>,<y>;
  LINE:<action>:<width>,<x>,<y>,<x>,<y>;
  PEN:<action>,<width>,<x>,<y>[,<x>,<y>]*;
\endverbatim
Where:
\verbatim
 <action> = DRAW | ERASE
\endverbatim
The circle command draws or erases a filled circle which is specified
by it's radius and centre parameters.
The line command draws a rectangle using the given width and end
coordinates of the mid-line.
The pen command draws a polyline which is composed of a series of
rectangular segments, with each segment ending in a semi-circular cap.
The parameters of the pen command are the width and line segment end
point coordinates.
All widths, radii and coordinates may be in any floating point format
recognised by scanf(3).
Other commands may be present provided they have the same syntax
described above, but they will be ignored.
All white space characters are ignored.
\par Examples
\verbatim
WlzDrawDomainObj -o out.wlz \\
                 -s 'WLZ_DRAW_DOMAIN:1; CIRCLE:DRAW,100,200,300; END:;'
\endverbatim
This creates a 3D domain object with a single plane at z = 0
which has a domain that is  a single circle (radius = 100,
centre = 200,300).
to the file out.num.
\par File
\ref WlzDrawDomainObj.c "WlzDrawDomainObj.c"
\par See Also
\ref WlzDomainOps WlzDrawDomainObj().
\ref BinWlzApp "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <float.h>
#include <Wlz.h>
#include <bibFile.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

static int			WlzDrawDomObjScanPair(
				  char *str,
				  double *d0,
				  double *d1);
static int			WlzDrawDomObjScanTriple(
				  char *str,
				  double *d0,
				  double *d1,
				  double *d2);
static char			*WlzDrawDomObjReadStr(
				  char *fStr,
				  WlzErrorNum *dstErr);
static  WlzErrorNum 		WlzDrawDomObjReadView(
				  char *fStr,
				  WlzThreeDViewStruct *view,
				  char **errMsg);

int             main(int argc, char **argv)
{
  int		tI,
		errIdx,
  		option,
		keep2D = 0,
		ok = 1,
		usage = 0;
  WlzDVertex2	org;
  WlzDVertex3	upVector;
  WlzThreeDViewStruct *view = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*dwnObj = NULL;
  FILE		*fP = NULL;
  char 		*errMsg0,
  		*cmdStr = NULL,
  		*outFileStr,
        	*inFileStr;
  const char	*errMsg;
  static char	optList[] = "a:b:d:f:g:m:o:s:u:2h",
		outFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  opterr = 0;
  errMsg = errMsg0 = "";
  org.vtX = org.vtY = 0.0;
  upVector.vtX = 0.0; upVector.vtY = 0.0; upVector.vtZ = 1.0;
  inFileStr = inFileStrDef;
  outFileStr = outFileStrDef;
  if((view = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum)) == NULL)
  {
    ok = 0;
    (void )WlzStringFromErrorNum(errNum, &errMsg);
    (void )fprintf(stderr,
		   "%s: failed to create view structure (%s)\n",
		   *argv, errMsg);
  }
  if(ok)
  {
    view->dist = 0.0;
    view->scale = 1.0;
    view->ref_obj = NULL;
    view->view_mode = WLZ_UP_IS_UP_MODE;
    view->phi = view->theta = view->zeta = 0.0;
    view->up.vtX = view->up.vtY = 0.0; view->up.vtZ = 1.0;
    view->fixed.vtX = view->fixed.vtY = view->fixed.vtZ = 0.0;
    while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
    {
      switch(option)
      {
	case '2':
	  keep2D = 1;
	  break;
	case 'a':
	  usage = WlzDrawDomObjScanTriple(optarg,
	                                  &(view->phi),
					  &(view->theta),
					  &(view->zeta)) < 1;
	  break;
	case 'f':
	  usage = WlzDrawDomObjScanTriple(optarg,
	                                  &(view->fixed.vtX),
					  &(view->fixed.vtY),
					  &(view->fixed.vtZ)) < 1;
	  break;
	case 'd':
	  usage = sscanf(optarg, "%lg", &(view->dist)) != 1;
	  break;
	case 'b':
	  usage = WlzDrawDomObjReadView(optarg, view,
					&errMsg0) != WLZ_ERR_NONE;
	  if(usage)
	  {
	    (void )fprintf(stderr,
			   "%s: failed to read view parameters from file "
			   "%s (%s)\n",
			   *argv, optarg, errMsg0);
	  }
	  break;
	case 'm':
	  if(WlzStringMatchValue(&tI, optarg,
				 "up-is-up", WLZ_UP_IS_UP_MODE,
				 "statue", WLZ_STATUE_MODE,
				 "absolute", WLZ_ZETA_MODE,
				 NULL))
	  {
	    view->view_mode = (WlzThreeDViewMode )tI;
	  }
	  else
	  {
	    usage = 1;
	  }
	  break;
	case 'u':
	  usage = WlzDrawDomObjScanTriple(optarg,
					  &(view->up.vtX),
					  &(view->up.vtY),
					  &(view->up.vtZ)) < 1;
	  break;
	case 'g':
	  usage = WlzDrawDomObjScanPair(optarg, &(org.vtX), &(org.vtY)) < 1;
	  break;
	case 'o':
	  outFileStr = optarg;
	  break;
	case 's':
	  if((cmdStr = AlcStrDup(optarg)) == NULL)
	  {
	    usage = 1;
	    errNum = WLZ_ERR_MEM_ALLOC;
	  }
	  break;
	case 'h':
	default:
	  usage = 1;
	  ok = 0;
	  break;
      }
    }
    if((usage == 0) && (cmdStr == NULL))
    {
      if((optind == argc))
      {
	inFileStr = inFileStrDef;
      }
      else if(optind + 1 == argc )
      {
	inFileStr = *(argv + optind);
      }
      else
      {
	usage = 1;
      }
      if(usage == 0)
      {
	cmdStr = WlzDrawDomObjReadStr(inFileStr, &errNum);
      }
    }
    ok = usage == 0;
  }
  if(ok)
  {
    dwnObj = WlzDrawDomainObj(org, view, keep2D, cmdStr, &errIdx, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: failed to draw object in command prior to string\n"
		     "position %d (%s)\n",
		     *argv, errIdx, errMsg);
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if((fP = (strcmp(outFileStr, "-")? fopen(outFileStr, "w"):
	      			       stdout)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to open output file (%s).\n",
		     *argv, errMsg);
    }
    else
    {
      if((errNum = WlzWriteObj(fP, dwnObj)) != WLZ_ERR_NONE)
      {
        ok = 0;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	               "%s: Failed to write output object (%s).\n",
		       argv[0], errMsg);
      }
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      fclose(fP);
    }
  }
  AlcFree(cmdStr);
  (void )WlzFree3DViewStruct(view);
  (void )WlzFreeObj(dwnObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%sExample: %s%s",
    *argv,
    " [-2] [-a<pitch,yaw,roll>] [-f <fx,fy,fz>]\n"
    "                 [-d <dist> [-b <view bib file>]\n"
    "                 [-m <mode>] [-u<ux,uy,uz>]\n"
    "                 [-g <ox,oy,oz>] [-h] [-o<out>]\n"
    "                 [-s <cmd str>] [<cmd str file>>]\n"
    "Options:\n"
    "  -2  Ignore the view struct and keep as a 2D object.\n"
    "  -a  Viewing angles: pitch (phi), yaw (theta) and roll (zeta),\n"
    "      default 0.0,0.0,0.0.\n"
    "  -f  Fixed point position, default 0.0,0.0,0.0.\n"
    "  -d  Distance parameter, default 0.0.\n"
    "  -b  Bib file with view parameters, e.g. saved from MAPaint.\n"
    "  -m  Viewing mode, one of: up-is-up, statue or absolute, default\n"
    "      is up-is-up.\n"
    "  -u  Up vector, default 0.0,0.0,1.0.\n"
    "  -g  Origin of the drawing with respect to the 2D Woolz object\n"
    "      cut using the view transform, default 0.0,0.0.\n"
    "  -s  Drawing command string.\n"
    "  -h  Help, prints this usage message.\n"
    "  -o  Output file name.\n"
    "Reads drawing commands from either a string given on the command line,\n"
    "a file given on the command line or the standard input (in this order\n"
    "of precidence). The drawing commands are used to create a drawn\n"
    "section which is placed into 3D using the view transform.\n"
    "\n"
    "The command string must have the following syntax:\n"
    "  <command string> = <init command>[<command>]+<end command>\n"
    "  <init command> = <ident>:<version>;\n"
    "  <end command> = END:;\n"
    "  <ident> = WLZ_DRAW_DOMAIN\n"
    "  <version> = 1\n"
    "  <command> = <command name>:[<parameter>[,<parameter>]*];\n"
    "  <command name> = PEN | LINE | CIRCLE\n"
    "In addition to the init and end commands, the following\n"
    "drawing commands are recognised:\n"
    "  CIRCLE:<action>:<radius>,<x>,<y>;\n"
    "  LINE:<action>:<width>,<x>,<y>,<x>,<y>;\n"
    "  PEN:<action>,<width>,<x>,<y>[,<x>,<y>]*;\n"
    "Where:\n"
    " <action> = DRAW | ERASE\n"
    "The circle command draws or erases a filled circle which is specified\n"
    "by it's radius and centre parameters.\n"
    "The line command draws a rectangle using the given width and end\n"
    "coordinates of the mid-line.\n"
    "The pen command draws a polyline which is composed of a series of\n"
    "rectangular segments, with each segment ending in a semi-circular cap.\n"
    "The parameters of the pen command are the width and line segment end\n"
    "point coordinates.\n"
    "All widths, radii and coordinates may be in any floating point format\n"
    "recognised by scanf(3).\n"
    "Other commands may be present provided they have the same syntax\n"
    "described above, but they will be ignored.\n"
    "All white space characters are ignored.\n",
    *argv,
    " -o out.wlz -s 'WLZ_DRAW_DOMAIN:1; CIRCLE:DRAW,100,200,300; END:;'\n"
    "This creates a 3D domain object with a single plane at z = 0\n"
    "which has a domain that is  a single circle (radius = 100,\n"
    "centre = 200,300).\n");
  }
  return(!ok);
}

/*!
* \return	Number of fields found, which will be -ve if there is a
* 		parsing error.
* \ingroup	BinWlzApp
* \brief	Parses the input string for two coma seperated double
* 		values, any of which may be missing, eg "4, " and
* 		",1.0 " are both valid (white space is ignored).
* \param	str			String to parse for double values.
* \param	d0			Destination pointer for first value,
* 					must not be NULL.
* \param	d1			Destination pointer for second value,
* 					must not be NULL.
* \param	d2			Destination pointer for third value,
* 					must not be NULL.
*/
static int	WlzDrawDomObjScanPair(char *str, double *d0, double *d1)
{
  int		idx,
  		cnt = 0;
  char		*subStr[2];
  double	*dP[2];

  if(str)
  {
    dP[0] = d0; dP[1] = d1;
    while(*str && isspace(*str))
    {
      ++str;
    }
    if(*str == ',')
    {
      subStr[0] = NULL;
      subStr[1] = strtok(str, ",");
    }
    else
    {
      subStr[0] = strtok(str, ",");
      subStr[1] = strtok(NULL, ",");
    }
    if((subStr[0] == NULL) && (subStr[1] == NULL))
    {
      cnt = -1;
    }
    else
    {
      for(idx = 0; idx < 2; ++idx)
      {
	if(subStr[idx] && (sscanf(subStr[idx], "%lg", dP[idx]) == 1))
	{
	  ++cnt;
	}
	else
	{
	  cnt = -1;
	  break;
	}
      }
    }
  }
  return(cnt);
}

/*!
* \return	Number of fields found, which will be -ve if there is a
* 		parsing error.
* \ingroup	BinWlzApp
* \brief	Parses the input string for three coma seperated double
* 		values, any of which may be missing, eg "4, ,2.367 " and
* 		",1.0 , " are both valid (white space is ignored).
* \param	str			String to parse for double values.
* \param	d0			Destination pointer for first value,
* 					must not be NULL.
* \param	d1			Destination pointer for second value,
* 					must not be NULL.
* \param	d2			Destination pointer for third value,
* 					must not be NULL.
*/
static int	WlzDrawDomObjScanTriple(char *str,
					double *d0, double *d1, double *d2)
{
  int		idx,
  		cnt = 0;
  char		*subStr[3];
  double	*dP[3];

  if(str)
  {
    dP[0] = d0; dP[1] = d1; dP[2] = d2; 
    while(*str && isspace(*str))
    {
      ++str;
    }
    if(*str == ',')
    {
      subStr[0] = NULL;
      subStr[1] = strtok(str, ",");
    }
    else
    {
      subStr[0] = strtok(str, ",");
      subStr[1] = strtok(NULL, ",");
    }
    if((subStr[0] == NULL) && (subStr[1] == NULL))
    {
      cnt = -1;
    }
    else
    {
      for(idx = 0; idx < 3; ++idx)
      {
	if(subStr[idx] && (sscanf(subStr[idx], "%lg", dP[idx]) == 1))
	{
	  ++cnt;
	}
	else
	{
	  cnt = -1;
	  break;
	}
      }
    }
  }
  return(cnt);
}

/*!
* \return	String or NULL on error. AlcFree() should be used to free
* 		the string.
* \ingroup	BinWlzApp
* \brief	Reads a string from a file with the given file name.
* \param	fStr			File name string, "-" is standard
* 					input.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static char	*WlzDrawDomObjReadStr(char *fStr, WlzErrorNum *dstErr)
{
  int		tI;
  size_t	len = 0;
  char		*c,
  		*str = NULL;
  AlcVector	*vec = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  FILE		*fP = NULL;
  const size_t	initStrMax = 4096,
  		blkSz = 4096;

  if((fP = (strcmp(fStr, "-")?
	   fopen(fStr, "rb"): stdin)) == NULL)
  {
    errNum = WLZ_ERR_READ_EOF;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((vec = AlcVectorNew(initStrMax, sizeof(char), blkSz, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    c = AlcVectorExtendAndGet(vec, len);
    while((c != NULL) && ((tI = getc(fP)) != EOF))
    {
      *c = tI;
      c = AlcVectorExtendAndGet(vec, ++len);
    }
    if(c == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      *c = '\0';
      ++len;
    }
  }
  if(fP && strcmp(fStr, "-"))
  {
    (void )fclose(fP);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if((str = (char *)AlcVectorToArray1D(vec, 0, len, NULL)) == NULL)
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  (void )AlcVectorFree(vec);
  return(str);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzBinApp
* \brief	Reads the parameters of a 3D view data structure from the
* 		given file which uses the same bib file format as MAPaint.
* 		Code for reading the bib file has been partly copied from
* 		Wlz3DViewTransformObj.c.
* \param	fStr			File name string.
* \param	view			3D view structure in which to set
* 					field values, must not be NULL.
* \param	errMsg			Destination pointer for error message,
* 					must not be NULL.
*/
static  WlzErrorNum WlzDrawDomObjReadView(char *fStr,
				WlzThreeDViewStruct *view,
				char **errMsg)
{
  int           tI;
  char          viewModeStr[64];
  FILE		*fP = NULL;
  BibFileRecord *bibRec;
  BibFileField  *bibFld;
  BibFileError  bibErr;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((fP = (strcmp(fStr, "-")?
	   fopen(fStr, "rb"): stdin)) == NULL)
  {
    errNum = WLZ_ERR_READ_EOF;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    /* Read the bib file until the first section view entry. */
    bibErr = BibFileRecordRead(&bibRec, errMsg, fP);
    while((bibErr == BIBFILE_ER_NONE) &&
	  (strncmp(bibRec->name, "Wlz3DSectionViewParams", 22)))
    {
      BibFileRecordFree(&bibRec);
      bibErr = BibFileRecordRead(&bibRec, errMsg, fP);
    }
    if(bibErr == BIBFILE_ER_NONE)
    {
      /* Having found the record, read and parse it. */
      (void )BibFileFieldParseFmt(bibRec->field,
	       (void *)&(view->fixed.vtX), "%lg ,%*lg ,%*lg", "FixedPoint",
	       (void *)&(view->fixed.vtY), "%*lg ,%lg ,%*lg", "FixedPoint",
	       (void *)&(view->fixed.vtZ), "%*lg ,%*lg ,%lg", "FixedPoint",
	       (void *)&(view->dist), "%lg", "Distance",
	       (void *)&(view->phi), "%lg", "Pitch",
	       (void *)&(view->theta), "%lg", "Yaw",
	       (void *)&(view->zeta), "%lg", "Roll",
	       (void *)&(view->scale), "%lg", "Scale",
	       (void *)&(view->up.vtX), "%lg ,%*lg ,%*lg", "UpVector",
	       (void *)&(view->up.vtY), "%*lg ,%lg ,%*lg", "UpVector",
	       (void *)&(view->up.vtZ), "%*lg ,%*lg ,%lg", "UpVector",
	       (void *)viewModeStr, "%s", "ViewMode",
	       NULL);
      /* Doesn't read the view mode correctly. */
      bibFld = bibRec->field;
      while(bibFld)
      {
	if(strncmp(bibFld->name, "ViewMode", 8) == 0)
	{
	  strcpy(viewModeStr, bibFld->value);
	  break;
	}
      }
      bibFld = bibFld->next;
      /* Convert angles to radians and get mode */
      view->phi *= WLZ_M_PI /180.0;
      view->theta *= WLZ_M_PI /180.0;
      view->zeta *= WLZ_M_PI /180.0;
      if(WlzStringMatchValue(&tI, optarg,
			     "up-is-up", WLZ_UP_IS_UP_MODE,
			     "statue", WLZ_STATUE_MODE,
			     "absolute", WLZ_ZETA_MODE,
			     NULL))
      {
	view->view_mode = (WlzThreeDViewMode )tI;
      }
    }
    BibFileRecordFree(&bibRec);
  }
  if(fP && strcmp(fStr, "-"))
  {
    (void )fclose(fP);
  }
  return(errNum);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */
