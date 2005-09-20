#pragma ident "MRC HGU $Id$"
/*!
* \file         binWlzApp/WlzReconstruct.c
* \author       Bill Hill
* \date         June 2005
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2005 Medical research Council, UK.
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
* \brief	Reconstructs a 3D object using the information a
* 		Reconstruct bibfile.
* \ingroup	BinWlzApp
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzreconstruct "WlzReconstruct"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzreconstruct WlzReconstruct
\par Name
WlzReconstruct - reconstructs a 3D object using the information in a
                Reconstruct bibfile.
\par Synopsis
\verbatim
WlzReconstruct  [-h] [-o <out obj file>] [-T <out obj fmt>] [-b <out bibfile>]
                [-f] [-g] [-h] [-i] [-l] [-L] [-p]
		[-s #] [-m #] [-c #,#,#] [-C #,#,#] <input bibfile>
                
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr>
    <td>o</td>
    <td>Output object file name.</td>
  </tr>
  <tr>
    <td>T</td>
    <td>Output object file format.</td>
  </tr>
  <tr>
    <td>b</td>
    <td>Output bib file, the default is the same file name as the output
      object file name but with a .bib file extension.</td>
  </tr>
  <tr>
    <td>o</td>
    <td>Output file name.</td>
  </tr>
  <tr>
    <td>f</td>
    <td>Force overwriting of files, without this option confirmation is
      needed.</td>
  </tr>
  <tr>
    <td>g</td>
    <td>Use resource greedy algorithms.</td>
  </tr>
  <tr>
    <td>h</td>
    <td>Help - show this usage information.</td>
  </tr>
  <tr>
    <td>i</td>
    <td>Use integer scaling.</td>
  </tr>
  <tr>
    <td>l</td>
    <td>Use slow but accurate filter when sampling output.</td>
  </tr>
  <tr>
    <td>p</td>
    <td>Use fast sampling.</td>
  </tr>
  <tr>
    <td>s</td>
    <td>Scaling factor (eg 0.5 implies half linear dimension of source).</td>
  </tr>
  <tr>
    <td>m</td>
    <td>Match all section histograms to the section with the given index.</td>
  </tr>
  <tr>
    <td>c</td>
    <td>Clip the source (2D sections) using the given bounding box.</td>
  </tr>
  <tr>
    <td>C</td>
    <td>Clip the destination (3D object) using the given bounding box.</td>
  </tr>
</table>
\par Description
Reads a reconstruction bibfile either from the command line or the
standard input and constructs a 3D object. The reconstruction process
is controlled by the values in the bibfile, many of which can be
overridden on the command line. Unless the force flag is set
overwriting files will require a confirmation.

By default ./WlzReconstruct reads from it's standard input and
writes to it's standard output.
Any bounding boxes specified on the command line should have the format
\verbatim
  xMin,yMin,zMin,xMax,yMax,zMax
\endverbatim
Both bounding box and scale components may be specified using standard
floating point formats.
\par Examples
\verbatim
WlzReconstruct -f -b new.bib recon.bib
\endverbatim
Reads the reconstruction bib file from recon.bib and constructs a 3D
object which is written to the object file name specified in the
bib file. a new bib file will be created for the reconstruction and
written to new.bib. Because of the -f flag any existing files may
be overwriten without confirmation.
\par File
\ref WlzReconstruct.c "WlzReconstruct.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref wlzconstruct3d "WlzConstruct3D(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <errno.h>
#include <ctype.h>
#include <Wlz.h>
#include <bibFile.h>
#include <Reconstruct.h>

typedef	struct _WlzRecCmdLnOptions
{
  int		fastSamFlg;
  int		filterFlg;
  int		forceFlg;
  int		greedyFlg;
  int		intScaleFlg;
  int		interpFlg;
  int		matchHFlg;
  int		matchHSec;
  char		*objFName;
  char		*inBibFName;
  char		*outBibFName;
  char		*dstBoxStr;
  char		*srcBoxStr;
  char		*scaleStr;
  WlzEffFormat	objFFmt;
} WlzRecCmdLnOptions;

static unsigned 		WlzRecParseBox3D(
    				  WlzDBox3 *box,
				  char *gStr);
static unsigned 		WlzRecParseVtx3D(
    				  WlzDVertex3 *vtx,
				  char *gStr);
static unsigned 		WlzRecParseCSD(
    				  int nVal,
				  double *val,
				  char *gStr);
static	int			WlzRecConfirmFileOverwrite(
    				  char *app,
    				  char *fName);
static RecError			WlzRecMakeFilePathAbs(
    				  char **path,
				  char *file);

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;

int		 main(int argc, char *argv[])
{
  int		option,
		ok = 1,
  		usage = 0,
		numSec = 0;
  FILE		*fP = NULL;
  RecError	recErr = REC_ERR_NONE;
  WlzRecCmdLnOptions recOptions;
  RecSectionList *secList = NULL;
  char		*cp0 = NULL,
		*cp1,
		*cp2,
		*errMsg;
  const char	pathSep = '/';
  static char	optList[] = "fghilLpb:c:C:m:o:s:T:",
  		defFile[] = "-",
		defErrMsg[] = "";

  opterr = 0;
  errMsg = defErrMsg;
  memset(&recOptions, 0, sizeof(WlzRecCmdLnOptions));
  recOptions.inBibFName = defFile;
  while((usage == 0) &&
        ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'f':
	recOptions.forceFlg = 1;
	break;
      case 'g':
	recOptions.greedyFlg = 1;
	break;
      case 'i':
	recOptions.intScaleFlg = 1;
	break;
      case 'l':
	recOptions.filterFlg = 1;
	break;
      case 'p':
	recOptions.fastSamFlg = 1;
	break;
      case 'b':
	recOptions.outBibFName = optarg;
	break;
      case 'c':
	recOptions.srcBoxStr = optarg;
	break;
      case 'C':
	recOptions.dstBoxStr = optarg;
	break;
      case 'L':
	recOptions.interpFlg = 1;
	break;
      case 'm':
	recOptions.matchHFlg = 1;
	usage = sscanf(optarg, "%d", &(recOptions.matchHSec)) != 1;
	break;
      case 'o':
	recOptions.objFName = optarg;
	break;
      case 's':
	recOptions.scaleStr = optarg;
	break;
      case 'T':
	recOptions.objFFmt = WlzEffStringToFormat(optarg);
	usage = recOptions.objFFmt == WLZEFF_FORMAT_NONE;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      recOptions.inBibFName = *(argv + optind);
    }
  }
  ok = !usage;
  /* Read input bibfile. */
  if(ok)
  {
    if((secList = RecSecNewSectionList(&recErr)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to memory for allocate section list.\n",
		     argv[0]);
    }
  }
  if(ok)
  {
    secList->reconstruction.srcBox.xMin = INT_MIN;
    secList->reconstruction.srcBox.yMin = INT_MIN;
    secList->reconstruction.srcBox.zMin = INT_MIN;
    secList->reconstruction.srcBox.xMax = INT_MAX;
    secList->reconstruction.srcBox.yMax = INT_MAX;
    secList->reconstruction.srcBox.zMax = INT_MAX;
    secList->reconstruction.dstBox = secList->reconstruction.srcBox;
    if((recOptions.inBibFName == NULL) ||
       (*recOptions.inBibFName == '\0') ||
       ((fP = (strcmp(recOptions.inBibFName, "-")?
	      fopen(recOptions.inBibFName, "r"): stdin)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open input file %s\n",
		     argv[0], (recOptions.inBibFName)?
		              recOptions.inBibFName: "(null)");
    }
  }
  if(ok)
  {
    if(RecFileSecListRead(secList, &numSec, fP, &errMsg) != REC_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read section list from file %s.\n%s%s",
		     argv[0],
		     recOptions.inBibFName,
		     (errMsg && *errMsg)? "Error message - ": "",
		     (errMsg && *errMsg)? errMsg: "");
    }
  }
  if(fP && (strcmp(recOptions.inBibFName, "-") == 0))
  {
    (void )fclose(fP);
  }
  fP = NULL;
  /* Over write reconstruction parameters specified on the command line. */
  if(ok)
  {

    if(recOptions.fastSamFlg)
    {
      secList->reconstruction.fastSam = recOptions.fastSamFlg;
    }
    if(recOptions.filterFlg)
    {
      secList->reconstruction.gaussFlt = recOptions.filterFlg;
    }
    if(recOptions.greedyFlg)
    {
      secList->reconstruction.greedy = recOptions.greedyFlg;
    }
    if(recOptions.intScaleFlg)
    {
      secList->reconstruction.intScale = recOptions.intScaleFlg;
    }
    if(recOptions.matchHFlg)
    {
      secList->reconstruction.matchHst = recOptions.matchHFlg;
      secList->reconstruction.matchSec = recOptions.matchHSec;
    }
    if(recOptions.objFFmt != WLZEFF_FORMAT_NONE)
    {
      secList->reconstruction.fileFormat = recOptions.objFFmt;
    }
    if(recOptions.scaleStr)
    {
      if(sscanf(recOptions.scaleStr, "%lg",
	        &(secList->reconstruction.scale.vtX)) != 1)
      {
	usage = 1;
        ok = 0;
      }
      else
      {
	secList->reconstruction.scale.vtY = secList->reconstruction.scale.vtX;
      }
    }
  }
  if(ok)
  {
    if(recOptions.dstBoxStr)
    {
      secList->reconstruction.clipDst = 1;
      if(WlzRecParseBox3D(&(secList->reconstruction.dstBox),
	    	          recOptions.dstBoxStr) == 0)
       {
	 usage = 1;
         ok = 0;
       }
    }
  }
  if(ok)
  {
    if(recOptions.srcBoxStr)
    {
      secList->reconstruction.clipSrc = 1;
      if(WlzRecParseBox3D(&(secList->reconstruction.srcBox),
	    	          recOptions.srcBoxStr) == 0)
       {
	 usage = 1;
         ok = 0;
       }
    }
  }
  if(ok)
  {
    recErr = RecConstruct3DObj(&(secList->reconstruction.obj), secList->list,
			0.0, secList->reconstruction.gaussFlt,
                        (recOptions.interpFlg)?
			WLZ_INTERPOLATION_LINEAR: WLZ_INTERPOLATION_NEAREST,
			secList->reconstruction.fastSam,
			secList->reconstruction.greedy,
			secList->reconstruction.intScale,
			secList->reconstruction.scale,
			secList->reconstruction.matchHst,
			secList->reconstruction.matchSec,
			&(secList->reconstruction.srcBox),
			&(secList->reconstruction.dstBox),
			&errMsg);
    if(recErr != REC_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to construct 3D object.\n%s%s",
		     argv[0],
		     (errMsg && *errMsg)? "Error message - ": "",
		     (errMsg && *errMsg)? errMsg: "");
    }
  }
  if(ok)
  {
    recErr = WlzRecMakeFilePathAbs(&(secList->reconstruction.fileName),
				   recOptions.objFName);
    if(recErr != REC_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to build output object file path.\n",
		     argv[0]);
    }
  }
  /* Write reconstructed object. */
  if(ok)
  {
    if(recOptions.forceFlg == 0)
    {
      ok = WlzRecConfirmFileOverwrite(argv[0],
	  			      secList->reconstruction.fileName);
    }
  }
  if(ok)
  {
    recErr = RecErrorFromWlz(
	     WlzEffWriteObj(NULL, secList->reconstruction.fileName,
	       		secList->reconstruction.obj,
			secList->reconstruction.fileFormat));
  }
  /* Write output bibfile. */
  if(ok)
  {
    cp0 = NULL;
    if(recOptions.outBibFName)
    {
      cp0 = AlcStrDup(recOptions.outBibFName);
    }
    else
    {
      cp0 = (char *)AlcMalloc(sizeof(char) *
	    		strlen(secList->reconstruction.fileName) + 5);
      if(cp0)
      {
	cp1 = strrchr(secList->reconstruction.fileName, '.');
	if(cp1)
	{
	  cp2 = strrchr(secList->reconstruction.fileName, '/');
	  if(cp2 && (cp1 < cp2))
	  {
	    (void )sprintf(cp0, "%s.bib", secList->reconstruction.fileName);
	  }
	  else
	  {
	    (void )strcpy(cp0, secList->reconstruction.fileName);
	    cp1 = strrchr(cp0, '.') + 1;
	    (void )strcpy(cp1, "bib");
	  }
	}
	else
	{
	  (void )sprintf(cp0, "%s.bib", secList->reconstruction.fileName);
	}
      }
    }
    if(cp0 == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: Failed to allocate memory.\n",
		     argv[0]);
    }
    else if(recOptions.forceFlg == 0)
    {
      ok = WlzRecConfirmFileOverwrite(argv[0], cp0);
    }
  }
  if(ok)
  {
    if((*cp0 == '\0') ||
       ((fP = (strcmp(cp0, "-")?  fopen(cp0, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s\n",
		     argv[0], cp0);
    }
  }
  if(ok)
  {
    if(RecFileSecListWrite(fP, secList, numSec, &errMsg) != REC_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to write section list to file %s.\n%s%s",
		     argv[0],
		     cp0,
		     (errMsg && *errMsg)? "Error message - ": "",
		     (errMsg && *errMsg)? errMsg: "");
      ok = 0;
    }
  }
  if(fP && (strcmp(cp0, "-") == 0))
  {
    (void )fclose(fP);
  }
  AlcFree(cp0);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage:\n"
    "  %s [-o <out obj file>] [-T <out obj fmt>] [-b <out bibfile>]\n"
    "                 [-f] [-g] [-h] [-i] [-l] [-L] [-p]\n"
    "                 [-s #] [-m #] [-c #,#,#] [-C #,#,#] <in bibfile>\n"
    "Reads a reconstruction bibfile either from the command line or the\n"
    "standard input and constructs a 3D object. The reconstruction process\n"
    "is controlled by the values in the bibfile, many of which can be\n"
    "overridden on the command line. Unless the force flag is set\n"
    "overwriting files will require a confirmation.\n"
    "Options are:\n"
    "  -o  Output object file name.\n"
    "  -T  Output object file format.\n"
    "  -b  Output bibfile, the default is the same file name as the output\n"
    "      object file name but with a .bib file extension.\n"
    "  -o  Output file name.\n"
    "  -f  Force overwriting of files, without this option confirmation is\n"
    "      needed.\n"
    "  -g  Use resource greedy algorithms.\n"
    "  -h  Help - show this usage information.\n"
    "  -i  Use integer scaling.\n"
    "  -l  Use slow but accurate filter when sampling output.\n"
    "  -p  Use fast sampling.\n"
    "  -s  Scaling factor (eg 0.5 implies half linear dimension of source).\n"
    "  -m  Match all section histograms to the section with the given index.\n"
    "  -c  Clip the source (2D sections) using the given bounding box.\n"
    "  -C  Clip the destination (3D object) using the given bounding box.\n"
    "By default %s reads from it's standard input and\n"
    "writes to it's standard output.\n"
    "Any bounding boxes specified on the command line should have the format\n"
    "  xMin,yMin,zMin,xMax,yMax,zMax\n"
    "Both bounding box and scale components may be specified using standard\n"
    "floating point formats.\n"
    "Example:\n"
    "  %s -f -b new.bib recon.bib\n"
    "Reads the reconstruction bibfile from recon.bib and constructs a 3D\n"
    "object which is written to the object file name specified in the\n"
    "bibfile. a new bibfile will be created for the reconstruction and\n"
    "written to new.bib. Because of the -f flag any existing files may\n"
    "be overwriten without confirmation.\n",
    argv[0], argv[0], argv[0]);
  }
  exit(!ok);
}

/*!
* \return	Integer with bits set to indicate which components of the
* 		box are in the string:
* 		<ul>
* 		<li> (1<<0) xMin
* 		<li> (1<<1) yMin
* 		<li> (1<<2) zMin
* 		<li> (1<<3) xMax
* 		<li> (1<<4) yMax
* 		<li> (1<<5) zMax
* 		</ul>
*
* \brief	Parses the given string for a floating point bounding
* 		box with the following format
* 		  xmin,ymin,zmin,xmax,ymax,zmax
*		Any of the components of the box definition may be missing
*		and the return value indicates which components
*		have been parsed.
* \param	box			Box in which values will be set.
* \param	gStr			Given string.
*/
static unsigned WlzRecParseBox3D(WlzDBox3 *box, char *gStr)
{
  unsigned int	cmpMsk;
  double	buf[6];

  cmpMsk = WlzRecParseCSD(6, buf, gStr);
  if(cmpMsk && (1 << 0))
  {
    box->xMin = buf[0];
  }
  if(cmpMsk && (1 << 1))
  {
    box->yMin = buf[1];
  }
  if(cmpMsk && (1 << 2))
  {
    box->zMin = buf[2];
  }
  if(cmpMsk && (1 << 3))
  {
    box->xMax = buf[3];
  }
  if(cmpMsk && (1 << 4))
  {
    box->yMax = buf[4];
  }
  if(cmpMsk && (1 << 5))
  {
    box->zMax = buf[5];
  }
  return(cmpMsk);
}

/*!
* \return	Integer with bits set to indicate which components of the
* 		vertex are in the string:
* 		<ul>
* 		<li> (1<<0) x
* 		<li> (1<<1) y
* 		<li> (1<<2) z
* 		</ul>
*
* \brief	Parses the given string for a floating point vertex
* 		with the following format
* 		  x,y,z
*		Any of the components of the vertex may be missing
*		and the return value indicates which components
*		have been parsed.
* \param	vtx			Vertex in which values will be set.
* \param	gStr			Given string.
*/
static unsigned WlzRecParseVtx3D(WlzDVertex3 *vtx, char *gStr)
{
  unsigned int	cmpMsk;
  double	buf[3];

  cmpMsk = WlzRecParseCSD(3, buf, gStr);
  if(cmpMsk && (1 << 0))
  {
    vtx->vtX = buf[0];
  }
  if(cmpMsk && (1 << 1))
  {
    vtx->vtY = buf[1];
  }
  if(cmpMsk && (1 << 2))
  {
    vtx->vtZ = buf[2];
  }
  return(cmpMsk);
}

/*!
* \return	Integer with bits set to indicate which elements of the
* 		buffer are in the string:
* 		<ul>
* 		<li> (1<<0) v0
* 		<li> (1<<1) v1
* 		<li> (1<<2) v2
* 		</ul>
*
* \brief	Parses the given string for coma seperated floating point
* 		values with the following format
* 		  v0,v1,v2,v3,v4...vn
*		Any of the values may be missing and the return value
*		indicates which components have been parsed.
* \param	nVal			Number of values to attempt parsing.
* \param	val			Buffer of values.
* \param	gStr			Given string.
*/
static unsigned WlzRecParseCSD(int nVal, double *val, char *gStr)
{
  int		idx,
		ok = 0;
  unsigned int	cmpMsk = 0;
  char		*cmpStr0,
		*cmpStr1;

  if(gStr)
  {
    cmpStr0 = gStr;
    while(*cmpStr0 && isspace(*cmpStr0))
    {
      ++cmpStr0;
    }
    ok = *cmpStr0 != 0;
  }
  idx = 0;
  while(ok && (idx < (nVal - 1)))
  {
    ok = ((cmpStr1 = strchr(cmpStr0, ',')) != NULL) && (cmpStr1 != '\0');
    if(ok)
    {
      *cmpStr1 = '\0';
      if(cmpStr0 != cmpStr1)
      {
        ok = sscanf(cmpStr0, "%lg", val + idx) == 1;
        if(ok)
	{
	  cmpMsk |= (1 << idx);
	}
      }
      cmpStr0 = cmpStr1 + 1;
    }
     ++idx;
  }
  if(ok)
  {
    if(*cmpStr0 != '\0')
    {
      ok = sscanf(cmpStr0, "%lg", val + idx) == 1;
      if(ok)
      {
	cmpMsk |= (1 << idx);
      }
    }
  }
  if(ok == 0)
  {
    cmpMsk = 0;
  }
  return(cmpMsk);
}

/*!
* \return	Non zero if the file doesn't existi or can be overwritten.
* \brief	Check if a file with the given file name exists and if
* 		it does, asks for confirmation to over write the file.
* \param	app			Application name.
* \param	fName			Given file name with path.
*/
static	int	WlzRecConfirmFileOverwrite(char *app, char *fName)
{
  int		res = 1;
  struct stat	buf;

  if(stat(fName, &buf) != ENOENT)
  {
    do
    {
      (void )fprintf(stderr, "\n%s: Overwite file %s (y/n)", app, fName);
      res = getchar();
      switch(res)
      {
	case 'Y':
	  res = 'y';
	  break;
	case 'N':
	  res = 'n';
	  break;
	default:
	  break;
      }
    } while((res != 'y') && (res != 'n'));
    res = res == 'y';
  }
  return(res);
}

/*!
* \return	Reconstruct error code.
* \brief	Creates an absolute path from the given file name.
* \param	path		Pointer to absolute path string which
* 				may be reallocated.
* \param	file		File string which may be either absolute or
* 				relative.
*/
static RecError	WlzRecMakeFilePathAbs(char **path, char *file)
{
  RecError	recErr = REC_ERR_NONE;
  char		cwd[PATH_MAX + 1];
  const char	pathSep = '/';

  if(file == NULL)
  {
    recErr = REC_ERR_WRITE;
  }
  else
  {
    /* Make sure the output object file name has an absolute path. */
    if(*file == pathSep)
    {
      if((*path = AlcStrDup(file)) == NULL)
      {
	recErr = REC_ERR_MALLOC;
      }
    }
    else
    {
      if(getcwd(cwd, PATH_MAX) == NULL)
      {
	recErr = REC_ERR_WRITE;
      }
      else
      {
	AlcFree(*path);
	*path = (char *)AlcMalloc(sizeof(char) *
		                  (strlen(cwd) + strlen(file) + 2));
	if(*path == NULL)
	{
	  recErr = REC_ERR_MALLOC;
	}
      }
      if(recErr == REC_ERR_NONE)
      {
	if(*file == '.')
	{
	  (void )sprintf(*path, "%s%s", cwd, file + 1);
	}
	else
	{
	  (void )sprintf(*path, "%s%c%s", cwd, pathSep, file);
	}
      }
    }
  }
  return(recErr);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
