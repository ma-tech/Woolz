#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _Wlz3DGetProjection_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlzApp/Wlz3DGetProjection.c
* \author       Richard Baldock
* \date         March 2005
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
* \brief        Command line binary to project a Woolz object
*		using a view transform.
*               
* \ingroup	BinWlzApp
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlz3dgetprojection "Wlz3DGetProjection"
*/

/*!
\ingroup BinWlzApp
\defgroup wlz3dgetprojection Wlz3DGetProjection
\par Name
Wlz3DGetProjection - projects a Woolz object using a view transform.
\par Synopsis
\verbatim
Wlz3DGetProjection [-h] [-v]
                   [-a <pitch,yaw[,roll]>] [-f <fx,fy,fz>] [-d <dist>]
		   [-b <parameter bibfile>] [-m <mode>] [-s <scale>]
		   [-o <output file>] [-u<ux,uy,uz>] [<3D object input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose operation.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Viewing angles in degrees.
        If roll is defined then the mode is "absolute"..</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Fixed point position, default (0,0,0).</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Distance parameter, default 0.0.</td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>Bibfile defining the view parameters e.g. from MAPaint
        or warp input I/O.
	Override all other parameter input.</td>
  </tr>
  <tr> 
    <td><b>-m</b></td>
    <td>Viewing mode, possible values:
    <table width="500" border="0">
      <tr> <td><b>Parameter value</b></td> <td><b>Viewing mode</b></td> </tr>
      <tr> <td>0</td> <td>up-is-up, default</td> </tr>
      <tr> <td>1</td> <td>statue</td> </tr>
      <tr> <td>2</td> <td>absolute</td> </tr>
    </table>
    </td>
  </tr>
  <tr> 
    <td><b>-s</b></td>
    <td>Scale factor, default 1.0.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output filename, default is stdout.</td>
  </tr>
  <tr> 
    <td><b>-u</b></td>
    <td>Up vector for up-is-up mode, default (0, 0, -1).</td>
  </tr>
</table>
\par Description
Gets an arbitrary slice projection from a 3D object,
writing the 2D object to standard output.
\par Examples
\verbatim
\endverbatim
\par File
\ref Wlz3DGetProjection.c "Wlz3DGetProjection.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref WlzGetProjectionFromObject "WlzGetProjectionFromObject(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <Wlz.h>
#include <WlzExtFF.h>
#include <bibFile.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s [-a <pitch,yaw[,roll]>] [-f <fx,fy,fz>] [-d <dist>]"
	  " [-b <parameter bibfile>] [-m <mode>] [-s <scale>]"
	  " [-o <output file>] [-u<ux,uy,uz>] [-h] [-v]"
	  " [<3D object input file>]\n"
	  "\tGet an arbitrary sliceprojection from a 3D object\n"
	  "\twriting the 2D object to standard output\n"
	  "\tOptions are:\n"
	  "\t  -a<pitch,yaw[,roll]> viewing angles in degrees. If roll\n"
	  "\t                       is defined then the mode is \"absolute\"\n"
	  "\t  -b<bibfile>        bibfile defining the view parameters e.g.\n"
	  "\t                     from MAPaint view or warp input I/O\n"
	  "\t                     Override all other parameter input\n"
	  "\t  -f<fx,fy,fz>       fixed point position, default - (0,0,0)\n"
	  "\t  -d<dist>           distance parameter, default - 0.0\n"
	  "\t  -m<mode>           viewing mode, possible values:\n"
	  "\t                       mode = 0 - up-is-up (default)\n"
	  "\t                       mode = 1 - statue\n"
	  "\t                       mode = 2 - absolute\n"
	  "\t  -o<output file>    Output filename, default to stdout\n"
	  "\t  -s<scale>          Scale factor, default - 1.0\n"
	  "\t  -u<ux,uy,uz>       Up vector for up-is-up mode.\n"
	  "\t			  Default: (0,0,-1)\n"
	  "\t  -h                 Help - prints this usage message\n"
	  "\t  -v                 verbose operation\n"
	  "",
	  proc_str);
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject	*obj, *nobj;
  FILE		*inFP, *outFP, *bibFP;
  char		*outFile, *bibFile;
  char 		optList[] = "a:b:d:f:m:o:s:u:hv";
  int		option;
  int		iVal;
  double	dist=0.0, pitch=0.0, yaw=0.0, roll=0.0;
  double	scale=1.0;
  WlzDVertex3	fixed={0.0,0.0,0.0};
  WlzDVertex3	up={0.0,0.0,-1.0};
  WlzThreeDViewStruct	*viewStr=NULL;
  WlzThreeDViewMode mode=WLZ_UP_IS_UP_MODE;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr;
  int		verboseFlg=0;
  char		*errMsg;
    
  /* additional defaults */
  outFile = "-";
  bibFile = NULL;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'a':
      switch( sscanf(optarg, "%lg,%lg,%lg", &pitch, &yaw, &roll) ){
      default:
	usage(argv[0]);
	return 1;
      case 2:
	break;
      case 3:
	mode = WLZ_ZETA_MODE;
	break;
      }
      break;

    case 'b':
      bibFile = optarg;
      break;

    case 'd':
      if( sscanf(optarg, "%lg", &dist) < 1 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'f':
      if( sscanf(optarg, "%lg,%lg,%lg", &(fixed.vtX), &(fixed.vtY),
		 &(fixed.vtZ)) < 3 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'm':
      if( sscanf(optarg, "%d", &iVal) < 1 ){
	usage(argv[0]);
	return 1;
      }
      else if( mode != WLZ_ZETA_MODE ){
	switch( iVal ){
	default:
	  usage(argv[0]);
	  return 1;
	case 0:
	  mode = WLZ_UP_IS_UP_MODE;
	  break;
	case 1:
	  mode = WLZ_STATUE_MODE;
	  break;
	case 2:
	  mode = WLZ_ZETA_MODE;
	  break;
	}
      }
      break;

    case 'o':
      outFile = optarg;
      break;

    case 's':
      if( sscanf(optarg, "%lg", &scale) < 1 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'u':
      if( sscanf(optarg, "%lg,%lg,%lg", &(up.vtX), &(up.vtY),
		 &(up.vtZ)) < 3 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'h':
    default:
      usage(argv[0]);
      return 0;

    case 'v':
      verboseFlg = 1;
      break;

    }
  }

  /* check input file/stream */
  inFP = stdin;
  if( optind < argc ){
    if( (inFP = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* check output file/stream */
  if(strcmp(outFile, "-"))
  {
    if((outFP = fopen(outFile, "w")) == NULL)
    {
      errNum = WLZ_ERR_WRITE_EOF;
    }
  }
  else
  {
    outFP = stdout;
  }

  /* create view structure */
  if((viewStr = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum)) != NULL){
    viewStr->theta = yaw * WLZ_M_PI / 180.0;
    viewStr->phi = pitch * WLZ_M_PI / 180.0;
    viewStr->zeta = roll * WLZ_M_PI / 180.0;
    viewStr->dist = dist;
    viewStr->fixed = fixed;
    viewStr->up = up;
    viewStr->view_mode = mode;
    viewStr->scale = scale;
  }

  /* check bibfile - select first section parameters in the file */
  if((errNum == WLZ_ERR_NONE) && (bibFile != NULL)){
    if((bibFP = fopen(bibFile, "r")) != NULL){
      bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, bibFP);
      while((bibFileErr == BIBFILE_ER_NONE) &&
	    (strncmp(bibfileRecord->name, "Wlz3DSectionViewParams", 22))){
	BibFileRecordFree(&bibfileRecord);
	bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, bibFP);
      }
      fclose( bibFP );
      if( bibFileErr != BIBFILE_ER_NONE ){
	fprintf(stderr, "%s: error reading bibfile: %s\n", argv[0], errMsg);
	AlcFree((void *) errMsg);
	return 1;
      }
      WlzEffBibParse3DSectionViewParamsRecord(bibfileRecord, viewStr);
      BibFileRecordFree(&bibfileRecord);
    }
    else {
      fprintf(stderr, "%s: can't open parameter bibfile %s\n", argv[0],
	      bibFile);
      return 1;
    }
  }
  

  /* read objects and section if possible */
  while((errNum == WLZ_ERR_NONE) &&
        ((obj = WlzReadObj(inFP, &errNum)) != NULL))
  {
    switch( obj->type )
    {
    case WLZ_3D_DOMAINOBJ:
      WlzInit3DViewStruct(viewStr, obj);
      nobj = WlzGetProjectionFromObject(obj, viewStr, NULL, NULL, &errNum);
      if( nobj != NULL){
	WlzWriteObj(outFP, nobj);
      }
      else {
	return errNum;
      }
      WlzFreeObj(nobj);
      break;

    default:
	WlzWriteObj(outFP, obj);
	break;
    }

    WlzFreeObj(obj);
  }
  if(errNum == WLZ_ERR_READ_EOF)
  {
    errNum = WLZ_ERR_NONE;
  }

  return errNum;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
