#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _Wlz3DViewTransformObj_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/Wlz3DViewTransformObj.c
* \author       Richard Baldock
* \date         November 1999
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
* \brief	Transforms a section view to a 3D object.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlz3dviewtransformobj "Wlz3DViewTransformObj"
*/

/*!
\ingroup BinWlzApp
\defgroup wlz3dviewtransformobj Wlz3DViewTransformObj
\par Name
Wlz3DViewTransformObj - transforms a section view to a 3D object.
\par Synopsis
\verbatim
Wlz3DViewTransformObj  [-h] [-v] [-D]
                 [-a <pitch,yaw[,roll]>] [-f <fx,fy,fz>]
                 [-d <dist>] [-b <view bibfile>] [-m <mode>]
		 [-u<ux,uy,uz>] [<3D object input file>]
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
    <td>Viewing angles in degrees, defined (0.0, 0.0,0.0).</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Bibfile defining the view parameters e.g. from MAPaint
        or warp input I/O.
	Override all other parameter input.</td>
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
    <td><b>-m</b></td>
    <td>Viewing mode, possible values:
    <table width="500" border="0">
      <tr> <td><b>Parameter value</b></td> <td><b>Viewing mode</b></td> </tr>
      <tr> <td>"up-is-up"</td> <td>up-is-up, default</td> </tr>
      <tr> <td>"statue"</td> <td>statue</td> </tr>
      <tr> <td>"absolute"</td> <td>absolute</td> </tr>
    </table>
    </td>
  </tr>
  <tr> 
    <td><b>-u</b></td>
    <td>Up vector for up-is-up mode, default (0, 0, -1).</td>
  </tr>
  <tr>
    <td><b>-D</b></td>
    <td>Command line viewing angles are degrees and not radians.</td>
  </tr>
</table>
\par Description
Transforms a section view to a 3D object
writing the 3D object to standard output.
\par Examples
\verbatim
\endverbatim
\par File
\ref Wlz3DViewTransformObj.c "Wlz3DViewTransformObj.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
\ref Wlz3DViewTransformObj "Wlz3DViewTransformObj(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <Wlz.h>
#include <bibFile.h>

/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  (void )fprintf(stderr,
	  "Usage:\t%s [-t<pitch,yaw[,roll]>] [-b<bibfile>]\n"
	  "\t[-f<fx,fy,fz>] [-d<dist>] [-D]"
	  " [-h] [-m<mode>] [<2D object input file>]\n"
	  "\tTransform an section view to a 3D object\n"
	  "\twriting the 3D object to standard output\n"
	  "Version: %s\n"
	  "Options:\n"
	  "\t  -a<pitch,yaw[,roll]> viewing angles - default 0,0,0\n"
	  "\t  -b<view-bibfile>   input parameters from the view\n"
	  "\t                     bibfile - e.g. saved from MAPaint\n"
	  "\t  -f<fx,fy,fz>       fixed point position - default zero\n"
	  "\t  -d<dist>           distance parameter - default zero\n"
	  "\t  -m<mode>           viewing mode, one of: up-is-up, statue,\n"
	  "\t                     absolute\n"
	  "\t  -u<ux,uy,uz>       up vector - default (0.0, 0.0, 1.0)\n"
	  "\t  -D                 Command line viewing angles are degrees\n"
	  "\t                     not radians.\n"
	  "\t  -h                 Help - prints this usage message\n"
	  "\t  -v                 Verbose operation\n",
	  proc_str,
	  WlzVersion());
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  WlzObject		*obj, *nobj;
  FILE			*inFile;
  char 			optList[] = "a:b:d:f:m:hvD";
  int			option,
  			useDegrees = 0,
			anglesSetOnCmdLine = 0;
  double		dist=0.0, theta=0.0, phi=0.0, zeta=0.0;
  WlzDVertex3		fixed={0.0,0.0,0.0};
  WlzThreeDViewStruct	*viewStr=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  int			verboseFlg=0;
  WlzThreeDViewMode	viewMode=WLZ_UP_IS_UP_MODE;
  double		scale=1.0;
  WlzDVertex3		upVectorVtx={0.0,0.0,1.0};
  char			*errMsg;
    
  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'a':
      anglesSetOnCmdLine = 1;
      if( sscanf(optarg, "%lg,%lg,%lg", &phi, &theta, &zeta) < 2 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'b':
      if((inFile = fopen(optarg, "r")) != NULL){
	char		viewModeStr[64];
	int		numParsedFields=0;
	BibFileRecord	*bibfileRecord;
	BibFileField	*bibfileField;
	BibFileError	bibFileErr;

	/* read the bibfile - get the first section view entry */
	bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
	while((bibFileErr == BIBFILE_ER_NONE) &&
	      (strncmp(bibfileRecord->name, "Wlz3DSectionViewParams", 22))){
	  BibFileRecordFree(&bibfileRecord);
	  bibFileErr = BibFileRecordRead(&bibfileRecord, &errMsg, inFile);
	}
	(void) fclose(inFile);
	if( bibFileErr != BIBFILE_ER_NONE ){
	  fprintf(stderr, "%s: can't read bibfile: %s\n", argv[0], optarg);
	  usage(argv[0]);
	  return 1;
	}

	/* parse the record */
        anglesSetOnCmdLine = 0;
	numParsedFields = BibFileFieldParseFmt
	  (bibfileRecord->field,
	   (void *) &(fixed.vtX), "%lg ,%*lg ,%*lg", "FixedPoint",
	   (void *) &(fixed.vtY), "%*lg ,%lg ,%*lg", "FixedPoint",
	   (void *) &(fixed.vtZ), "%*lg ,%*lg ,%lg", "FixedPoint",
	   (void *) &dist, "%lg", "Distance",
	   (void *) &phi, "%lg", "Pitch",
	   (void *) &theta, "%lg", "Yaw",
	   (void *) &zeta, "%lg", "Roll",
	   (void *) &scale, "%lg", "Scale",
	   (void *) &(upVectorVtx.vtX), "%lg ,%*lg ,%*lg", "UpVector",
	   (void *) &(upVectorVtx.vtY), "%*lg ,%lg ,%*lg", "UpVector",
	   (void *) &(upVectorVtx.vtZ), "%*lg ,%*lg ,%lg", "UpVector",
	   (void *) &(viewModeStr[0]), "%s", "ViewMode",
	   NULL);
	if(numParsedFields < 1) {
	  (void )fprintf(stderr,
	                 "%s: can't read section transform from bibfile: %s\n",
	                 argv[0], optarg);
	  usage(argv[0]);
	  return 1;
	}
	/* doesn't read the view mode correctly - ask Bill */
	bibfileField = bibfileRecord->field;
	while( bibfileField ){
	  if( strncmp(bibfileField->name, "ViewMode", 8) == 0 ){
	    strcpy(viewModeStr, bibfileField->value);
	    break;
	  }
	  bibfileField = bibfileField->next;
	}
	BibFileRecordFree(&bibfileRecord);

	/* convert angles to radians and get mode */
	phi *= (WLZ_M_PI/180.0);
	theta *= (WLZ_M_PI/180.0);
	zeta *= (WLZ_M_PI/180.0);
	if( !strncmp(viewModeStr, "up-is-up", 8) ){
	  viewMode = WLZ_UP_IS_UP_MODE;
	}
	else if( !strncmp(viewModeStr, "statue", 8) ){
	  viewMode = WLZ_STATUE_MODE;
	}
	else {
	  viewMode = WLZ_ZETA_MODE;
	}


      }
      else {
	fprintf(stderr, "%s: can't open bibfile: %s\n", argv[0], optarg);
	usage(argv[0]);
	return 1;
      }
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
      if( strcmp(optarg, "up-is-up") == 0 ){
	viewMode = WLZ_UP_IS_UP_MODE;
      }
      else if( strcmp(optarg, "statue") == 0 ){
	viewMode = WLZ_STATUE_MODE;
      }
      else if( strcmp(optarg, "absolute") == 0 ){
	viewMode = WLZ_ZETA_MODE;
      }
      else {
	fprintf(stderr, "%s: invalid view mode: %s\n", argv[0], optarg);
	usage(argv[0]);
	return 1;
      }
      break;

    case 'u':
      if( sscanf(optarg, "%lg,%lg,%lg", &(upVectorVtx.vtX), &(upVectorVtx.vtY),
		 &(upVectorVtx.vtZ)) < 3 ){
	usage(argv[0]);
	return 1;
      }
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'D':
      useDegrees = 1;
      break;
    case 'h':
    default:
      usage(argv[0]);
      return 1;

    }
  }

  inFile = stdin;
  if( optind < argc ){
    if( (inFile = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      usage(argv[0]);
      return 1;
    }
  }

  /* create the view structure */
  if( errNum == WLZ_ERR_NONE ){
    if(useDegrees && anglesSetOnCmdLine) {
      const double	dtr = WLZ_M_PI / 180.0;

      theta *= dtr;
      phi *= dtr;
      zeta *= dtr;
    }
    if((viewStr = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum)) != NULL){
      viewStr->fixed 	= fixed;
      viewStr->theta 	= theta;
      viewStr->phi 	= phi;
      viewStr->zeta 	= zeta;
      viewStr->dist 	= dist;
      viewStr->scale 	= 1.0;
      viewStr->view_mode = viewMode;
      viewStr->up 	= upVectorVtx;
    }
    else {
      fprintf(stderr, "%s: can't create view structure\n", argv[0]);
      usage(argv[0]);
      return 1;
    }
  }

  if(verboseFlg) {
    (void )fprintf(stderr,
		   "%s: Section parameters are\n"
		   "  fixed = %lg,%lg,%lg\n"
		   "  theta = %lg\n"
		   "  phi = %lg\n"
		   "  zeta = %lg\n"
		   "  dist = %lg\n"
		   "  scale = %lg\n"
		   "  viewMode = %s\n"
		   "  up = %lg,%lg,%lg\n",
		   argv[0],
		   fixed.vtX, fixed.vtY, fixed.vtZ,
		   theta,
		   phi,
		   zeta,
		   dist,
		   scale,
		   WlzStringFromThreeDViewMode(viewMode, NULL),
		   upVectorVtx.vtX, upVectorVtx.vtY, upVectorVtx.vtZ);
  }

  /* read objects and section if possible */
  while((obj = WlzReadObj(inFile, &errNum)) != NULL) 
  {
    obj = WlzAssignObject(obj, &errNum);
    switch( obj->type )
    {
      case WLZ_2D_DOMAINOBJ:
	nobj = NULL;
	/* initialise the view structure and transform the section */
	WlzInit3DViewStruct(viewStr, obj);
	nobj = Wlz3DViewTransformObj(obj, viewStr, &errNum);

	if( nobj != NULL){
	  WlzWriteObj(stdout, nobj);
	}
	else {
	  return errNum;
	}
	WlzFreeObj(nobj);
	break;

      default:
	WlzWriteObj(stdout, obj);
	break;
    }

    WlzFreeObj(obj);
  }
  WlzFree3DViewStruct(viewStr);

  return WLZ_ERR_NONE;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
