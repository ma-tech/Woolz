#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _Wlz3DGetProjection_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/Wlz3DGetProjection.c
* \author       Richard Baldock, Bill Hill
* \date         March 2005
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
* \ingroup	BinWlzApp
* \brief        Command line binary to project a Woolz object
*		using a view transform.
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
Wlz3DGetProjection [-h]
                   [-a <pitch,yaw[,roll]>] [-b <parameter bibfile>]
		   [-d <dist>] [-f <fx,fy,fz>] [-i <int mod>]
		   [-m <mode>] [-o <output file>] [-r <vox rescale>]
		   [-s <scale>] [-u<ux,uy,uz>] [-D <den>] [-V <val lut>]
		   [<3D object input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
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
    <td><b>-r</b></td>
    <td>Voxel size rescaling mode flags:
         bit 1 set - use voxel-size rescaling,
         bit 2 set - enable global scaling,
        default 1.0.</td>
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
  <tr>
    <td><b>-t</b></td>
    <td>Depth or thickness of the projection,
        Default: 0.0 (implies whole volume)</td>
  </tr>
  <tr>
    <td><b>-i</b></td>
    <td>Integration mode, possible values:
    <table width="500" border="0">
    <tr> <td>n</td><td>none</td><td>shadow domain</td> </tr>
    <tr> <td>d</td><td>domain</td><td>uniform domain (default)</td> </tr>
    <tr> <td>v</td><td>values</td><td>value integration</td> </tr>
    </table>
    </td>
  </tr>
  <tr>
    <td><b>-D</b></td>
    <td>Domain density for use with uniform domain density integration
        mode, range [0-255] (default 255).</td>
  </tr>
  <tr>
    <td><b>-V</b></td>
    <td>Value to density look up table for use with value integration
        mode, all 256 entries to have range [0-255] (default identity).</td>
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
\ref WlzProjectObjToPlane "WlzProjectObjToPlane(3)"
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

static void ShowUsage(char *proc_str)
{
  (void )fprintf(stderr,
  "Usage:\t%s [-a <pitch,yaw[,roll]>] [-f <fx,fy,fz>] [-d <dist>]\n"
  "\t[-b <parameter bibfile>] [-m <mode>] [-s <scale>]\n"
  "\t[-o <output file>] [-u<ux,uy,uz>] [-r<vox rescale>]\n"
  "\t[-i<int mod> [-D <den>] [-V <val lut>] [-t<depth>] [-h]\n"
  "\t[<3D object input file>]\n"
  "\tGet an arbitrary sliceprojection from a 3D object\n"
  "\twriting the 2D object to standard output\n"
  "Version: %s\n"
  "Options:\n"
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
      "\t  -r<vox rescale>    Voxel size rescaling mode flags:\n"
      "\t                       bit 1 set - use voxel-size rescaling\n"
      "\t                       bit 2 set - enable global scaling\n"
      "\t                     default - 0.\n"
      "\t  -s<scale>          Scale factor, default - 1.0\n"
      "\t  -t<depth>          Depth or thickness of the projection,\n"
      "\t                         Default: 0.0 (implies whole volume)\n"
      "\t  -u<ux,uy,uz>       Up vector for up-is-up mode.\n"
      "\t			  Default: (0,0,-1)\n"
      "\t  -i<int mode>	      Integration mode, possible values:\n"
      "\t                       mode = n - none, shadow domain\n"
      "\t                       mode = d - domain, uniform domain (default)\n"
      "\t                       mode = v - value, value integration\n"
      "\t  -D                 Domain density for use with uniform domain\n"
      "\t                     density integration mode, range [0-255]\n"
      "\t  -V                 Look up table for use with value integration\n"
      "\t                     mode, default is identity\n"
      "\t  -h                 Help - prints this usage message\n",
  proc_str,
  WlzVersion());
  return;
}
 
int main(int	argc,
	 char	**argv)
{

  int		tI,
  		usage = 0,
		voxRescale = 0;
  WlzObject	*obj = NULL, *nobj = NULL;
  WlzProjectIntMode intMod = WLZ_PROJECT_INT_MODE_DOMAIN;
  WlzUByte	denDom = 255;
  WlzUByte	denVal[256];
  FILE		*inFP, *outFP, *bibFP;
  char		*outFile = NULL,
                *bibFile = NULL,
  		*lutFile = NULL;
  char 		optList[] = "a:b:d:D:f:i:m:o:r:s:t:u:V:h";
  int		option;
  int		iVal;
  double	depth=0.0, dist=0.0, pitch=0.0, yaw=0.0, roll=0.0;
  double	scale=1.0;
  WlzDVertex3	fixed={0.0,0.0,0.0};
  WlzDVertex3	up={0.0,0.0,-1.0};
  WlzThreeDViewStruct	*viewStr=NULL;
  WlzThreeDViewMode mode=WLZ_UP_IS_UP_MODE;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  BibFileRecord	*bibfileRecord;
  BibFileError	bibFileErr;
  char		*errMsg;
    
  /* additional defaults */
  outFile = "-";
  bibFile = NULL;

  /* read the argument list and check for an input file */
  opterr = 0;
  while(((option = getopt(argc, argv, optList)) != EOF) &&
        (usage == 0)){
    switch( option ){

    case 'a':
      switch( sscanf(optarg, "%lg,%lg,%lg", &pitch, &yaw, &roll) ){
      default:
	break;
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
	break;
      }
      break;

    case 'f':
      if( sscanf(optarg, "%lg,%lg,%lg", &(fixed.vtX), &(fixed.vtY),
		 &(fixed.vtZ)) < 3 ){
	break;
      }
      break;

    case 'm':
      if( sscanf(optarg, "%d", &iVal) < 1 ){
	break;
      }
      else if( mode != WLZ_ZETA_MODE ){
	switch( iVal ){
	default:
	  usage = 2;
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

    case 'r':
      if( sscanf(optarg, "%d", &voxRescale) < 1 ){
        usage = 2;
	}
	break;

    case 's':
      if( sscanf(optarg, "%lg", &scale) < 1 ){
	usage = 2;
      }
      break;

    case 't':
      if( sscanf(optarg, "%lg", &depth) < 1 ){
	break;
      }
      break;

    case 'u':
      if( sscanf(optarg, "%lg,%lg,%lg", &(up.vtX), &(up.vtY),
		 &(up.vtZ)) < 3 ){
	usage = 2;
      }
      break;

    case 'i':
      switch(*optarg)
      {
        case 'n':
	  intMod = WLZ_PROJECT_INT_MODE_NONE;
	  break;
	case 'd':
	  intMod = WLZ_PROJECT_INT_MODE_DOMAIN;
	  break;
	case 'v':
	  intMod = WLZ_PROJECT_INT_MODE_VALUES;
	  break;
	default:
	  usage = 2;
	  break;
      }
      break;

    case 'D':
      if((sscanf(optarg, "%d", &tI) < 1) || (tI < 0) || (tI > 255))
      {
        usage = 2;
      }
      else
      {
        denDom = tI;
      }
      break;

    case 'V':
      lutFile = optarg;
      break;

    case 'h':
    default:
      usage = 1;
    }
  }
  if(usage)
  {
    ShowUsage(argv[0]);
    return(usage - 1);
  }

  /* check input file/stream */
  inFP = stdin;
  if( optind < argc ){
    if( (inFP = fopen(*(argv+optind), "r")) == NULL ){
      fprintf(stderr, "%s: can't open file %s\n", argv[0], *(argv+optind));
      ShowUsage(argv[0]);
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

  /* check for values density look up table */
  if(intMod == WLZ_PROJECT_INT_MODE_VALUES)
  {
    if(lutFile)
    {
      WlzObject *lutObj = NULL;
      FILE 	*lutFP = NULL;

      if(((lutFP = (strcmp(lutFile, "-")?
		    fopen(lutFile, "r"): stdin)) == NULL) ||
		   ((lutObj = WlzAssignObject(
			      WlzReadObj(lutFP, &errNum), NULL)) == NULL))
      {
	(void )fprintf(stderr,
	    "%s: failed to read look up table object from file %s\n",
	    *argv, lutFile);
        ShowUsage(argv[0]);
	return(1);
      }
      if(lutFP && strcmp(lutFile, "-"))
      {
	fclose(lutFP);
      }
      if(lutObj)
      {
	int	lutOK = 1;

        if((lutObj->type != WLZ_LUT) ||
           (lutObj->domain.core == NULL) ||
           (lutObj->values.core == NULL) ||
           (lutObj->values.lut->vType != WLZ_GREY_INT) ||
	   (lutObj->domain.lut->bin1 != 0) ||
	   (lutObj->domain.lut->lastbin != 255))
	{
	  lutOK = 0;
	}
	else
	{
	  int	i;
	  int	*lut;

	  lut = lutObj->values.lut->val.inp;
	  for(i = 0; i < 256; ++i)
	  {
	    if((lut[i] < 0) || (lut[i] > 255))
	    {
	      lutOK = 0;
	      break;
	    }
	    else
	    {
	      denVal[i] = lut[i];
	    }
	  }
	}
	if(lutOK == 0)
        {
	  (void )fprintf(stderr,
	      "%s: values density look up table must have 256 integer\n"
	      "entries with each entry in the range [0-255]\n",
	      *argv);
	}
        (void )WlzFreeObj(lutObj);
      }
    }
    else
    {
      int	i;

      for(i = 0; i < 256; ++i)
      {
        denVal[i] = i;
      }
    }
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
    viewStr->voxelRescaleFlg = voxRescale;
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
        ((obj = WlzAssignObject(
	        WlzReadObj(inFP, &errNum), NULL)) != NULL))
  {
    switch( obj->type )
    {
    case WLZ_3D_DOMAINOBJ:
      if(voxRescale && obj->domain.core)
      {
        viewStr->voxelSize[0] = obj->domain.p->voxel_size[0];
        viewStr->voxelSize[1] = obj->domain.p->voxel_size[1];
        viewStr->voxelSize[2] = obj->domain.p->voxel_size[2];
      }
      WlzInit3DViewStruct(viewStr, obj);
      nobj = WlzProjectObjToPlane(obj, viewStr, intMod, denDom, denVal,
      				  depth, &errNum);
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
  if(inFP && (inFP != stdin)) {
    (void )fclose(inFP);
  }
  if(outFP && (outFP != stdout)) {
    (void )fclose(outFP);
  }
  (void )WlzFree3DViewStruct(viewStr);
  if(errNum == WLZ_ERR_READ_EOF)
  {
    errNum = WLZ_ERR_NONE;
  }
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
