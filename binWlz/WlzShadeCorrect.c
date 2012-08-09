#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzShadeCorrect_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzShadeCorrect.c
* \author       Richard Baldock
* \date         June 2004
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
* \brief	Shade corrects a domain object with values.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzshadecorrect "WlzShadeCorrect"
*/

/*!
\ingroup BinWlz
\defgroup wlzshadecorrect WlzShadeCorrect
\par Name
WlzShadeCorrect - shade corrects a domain object with values.
\par Synopsis
\verbatim
WlzShadeCorrect [-h] [-v] [-p] [-d <dark-field image>]
                -b <bright-field image> [<input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-v</b></td>
    <td>Verbose operation.</td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Patch image, shade images shifted to correct each patch.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Dark-field image (optional)..</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Bright-field image, <b>required</b>.</td>
  </tr>
</table>
\par Description
Shade corrects a domain object with values.
The shade using the bright-field and dark-field images
assumes intensity images and calculates a normalised
transmission coefficient.
\par Examples
\verbatim
WlzShadeCorrect bright.wlz slide.wlz >side_sc.wlz
\endverbatim
Reads objects from files bright.wlz and slide.wlz,
the bright-field and source objects.
Shade corrects the source object and then writes the shade corrected object
to the file side_sc.wlz.
\par File
\ref WlzShadeCorrect.c "WlzShadeCorrect.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzShadeCorrect "WlzShadeCorrect(3)"
\ref WlzShadeCorrectBFDF "WlzShadeCorrectBFDF(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>

#include <Wlz.h>


/* externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern int 	optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static void usage(char *proc_str)
{
  fprintf(stderr,
	  "Usage:\t%s -b <bright-field image> -d <dark-field image> "
	  "-p -h -v [<input file>]\n"
	  "\tApply a shade correction to the input objects. The shade\n"
	  "\tcorrection using the bright-field and dark-field images\n"
	  "\tassumes intensity images and calculates a normalised\n"
	  "\ttransmission coefficient.\n"
	  "Version: %s\n"
	  "Options:\n" 
	  "\t  -b <file> Bright-field image - REQUIRED\n"
	  "\t  -d <file> Dark-field image - optional\n"
	  "\t  -p        Patch image, shade images shifted to correct each patch\n"
	  "\t  -h        Help - prints this usage message\n"
	  "\t  -v        Verbose operation\n",
	  proc_str,
	  WlzVersion());
  return;
}

int main(int	argc,
	 char	**argv)
{
  char 		optList[] = "b:d:phv";
  int		option;
  int		verboseFlg=0;
  int		patchFlg=0;
  int		i, xShift, yShift, zShift;
  FILE		*fp;
  WlzObject	*inObj=NULL, *outObj, *bfObj=NULL, *dfObj=NULL;
  WlzObject	**outObjs, *obj1, *obj2;
  WlzCompoundArray	*cObj;
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* read the argument list and check for an input file */
  opterr = 0;
  while( (option = getopt(argc, argv, optList)) != EOF ){
    switch( option ){

    case 'b':
      if( (fp = fopen(optarg, "rb")) == NULL ){
	fprintf(stderr, "%s: can't open bright-field file %s\n",
		argv[0], argv[optind]);
	usage(argv[0]);
	return WLZ_ERR_UNSPECIFIED;
      }
      else {
	bfObj = WlzReadObj(fp, &errNum);
	fclose(fp);
      }
      break;

    case 'd':
      if( (fp = fopen(optarg, "rb")) == NULL ){
	fprintf(stderr, "%s: can't open dark-field file %s\n",
		argv[0], argv[optind]);
	usage(argv[0]);
	return WLZ_ERR_UNSPECIFIED;
      }
      else {
	dfObj = WlzReadObj(fp, &errNum);
	fclose(fp);
	break;
      }
      break;

    case 'p':
      patchFlg = 1;
      break;

    case 'v':
      verboseFlg = 1;
      break;

    case 'h':
    default:
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }

  /* bright field must have been defined */
  if( bfObj == NULL ){
    fprintf(stderr, "%s: bright field image required\n",
	    argv[0]);
    usage(argv[0]);
    return WLZ_ERR_UNSPECIFIED;
  }

  /* get the input stream */
  if( optind < argc ){
    if( (fp = fopen(*(argv+optind), "rb")) == NULL ){
      fprintf(stderr, "%s: can't open input file %s\n",
	      argv[0], argv[optind]);
      usage(argv[0]);
      return WLZ_ERR_UNSPECIFIED;
    }
  }
  else {
    fp = stdin;
  }

  /* now read the objects and shade correct */
  while((inObj = WlzReadObj(fp, &errNum)) != NULL){
    switch( inObj->type ){

    case WLZ_2D_DOMAINOBJ:
    case WLZ_3D_DOMAINOBJ:
      if((outObj = WlzShadeCorrectBFDF(inObj, bfObj, dfObj,
				       255.0, 0, &errNum)) != NULL){
	WlzWriteObj(stdout, outObj);
	WlzFreeObj(outObj);
      }
      else {
	fprintf(stderr, "%s: shade correction failed\n", argv[0]);
	return errNum;
      }
      break;

    case WLZ_COMPOUND_ARR_1:
    case WLZ_COMPOUND_ARR_2:
      if( patchFlg ){
	/* shade correct each patch, shifting the shade image
	   as required */
	cObj = (WlzCompoundArray *) inObj;
	if((outObjs = AlcMalloc(sizeof(WlzObject *) * cObj->n)) != NULL){
	  for(i=0; (i < cObj->n) && (errNum == WLZ_ERR_NONE); i++){
	    /* could test for 3D here */
	    xShift = cObj->o[i]->domain.i->kol1 - bfObj->domain.i->kol1;
	    yShift = cObj->o[i]->domain.i->line1 - bfObj->domain.i->line1;
	    zShift = 0;
	    obj1 = WlzShiftObject(bfObj, xShift, yShift, zShift,
				  &errNum);
	    if( dfObj ){
	      obj2 = WlzShiftObject(dfObj, xShift, yShift, zShift,
				    &errNum);
	    }
	    else {
	      obj2 = NULL;
	    }
	    outObjs[i] = WlzShadeCorrectBFDF(cObj->o[i], obj1, obj2,
					     255.0, 0, &errNum);
	    WlzFreeObj(obj1);
	    if( obj2 ){
	      WlzFreeObj(obj2);
	    }
	  }
	  if((outObj = (WlzObject *)
	     WlzMakeCompoundArray(cObj->type, 3, cObj->n,
				  &(outObjs[0]), cObj->o[0]->type,
				  &errNum)) != NULL){
	    WlzWriteObj(stdout, outObj);
	    WlzFreeObj(outObj);
	    AlcFree(outObjs);
	  }
	}
	else {
	  errNum = WLZ_ERR_MEM_ALLOC;
	}
      }
      else if((outObj = WlzShadeCorrectBFDF(inObj, bfObj, dfObj,
				       255.0, 0, &errNum)) != NULL){
	WlzWriteObj(stdout, outObj);
	WlzFreeObj(outObj);
      }
      else {
	fprintf(stderr, "%s: shade correction failed\n", argv[0]);
	return errNum;
      }
      break;

    default:
      WlzWriteObj(stdout, inObj);
      break;
    }
    WlzFreeObj(inObj);
  }

  if( bfObj ){
    WlzFreeObj(bfObj);
  }
  if( dfObj ){
    WlzFreeObj(dfObj);
  }

  return 0;
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
