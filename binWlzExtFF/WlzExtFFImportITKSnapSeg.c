#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzExtFFImportITKSnapSeg_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         WlzExtFFImportITKSnapSeg.c
* \author       Bill Hill
* \date         July 2015
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               MRC Institute of Genetics and Molecular Medicine,
*               University of Edinburgh,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C), [2015],
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
* \brief	Imports an ITKSnap segmentation to Woolz.
* \ingroup	BinWlzExtFF
*
* \par Binary
* \ref wlzextffimportitksnapseg "WlzExtFFImportITKSnapSeg"
*/

/*!
\ingroup BinWlzExtFF
\defgroup wlzextffimportitksnapseg WlzExtFFImportITKSnapSeg
\par Name
WlzExtFFImportITKSnapSeg - imports an ITKSnap segmentation to Woolz.
\par Synopsis
\verbatim
WlzExtFFImportITKSnapSeg [-h] [-o<output file>] [-f<image format>]
			 [-l<label file>] [-b] [<segmentation image>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
  <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr>
    <td><b>-f</b></td>
    <td>Source file format. With simple file names, the format can be
        inferred from the file extension.</td>
  </tr>
  <tr>
    <td><b>-l</b></td>
    <td>ITKSnap segmentation label description file (for stdin use -).</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr>
    <td><b>-b</b></td>
    <td>Scub out non-alphanumeric characters replacing them with an
	underscore in the label name property.</td>
  </tr>
</table>
Reads a segmentation image and an optional ITKSnap segmentation label
description file, then uses these to create a compound object with
domains indexed as for the input image labels. If the label description file
is provided then the objects will have their name and grey (for colour)
properties set.
By default files are read from the standard input (image before labels
if the standard input is specified for the labels as well as the image).
\par Examples
\verbatim
WlzExtFFImportITKSnapSeg -l seg.txt seg.nii | WlzExplode -p
\endverbatim
Reads the segmentation image from seg.nii using the NIfTI format
(inferred from file extension) and the label descriptions in
ITKSnap text format from seg.txt.
These are then used to create a compound object which is written via a pipe
to WlzExplode. WlzExplode creates seperate Woolz domain files,
each with the file body of the label description text.
\par File
WlzExtFFImportITKSnapSeg.c "WlzExtFFImportITKSnapSeg.c"
\par See Also
\ref BinWlzExtFF "WlzIntro(1)"
\ref wlzextffconvert "WlzExtFFConvert(1)"
\ref wlzexplode "WlzExplode(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <Wlz.h>
#include <WlzExtFF.h>

#define WLZEXTFF_ITK_LN_MAX (256)

static WlzPropertyList 		*WlzEffITKSnapPList(
				  AlcVector *vec,
				  int *idS,
				  int idL,
				  WlzErrorNum *dstErr);

typedef struct _WlzEffITKSnapSegRecord
{
  int		idx;
  WlzUInt	colour;
  int		vis;
  int		lmv;
  char 		*label;
} WlzEffITKSnapSegRecord;

extern int      		getopt(
				  int argc,
				  char * const *argv,
				  const char *optstring);

extern char     *optarg;
extern int      optind,
                opterr,
                optopt;

int             		main(
				  int argc,
				  char **argv)
{
  int           option,
                ok = 1,
                usage = 0,
		dim = 0,
		scrub = 0,
		maxLabelIdx = 0;
  AlcVector	*textVec = NULL,
  		*labelVec = NULL;
  WlzEffFormat  inFmt = {WLZEFF_FORMAT_NIFTI};
  WlzObject	*inObj = NULL,
  	        *outObj = NULL;
  char		*outFileStr = NULL,
  		*inImgFileStr = NULL,
  		*inLabFileStr = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;
  const char    *errMsg;
  const char    * const itkSnapMagic = "#ITK-SnAPLabelDescriptionFile";
  static char   optList[] = "bhf:l:o:",
                fileStrDef[] = "-";

  opterr = 0;
  inImgFileStr = fileStrDef;
  outFileStr = fileStrDef;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
	scrub = 1;
	break;
      case 'f':
	if((inFmt = WlzEffStringExtToFormat(optarg)) == 0)
	{
	  usage = 1;
	}
	break;
      case 'l':
	inLabFileStr = optarg;
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'h': /* FALLTHROUGH */
      default:
	usage = 1;
	break;
    }
  }
  if(!usage)
  {
    usage = (optind + 1 > argc);
  }
  if(!usage && (argc > optind))
  {
    inImgFileStr = *(argv + optind);
  }
  ok = !usage;
  /* Read input label image. */
  if(ok)
  {
    char *fStr = NULL;
    FILE *fP = stdin;

    if(strcmp(inImgFileStr, "-"))
    {
      fP = NULL;
      fStr = inImgFileStr;
    }
    inObj = WlzAssignObject(
    	    WlzEffReadObj(fP, fStr, inFmt, 0, 0, 0, &errNum), NULL);
    if(errNum == WLZ_ERR_NONE)
    {
      if(inObj == NULL)
      {
        errNum = WLZ_ERR_OBJECT_NULL;
      }
      else
      {
        switch(inObj->type)
	{
	  case WLZ_2D_DOMAINOBJ:
	    dim = 2;
	    break;
	  case WLZ_3D_DOMAINOBJ:
	    dim = 3;
	    break;
	  default:
	    errNum = WLZ_ERR_OBJECT_TYPE;
	    break;
	}
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read image object from file(s) %s (%s)\n",
	             *argv, inImgFileStr, errMsg);

    }
  }
  if(ok)
  {
    if(((labelVec = AlcVectorNew(256, sizeof(WlzEffITKSnapSegRecord), 256,
    				 NULL)) == NULL) ||
       ((textVec = AlcVectorNew(256, WLZEXTFF_ITK_LN_MAX, 256,
    				NULL)) == NULL))
    {
      ok = 0;
      errNum = WLZ_ERR_MEM_ALLOC;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to allocate label buffer (%s).\n",
		     *argv, errMsg);
    }
  }
  /* Optionally read input label description text. */
  if(ok && inLabFileStr)
  {
    FILE *fP = stdin;

    if((fP = (strcmp(inLabFileStr, "-"))?
       fopen(inLabFileStr, "r"): stdin) == NULL)
    {
      errNum = WLZ_ERR_FILE_OPEN;
    }
    else
    {
      int	idl = 0,
      		idr = 0,
		valid = 0,
		lastIdx = -1;
      char 	lnBuf[WLZEXTFF_ITK_LN_MAX];

      while((errNum == WLZ_ERR_NONE) &&
            (fgets(lnBuf, WLZEXTFF_ITK_LN_MAX, fP) != NULL))
      {
        lnBuf[WLZEXTFF_ITK_LN_MAX - 1] = '\0';
	if(lnBuf[0] == '#')
	{
	  if(!valid)
	  {
	    (void )WlzStringWhiteSpSkip(lnBuf);
	    valid = !strcmp(itkSnapMagic, lnBuf);
	  }
	}
	else if(valid)
	{
	  char		*txt;
	  float		alpha;
	  WlzUInt	col[4];
	  WlzEffITKSnapSegRecord *lab;

	  if(((lab = (WlzEffITKSnapSegRecord *)
	             AlcVectorExtendAndGet(labelVec, idr)) == NULL) ||
	     ((txt = (char *)
	             AlcVectorExtendAndGet(textVec, idr)) == NULL))
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  else if(sscanf(lnBuf, "%d %d %d %d %f %d %d",
	                 &(lab->idx),
			 &(col[0]), &(col[1]), &(col[2]), &alpha,
			 &(lab->vis), &(lab->lmv)) != 7)
	  {
	    errNum = WLZ_ERR_READ_INCOMPLETE;
	  }
	  else
	  {
	    char	*s0,
	    		*s1;

	    col[3] = WLZ_NINT(alpha * 255.0);
	    WLZ_RGBA_RGBA_SET(lab->colour,
	                      col[0], col[1], col[2], col[3]);
	    if(((s0 = strchr(lnBuf, '"')) == NULL) ||
	       ((s1 = strchr(s0 + 1, '"')) == NULL))
	    {
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      char *s,
	      	   *t;

	      t = txt;
	      for(s = s0 + 1; s < s1; ++s)
	      {
		if(!scrub || isalnum(*s))
		{
	          *t++ = *s;
		}
		else
		{
		  *t++ = '_';
		}
	      }
	      *t = '\0';
	      lab->label = txt;
	    }
	    if(((idr == 0) && (lab->idx != 0)) || (lastIdx >= lab->idx))
	    {
	      /* Indices must increment from zero. */
	      errNum = WLZ_ERR_READ_INCOMPLETE;
	    }
	    else
	    {
	      ++idr;
	    }
	  }
	  lastIdx = lab->idx;
	}
        ++idl;
	if(errNum != WLZ_ERR_NONE)
	{
	  ok = 0;
	  (void )WlzStringFromErrorNum(errNum, &errMsg);
	  (void )fprintf(stderr,
	      "%s: Failed to read label description file at line %d (%s).\n",
	      *argv, idl, errMsg);
	}
      }
    }
    if(fP && strcmp(inLabFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  /* Find upper limit on the number of domains. */
  if(ok)
  {
    WlzPixelV	minV = {0},
		maxV = {0};

    errNum = WlzGreyRange(inObj, &minV, &maxV);
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzValueConvertPixel(&maxV, maxV, WLZ_GREY_INT);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      maxLabelIdx = maxV.v.inv;
    }
    else
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
          "%s: Failed to get maximum index image value (%s).\n",
	  *argv, errMsg);
    }
  }
  /* Create a compound object with domains corresponding to the values
   * in the input image. */
  if(ok)
  {
    int		ids,
		idx,
    		emptyObj = 0,
    		allEmptyObj = 0;
    WlzPixelV	thrV;
    WlzObjectType objType;
    WlzPropertyList *pList = NULL;
    WlzObject	*tObj0 = NULL;
    WlzCompoundArray *aObj;

    ids = 0;
    thrV.type = WLZ_GREY_INT;
    objType = (dim == 3) ? WLZ_3D_DOMAINOBJ: WLZ_2D_DOMAINOBJ;
    aObj = WlzMakeCompoundArray(WLZ_COMPOUND_ARR_1, 1, maxLabelIdx,
	       			NULL, objType, &errNum);
    for(idx = 0; (errNum == WLZ_ERR_NONE) && (idx < maxLabelIdx); ++idx)
    {
      int	idv;
      WlzObject	*tObj1 = NULL;

      idv = idx + 1;
      thrV.v.inv = idv;
      if(thrV.v.inv == 1)
      {
        tObj0 = WlzAssignObject(
	        WlzThreshold(inObj, thrV, WLZ_THRESH_HIGH, &errNum), NULL);
	if(errNum == WLZ_ERR_NONE)
	{
          allEmptyObj = WlzIsEmpty(tObj0, &errNum);
        }
      }
      else if(!allEmptyObj)
      {
	WlzObject *tObj2;

        tObj2 = WlzAssignObject(
	        WlzThreshold(tObj0, thrV, WLZ_THRESH_HIGH, &errNum), NULL);
        (void )WlzFreeObj(tObj0);
	tObj0 = tObj2;
      }
      if(errNum == WLZ_ERR_NONE)
      {
        if(allEmptyObj)
	{
	  emptyObj = 1;
	}
	else
	{
	  tObj1 = WlzThreshold(tObj0, thrV, WLZ_THRESH_EQUAL, &errNum);
	  if(errNum == WLZ_ERR_NONE)
	  {
	    emptyObj = WlzIsEmpty(tObj1, &errNum);
	  }
	}
      }
      if(errNum == WLZ_ERR_NONE)
      {
        pList = WlzEffITKSnapPList(labelVec, &ids, idv, &errNum);
      }
      if(errNum == WLZ_ERR_NONE)
      {
	WlzObject *tObj2;
	WlzDomain nullDom = {0};
	WlzValues nullVal = {0};

        if(emptyObj)
	{
	  tObj2 = WlzMakeMain(WLZ_EMPTY_OBJ, nullDom, nullVal, pList,
	  		      NULL, &errNum);
	}
	else
	{
	  tObj2 = WlzMakeMain(objType, tObj1->domain, nullVal, pList,
	                      NULL, &errNum);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  aObj->o[idx] = WlzAssignObject(tObj2, NULL);
	}
	else
	{
	  (void )WlzFreePropertyList(pList);
	}
      }
      (void )WlzFreeObj(tObj1);
    }
    (void )WlzFreeObj(tObj0);
    if(errNum == WLZ_ERR_NONE)
    {
      outObj = (WlzObject *)aObj;
    }
    else
    {
      ok = 0;
      (void )WlzFreeObj((WlzObject *)aObj);
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to create compound object (%s).\n",
		     *argv, errMsg);
    }
  }
  /* Write the output object. */
  if(ok)
  {
    FILE 	*fP;

    errNum = WLZ_ERR_FILE_OPEN;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output object to file %s (%s)\n",
		     *argv, outFileStr, errMsg);
      
    }
    if(fP && strcmp(outFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  (void )AlcVectorFree(textVec);
  (void )AlcVectorFree(labelVec);
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(usage)
  {
    char *fmtStr = NULL;

    fmtStr = WlzEffFormatTable(2, 50, 10, NULL);
    (void )fprintf(
        stderr,
	"Usage: %s%s%s %s%s%s%s%s\n",
	*argv,
	" [-h] [-o<output file>] [-f<image format>]\n"
	"\t\t[-l<label file>] [-b] [<segmentation image>]\n"
        "Imports an ITKSnap segmentation to Woolz.  Reads a segmentation\n"
	"image and an optional ITKSnap segmentation label description file,\n"
	"then uses these to create a compound object with domains indexed\n"
	"as for the input image labels. If the label description file is\n"
	"provided then the objects will have their name and grey (for\n"
	"colour) properties set. By default files are read from the\n"
	"standard input (image before labels if the standard input is\n"
	"specified for the labels as well as the image).\n"
	"Version: ",
	WlzVersion(),
	"\n"
	"  -h Help, prints this usage information.\n"
	"  -o Output file.\n"
	"  -f Source file format. With simple file names, the format can be\n"
	"     inferred from the file extension.\n"
	"  -l ITKSnap segmentation label description file (for stdin use -).\n"
	"  -b Scub out non-alphanumeric characters replacing them with an\n"
	"     underscore in the label name property.\n"
        "The known file formats (only image formats are appropriate):\n"
	"  Description                                       Extension\n"
	"  ***********                                       *********\n",
	fmtStr,
	"Example:\n  ",
	*argv,
	"-l seg.txt seg.nii | WlzExplode -p\n"
	"Reads the segmentation image from seg.nii using the NIfTI format\n"
	"(inferred from file extension) and the label descriptions in\n"
	"ITKSnap text format from seg.txt. These are then used to create\n"
	"a compound object which is written via a pipe to WlzExplode.\n"
	"WlzExplode creates seperate Woolz domain files, each with the file\n"
	"body of the label description text.\n");
    AlcFree(fmtStr);
  }
  return(!ok);
}

/*!
* \return	New property list or NULL on error.
* \ingroup	BinWlzExtFF
* \brief	Creates a new property list for the ITKSnap label with
* 		the givven index value.
* \param	vec		Vector of ITKSnap labels each a
* 				WlzEffITKSnapSegRecord structure.
* \param	idS		Index into the vector from which to start
* 				searching for a matching label index, this
* 				is updated on a successful match.
* \param	idL		Required label index.
* \param	dstErr		Destination error pointer, may be NULL.
*/
static WlzPropertyList 		*WlzEffITKSnapPList(
				  AlcVector *vec,
				  int *idS,
				  int idL,
				  WlzErrorNum *dstErr)
{
  int		idx;
  WlzUInt	colour;
  char		*label;
  WlzProperty	prop = {0};
  WlzPropertyList *pList;
  WlzEffITKSnapSegRecord *lab;
  WlzErrorNum errNum = WLZ_ERR_NONE;
  const WlzUInt	defColour = 0xffffffff;
  char	* const defLabel = "Unknown";

  /* Find element in label vector corresponding to the given index. */
  idx = *idS;
  do
  {
    lab = (WlzEffITKSnapSegRecord *)AlcVectorItemGet(vec, idx++);
  } while(lab && (lab->idx < idL));
  /* Get property values. */
  if(lab && (lab->idx == idL))
  {
    *idS = idx;
    label = lab->label;
    colour = lab->colour;
  }
  else
  {
    label = defLabel;
    colour = defColour;
  }
  /* Create property list. */
  if((pList = WlzMakePropertyList(NULL)) == NULL)
  {
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    prop.name = WlzMakeNameProperty(label, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcDLPListEntryAppend(pList->list, NULL, (void *)(prop.core),
			     WlzFreePropertyListEntry) != ALC_ER_NONE)
    {
      AlcFree((void *)prop.core);
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    WlzPixelV		pix;

    pix.type = WLZ_GREY_RGBA;
    pix.v.rgbv = colour;
    prop.greyV = WlzMakeGreyProperty("Colour", pix, &errNum);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(AlcDLPListEntryAppend(pList->list, NULL, (void *)(prop.core),
			     WlzFreePropertyListEntry) != ALC_ER_NONE)
    {
      AlcFree((void *)prop.core);
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }
  if(errNum != WLZ_ERR_NONE)
  {
    (void )WlzFreePropertyList(pList);
    pList = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(pList);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
