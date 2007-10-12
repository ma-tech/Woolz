#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzBibTransAbsRel_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlzApp/WlzBibTransAbsRel.c
* \author       Bill Hill
* \date         February 2005
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
* \brief	Reads a reconstruction bibfile and writes it back out
*		with either absolute or relative transforms.
* \ingroup	BinWlzApp
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzbibtransabsrel "WlzBibTransAbsRel"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzbibtransabsrel WlzBibTransAbsRel
\par Name
WlzBibTransAbsRel - converts reconstruction bibfile between absolute and
                    relative transforms.
\par Synopsis
\verbatim
WlzBibTransAbsRel [-a] [-h] [-r] [-o outfile] [<infile>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr> 
    <td><b>-a</b></td>
    <td>Output bibfile with absolute transforms, default.</td>
  </tr>
  <tr> 
    <td><b>-r</b></td>
    <td>Output bibfile with relative transforms.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
</table>
\par Description
WlzBibTransAbsRel reads a reconstruction bibfile
and writes it to the output file
with either absolute or relative transforms.
\par Examples
\verbatim
WlzBibTransAbsRel <recon_rel.bib >recon_abs.bib
\endverbatim
Reads a releative transform reconstruction bibfile from recon_rel
computes the absolute transforms and then writes the reconstruction
bibfile with absolute transforms to the file recon_abs.
\par File
\ref WlzFacts.c "WlzFacts.c"
\par See Also
\ref BinWlzApp "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <Wlz.h>
#include <bibFile.h>
#include <Reconstruct.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;

/* TODO move these functions to libReconstruct and libWlz */
RecError			RecSecListAbsFromRel(
				  RecSectionList *secList);
RecError        		RecSecListRelFromAbs(
				  RecSectionList *secList);
WlzErrorNum			WlzTransformNearIdentity2D(
				  WlzAffineTransform *tr,
				  double tol);


int		 main(int argc, char *argv[])
{
  int		option,
  		numSec,
  		absMode = 1,
		ok = 1,
  		usage = 0;
  char		*inFile,
  		*outFile;
  FILE		*fP = NULL;
  RecError	recErr = REC_ERR_NONE;
  RecSectionList *secList = NULL;
  char		*errMsg;
  static char	optList[] = "ahro:",
  		defFile[] = "-",
		defErrMsg[] = "none.";

  opterr = 0;
  errMsg = defErrMsg;
  inFile = outFile = defFile;
  while((usage == 0) &&
        ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
	absMode = 1;
	break;
      case 'r':
	absMode = 0;
	break;
      case 'o':
	outFile = optarg;
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
      inFile = *(argv + optind);
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
		     *argv);
    }
  }
  if(ok)
  {
    if((inFile == NULL) ||
       (*inFile == '\0') ||
       ((fP = (strcmp(inFile, "-")?  fopen(inFile, "r"): stdin)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open input file %s\n",
		     *argv, (inFile)? inFile: "(null)");
    }
  }
  if(ok)
  {
    if(RecFileSecListRead(secList, &numSec, fP, &errMsg) != REC_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read section list from file %s.\n"
		     "Error message - %s\n",
		     *argv, inFile, errMsg);
    }
  }
  if(fP && (strcmp(inFile, "-") == 0))
  {
    (void )fclose(fP);
  }
  /* Convert transforms. */
  if(ok)
  {
    if(absMode)
    {
      recErr = RecSecListAbsFromRel(secList);
    }
    else
    {
      recErr = RecSecListRelFromAbs(secList);
    }
    if(recErr != REC_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to convert the section list from %s to\n"
		     "%s transforms.\n",
		     *argv,
		     (absMode)? "relative": "absolute",
		     (absMode)? "absolute": "relative");

    }
  }
  /* Write output bibfile. */
  if(ok)
  {
    if((outFile == NULL) ||
       (*outFile == '\0') ||
       ((fP = (strcmp(outFile, "-")?  fopen(outFile, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s\n",
		     *argv, (inFile)? inFile: "(null)");
    }
  }
  if(ok)
  {
    if(RecFileSecListWrite(fP, secList, numSec, &errMsg) != REC_ERR_NONE)
    {
      (void )fprintf(stderr,
      		     "%s: Failed to write section list to file %s.\n"
		     "Error message - %s\n",
		     *argv, inFile, errMsg);
      ok = 0;
    }
  }
  if(fP && (strcmp(outFile, "-") == 0))
  {
    (void )fclose(fP);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-a] [-h] [-r] [-o outfile] [<infile>]\n"
    "Reads a reconstruction bibfile and writes it to the output file\n"
    "with either absolute or relative transforms.\n"
    "Options are:\n"
    "  -a  Output bibfile with absolute transforms (%s).\n"
    "  -h  Help - print this usage information.\n"
    "  -r  Output bibfile with relative transforms (%s).\n"
    "  -o  Output file name.\n"
    "By default %s reads from it's standard input and\n"
    "writes to it's standard output.\n",
    *argv,
    (absMode != 0)? "set": "not set",
    (absMode == 0)? "set": "not set",
    *argv);
  }
  exit(!ok);
}

/* TODO move these functions to libReconstruct. */

RecError	RecSecListAbsFromRel(RecSectionList *secList)
{
  RecSection	*sec;
  HGUDlpListItem *item;
  WlzAffineTransform *tr0 = NULL,
  		*tr1;
  RecError	recErr = REC_ERR_NONE;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  if(secList == NULL)
  {
    recErr = REC_ERR_FUNC;
  }
  else if(secList->attributes.trMode == REC_TRMODE_REL)
  {
    item = HGUDlpListHead(secList->list);
    do
    {
      sec = (RecSection *)HGUDlpListEntryGet(secList->list, item);
      if(sec && strcmp(sec->imageFile, "empty"))
      {
	if(tr0)
	{
	  tr1 = WlzAffineTransformProduct(sec->transform, tr0, &wlzErr);
	  if(wlzErr == WLZ_ERR_NONE)
	  {
	    (void )WlzFreeAffineTransform(sec->transform);
	    sec->transform = WlzAssignAffineTransform(tr1, NULL);
	  }
	  tr0 = tr1;
	}
	else
	{
	  tr0 = sec->transform;
	}
      }
      item = HGUDlpListNext(secList->list, item);
    } while((wlzErr == WLZ_ERR_NONE) && item &&
            (item != HGUDlpListHead(secList->list)));
    secList->attributes.trMode = REC_TRMODE_ABS;
    recErr = RecErrorFromWlz(wlzErr);
  }
  return(recErr);
}

RecError        RecSecListRelFromAbs(RecSectionList *secList)
{
  RecSection	*sec;
  HGUDlpListItem *item;
  WlzAffineTransform *tr0 = NULL,
  		*tr1,
		*tr2;
  RecError	recErr = REC_ERR_NONE;
  WlzErrorNum	wlzErr = WLZ_ERR_NONE;

  if(secList == NULL)
  {
    recErr = REC_ERR_FUNC;
  }
  else if(secList->attributes.trMode == REC_TRMODE_ABS)
  {
    item = HGUDlpListHead(secList->list);
    do
    {
      sec = (RecSection *)HGUDlpListEntryGet(secList->list, item);
      if(sec && strcmp(sec->imageFile, "empty"))
      {
        if(tr0)
	{
	  tr1 = WlzAffineTransformInverse(tr0, &wlzErr);
	  if(wlzErr == WLZ_ERR_NONE)
	  {
	    (void )WlzFreeAffineTransform(tr0);
	    tr0 = WlzAssignAffineTransform(sec->transform, NULL);
	    tr2 = WlzAffineTransformProduct(tr0, tr1, &wlzErr);
	  }
	  if(wlzErr == WLZ_ERR_NONE)
	  {
	    (void )WlzTransformNearIdentity2D(tr2, 1.0e-6);
	    (void )WlzFreeAffineTransform(sec->transform);
	    sec->transform = WlzAssignAffineTransform(tr2, NULL);
	  }
	}
	else
	{
	  tr0 = WlzAssignAffineTransform(sec->transform, NULL);
	}
      }
      item = HGUDlpListNext(secList->list, item);
    } while((wlzErr == WLZ_ERR_NONE) && item &&
            (item != HGUDlpListHead(secList->list)));
    if(tr0)
    {
      (void )WlzFreeAffineTransform(tr0);
    }
    secList->attributes.trMode = REC_TRMODE_REL;
    recErr = RecErrorFromWlz(wlzErr);
  }
  return(recErr);
}


/* TODO move these functions to libWlz */

/*!
* \return	Woolz error code.
* \ingroup	WlzTransform
* \brief	Checks whether the given 2D affine transform is near
*		to an identity transform. If it's near it is set to
*		an identity transform.
* \param	tr			Given 2D affine transform.
* \param	tol			Tolerance value.
*/
WlzErrorNum	WlzTransformNearIdentity2D(WlzAffineTransform *tr, double tol)
{
  WlzAffineTransformPrim prim;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(tr)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzAffineTransformPrimGet(tr, &prim);
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(!finite(prim.tx) || (fabs(prim.tx) <= tol))
    {
      prim.tx = 0.0;
    }
    if(!finite(prim.ty) || (fabs(prim.ty) <= tol))
    {
      prim.ty = 0.0;
    }
    if(!finite(prim.tz) || (fabs(prim.tz) <= tol))
    {
      prim.tz = 0.0;
    }
    if(!finite(prim.scale) || (fabs(prim.scale - 1.0) <= tol))
    {
      prim.scale = 1.0;
    }
    if(!finite(prim.theta) || (fabs(prim.theta) <= tol))
    {
      prim.theta = 0.0;
    }
    if(!finite(prim.phi) || (fabs(prim.phi) <= tol))
    {
      prim.phi = 0.0;
    }
    if(!finite(prim.alpha) || (fabs(prim.alpha) <= tol))
    {
      prim.alpha = 0.0;
    }
    if(!finite(prim.psi) || (fabs(prim.psi) <= tol))
    {
      prim.psi = 0.0;
    }
    if(!finite(prim.xsi) || (fabs(prim.xsi) <= tol))
    {
      prim.xsi = 0.0;
    }
    errNum = WlzAffineTransformPrimSet(tr, prim);
  }
  return(errNum);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
