#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzInteriority_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzInteriority.c
* \author       Bill Hill
* \date         January 2015
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
* \brief 	Computes interiority scores for object domains.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzinteriority "WlzInteriority"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzinteriority WlzInteriority
\par Name
WlzInteriority - Computes interiority scores for object domains.
\par Synopsis
\verbatim
WlzInteriority [-a] [-h] -r<ref obj> [-o<out file>] [-v] [<test object>]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-a</b></td>
    <td>Read all objects before computing scores (default is read each
        object in turn).</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output file name.</td>
  </tr>
  <tr>
    <td><b>-r</b></td>
    <td>Reference object.</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose output; prefixes scores with index and filename.</td>
  </tr>
</table>
The input objects are read from stdin (first the reference object then
the test objects from the command line followed by test objects from the
standard input) and scores are written to stdout unless filenames are given.
Use - to indicate the standard input and output.

\par Description
Computes interiority scores for test domain objects with respect to
the reference object. The interiority is the mean distance that the
(intersection of the reference and the) test object's domain is
inside the reference object.

\par Examples
\verbatim
 
# Make a pair of rectangles and then compute the interiority of the
# smaller with respect to the larger:
 
WlzMakeRect -2 -x 0,99 -y 0,99 >ref.wlz
WlzMakeRect -2 -x 0,49 -y 0,49 >tst.wlz
WlzInteriority -r ref.wlz tst.wlz
18.938
\endverbatim

\par File
\ref WlzInteriority.c "WlzInteriority.c"
\par See Also
\ref WlzInteriority  "WlzInteriority(3)"
\ref WlzInteriorityN "WlzInteriorityN(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <string.h>
#include <Wlz.h>

static WlzErrorNum		WlzInteriorityRealloc(
				  WlzObject ***objs,
				  char ***fnames,
				  double **scores,
				  int *max,
				  int n);
static WlzObject		*WlzInteriorityReadObj(
				  FILE **fP,
				  const char *fileStr,
				  WlzErrorNum *dstErr);


int             main(int argc, char **argv)
{
  int		option,
		all = 0,
		ok = 1,
		usage = 0,
		nFile = 0,
		nScore = 0,
		nTstObj = 0,
		verbose = 0,
		maxTstObj = 0;
  FILE		*fP = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  WlzObject	*refObj = NULL;
  WlzObject	**tstObjs = NULL;
  double	*scores = NULL;
  char 		*outFileStr,
  		*tstObjFileStr,
		*refObjFileStr;
  char		**filenames = NULL;
  const char	*errMsg;
  static char	optList[] = "ahvo:r:",
		stdStr[] = "-";

  outFileStr = stdStr;
  tstObjFileStr = stdStr;
  refObjFileStr = NULL;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'a':
        all = 1;
	break;
      case 'o':
        outFileStr = optarg;
	break;
      case 'r':
        refObjFileStr = optarg;
	break;
      case 'v':
        verbose = 1;
	break;
      case 'h':
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if((tstObjFileStr == NULL) || (*tstObjFileStr == '\0') ||
     (outFileStr == NULL) || (*outFileStr == '\0') ||
     (refObjFileStr == NULL))
  {
    ok = 0;
    usage = 1;
  }
  if(ok)
  {
    if((refObj = WlzAssignObject(
                 WlzInteriorityReadObj(&fP, refObjFileStr,
		                       &errNum), NULL)) == NULL)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to read reference object from file %s (%s)\n",
		     *argv, refObjFileStr, errMsg);
    }
    if(fP)
    {
      (void )fclose(fP);
      fP = NULL;
    }
  }
  if(ok)
  {
    if(optind >= argc)
    {
      ok = 0;
      usage = 1;
    }
    else
    {
      int 	i;

      for(i = optind; (errNum == WLZ_ERR_NONE) && (i < argc); ++i)
      {
	int	readOne = 0;

	++nFile;
	tstObjFileStr = argv[i];
	do
	{
	  errNum = WlzInteriorityRealloc(&tstObjs, &filenames, &scores,
	  				 &maxTstObj, ALG_MAX(nScore, nTstObj));
	  if(errNum == WLZ_ERR_NONE)
	  {
	    tstObjs[nTstObj] = WlzAssignObject(
	                       WlzInteriorityReadObj(&fP, tstObjFileStr,
	                                             &errNum), NULL);
	  }
	  if(errNum == WLZ_ERR_NONE)
	  {
	    readOne = 1;
	    if(all)
	    {
	      filenames[nTstObj] = tstObjFileStr;
	      ++nTstObj;
	    }
	    else
	    {
	      filenames[nScore] = tstObjFileStr;
	      scores[nScore] = WlzInteriority(refObj, tstObjs[0], &errNum);
	      (void )WlzFreeObj(tstObjs[0]);
	      tstObjs[0] = NULL;
	      ++nScore;
	    }
	  }
	} while(errNum == WLZ_ERR_NONE);
	if(errNum != WLZ_ERR_NONE)
	{
	  if(fP)
	  {
	    (void )fclose(fP);
	    fP = NULL;
	  }
	  if((errNum == WLZ_ERR_READ_EOF) && readOne)
	  {
	    errNum = WLZ_ERR_NONE;
	  }
	}
      }
    }
    if((errNum == WLZ_ERR_NONE) && all)
    {
      double	*s = NULL;

      nScore = nTstObj;
      s = WlzInteriorityN(refObj, nTstObj, tstObjs, &errNum);
      AlcFree(scores);
      scores = s;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to computes interiority score at file %d\n"
		     "\t\t(%s)\n",
		     *argv, nFile, errMsg);
    }
  }
  (void )WlzFreeObj(refObj);
  if(ok)
  {
    if(errNum == WLZ_ERR_NONE)
    {
      if((fP = (strcmp(outFileStr, "-")?
               fopen(outFileStr, "w"): stdout)) == NULL)
      {
        ok = 0;
	errNum = WLZ_ERR_FILE_OPEN;
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	               "%s: failed to open output file %s (%s)\n",
		       *argv, outFileStr, errMsg);
      }
      else
      {
	int	i;

	for(i = 0; i < nScore; ++i)
	{
	  if(verbose)
	  {
	    if(fprintf(fP, "%d %s %g\n", i, filenames[i], scores[i]) < 0)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	  }
	  else
	  {
	    if(fprintf(fP, "%g\n", scores[i]) < 0)
	    {
	      errNum = WLZ_ERR_WRITE_INCOMPLETE;
	    }
	  }
	}
        if(strcmp(outFileStr, "-"))
	{
	  (void )fclose(fP);
	}
      }
      if(errNum != WLZ_ERR_NONE)
      {
        ok = 0;   
	(void )WlzStringFromErrorNum(errNum, &errMsg);
	(void )fprintf(stderr,
	               "%s: failed to write to file %s (%s)\n",
		       *argv, outFileStr, errMsg);
      }
    }
  }
  if(tstObjs)
  {
    int	i;

    for(i = 0; i < nTstObj; ++i)
    {
      (void )WlzFreeObj(tstObjs[i]);
    }
    AlcFree(tstObjs);
  }
  AlcFree(filenames);
  AlcFree(scores);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s%s%sExample: %s%s",
    *argv,
    " [-a] [-h] [-o<out file>] -r<ref obj> [-v] [<test objects>]\n"
    "Version: ",
    WlzVersion(),
    "\n"
    "Options:\n"
    "  -a  Read all objects before computing scores (default is read each\n"
    "      object in turn).\n"
    "  -h  Help, prints this usage message.\n"
    "  -r  Reference object, must be supplied.\n"
    "  -o  Output file name.\n"
    "  -v  Verbose output; prefixes scores with index and filename.\n"
    "Computes interiority scores for test object domains with respect to\n"
    "a reference object's domain.\n"
    "The input objects are read from stdin (first the reference object then\n"
    "the test objects from the command line followed by test objects from\n"
    "the standard input) and scores are written to stdout unless filenames\n"
    "are given. Use - to indicate the standard input and output.\n",
    *argv,
    "WlzMakeRect -2 -x 0,99 -y 0,99 >ref.wlz\n"
    "WlzMakeRect -2 -x 0,49 -y 0,49 >tst.wlz\n"
    "%s -r ref.wlz tst.wlz\n"
    "18.938\n"
    "A pair of rectangles are created then the interiority of the smaller\n"
    "with respect to the larger is computed.\n");
  }
  return(!ok);
}

/*!
* \return	Woolz object read from file or NULL if either the file was
* 		empty or there was an error.
* \ingroup	BinWlz
* \brief	Reads object from file given file name.
* \param	fP			Source and destination file pointer
* 					which must point to a NULL file
* 					pointer when a new file is to be
* 					opened.
* \param	fileStr			Given file string or "-" for the
* 					standard input.
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject		*WlzInteriorityReadObj(
				  FILE **fP,
				  const char *fileStr,
				  WlzErrorNum *dstErr)
{
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(*fP == NULL)
  {
    if((*fP = strcmp(fileStr, "-")?  fopen(fileStr, "r"): stdin) == NULL)
    {
      errNum = WLZ_ERR_FILE_OPEN;
    }
  }
  if(*fP)
  {
    obj = WlzReadObj(*fP, &errNum);
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

			
/*!
* \return	Woolz error code.
* \ingroup	BinWlz
* \brief
* \param	objs		Destination pointer for the array of objects.
* \param	fnames		Destination pointer for the array of file
* 				name pointers.
* \param	scores		Destination pointer for the array of scores.
* \param	max		Source and destination pointer for the
* 				current maximum number of objects and
* 				scores.
* \param	n		Number of objects and scores required.
*/
static WlzErrorNum		WlzInteriorityRealloc(
				  WlzObject ***objs,
				  char ***fnames,
				  double **scores,
				  int *max,
				  int n)
{
  size_t	newMax;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  newMax = *max;
  if(n >= *max)
  {
    char	**f = NULL;
    double	*s = NULL;
    WlzObject	**o = NULL;

    do
    {
      newMax = (newMax <= 0)? 8: newMax * 2;
    } while(n >= newMax);
    if(((o = (WlzObject **)
	     AlcRealloc(*objs, sizeof(WlzObject *) * newMax)) == NULL) ||
       ((f = (char **)
	     AlcRealloc(*fnames, sizeof(char *) * newMax)) == NULL) ||
       ((s = (double *)
	     AlcRealloc(*scores, sizeof(double) * newMax)) == NULL))
    {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else
    {
      *objs = o;
      *fnames = f;
      *scores = s;
      *max = newMax;
    }
  }
  return(errNum);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */

