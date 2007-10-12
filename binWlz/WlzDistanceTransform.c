#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzDistanceTransform_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlz/WlzDistanceTransform.c
* \author       Konstantinos Liakos, Bill Hill
* \date         October 2004
* \version      $Id$
* \note
*               Copyright
*               2004 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Computes distance transform objects.
* \ingroup	BinWlz
* \todo         -
* \bug          None known.
*
* \par Binary
* \ref wlzdistancetransform "WlzDistanceTransform"
*/

/*!
\ingroup BinWlz
\defgroup wlzdistancetransform WlzDistanceTransform
\par Name
WlzDistanceTransform - computes distance transform objects.
\par Synopsis
\verbatim
WlzDistanceTransform [-b] [-d#] [-f#] [-p#] [-o#] [-h] [<Reference input file>]
\endverbatim
\par Options
<table width="500" border="0">
  <tr> 
    <td><b>-h</b></td>
    <td>Help, prints usage message.</td>
  </tr>
  <tr> 
    <td><b>-o</b></td>
    <td>Output file.</td>
  </tr>
  <tr> 
    <td><b>-b</b></td>
    <td>Use the boundary of the reference object.</td>
  </tr>
  <tr> 
    <td><b>-d</b></td>
    <td>Distance function:
      <table width="500" border="0">
      <tr> <td>0</td> <td>Euclidean (2D and 3D, unimplemented)</td></tr>
      <tr> <td>1</td> <td>octagonal (2D and 3D) - default</td></tr>
      <tr> <td>2</td> <td>approximate Euclidean (2D and 3D)</td></tr>
      <tr> <td>4</td> <td>4-connected (2D)</td></tr>
      <tr> <td>8</td> <td>8-connected (2D)</td></tr>
      <tr> <td>6</td> <td>6-connected (3D)</td></tr>
      <tr> <td>18</td> <td>18-connected (3D)</td></tr>
      <tr> <td>26</td> <td>26-connected (3D)</td></tr>
    </td>
  </tr>
  <tr> 
    <td><b>-f</b></td>
    <td>The foreground object file.</td>
  </tr>
  <tr> 
    <td><b>-p</b></td>
    <td>Distance function parameter, which is curently only used for
        approximate Euclidean distance transforms. The value should
	be \f$geq\f$ 1, with larger values giving greater accuracy
	at the cost of increased time. Default value - 10.0.</td>
  </tr>
</table>
\par Description
Computes a distance transform object which has the domain of the
foreground object and values which are the distance from the reference
object.
\par Examples
\verbatim
WlzDistanceTransform -b -f fish3.wlz fish4.wlz > out.wlz
\endverbatim
Creates a new domain object with values, in which the domain is the same as 
that of fish3.wlz and the values are the minimum octagonal distance from their
location to the boundary of the reference object fish4.wlz.
\par File
\ref WlzDistanceTransform.c "WlzDistanceTransform.c"
\par See Also
\ref BinWlz "WlzIntro(1)"
\ref WlzDistanceTransform "WlzDistanceTransform(3)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <Wlz.h>

extern int      getopt(int argc, char * const *argv, const char *optstring);
 
extern char     *optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		ok,
		bnd = 0,
  		con,
  		option,
		usage = 0;
  char		*inFileStr,
  		*forFileStr,
  		*refFileStr,
		*outFileStr;
  double	dParam = 10.0;
  FILE		*fP = NULL;
  WlzObject	*tObj0,
  		*tObj1,
		*forObj = NULL,
		*refObj = NULL,
		*dstObj = NULL;
  WlzDistanceType dFn;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsgStr;
  static char	optList[] = "bhd:f:p:o:";
  const char    fileStrDef[] = "-";
  const WlzConnectType defDFn = WLZ_OCTAGONAL_DISTANCE;

  /* Parse the argument list and check for input files. */
  opterr = 0;
  inFileStr = (char *)fileStrDef;
  forFileStr = (char *)fileStrDef;
  outFileStr = (char *)fileStrDef;
  dFn = defDFn;
  while((option = getopt(argc, argv, optList)) != EOF)
  {
    switch(option)
    {
      case 'b':
        bnd = 1;
	break;
      case 'd':
	if(sscanf(optarg, "%d", &con) != 1)
	{
	  usage = 1;
	}
	else
	{
	  switch(con)
	  {
	    case 0:
	      dFn = WLZ_EUCLIDEAN_DISTANCE;
	      break;
	    case 1:
	      dFn = WLZ_OCTAGONAL_DISTANCE;
	      break;
	    case 2:
	      dFn = WLZ_APX_EUCLIDEAN_DISTANCE;
	      break;
	    case 4:
	      dFn = WLZ_4_DISTANCE;
	      break;
	    case 6:
	      dFn = WLZ_6_DISTANCE;
	      break;
	    case 8:
	      dFn = WLZ_8_DISTANCE;
	      break;
	    case 18:
	      dFn = WLZ_18_DISTANCE;
	      break;
	    case 26:
	      dFn = WLZ_26_DISTANCE;
	      break;
	    default:
	      usage = 1;
	      break;
	  }
	}
	break;
      case 'f':
	forFileStr = optarg;
	break;
      case 'p':
	if(sscanf(optarg, "%lg", &dParam) != 1)
	{
	  usage = 1;
	}
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = !usage;
  if(ok && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
      ok = 0;
    }
    else
    {
      refFileStr = argv[optind];
    }
  }
  /* Read the reference object. */
  if(ok)
  {
    if((refFileStr == NULL) ||
       (*refFileStr == '\0') ||
       ((fP = (strcmp(refFileStr, "-")?
              fopen(refFileStr, "r"): stdin)) == NULL) ||
       ((refObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read reference object from file %s.\n",
                     argv[0], refFileStr);
    }
    if(fP && strcmp(refFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Read the foreground object. */
  if(ok)
  {
    if((forFileStr == NULL) ||
       (*forFileStr == '\0') ||
       ((fP = (strcmp(forFileStr, "-")?
              fopen(forFileStr, "r"): stdin)) == NULL) ||
       ((forObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to read foreground object from file %s.\n",
                     argv[0], forFileStr);
    }
    if(fP && strcmp(forFileStr, "-"))
    {
      (void )fclose(fP); fP = NULL;
    }
  }
  /* Check object types and parse the promote the default distance function
   * if required. */
  if(ok)
  {
    if(refObj->type != forObj->type)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Foreground and reference object types differ.\n",
		     argv[0]);
    }
  }
  if((errNum == WLZ_ERR_NONE) && bnd)
  {
    tObj0 = tObj1 = NULL;
    tObj0 = WlzObjToBoundary(refObj, 1, &errNum);
    if(errNum == WLZ_ERR_NONE)
    {
      tObj1 = WlzBoundaryToObj(tObj0, WLZ_VERTEX_FILL, &errNum);
    }
    if(errNum == WLZ_ERR_NONE)
    {
      (void )WlzFreeObj(refObj);
      refObj = tObj1;
      tObj1 = NULL;
    }
    (void )WlzFreeObj(tObj0);
    (void )WlzFreeObj(tObj1);
  }
  /* Compute the distance transform object. */
  if(ok)
  {
    dstObj = WlzAssignObject(
    	     WlzDistanceTransform(forObj, refObj, dFn, dParam, &errNum), NULL);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to compute distance object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  /* Output the distance transform object. */
  if(ok)
  {
    if((fP = (strcmp(outFileStr, "-")?
             fopen(outFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, dstObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
                     argv[0], errMsgStr);
    }
  }
  (void )WlzFreeObj(forObj);
  (void )WlzFreeObj(refObj);
  (void )WlzFreeObj(dstObj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-b] [-d#] [-f#] [-p#] [-o#] [-h] [<Reference input file>]\n"
    "Computes a distance transform object which has the domain of the\n"
    "foreground object and values which are the distance from the reference\n"
    "object.\n"
    "Options are:\n"
    "  -b  Use the boundary of the reference object.\n"
    "  -d  Distance function:\n"
    "              0: Euclidean (2D and 3D, unimplemented)\n"
    "              1: octagonal (2D and 3D) - default\n"
    "              2: approximate Euclidean (2D and 3D)\n"
    "              4: 4-connected (2D)\n"
    "              8: 8-connected (2D)\n"
    "              6: 6-connected (3D)\n"
    "             18: 18-connected (3D)\n"
    "             26: 26-connected (3D)\n"
    "  -f  The foreground object file.\n" 
    "  -p  Distance function parameter, which is curently only used for\n"
    "      approximate Euclidean distance transforms. The value should\n"
    "      be >= 1, with larger values giving greater accuracy at the\n"
    "      cost of increased time. Default value - 10.0.\n"
    "  -o  Output object file.\n"
    "  -h  Help - prints this usage message\n",
    argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
