#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzGeoModelFilters.c
* \author       Bill Hill
* \date         November 2000
* \version      $Id$
* \note
*               Copyright
*               2001 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Filters for Woolz geometric models (GM's).
* \ingroup      WlzGeoModel
* \todo         -
* \bug          None known.
*/
#include <Wlz.h>
#include <limits.h>

static WlzErrorNum		WlzGMFilterRmSmShells2(
				  WlzGMModel *model,
				  int maxElm);
static WlzErrorNum		WlzGMFilterRmSmShells3(
				  WlzGMModel *model,
				  int maxElm);

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Removes small shells from the given geometric model.
* \param	model			Given model.
* \param	minSpx			Minimum number of simplicies
*                                       (edges in a 2D model or faces
*                                       in a 3D model) in a shell
*                                       for it to stay in the model.
*/
WlzErrorNum	WlzGMFilterRmSmShells(WlzGMModel *model, int minSpx)
{
  int		sSz;
  WlzGMShell	*cS,
  		*fS,
		*nS;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if((fS = model->child) != NULL)
  {
    cS = fS->next;
    while(cS != fS)
    {
       nS = cS->next;
       if( WlzGMShellSimplexCnt(cS) < minSpx)
      {
         WlzGMModelDeleteS(model, cS);
       }
      cS = nS;
    }
    if( WlzGMShellSimplexCnt(fS) < minSpx)
    {
      WlzGMModelDeleteS(model, fS);
    }
  }
  return(errNum);
}

/* #define WLZ_GMFILTER_TEST */
#ifdef WLZ_GMFILTER_TEST

/* Test main() for WlzGMFilterRmSmShells(). */

extern int	getopt(int argc, char * const *argv, const char *optstring);

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int             main(int argc, char *argv[])
{
  int           idx,
  		option,
  		ok = 1,
		usage = 0,
		smShellSz = -1;
  FILE		*fP = NULL;
  char		*outObjFileStr;
  char		*inObjFileStr;
  WlzObject	*inObj;
  WlzDomain	ctrDom;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "o:s:h:";
  const char	inObjFileStrDef[] = "-",
  		outObjFileStrDef[] = "-";

  inObj = NULL;
  ctrDom.core = NULL;
  inObjFileStr = (char *)inObjFileStrDef;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 's':
        if(sscanf(optarg, "%d", &smShellSz) != 1)
	{
	  usage = 1;
	  ok = 0;
	}
	break;
      case 'h':
        usage = 1;
	ok = 0;
	break;
      default:
        usage = 1;
	ok = 0;
	break;
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
    if(ok && (optind < argc))
    {
      if((optind + 1) != argc)
      {
        usage = 1;
	ok = 0;
      }
      else
      {
        inObjFileStr = *(argv + optind);
      }
    }
  }
  if(ok)
  {
    if((inObjFileStr == NULL) || (*inObjFileStr == '\0') ||
	((fP = (strcmp(inObjFileStr, "-")?
		fopen(inObjFileStr, "r"): stdin)) == NULL) ||
	((inObj = WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL) ||
	(errNum != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to read object from file %s\n",
		     *argv, inObjFileStr);
    }
    if(fP && strcmp(inObjFileStr, "-"))
    {
      (void )fclose(fP);
      fP = NULL;
    }
  }
  if(ok)
  {
    if(inObj->type != WLZ_CONTOUR)
    {
      ok = 0;
      (void )fprintf(stderr, "%s Input object not a contour.\n",
		     argv[0]);
    }
  }
  if(ok && (smShellSz > 0))
  {
    errNum = WlzGMFilterRmSmShells(inObj->domain.ctr->model, smShellSz);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
		     "%s: failed to filter small shells from contour (%s).\n",
		     *argv, errMsg);
    }
  }
  if(ok)
  {
    if((fP = fopen("tstOut1", "w")) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to open test output file (%s)\n",
		     *argv, "tstOut");
    }
  }
  if(ok)
  {
  if((fP = (strcmp(outObjFileStr, "-")?
           fopen(outObjFileStr, "w"): stdout)) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to open output file %s.\n",
                     argv[0], outObjFileStr);
    }
  }
  if(ok)
  {
    errNum = WlzWriteObj(fP, inObj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr,
                     "%s: Failed to write output object.\n",
                     argv[0]);
    }
  }
  if(fP)
  {
    (void )fclose(fP);
    fP = NULL;
  }
  if(inObj)
  {
    (void )WlzFreeObj(inObj);
  }
  if(usage)
  {
      (void )fprintf(stderr,
      "Usage: %s %s",
      *argv,
      " [-h] [<input object>]\n"
      "Options:\n"
      "  -o  Output file.\n"
      "  -s  Threshold shell size (edges in 2D, loops in 3D) for small\n"
      "      shells to be removed from the model.\n"
      "  -h  Prints this usage information.\n"
      "Test for various Woolz Geometric Model filters.\n");
  }
  return(!ok);
}

#endif /* WLZ_GMFILTER_TEST */
