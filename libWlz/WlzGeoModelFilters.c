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
* \param	maxElm			Maximum number of elements
*                                       (edges in a 2D model or loops
*                                       in a 3D model) in a shell
*                                       for it to be a small shell.
*/
WlzErrorNum	WlzGMFilterRmSmShells(WlzGMModel *model, int maxElm)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(model == NULL)
  {
    errNum = WLZ_ERR_DOMAIN_NULL;
  }
  else if(model->child)
  {
    switch(model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
        errNum = WlzGMFilterRmSmShells2(model, maxElm);
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
        errNum = WlzGMFilterRmSmShells3(model, maxElm);
	break;
      default:
        errNum = WLZ_ERR_DOMAIN_TYPE;
	break;
    }
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Removes small shells from the given 2D geometric model.
* \param	model			Given model.
* \param	maxElm			Maximum number of elements
*                                       (edges) in a shell for it to
*                                       be a small shell.
*/
static WlzErrorNum WlzGMFilterRmSmShells2(WlzGMModel *model, int maxElm)
{
  int		fETIdx,
  		fLTIdx,
  		eCnt;
  WlzGMEdgeT	*tET;
  WlzGMLoopT	*tLT;
  WlzGMShell	*lS,
		*nS,
		*tS;
  char		*eFlg = NULL,
  		*lFlg = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if(((eFlg = (char *)AlcCalloc(model->res.edge.numIdx,
  			        sizeof(char))) == NULL) ||
     ((lFlg = (char *)AlcCalloc(model->res.loop.numIdx,
  			        sizeof(char))) == NULL))
  {
    errNum == WLZ_ERR_MEM_ALLOC;
  }
  if((errNum == WLZ_ERR_NONE) && ((tS = model->child) != NULL))
  {
    /* For each shell. */
    lS = tS;
    nS = tS->next;
    do
    {
      eCnt = 0;
      tS = nS;
      nS = nS->next;
      tLT = tS->child;
      fLTIdx = tLT->idx;
      /* For each loopT. */
      do
      {
	if(*(lFlg + tLT->loop->idx) == 0)
	{
	  /* For each loop, there are two loopT's per loop, count the
	   * number of edges. */
	  tET = tLT->edgeT;
          fETIdx = tET->idx;
	  do
	  {
	    if(*(eFlg + tET->edge->idx) == 0)
	    {
	      ++eCnt;
	      *(eFlg + tET->edge->idx) = 1; 		/* Mark edge visited */
	    }
	    tET = tET->next;
	  } while((tET->idx != fETIdx) && (eCnt <= maxElm));
	  *(lFlg + tLT->loop->idx) = 1; 		/* Mark loop visited */
	}
	tLT = tLT->next;
      } while((tLT->idx != fLTIdx) && (eCnt <= maxElm));
      if(eCnt < maxElm)
      {
        /* This is a small shell, so delete it from the model. */
	errNum = WlzGMModelDeleteS(model, tS);
      }
    } while((errNum == WLZ_ERR_NONE) && model->child && (tS->idx != lS->idx));
  }
  if(eFlg)
  {
    AlcFree(eFlg);
  }
  if(lFlg)
  {
    AlcFree(lFlg);
  }
  return(errNum);
}

/*!
* \return				Woolz error code.
* \ingroup      WlzGeoModel
* \brief	Removes small shells from the given 3D geometric model.
* \param	model			Given model.
* \param	maxElm			Maximum number of elements
*                                       (loops) in a shell for it to
*                                       be a small shell.
*/
static WlzErrorNum WlzGMFilterRmSmShells3(WlzGMModel *model, int maxElm)
{
  int		fLTIdx,
  		lCnt;
  WlzGMLoopT	*tLT;
  WlzGMShell	*lS,
		*nS,
		*tS;
  char		*lFlg = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  if((lFlg = (char *)AlcCalloc(model->res.loop.numIdx, sizeof(char))) == NULL)
  {
    errNum == WLZ_ERR_MEM_ALLOC;
  }
  if((errNum == WLZ_ERR_NONE) && ((tS = model->child) != NULL))
  {
    /* For each shell. */
    lS = tS;
    nS = tS->next;
    do
    {
      lCnt = 0;
      tS = nS;
      nS = nS->next;
      tLT = tS->child;
      fLTIdx = tLT->idx;
      /* For each loopT. */
      do
      {
	if(*(lFlg + tLT->loop->idx) == 0)
	{
	  ++lCnt;
	  *(lFlg + tLT->loop->idx) = 1;
	}
	tLT = tLT->next;
      } while((tLT->idx != fLTIdx) && (lCnt <= maxElm));
      if(lCnt < maxElm)
      {
        /* This is a small shell, so delete it from the model. */
	errNum = WlzGMModelDeleteS(model, tS);
      }
    } while((errNum == WLZ_ERR_NONE) && model->child && (tS->idx != lS->idx));
  }
  if(lFlg)
  {
    AlcFree(lFlg);
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
  WlzDVertex2	tstOutOffset2D,
  		tstOutScale2D;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "s:h:";
  const char	inObjFileStrDef[] = "-";

  inObj = NULL;
  ctrDom.core = NULL;
  tstOutOffset2D.vtX = 100.0;
  tstOutOffset2D.vtY = 100.0;
  tstOutScale2D.vtX = 1.0;
  tstOutScale2D.vtY = 1.0;
  inObjFileStr = (char *)inObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
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
  if(ok)
  {
    if((fP = fopen("tstOut0", "w")) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
		     "%s: failed to open test output file (%s)\n",
		     *argv, "tstOut");
    }
  }
  if(ok)
  {
    switch(inObj->domain.ctr->model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
	(void )fprintf(fP, "%%!\n"
			   "0.001 setlinewidth\n"
			   "255 0 0 setrgbcolor\n");
	errNum = WlzGMModelTestOutPS(inObj->domain.ctr->model, fP,
				     tstOutOffset2D, tstOutScale2D, 1);
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	errNum = WlzGMModelTestOutVTK(inObj->domain.ctr->model, fP);
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Failed to write to test file.\n",
		     *argv);
    }
  }
  if(fP)
  {
    (void )fclose(fP);
    fP = NULL;
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
    switch(inObj->domain.ctr->model->type)
    {
      case WLZ_GMMOD_2I: /* FALLTHROUGH */
      case WLZ_GMMOD_2D:
        (void )fprintf(fP, "%%!\n"
			   "0.001 setlinewidth\n"
			   "0 64 64 setrgbcolor\n");
	errNum = WlzGMModelTestOutPS(inObj->domain.ctr->model, fP,
				     tstOutOffset2D, tstOutScale2D, 1);
	break;
      case WLZ_GMMOD_3I: /* FALLTHROUGH */
      case WLZ_GMMOD_3D:
	errNum = WlzGMModelTestOutVTK(inObj->domain.ctr->model, fP);
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )fprintf(stderr, "%s: Failed to write to test file.\n",
		     *argv);
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
      "  -s  Threshold shell size (edges in 2D, loops in 3D) for small\n"
      "      shells to be removed from the model.\n"
      "  -h  Prints this usage information.\n"
      "Test for various Woolz Geometric Model filters.\n");
  }
  return(!ok);
}

#endif /* WLZ_GMFILTER_TEST */
