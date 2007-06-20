#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzLBTDomain_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         binWlzTst/WlzTstLBTDomain.c
* \author       Bill Hill
* \date         June 2007
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
* \brief	Test for creating and manipulating linear binary
* 		tree domains.
* \ingroup	Tst
* \todo         -
* \bug          None known.
*/

#include <Wlz.h>
#include <limits.h>
#include <float.h>

typedef enum _WlzLBTTestOutType
{
  WLZ_LBTDOMAIN_TEST_OUT_IDOM,
  WLZ_LBTDOMAIN_TEST_OUT_IDX,
  WLZ_LBTDOMAIN_TEST_OUT_TXT,
  WLZ_LBTDOMAIN_TEST_OUT_VTK
} WlzLBTTestOutType;

extern char	*optarg;
extern int	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		dim = 0,
  		balance = 0,
  		option,
  		ok = 1,
		maxNodSz = INT_MAX,
		maxBndNodSz = INT_MAX,
		usage = 0;
  WlzLBTTestOutType outType = WLZ_LBTDOMAIN_TEST_OUT_VTK;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  FILE		*fP = NULL;
  WlzObject	*idxObj = NULL,
  		*inObj = NULL,
  		*outObj = NULL;
  WlzDomain	lDom,
  		outDom;
  WlzValues	nullVal;
  char		*inFileStr,
  		*outFileStr;
  static char	optList[] = "bditvho:",
  		outFileStrDef[] = "-",
  		inFileStrDef[] = "-";

  lDom.core = NULL;
  opterr = 0;
  nullVal.core = NULL;
  outFileStr = outFileStrDef;
  inFileStr = inFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'b':
        balance = 1;
	break;
      case 'd':
        outType = WLZ_LBTDOMAIN_TEST_OUT_IDOM;
	break;
      case 'i':
        outType = WLZ_LBTDOMAIN_TEST_OUT_IDX;
        break;
      case 't':
        outType = WLZ_LBTDOMAIN_TEST_OUT_TXT;
        break;
      case 'v':
	outType = WLZ_LBTDOMAIN_TEST_OUT_VTK;
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
  if((inFileStr == NULL) || (*inFileStr == '\0') ||
     (outFileStr == NULL) || (*outFileStr == '\0'))
  {
    ok = 0;
    usage = 1;
  }
  if(ok && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      inFileStr = *(argv + optind);
    }
  }
  ok = !usage;
  if(ok)
  {
    errNum = WLZ_ERR_READ_EOF;
    if((inFileStr == NULL) || (*inFileStr == '\0') ||
       ((fP = (strcmp(inFileStr, "-")? fopen(inFileStr, "r"):
                                       stdin)) == NULL) ||
       ((inObj= WlzAssignObject(WlzReadObj(fP, &errNum), NULL)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to read object from file %s (%s).\n",
		     *argv, inFileStr, errMsg);
    }
    if(fP && strcmp(inFileStr, "-"))
    {
      (void )fclose(fP);
    }
  }
  if(ok)
  {
    if(inObj == NULL)
    {
      errNum = WLZ_ERR_OBJECT_NULL;
    }
    else if(inObj->domain.core == NULL)
    {
      errNum = WLZ_ERR_DOMAIN_NULL;
    }
    else
    {
      switch(inObj->type)
      {
        case WLZ_2D_DOMAINOBJ:
	  dim = 2;
	  break;
        case WLZ_3D_DOMAINOBJ:
	  dim =3;
	  break;
	default:
	  errNum = WLZ_ERR_OBJECT_TYPE;
	  break;
      }
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: invalid object read from file %s (%s).\n",
		     *argv, inFileStr, errMsg);
    }
  }
  if(ok)
  {
    lDom = WlzLBTDomainFromObj(inObj, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to compute LBT domain (%d).\n",
		     *argv, errMsg);
    }
  }
  if(ok && (balance || (outType == WLZ_LBTDOMAIN_TEST_OUT_IDX)))
  {
    switch(lDom.core->type)
    {
      case WLZ_LBTDOMAIN_2D:
        idxObj = WlzLBTMakeNodeIndexObj2D(lDom.l2, inObj->domain.i, &errNum);
        if((errNum == WLZ_ERR_NONE) && balance)
        {
	   errNum = WlzLBTBalanceDomain2D(lDom.l2, idxObj, maxNodSz,
	                                  maxBndNodSz);
        }
	break;
      case WLZ_LBTDOMAIN_3D:
        idxObj = WlzLBTMakeNodeIndexObj3D(lDom.l3, inObj->domain.p, &errNum);
        if((errNum == WLZ_ERR_NONE) && balance)
        {
	   errNum = WlzLBTBalanceDomain3D(lDom.l3, idxObj, maxNodSz,
	                                  maxBndNodSz);
        }
	break;
    }
  }
  if(ok)
  {
    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileStr, "-")?
              fopen(outFileStr, "w"): stdout)) == NULL))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
  }
  if(ok)
  {
    switch(outType)
    {
      case WLZ_LBTDOMAIN_TEST_OUT_IDX:
        errNum = WlzWriteObj(fP, idxObj);
        break;
      case WLZ_LBTDOMAIN_TEST_OUT_IDOM:
        switch(lDom.core->type)
	{
          case WLZ_LBTDOMAIN_2D:
            outDom.i = WlzLBTDomainToIDomain(lDom.l2, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      outObj = WlzMakeMain(WLZ_2D_DOMAINOBJ, outDom, nullVal,
				   NULL, NULL, &errNum);
	    }
	    break;
          case WLZ_LBTDOMAIN_3D:
            outDom.p = WlzLBTDomainToPDomain(lDom.l3, &errNum);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      outObj = WlzMakeMain(WLZ_3D_DOMAINOBJ, outDom, nullVal,
				   NULL, NULL, &errNum);
	    }
	    break;
	  default:
	    errNum = WLZ_ERR_DOMAIN_TYPE;
	    break;
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzWriteObj(fP, outObj);
	}
	if(outObj)
	{
	  (void )WlzFreeObj(outObj);
	}
	else if(outDom.i)
	{
	  (void )WlzFreeIntervalDomain(outDom.i);
	}
	break;
      case WLZ_LBTDOMAIN_TEST_OUT_TXT:
        errNum = WlzLBTTestOutputNodesTxt(fP, lDom);
	break;
      case WLZ_LBTDOMAIN_TEST_OUT_VTK:
        errNum = WlzLBTTestOutputNodesVtk(fP, lDom);
	break;
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write to file %s (%s).\n",
		     *argv, outFileStr, errMsg);
    }
  }
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(idxObj);
  (void )WlzFreeLBTDomain2D(lDom.l2);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s%s",
    *argv,
    " [-o<output object>] [-h] [-b] [-i|t|v] [<input object>]\n"
    "Options:\n"
    "  -h  Prints this usage information.\n"
    "  -o  Output object file name.\n"
    "  -b  Produce a balanced tree.\n"
    "  -d  Woolz interval or plane domain output.\n"
    "  -i  Index object output.\n"
    "  -t  Text output.\n"
    "  -v  VTK output.\n");
  }
  return(!ok);
}
