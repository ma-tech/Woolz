#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzConstruct3D.c
* \author       Bill Hill
* \date         April 2003
* \version      $Id$
* \note
*               Copyright
*               2003 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief	Constructs a 3D Woolz domain object from the
*		given 2D Woolz domain objects.
* \todo         -
* \bug          None known.
*/

#include <stdio.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <Wlz.h>

extern char 	*optarg;
extern int 	optind,
		opterr,
		optopt;

int		main(int argc, char *argv[])
{
  int		idx,
  		ok = 1,
  		option,
  		usage = 0,
		nFiles,
		plane1 = 0;
  FILE		*fP = NULL;
  char		*fStr,
  		*outObjFileStr;
  char		**inFileStr = NULL;
  const char	*errMsgStr;
  WlzDVertex3	voxSz;
  WlzObject	*obj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  static char   optList[] = "ho:p:s:";
  const char    outObjFileStrDef[] = "-";

  opterr = 0;
  voxSz.vtX = voxSz.vtY = voxSz.vtZ = 1.0;
  outObjFileStr = (char *)outObjFileStrDef;
  while(ok && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch(option)
    {
      case 'o':
        outObjFileStr = optarg;
	break;
      case 'p':
	if(sscanf(optarg, "%d", &plane1) != 1)
	{
	  usage = 1;
	}
        break;
      case 's':
	if(sscanf(optarg, "%lg,%lg,%lg",
	          &(voxSz.vtX), &(voxSz.vtY), &(voxSz.vtX)) != 3)
	{
	  usage = 1;
	}
        break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = !usage;
  if(ok)
  {
    if((outObjFileStr == NULL) || (*outObjFileStr == '\0'))
    {
      ok = 0;
      usage = 1;
    }
  }
  if(ok)
  {
    if((nFiles = argc - optind) <= 0)
    {
      ok = 0;
      usage = 1;
    }
  }
  if(ok)
  {
    if((inFileStr = (char **)AlcCalloc(nFiles, sizeof(char *))) == NULL)
    {
      ok = 0;
      (void )fprintf(stderr,
      		     "%s: failed to allocate list of files\n.",
		     *argv);
    }
  }
  if(ok)
  {
    for(idx = 0; idx < nFiles; ++idx)
    {
      fStr = *(argv + optind + idx);
      if(fStr && strcmp(fStr, "NULL"))
      {
        *(inFileStr + idx) = fStr;
      }
    }
  }
  if(ok)
  {
    obj = WlzConstruct3DObjFromFile(nFiles, inFileStr, plane1,
    				    voxSz.vtX, voxSz.vtY, voxSz.vtZ,
				    &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
      		     "%s: Failed to construct 3D object (%s).\n",
      		     argv[0],
		     errMsgStr);
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
    errNum = WlzWriteObj(fP, obj);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsgStr);
      (void )fprintf(stderr,
                     "%s: Failed to write output object, %s.\n",
		     argv[0], errMsgStr);
    }
  }
  (void )WlzFreeObj(obj);
  if(usage)
  {
    (void )fprintf(stderr,
    "Usage: %s [-h] [-o<output file>] [-p #] [-s #,#,#] [<input file>]\n"
    "  -h  Output this usage message.\n"
    "  -o  Output file name, default is the standard output.\n"
    "  -p  Coordinate of the first plane.\n"
    "  -s  Voxel size (x,y,z).\n"
    "Constructs a 3D Woolz domain object from the given 2D Woolz domain\n"
    "objects.\n"
    "Each file is read in turn and added to the 3D object. An empty plane\n"
    "can be specified by using the string NULL. Either all or none of the\n"
    "2D objects must have values. When the 2D objects have values then the\n"
    "background value of the first 2D object is set to be the background\n"
    "value of the 3D object.\n",
    argv[0]);
  }
  return(!ok);
}
