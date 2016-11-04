#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzCompound_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlz/WlzCompound.c
* \author       Richard Baldock, Bill Hill
* \date         March 1999
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
* \brief        Generates a Woolz compound object from objects.
* \ingroup	BinWlz
*
* \par Binary
* \ref wlzcompound "WlzCompound"
*/

/*!
\ingroup      BinWlz
\defgroup     wlzcompound WlzCompound
\par Name
WlzCompound - Create a Woolx compounf object.
\par Synopsis
\verbatim
WlzCompound [-t #] [-o <output file>] [-n <num objects>] [-h] [-v]

\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-t</b></td>
    <td>Type of compound array, either 1 or 2 (default) for
        WLZ_COMPOUND_ARR_1 (all objects same type) or
	WLZ_COMPOUND_ARR_2 (object may be different type)</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output object.</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>Maximum number of objects, default 1024.</td>
  </tr>
  <tr>
    <td><b>-h</b></td>
    <td>Help - print help message</td>
  </tr>
  <tr>
    <td><b>-v</b></td>
    <td>Verbose operation</td>
  </tr>
</table>
Reads woolz objects from standard input and writes a compound
 object to standard out.

\par Description
Create a woolz compound object from a given input stream.
Only objects of the same type as the first object read in are included
and the compound object will be of type WLZ_COMPOUND_ARR_2 unless
WLZ_COMPOUND_ARR_1 is specified.

\par Examples
\verbatim
#
# read the object tiled.wlz then up to 149 woolz objects matching the
# pattern (dom*.wlz) in the current directory and creates a compound
# object (compound.wlz)
#
cat dom*.wlz | WlzCompound -o compound.wlz -n 150 tiled.wlz -

\endverbatim

\par File
\ref WlzCompound.c "WlzCompound.c"
\par See Also
\ref wlzexplode "WlzExplode(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <Wlz.h>

/* Externals required by getopt  - not in ANSI C standard */
#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static WlzErrorNum ReadFile(int nmax, int *n, WlzObjectType *type,
			    WlzObject **objs, FILE *fP)
{
  int		nObj = 0;
  WlzObject	*o;
  WlzErrorNum 	errNum = WLZ_ERR_NONE;

  while((errNum == WLZ_ERR_NONE) &&
	(*n < nmax) &&
	((o = WlzAssignObject(
	      WlzReadObj(fP, &errNum), NULL)) != NULL))
  {
    ++nObj;
    if((*type == -1) &&
       (o->type == WLZ_2D_DOMAINOBJ || o->type == WLZ_3D_DOMAINOBJ))
    {
      *type = o->type;
    }
    if((o->type == *type) || (o->type == WLZ_EMPTY_OBJ))
    {
      objs[*n] = WlzAssignObject(o, NULL);
      ++*n;
    }
    (void )WlzFreeObj(o);
  }
  /* Demand that files have at least on object in them. */
  if((nObj > 0) && ((errNum == WLZ_ERR_EOO) || (errNum == WLZ_ERR_READ_EOF)))
  {
    errNum = WLZ_ERR_NONE;
  }
  return(errNum);
}

int main(int	argc,
	 char	**argv)
{

  int 		n = 0,
		ok = 1,
		usage = 0,
  		nmax = 1024,
		verbose = 0;
  WlzObject	**objlist = NULL;
  WlzCompoundArray *cobj = NULL;
  WlzObjectType	type = (WlzObjectType) -1;
  WlzObjectType	cpdType = WLZ_COMPOUND_ARR_2;
  char		*outFileStr;
  int		option;
  const char	*errMsg;
  static char	defFileStr[] = "-";
  static char	optList[] = "n:o:t:hv";
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  opterr = 0;
  outFileStr = defFileStr;
  /* read the argument list and check for an input file */
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != EOF))
  {
    switch( option )
    {
      case 'n':
	if((sscanf(optarg, "%i", &nmax) != 1) || (nmax < 1)) 
	{
	  (void )fprintf(stderr, "%s: nmax = %d is invalid\n", argv[0], nmax);
	  usage = 1;
	}
	break;
      case 'o':
	outFileStr = optarg;
	break;
      case 't':
	{
	  int	t;

	  if((sscanf(optarg, "%d", &t) != 1) || ((t != 1) && (t != 2)))
	  {
	    usage = 1;
	  }
	  else
	  {
	    cpdType = (t == 1)? WLZ_COMPOUND_ARR_1: WLZ_COMPOUND_ARR_2;
	  }
	}
	break;
      case 'v':
	verbose = 1;
	break;
      case 'h':
      default:
	usage = 1;
	break;
    }
  }
  ok = !usage;
  /* allocate space for the object pointers */
  if(ok) 
  {
    if((objlist = (WlzObject **)
	          AlcMalloc(sizeof(WlzObject *) * nmax)) == NULL )
    {
      (void )fprintf(stderr, "%s: memory allocation failed.\n", argv[0]);
      ok = 0;
    }
  }
  /* read objects accumulating compatible types */
  if(ok)
  {
    int		idx,
    		haveRead = 0;

    idx = optind;
    for(idx = optind; ok && (idx < argc); ++idx)
    {
      char	*fName;

      if(((fName = argv[idx]) != NULL) && (*fName != '\0'))
      {
	int	isNamedFile;
	FILE	*fP = NULL;

	haveRead = 1;
	isNamedFile = strcmp(fName, "-");
	if(verbose)
	{
	  (void )fprintf(stderr, "%s: Reading object(s) from %s\n",
	                 argv[0], fName);
	}
        if(isNamedFile)
	{
	  if((fP = fopen(fName, "r")) == NULL)
	  {
	    ok = 0;
	    (void )fprintf(stderr,
	                   "%s: failed to read object from file %s.\n",
			   argv[0], (isNamedFile)? fName: "stdin");
	  }
	}
	else
	{
	  fP = stdin;
	}
	if(ok)
	{
	  if((errNum = ReadFile(nmax, &n, &type, objlist, fP)) != WLZ_ERR_NONE)
	  {
	    (void )fprintf(stderr,
	                   "%s: failed to read object from file %s\n",
			   argv[0], fName);
	  }
	}
	if(isNamedFile)
	{
	  (void )fclose(fP);
	}
      }
    }
    if(ok && (haveRead == 0) )
    {
      if((errNum = ReadFile(nmax, &n, &type, objlist, stdin)) != WLZ_ERR_NONE)
      {
        (void )fprintf(stderr,
	               "%s: Failed to read object from stdin.\n",
		       argv[0]);
      }
    }
  }
  if(ok)
  {
    if(n < 1)
    {
      (void )fprintf(stderr,
		 "%s: Refusing to create a compound object with no objects.\n",
		 argv[0]);
    }
  }
  /* Create compound object. */
  if(ok)
  {
    if(verbose)
    {
      (void )fprintf(stderr, "%s: Creating compound object.\n", argv[0]);
    }
    if((cobj = WlzMakeCompoundArray(cpdType, 3, n, objlist,
		     (cpdType == WLZ_COMPOUND_ARR_1)? type: WLZ_NULL,
		     &errNum)) == NULL )
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: failed to create compound object (%s).\n",
		     argv[0], errMsg);
    }
  }
  /* Write compound object. */
  if(ok)
  {
    int		isNamedFile;
    FILE	*fP = NULL;

    if(verbose)
    {
      (void )fprintf(stderr, "%s: Writing compound object to file %s.\n",
      		     argv[0], outFileStr);
    }
    isNamedFile = strcmp(outFileStr, "-");
    if(isNamedFile)
    {
      if((fP = fopen(outFileStr, "w")) == NULL)
      {
        errNum = WLZ_ERR_FILE_OPEN;
      }
    }
    else
    {
      fP = stdout;
    }
    if(errNum == WLZ_ERR_NONE)
    {
      errNum = WlzWriteObj(fP, (WlzObject *)cobj);
    }
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr, "%s: Failed to write compound object (%s).\n",
		     argv[0], errMsg);
    }
    if(isNamedFile && fP)
    {
      (void )fclose(fP);
    }
  }
  if(objlist)
  {
    int	i;

    for(i = 0; i < n; ++i)
    {
      (void )WlzFreeObj(objlist[i]);
    }
    AlcFree(objlist);
  }
  (void )WlzFreeObj((WlzObject *)cobj);
  if(usage)
  {
    (void )fprintf(stderr,
	"Usage:\t%s [-h] [-n#]\n"
	"Generates a woolz compound object from the given input objects.\n"
	"Only objects of the same type as the first object read in are\n"
	"included and the compound object will be of type\n"
	"WLZ_COMPOUND_ARR_2 unless WLZ_COMPOUND_ARR_1 is specified.\n"
	"Objects are read from stdin and written to stdout unless\n"
	"filenames are given on the command line. Tiled objects must\n"
	"always be read from a named file.\n"
	"Version: %s\n"
	"Example:\n"
	"  cat dom*.wlz | %s -o compound.wlz -n 150 tiled.wlz -\n"
	"Reads the object tiled.wlz then up to 149 woolz objects matching\n"
	"the pattern (dom*.wlz) in the current directory and creates a\n"
	"compound object (compound.wlz).\n"
	"Options:\n"
	"\t  -h        Help - prints this usage message\n"
	"\t  -o        Output object.\n"
	"\t  -t        Type of compound array, either 1 or 2 (default) for\n"
	"\t            WLZ_COMPOUND_ARR_1 (all objects same type) or\n"
	"\t            WLZ_COMPOUND_ARR_2 (object may be different type).\n"
	"\t  -n#       Maximum number of objects -default=1024\n"
	"\t  -v        Verbose output to stderr.\n",
	argv[0],
	WlzVersion(),
	argv[0]);
  }
  return(!ok);
}
#endif /* DOXYGEN_SHOULD_SKIP_THIS */
