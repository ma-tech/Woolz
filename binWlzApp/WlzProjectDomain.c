#if defined(__GNUC__)
#ident "University of Edinburgh $Id$"
#else
static char _WlzProjectDomain_c[] = "University of Edinburgh $Id$";
#endif
/*!
* \file         binWlzApp/WlzProjectDomain.c
* \author       Bill Hill
* \date         June 2011
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
* \brief	Projects a 3D domain onto a single section 3D projection
* 		as used for EMAGE wholemounts.
* \ingroup	BinWlzApp
*
* \par Binary
* \ref wlzprojectdomain "WlzProjectDomain"
*/

/*!
\ingroup BinWlzApp
\defgroup wlzprojectdomain WlzProjectDomain
\par Name
WlzProjectDomain - Projects a given 3D domain into a section.
\par Synopsis
\verbatim
WlzProjectDomain [-h] [-o<out object>] [-p<proj>] [-s[<sec>]
		 [-n#] [-m<mask1>[,<mask2>,...,<maskn>]]
		 [-t<tr1>[,<tr2>,...,<trn>]] <input object>
\endverbatim
\par Options
<table width="500" border="0">
  <tr>
    <td><b>-h</b></td>
    <td>Prints usage information.</td>
  </tr>
  <tr>
    <td><b>-o</b></td>
    <td>Output file for the projected section view.</td>
  </tr>
  <tr>
    <td><b>-p</b></td>
    <td>Projection from the 3D space of the input object and mask
        domains onto a plane. This must be defined by a bibfile (as
        saved by MAPaint).</td>
  </tr>
  <tr>
    <td><b>-s</b></td>
    <td>Section transform which maps the plane back into a 3D space.
        This must be defined by a bibfile (as saved by MAPaint).</td>
  </tr>
  <tr>
    <td><b>-n</b></td>
    <td>The number of parts to be used for the projection.
        If the number of parts is one (which is the default) then
        no mask is required. In this case the default in-plane
        transform is the identity transform.</td>
  </tr>
  <tr>
    <td><b>-m</b></td>
    <td>Mask domains for each of the parts. If the mask is to cover
        the whole of the input object's domain then the string
        "null" (without the quotes) may be used.</td>
  </tr>
  <tr>
    <td><b>-t</b></td>
    <td>Transforms withing the plane for each part. If the identity
        transform is required then the string "null" (without
        the quotes) may be used.</td>
  </tr>
</table>
\par Description
WlzProjectDomain
projects the given input 3D (spatial domain) object into a section,
as used for EMAGE wholemount views.
The input domain may be projected and transformed in parts to allow
the separation of regions that may otherwise be obscured in a
projected view.
The list of masks and in-plane transforms should be comma separated
and a missing string represents either the universal domain or an
identity transform.
\par Examples
\verbatim
./WlzProjectDomain -o out.wlz -p prj.bib -s sec.bib domain.wlz
\endverbatim
The 3D domain in the file domain.wlz is projected onto a plane
using the projecttion transform defined in the bibfile prj.bib
this is then made a section defined by the section transform
read from the bibfile sec.bib.

\verbatim
WlzProjectDomain -o out.wlz -p prj.bib -s sec.bib -n 2 \
        -m body-msk.wlz,tail-msk.wlz \
        -t body-tr.wlz,tail-tr.wlz domain.wlz
\endverbatim
The 3D domain in the file domain.wlz is masked using the 3D domains
read from the files body-msk.wlz and tail-msk.wlz. The resulting
domains are projected onto a plane using the projection defined in
the bibfile  prj.bib. Each of these projected domains is then
transformed independently using the transforms read from the files
body-tr.wlz and tail-tr.wlz. The union of these transformed 2D
domains is then made a section defined by the section transform
read from the bibfile sec.bib.
\par File
\par See Also
\ref BinWlzApp "WlzIntro(1)"
*/

#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include <WlzExtFF.h>

#ifdef __STDC__ /* [ */
extern int      getopt(int argc, char * const *argv, const char *optstring);

extern int      optind, opterr, optopt;
extern char     *optarg;
#endif /* __STDC__ ] */

static int			WlzProjDomStripSpaces(
				  char *s);
static int 			WlzProjDomParseNameList(
				  char ***fNames,
				  char *fList,
				  int newCnt);
static WlzObject 		*WlzProjDomReadObjFromNamedFile(
				  const char *fName,
				  WlzErrorNum *dstErr);
static WlzThreeDViewStruct 	*WlzProjDomReadViewFromNamedFile(
				  const char *fName,
                                  WlzErrorNum *dstErr);

int		main(int argc, char *argv[])
{
  int		ok,
		idx,
		option,
		nMsk = 0,
		nTr = 0,
		nPart = 1,
  		usage = 0;
  char		*inFileName,
		*outFileName,
  		*prjFileName,
		*secFileName;
  char		**mskFileNames = NULL,
  		**trFileNames = NULL;
  WlzObject	*inObj = NULL,
  		*outObj= NULL;
  WlzObject	**mskObjs = NULL,
  		**trObjs = NULL;
  WlzThreeDViewStruct *prjView = NULL,
  		*secView = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const char	*errMsg;
  static char	optList[] = "ho:p:s:n:m:t:",
  		defFileName[] = "-";

  opterr = 0;
  inFileName = defFileName;
  outFileName = defFileName;
  prjFileName = defFileName;
  secFileName = defFileName;
  while((usage == 0) && ((option = getopt(argc, argv, optList)) != -1))
  {
    switch(option)
    {
      case 'o':
        outFileName = optarg;
	break;
      case 'p':
        prjFileName = optarg;
	break;
      case 's':
        secFileName = optarg;
	break;
      case 'n':
        if((sscanf(optarg, "%d\n", &nPart) != 1) || (nPart < 1))
	{
	  usage = 1;
	}
        break;
      case 'm':
	if((nMsk = WlzProjDomParseNameList(&mskFileNames, optarg, nPart)) < 1)
	{
	  usage = 1;
	}
        break;
      case 't':
	if((nTr = WlzProjDomParseNameList(&trFileNames, optarg, nPart)) < 1)
	{
	  usage = 1;
	}
        break;
      case 'h': /* FALLTHROUGH */
      default:
        usage = 1;
	break;
    }
  }
  if(nMsk != nTr)
  {
    usage = 0;
  }
  if((usage == 0) && (optind < argc))
  {
    if((optind + 1) != argc)
    {
      usage = 1;
    }
    else
    {
      inFileName = *(argv + optind);
    }
  }
  ok = (usage == 0);
  if(ok)
  {
    if(((prjView = WlzAssign3DViewStruct(
                   WlzProjDomReadViewFromNamedFile(prjFileName, &errNum),
	           NULL)) == NULL) ||
       ((secView = WlzAssign3DViewStruct(
                   WlzProjDomReadViewFromNamedFile(secFileName, &errNum),
	           NULL)) == NULL) ||
       ((inObj = WlzAssignObject(
                 WlzProjDomReadObjFromNamedFile(inFileName, &errNum),
		 NULL)) == NULL))
    {
      ok = 0;
    }
    if(ok)
    {
      if(((trObjs = (WlzObject **)
                    AlcCalloc(nPart, sizeof(WlzObject *))) == NULL) ||
         ((mskObjs = (WlzObject **)
                     AlcCalloc(nPart, sizeof(WlzObject *))) == NULL))
      {
        ok = 0;
      }
    }
    if(ok)
    {
      for(idx = 0; idx < nPart; ++idx)
      {
        if(strcmp(trFileNames[idx], "null"))
	{
	  if((trObjs[idx] = WlzAssignObject(
	               WlzProjDomReadObjFromNamedFile(trFileNames[idx],
		                                      &errNum), NULL)) == NULL)
	  {
	    ok = 0;
	    break;
	  }
	}
      }
    }
    if(ok)
    {
      for(idx = 0; idx < nPart; ++idx)
      {
        if(strcmp(mskFileNames[idx], "null"))
	{
	  if((mskObjs[idx] = WlzAssignObject(
	               WlzProjDomReadObjFromNamedFile(mskFileNames[idx],
		                                      &errNum), NULL)) == NULL)
	  {
	    ok = 0;
	    break;
	  }
	}
      }
    }
    if(ok == 0)
    {
      (void )fprintf(stderr,
                     "%s: Failed to read input view or object.\n",
		     *argv);
    }
  }
  AlcFree(trFileNames);
  AlcFree(mskFileNames);
  if(ok)
  {
    outObj = WlzProj3DToSection(inObj, nPart, mskObjs, prjView, nPart, trObjs,
                                secView, &errNum);
    if(errNum != WLZ_ERR_NONE)
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
      "%s: Failed to create projected section (%s).\n",
      *argv, errMsg);
    }
  }
  if(ok)
  {
    FILE 	*fP = NULL;

    errNum = WLZ_ERR_WRITE_EOF;
    if(((fP = (strcmp(outFileName, "-")?
              fopen(outFileName, "w"): stdout)) == NULL) ||
       ((errNum = WlzWriteObj(fP, outObj)) != WLZ_ERR_NONE))
    {
      ok = 0;
      (void )WlzStringFromErrorNum(errNum, &errMsg);
      (void )fprintf(stderr,
                     "%s: Failed to write output object (%s).\n",
                     *argv, errMsg);
    }
    if(fP && strcmp(outFileName, "-"))
    {
      (void )fclose(fP);
    }
  }
  (void )WlzFree3DViewStruct(prjView);
  (void )WlzFree3DViewStruct(secView);
  (void )WlzFreeObj(inObj);
  (void )WlzFreeObj(outObj);
  if(trObjs)
  {
    for(idx = 0; idx < nPart; ++idx)
    {
      WlzFreeObj(trObjs[idx]);
    }
    AlcFree(trObjs);
  }
  if(mskObjs)
  {
    for(idx = 0; idx < nPart; ++idx)
    {
      WlzFreeObj(mskObjs[idx]);
    }
    AlcFree(mskObjs);
  }
  if(usage)
  {
    (void )fprintf(stderr,
    	"Usage: %s%s%s%s%s%s%s%s\n",
	*argv,
	" [-h] [-o<out object>] [-p<proj>] [-s[<sec>]\n"
	"        [-n#] [-m<mask1>[,<mask2>,...,<maskn>]]\n"
	"        [-t<tr1>[,<tr2>,...,<trn>]] <input object>\n"
	"Projects the given input 3D (spatial domain) object into a section,\n"
	"as used for EMAGE wholemount views.\n"
	"The input domain may be projected and transformed in parts to allow\n"
	"the separation of regions that may otherwise be obscured in a\n"
	"projected view.\n"
	"Version: ",
	WlzVersion(),
	"\n"
	"Options:\n"
	"  -h    Help - prints this usage message\n"
	"  -o    Output file for the projected section view.\n"
	"  -p    Projection from the 3D space of the input object and mask\n"
	"        domains onto a plane. This must be defined by a bibfile (as\n"
	"        saved by MAPaint).\n"
	"  -s    Section transform which maps the plane back into a 3D space.\n"
	"        This must be defined by a bibfile (as saved by MAPaint).\n"
	"  -n    The number of parts to be used for the projection.\n"
	"        If the number of parts is one (which is the default) then\n"
	"        no mask is required. In this case the default in-plane\n"
	"        transform is the identity transform.\n"
	"  -m    Mask domains for each of the parts. If the mask is to cover\n"
	"        the whole of the input object's domain then the string\n"
	"        \"null\" (without the quotes) may be used.\n"
	"  -t    Transforms withing the plane for each part. If the identity\n"
	"        transform is required then the string \"null\" (without\n"
	"        the quotes) may be used.\n"
	"The list of masks and in-plane transforms should be comma separated\n"
	"and a missing string represents either the universal domain or an\n"
	"identity transform.\n"
	"Examples:\n",
	*argv,
	" -o out.wlz -p prj.bib -s sec.bib -n 2 \\\n"
	"        -m body-msk.wlz,tail-msk.wlz \\\n"
	"        -t body-tr.wlz,tail-tr.wlz domain.wlz\n"
	"The 3D domain in the file domain.wlz is masked using the 3D domains\n"
	"read from the files body-msk.wlz and tail-msk.wlz. The resulting\n"
	"domains are projected onto a plane using the projection defined in\n"
	"the bibfile  prj.bib. Each of these projected domains is then\n"
	"transformed independently using the transforms read from the files\n"
	"body-tr.wlz and tail-tr.wlz. The union of these transformed 2D\n"
	"domains is then made a section defined by the section transform\n"
	"read from the bibfile sec.bib.\n",
        *argv,
	" -o out.wlz -p prj.bib -s sec.bib domain.wlz\n"
	"The 3D domain in the file domain.wlz is projected onto a plane\n"
	"using the projecttion transform defined in the bibfile prj.bib\n"
	"this is then made a section defined by the section transform\n"
	"read from the bibfile sec.bib.\n");
  }
  exit(!ok);
}

/*!
* \return	New object or NULL on error.
* \ingroup	BinWlzApp
* \brief	Reads an object from the named file
* \param	fName			File name ("-" implies stdin).
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzObject *WlzProjDomReadObjFromNamedFile(const char *fName,
				WlzErrorNum *dstErr)
{
  WlzObject     *obj = NULL;
  FILE          *fP = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((((fP = (strcmp(fName, "-")?
              fopen(fName, "r"): stdin)) == NULL) ||
      ((obj = WlzReadObj(fP, &errNum)) == NULL)) &&
      (errNum == WLZ_ERR_NONE))
  {
    errNum = WLZ_ERR_READ_EOF;
  }
  if(fP)
  {
    if(strcmp(fName, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(obj);
}

/*!
* \return	New view or NULL on error.
* \ingroup	BinWlzApp
* \brief	Reads a 3D view from the named file.
* \param	fName			File name ("-" implies stdin).
* \param	dstErr			Destination error pointer, may be NULL.
*/
static WlzThreeDViewStruct *WlzProjDomReadViewFromNamedFile(const char *fName,
                                WlzErrorNum *dstErr)
{
  WlzThreeDViewStruct *view = NULL;
  FILE          *fP = NULL;
  WlzErrorNum   errNum = WLZ_ERR_NONE;

  if((((fP = (strcmp(fName, "-")?
              fopen(fName, "r"): stdin)) == NULL) ||
      ((view = WlzEffBibRead3DView(fP, NULL, &errNum)) == NULL)) &&
      (errNum == WLZ_ERR_NONE))
  {
    errNum = WLZ_ERR_READ_EOF;
  }
  if(fP)
  {
    if(strcmp(fName, "-"))
    {
      (void )fclose(fP);
    }
    fP = NULL;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(view);
}

/*!
* \return	Number of file names in given string or NULL on error.
* \ingroup	BinWlzApp
* \brief	Tries to parse file names from the given comma separated
* 		list. It's an error if there are insufficient file names.
* \param	fNames			Source and destination pointer for
* 					the array of file names.
* \param	fList			Comma separated list of file names.
* \param	newCnt			Number of file names expected.
*/
static int 	WlzProjDomParseNameList(char ***fNames, char *fList,
				int newCnt)
{
  if((fNames == NULL) || (newCnt < 1))
  {
    newCnt = 0;
  }
  else
  {
    if((*fNames = AlcRealloc(*fNames, sizeof(char *) * newCnt)) == NULL)
    {
      newCnt = 0;
    }
  }
  if(newCnt > 0)
  {
    int		idx;
    char	*tok;
    
    tok = strtok(fList, ",");
    for(idx = 0; (tok != NULL) && (idx < newCnt); ++idx)
    {
      if(WlzProjDomStripSpaces(tok) < 1)
      {
        tok = NULL;
      }
      else
      {
        (*fNames)[idx] = tok;
	tok = strtok(NULL, ",");
      }
    }
    if(idx != newCnt)
    {
      newCnt = 0;
    }
  }
  return(newCnt);
}

/*!
* \return	Length of striped string.
* \ingroup	BinWlzApp
* \brief	Strips leading and trailing space characters from the given
* 		string.
* \param	s			String to strip space characters from.
*/
static int	WlzProjDomStripSpaces(char *s)
{
  int		l = 0;

  if(s)
  {
    char	*s0,
    		*s1;

    s1 = s0 = s;
    while(*s1 && isspace(*s1))
    {
      ++s1;
    }
    while(*s1 && !isspace(*s1))
    {
      ++l;
      *s0++ = *s1++;
    }
    s0 = '\0';
  }
  return(l);
}

#endif /* DOXYGEN_SHOULD_SKIP_THIS */



