#if defined(__GNUC__)
#ident "MRC HGU $Id$"
#else
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC)
#pragma ident "MRC HGU $Id$"
#else
static char _WlzDrawDomain_c[] = "MRC HGU $Id$";
#endif
#endif
/*!
* \file         WlzDrawDomain.c
* \author       Bill Hill
* \date         August 2008
* \version      $Id$
* \par
* Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \par
* Copyright (C) 2008 Medical research Council, UK.
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
* \brief	Functions for composing a domain from simple drawing
* 		commands.
* \ingroup	WlzDomainOps
* \todo         -
* \bug          None known.
*/

#include <ctype.h>
#include <float.h>
#include <string.h>
#include <Wlz.h>

/*!
* \enum		_WlzDrawDomAct
* \ingroup      WlzDrawDomAct
* \brief	Drawing actions.
*		Typedef: ::WlzDrawDomAct.
*/
typedef enum _WlzDrawDomAct
{
  WLZ_DRAWDOM_ACT_NONE,			/*!< No action. */
  WLZ_DRAWDOM_ACT_DRAW,			/*!< Add drawn region. */
  WLZ_DRAWDOM_ACT_ERASE			/*!< Remove drawn region. */
} WlzDrawDomAct;

/*!
* \enum		_WlzDrawDomCmd
* \ingroup      WlzDrawDomCmd
* \brief	Commands in the drawing string which are understood.
*		Typedef: ::WlzDrawDomCmd.
*/
typedef enum _WlzDrawDomCmd
{
  WLZ_DRAWDOM_CMD_CIRCLE,		/*!< Draw a circle. */
  WLZ_DRAWDOM_CMD_END,			/*!< End of command string. */
  WLZ_DRAWDOM_CMD_ERROR,		/*!< Not a command but an error. */
  WLZ_DRAWDOM_CMD_LINE,			/*!< Draw a line. */
  WLZ_DRAWDOM_CMD_INIT,			/*!< Initial command of the string
                                             which identifies the string
					     as a drawing command string and
					     it's version. */
  WLZ_DRAWDOM_CMD_PEN,			/*!< Draws a polyline with rounded
                                             ends and joins. */
  WLZ_DRAWDOM_CMD_UNKNOWN		/*!< Unknown command. */
} WlzDrawDomCmd;

/*! \struct	_WlzDrawDomWSp
 *  \ingroup	WlzDrawDomWSp
 *  \brief	Workspace in which drawing commands are parsed.
 */
typedef struct _WlzDrawDomWSp
{
  char		*cmdStr;		/*!< Given command string, which is
  					     not modified. */
  int		cmdStrIdx;		/*!< Index into the command string
            				     for the start of the next
					     command. */
  char		*buf;			/*!< Dynamically sized buffer holding
  				             a copy of the current command
					     being parsed. */
  int		bufMax;			/*!< Current buffer size. */
#if defined _SVID_SOURCE || _BSD_SOURCE || _POSIX_C_SOURCE || _XOPEN_SOURCE
  char		*savPtr;		/*!< Used by strtok_r() if it's
                                             available, which makes the
					     code reentrant. */
#endif
} WlzDrawDomWSp;

static WlzErrorNum  		WlzDrawDomStrGetAct(
				  WlzDrawDomWSp *wSp,
				  WlzDrawDomAct *dstAct);
static WlzErrorNum 		WlzDrawDomStrGetVal1D(
				  WlzDrawDomWSp *wSp,
				  double *dstVal);
static WlzErrorNum 		WlzDrawDomStrGetVal2D(
				  WlzDrawDomWSp *wSp,
				  double *dstVal0,
				  double *dstVal1);
static WlzErrorNum 		WlzDrawDomStrGetValND(
				  WlzDrawDomWSp *wSp,
				  int nVal,
				  double *val);
static WlzErrorNum 		WlzDrawDomAddCircle(
				  WlzObject **curObj,
				  WlzDrawDomAct act,
				  double radius,
				  WlzDVertex2 centre);
static WlzErrorNum 		WlzDrawDomAddLine(
				  WlzObject **curObj,
				  WlzDrawDomAct act,
				  double width,
				  WlzDVertex2 pos0,
				  WlzDVertex2 pos1);
static WlzErrorNum 		WlzDrawDomAddObj(
				  WlzObject **curObj,
				  WlzDrawDomAct act,
				  WlzObject *obj);
static WlzDrawDomCmd		WlzDrawDomStrGetCmd(
				  WlzDrawDomWSp *wSp);
static WlzObject 		*WlzDrawDomain2D(
				  WlzDrawDomWSp *wSp,
				  WlzErrorNum *dstErr);

/*!
* \return	New 3D spatial domain object, an empty object or NULL on error.
* \ingroup	WlzDomainOps
* \brief	Constructs a 3D spatial domain object (without values) by
* 		drawing on a section in 3D space, using commands in the
* 		given command string.
*
* 		The command string must have the following syntax:
* 		\verbatim
 		<command string> = <init command>[<command>]+<end command>
 		<init command> = <ident>:<version>;
 		<end command> = END:;
 		<ident> = WLZ_DRAW_DOMAIN
 		<version> = 1
 		<command> = <command name>:[<parameter>[,<parameter>]*];
 		<command name> = PEN | LINE | CIRCLE
 		\endverbatim
* 		In addition to the init and end commands, the following
* 		drawing commands are recognised:
* 		\verbatim
		CIRCLE:<action>:<radius>,<x>,<y>;
		LINE:<action>:<width>,<x>,<y>,<x>,<y>;
                PEN:<action>,<width>,<x>,<y>[,<x>,<y>]*;
 		\endverbatim
* 		Where
* 		\verbatim
		<action> = DRAW | ERASE
 		\endverbatim
*		A simple example is:
* 		\verbatim
		WLZ_DRAW_DOMAIN:1;
		CIRCLE:DRAW,100,200,200;
		LINE:ERASE,10,0,0,200,200;
		PEN,DRAW,4,0,0,200,0,200,200,0,200,0,0;
		END:;
 		\endverbatim
*		The circle command draws or erases a filled circle which is
*		specified by it's radius and centre parameters.
*		The line command draws a rectangle using the given width
*		and end coordinates of the mid-line.
*		The pen command draws a polyline which is composed of a
*		series of rectangular segments, with each segment ending
*		in a semi-circular cap. The parameters of the pen command
*		are the width and line segment end point coordinates.
*		All widths, radii and coordinates may be in any floating
*		point format recognised by scanf().
*
*               Other commands may be present provided they have the
*               same syntax described above, but they will be ignored.
*
*               All white space characters are ignored.
*
* 		The drawn domain is offset by the given origin before being
* 		transformed to the section plane given by the view structure.
* \param	org			Origin of the 2D drawing frame with
* 					respect to the 2D Woolz object cut
* 					using the given view transform.
* \param	view			The section plane in 3D space on
* 					which the domain is drawn. May be
*					NULL, in which case the view transform
*					is the default transform created by
*					WlzMake3DViewStruct().
* \param	keep2D			Ignore the view struct and keep the
* 					object 2D if this flag is non-zero.
* \param	cmdStr			String with drawing commands using the
* 					syntax described above.
* \param	dstErrIdx		Used to return the index of the last
* 					character parsed in the command string.
* 					May be NULL.
* \param	dstErr			Destination error pointer, may be NULL.
*/
WlzObject	*WlzDrawDomainObj(WlzDVertex2 org, WlzThreeDViewStruct *view,
				  int keep2D, char *cmdStr, int *dstErrIdx,
				  WlzErrorNum *dstErr)
{
  double	version = 0.0;
  WlzDrawDomCmd	cmd;
  WlzDrawDomWSp	wSp;
  WlzObject	*obj0 = NULL,
		*obj1 = NULL,
  		*dwnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;
  const	double	versions[1] =
  {
    1.0
  };

  if(cmdStr == NULL)
  {
    errNum = WLZ_ERR_PARAM_NULL;
  }
  else
  {
    wSp.cmdStr = cmdStr;
    wSp.cmdStrIdx = 0;
    wSp.buf = NULL;
    wSp.bufMax = 0;
    cmd = WlzDrawDomStrGetCmd(&wSp);
  }
  if(cmd == WLZ_DRAWDOM_CMD_INIT)
  {
    errNum = WlzDrawDomStrGetVal1D(&wSp, &version);
  }
  else
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  if(errNum == WLZ_ERR_NONE)
  {
    if(fabs(version - versions[0]) < DBL_EPSILON)
    {
      obj0 = WlzDrawDomain2D(&wSp, &errNum);
    }
    else
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  if(errNum == WLZ_ERR_NONE)
  {
    obj1 = WlzShiftObject(obj0, org.vtX, org.vtY, 0.0, &errNum);
  }
  (void )WlzFreeObj(obj0);
  if(errNum == WLZ_ERR_NONE)
  {
    if(keep2D)
    {
      dwnObj = obj1;
      obj1 = NULL;
    }
    else
    {
      (void )WlzAssignObject(obj1, NULL);
      if(view == NULL)
      {
	view = WlzMake3DViewStruct(WLZ_3D_VIEW_STRUCT, &errNum);
	if(errNum == WLZ_ERR_NONE)
	{
	  view->view_mode = WLZ_UP_IS_UP_MODE;
	  errNum = WlzInit3DViewStruct(view, obj1);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  dwnObj = Wlz3DViewTransformObj(obj1, view, &errNum);
	}
        (void )WlzFree3DViewStruct(view); view = NULL;
      }
      else
      {
	errNum = WlzInit3DViewStruct(view, obj1);
	if(errNum == WLZ_ERR_NONE)
	{
	  dwnObj = Wlz3DViewTransformObj(obj1, view, &errNum);
	}
      }
    }
  }
  AlcFree(wSp.buf);
  (void )WlzFreeObj(obj1);
  if(dstErr)
  {
    *dstErrIdx = wSp.cmdStrIdx;
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dwnObj);
}

/*!
* \return	New Woolz spatial domain object drawn using the command
* 		string. 
* \ingroup	WlzDomainOps
* \brief	Creates a new Woolz spatial domain object, in which the
* 		domain is the result of the ordered drawing operations
* 		in the drawing command string.
* \param	wSp			Workspace data structure with
* 					a pointer to the drawing command
* 					string and an initialised command
* 					buffer (may be NULL).
* \param	dstErr			Destination error pointer, may be
* 					NULL.
*/
static WlzObject *WlzDrawDomain2D(WlzDrawDomWSp *wSp, WlzErrorNum *dstErr)
{
  double	radius,
  		width;
  WlzDVertex2	pos0,
  		pos1;
  WlzDrawDomCmd	cmd;
  WlzDrawDomAct act;
  WlzObject	*dwnObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Read and execute commands. */
  do
  {
    cmd = WlzDrawDomStrGetCmd(wSp);
    switch(cmd)
    {
      case WLZ_DRAWDOM_CMD_CIRCLE:
	if((errNum = WlzDrawDomStrGetAct(wSp, &act)) == WLZ_ERR_NONE)
	{
	  errNum = WlzDrawDomStrGetVal1D(wSp, &radius);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzDrawDomStrGetVal2D(wSp, &(pos0.vtX), &(pos0.vtY));
	  if(errNum == WLZ_ERR_EOO)
	  {
	    errNum = WLZ_ERR_NONE;
	  }
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzDrawDomAddCircle(&dwnObj, act, radius, pos0);
	}
	break;
      case WLZ_DRAWDOM_CMD_END:
	break;
      case WLZ_DRAWDOM_CMD_ERROR:
	errNum = WLZ_ERR_PARAM_DATA;
	break;
      case WLZ_DRAWDOM_CMD_LINE:
	if((errNum = WlzDrawDomStrGetAct(wSp, &act)) == WLZ_ERR_NONE)
	{
	  errNum = WlzDrawDomStrGetVal1D(wSp, &width);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzDrawDomStrGetVal2D(wSp, &(pos0.vtX), &(pos0.vtY));
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzDrawDomStrGetVal2D(wSp, &(pos1.vtX), &(pos1.vtY));
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzDrawDomAddLine(&dwnObj, act, width, pos0, pos1);
	}
	break;
      case WLZ_DRAWDOM_CMD_INIT:
	break;
      case WLZ_DRAWDOM_CMD_PEN:
	if((errNum = WlzDrawDomStrGetAct(wSp, &act)) == WLZ_ERR_NONE)
	{
	  errNum = WlzDrawDomStrGetVal1D(wSp, &width);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  errNum = WlzDrawDomStrGetVal2D(wSp, &(pos0.vtX), &(pos0.vtY));
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  radius = 0.5 * width;
	  errNum = WlzDrawDomAddCircle(&dwnObj, act, radius, pos0);
	}
	if(errNum == WLZ_ERR_NONE)
	{
	  while((errNum == WLZ_ERR_NONE) &&
		((errNum = WlzDrawDomStrGetVal2D(wSp,
				 &(pos1.vtX), &(pos1.vtY))) == WLZ_ERR_NONE))
	  {
	    errNum = WlzDrawDomAddLine(&dwnObj, act, width, pos0, pos1);
	    if(errNum == WLZ_ERR_NONE)
	    {
	      errNum = WlzDrawDomAddCircle(&dwnObj, act, radius, pos1);
	    }
	    pos0 = pos1;
	  }
	  if(errNum == WLZ_ERR_EOO)
	  {
	    errNum = WLZ_ERR_NONE;
	  }
	}
	break;
      case WLZ_DRAWDOM_CMD_UNKNOWN:
	break;
    }
  }
  while((errNum == WLZ_ERR_NONE) && (cmd != WLZ_DRAWDOM_CMD_END));
  if(errNum == WLZ_ERR_NONE)
  {
    if(dwnObj == NULL)
    {
      dwnObj = WlzMakeEmpty(&errNum);
    }
  }
  if(dstErr)
  {
    *dstErr = errNum;
  }
  return(dwnObj);
}

/*!
* \return	The parsed drawing command.
* \ingroup	WlzDomainOps
* \brief	Copies the next drawing command string from the string
* 		of drawing commands. The copied command is partial parsed
* 		using strtok_r (or strtok if strtok_r is not available)
* 		for the command name.
* \param	wSp			Drawing command workspace.
*/
static WlzDrawDomCmd	WlzDrawDomStrGetCmd(WlzDrawDomWSp *wSp)
{
  int		id0,
  		errFlg = 0;
  char		*str0,
  		*str1;
  WlzDrawDomCmd	cmd = WLZ_DRAWDOM_CMD_ERROR;
  const int	bufInc = 1024;

  /* Copy command to buffer. */
  id0 = 0;
  while((*(wSp->cmdStr + wSp->cmdStrIdx) != '\0') &&
        (*(wSp->cmdStr + wSp->cmdStrIdx) != ';'))
  {
    if(id0 >= wSp->bufMax)
    {
      wSp->bufMax += bufInc;
      if((wSp->buf = (char *)AlcRealloc(wSp->buf,
                                        wSp->bufMax * sizeof(char))) == NULL)
      {
	errFlg = 1;
        break;
      }
    }
    *(wSp->buf + id0++) = *(wSp->cmdStr + (wSp->cmdStrIdx)++);
  }
  /* Make sure buffer copy of command is NULL termainated. */ 
  if((errFlg == 0) && (id0 >= wSp->bufMax))
  {
    wSp->bufMax += bufInc;
    if((wSp->buf = (char *)AlcRealloc(wSp->buf,
                                     wSp->bufMax * sizeof(char))) == NULL)
    {
      errFlg = 1;
    }
  }
  if(errFlg != 0)
  {
    cmd = WLZ_DRAWDOM_CMD_ERROR;
  }
  else
  {
    ++(wSp->cmdStrIdx);
    *(wSp->buf + id0) = '\0';
    /* Strip out white space from the buffered command and make all
     * letters upper case. */
    str0 = str1 = wSp->buf;
    while((*str1 = toupper(*str0++)) != '\0')
    {
      if(isspace(*str1) == 0)
      {
        ++str1;
      }
    }
    /* Tokenize the buffer command string and get command name. */
#if defined _SVID_SOURCE || _BSD_SOURCE || _POSIX_C_SOURCE || _XOPEN_SOURCE
    str0 = strtok_r(wSp->buf, ":", &(wSp->savPtr));
#else
    str0 = strtok(wSp->buf, ":");
#endif
    if(str0 != NULL)
    {
      if(WlzStringMatchValue(&id0, str0,
			     "CIRCLE", WLZ_DRAWDOM_CMD_CIRCLE,
			     "END", WLZ_DRAWDOM_CMD_END,
			     "LINE", WLZ_DRAWDOM_CMD_LINE,
			     "WLZ_DRAW_DOMAIN", WLZ_DRAWDOM_CMD_INIT,
			     "PEN", WLZ_DRAWDOM_CMD_PEN,
			     NULL) != 0)
      {
	cmd = (WlzDrawDomCmd )id0;
      }
      else
      {
	cmd = WLZ_DRAWDOM_CMD_UNKNOWN;
      }
    }
  }
  return(cmd);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Continues parsing the copied drawing command using
* 		strtok_r (or strtok if strtok_r is not available)
* 		for the command action.
* \param	wSp			Drawing command workspace.
*/
static WlzErrorNum  WlzDrawDomStrGetAct(WlzDrawDomWSp *wSp,
				        WlzDrawDomAct *dstAct)
{
  int		id0;
  char		*str0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  /* Tokenize the buffer command string and get action. */
#if defined _SVID_SOURCE || _BSD_SOURCE || _POSIX_C_SOURCE || _XOPEN_SOURCE
  str0 = strtok_r(NULL, ",", &(wSp->savPtr));
#else
  str0 = strtok(NULL, ",");
#endif
  if(str0 == NULL)
  {
    errNum = WLZ_ERR_PARAM_DATA;
  }
  else
  {
    if(WlzStringMatchValue(&id0, str0,
			   "NONE", WLZ_DRAWDOM_ACT_NONE,
			   "DRAW", WLZ_DRAWDOM_ACT_DRAW,
			   "ERASE", WLZ_DRAWDOM_ACT_ERASE,
			   NULL) != 0)
    {
      *dstAct = (WlzDrawDomAct )id0;
    }
    else
    {
      errNum = WLZ_ERR_PARAM_DATA;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Continues parsing the copied drawing command using
* 		strtok_r (or strtok if strtok_r is not available)
* 		for a single double precision floating point value.
* \param	wSp			Drawing command workspace.
* \param	dstVal			Destination pointer for the
* 					value.
*/
static WlzErrorNum WlzDrawDomStrGetVal1D(WlzDrawDomWSp *wSp, double *dstVal)
{
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzDrawDomStrGetValND(wSp, 1, dstVal);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Continues parsing the copied drawing command using
* 		strtok_r (or strtok if strtok_r is not available)
* 		for a pair of double precision floating point values.
* \param	wSp			Drawing command workspace.
* \param	dstVal0			Destination pointer for the first
* 					value.
* \param	dstVal1			Destination pointer for the second
* 					value.
*/
static WlzErrorNum WlzDrawDomStrGetVal2D(WlzDrawDomWSp *wSp,
                                        double *dstVal0, double *dstVal1)
{
  double	val[2];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  errNum = WlzDrawDomStrGetValND(wSp, 2, val);
  if(errNum == WLZ_ERR_NONE)
  {
    *dstVal0 = val[0];
    *dstVal1 = val[1];
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Continues parsing the copied drawing command using
* 		strtok_r (or strtok if strtok_r is not available)
* 		for the specified number of double precision floating
* 		point values.
* \param	wSp			Drawing command workspace.
* \param	nVal			Number of values.
* \param	val			Destination pointer for the
* 					contiguious values.
*/
static WlzErrorNum WlzDrawDomStrGetValND(WlzDrawDomWSp *wSp,
					 int nVal, double *val)
{
  char		*str0;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  while(nVal-- > 0)
  {
    /* Tokenize the buffer command string and get action. */
#if defined _SVID_SOURCE || _BSD_SOURCE || _POSIX_C_SOURCE || _XOPEN_SOURCE
    str0 = strtok_r(NULL, ",", &(wSp->savPtr));
#else
    str0 = strtok(NULL, ",");
#endif
    if(str0 == NULL)
    {
      errNum = WLZ_ERR_EOO;
      break;
    }
    else if(sscanf(str0, "%lg", val++) != 1)
    {
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Adds or removes a circular region with the given centre
* 		and radius.
* \param	curObj			Pointer to current object pointer
* 					which will be replaced.
* \param	act			Drawing action.
* \param	radius			Circle radius.
* \param	centre			Circle centre.
*/
static WlzErrorNum WlzDrawDomAddCircle(WlzObject **curObj,
				       WlzDrawDomAct act,
				       double radius, WlzDVertex2 centre)
{
  WlzObject	*cObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  cObj = WlzMakeCircleObject(radius, centre.vtX, centre.vtY, &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzDrawDomAddObj(curObj, act, cObj);
  }
  (void )WlzFreeObj(cObj);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Adds or removes a rectangular region with the given
* 		width and centres for the opposing edges (of given width).
* \param	curObj			Pointer to current object pointer
* 					which will be replaced.
* \param	act			Drawing action.
* \param	width			rectangle width.
* \param	pos0			First edge centre.
* \param	pos1			Second edge centre.
*/
static WlzErrorNum WlzDrawDomAddLine(WlzObject **curObj, 
				       WlzDrawDomAct act,
				       double width,
				       WlzDVertex2 pos0, WlzDVertex2 pos1)
{
  WlzObject	*qObj = NULL;
  WlzDVertex2	vtx[4];
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  (void )WlzGeomRectFromWideLine(pos0, pos1, width, vtx);
  qObj = WlzMakeQuadrilateral(vtx[0].vtX, vtx[0].vtY, vtx[1].vtX, vtx[1].vtY,
                            vtx[2].vtX, vtx[2].vtY, vtx[3].vtX, vtx[3].vtY,
			    &errNum);
  if(errNum == WLZ_ERR_NONE)
  {
    errNum = WlzDrawDomAddObj(curObj, act, qObj);
  }
  (void )WlzFreeObj(qObj);
  return(errNum);
}

/*!
* \return	Woolz error code.
* \ingroup	WlzDomainOps
* \brief	Adds or removes a given domain from the current drawn
* 		domain.
* \param	curObj			Pointer to current object pointer
* 					which will be replaced, must not
* 					be NULL although *curObj may be.
* \param	act			Drawing action.
* \param	obj			Object to be added to or removed
* 					from the current domain object.
*/
static WlzErrorNum WlzDrawDomAddObj(WlzObject **curObj, WlzDrawDomAct act,
				    WlzObject *obj)
{
  WlzObject	*nObj = NULL;
  WlzErrorNum	errNum = WLZ_ERR_NONE;

  switch(act)
  {
    case WLZ_DRAWDOM_ACT_DRAW:
      if(*curObj == NULL)
      {
        nObj = WlzCopyObject(obj, &errNum); 
      }
      else
      {
        nObj = WlzUnion2(*curObj, obj, &errNum);
      }
      break;
    case WLZ_DRAWDOM_ACT_ERASE:
      if(*curObj != NULL)
      {
        nObj = WlzDiffDomain(*curObj, obj, &errNum);
      }
      break;
    default:
      break;
  }
  if(nObj)
  {
    (void )WlzFreeObj(*curObj);
    *curObj = nObj;
  }
  return(errNum);
}
