#pragma ident "MRC HGU $Id$"
/************************************************************************
* Project:	Mouse Atlas
* Title:        ReconstructMisc.c				
* Date:         April 1999
* Author:       Bill Hill                                              
* Copyright:    1999 Medical Research Council, UK.
*		All rights reserved.				
* Address:	MRC Human Genetics Unit,			
*		Western General Hospital,			
*		Edinburgh, EH4 2XU, UK.				
* Purpose:      Provides miscellaneous small functions for the the
*		MRC Human Genetics Unit reconstruction library.	
* $Revision$
* Maintenance:  Log changes below, with most recent at top of list.    
************************************************************************/
#include <Reconstruct.h>

int		RecMinOf4I(int val1, int val2, int val3, int val4),
		RecMaxOf4I(int val1, int val2, int val3, int val4);
RecError	RecErrorFromWlz(WlzErrorNum wlzErr);

/************************************************************************
* Function:	RecPowerOfTwoC2I				
* Returns:	WlzIVertex2:		Calculated powers of two.
* Purpose:	Uses RecPowerOfTwo() for x,y coordinate pairs.	
*		Calculates the pair of integers that are an integer
*		power of two and greater than or equal to the given
*		pair.						
* Global refs:	-						
* Parameters:	WlzIVertex2 *ipLarger:	Destination for pair that
*					has each component an integer
*					power of two that is greater
*					than or equal to the given
*					component.		
*		WlzIVertex2 given:	Given coordinate pair.	
************************************************************************/
WlzIVertex2	RecPowerOfTwoC2I(WlzIVertex2 *ipLarger, WlzIVertex2 given)
{
  WlzIVertex2	powOfTwo;

  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecPowerOfTwoC2I FE 0x%lx {%d %d}\n",
	   (unsigned long )ipLarger, given.vtX, given.vtY));
  if(ipLarger)
  {
    powOfTwo.vtX = RecPowerOfTwo(&(ipLarger->vtX), given.vtX);
    powOfTwo.vtY = RecPowerOfTwo(&(ipLarger->vtY), given.vtY);
  }
  else
  {
    powOfTwo.vtX = RecPowerOfTwo(NULL, given.vtX);
    powOfTwo.vtY = RecPowerOfTwo(NULL, given.vtY);
  }
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecPowerOfTwoC2I FX {%d %d}\n",
	   powOfTwo.vtX, powOfTwo.vtY));
  return(powOfTwo);
}
/************************************************************************
* Function:	RecPowerOfTwo					
* Returns:	int:			Calculated power of two.
* Purpose:	Calculates an integer that is an integer power of two
*		and greater than or equal to the given integer.	
* Global refs:	-						
* Parameters:	int *ipLarger:		Destination for integer that is
*					an integer power of two and is
*					for greater than or equal to
*					the given integer.	
*		unsigned int given:	Given integer.		
************************************************************************/
int		RecPowerOfTwo(int *ipLarger, unsigned int given)
{
  int		n = 1,
		powOfTwo = 0;

  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecPowerOfTwo FE 0x%lx %d\n",
	   (unsigned long )ipLarger, given));
  while(n < given)
  {
    n <<= 1;
    ++powOfTwo;
  }
  if(ipLarger)
  {
    *ipLarger = n;
  }
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecPowerOfTwo FX %d\n",
	   powOfTwo));
  return(powOfTwo);
}

/************************************************************************
* Function:	RecMinOf4I					
* Returns:	int:			Minimum value.		
* Purpose:	Given four values, finds and returns the minimum.
* Global refs:	-						
* Parameters:	int val1:		First value.		
*		int val2:		Second value.		
*		int val3:		Third value.		
*		int val4:		Forth value.		
************************************************************************/
int		RecMinOf4I(int val1, int val2, int val3, int val4)
{
  int		minValue;

  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecMinOf4I FE %d %d %d %d\n",
	   val1, val2, val3, val4));
  if((val1 < val2) && (val1 < val3) && (val1 < val4))
  {
    minValue = val1;
  }
  else if((val2 < val3) && (val2 < val4))
  {
    minValue = val2;
  }
  else if(val3 < val4)
  {
    minValue = val3;
  }
  else
  {
    minValue = val4;
  }
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecMinOf4I FX %d\n",
	   minValue));
  return(minValue);
}

/************************************************************************
* Function:	RecMaxOf4I					
* Returns:	int:			Maximum value.		
* Purpose:	Given four values, finds and returns the maximum.
* Global refs:	-						
* Parameters:	int val1:		First value.		
*		int val2:		Second value.		
*		int val3:		Third value.		
*		int val4:		Forth value.		
************************************************************************/
int		RecMaxOf4I(int val1, int val2, int val3, int val4)
{
  int		maxValue;

  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecMaxOf4I FE %d %d %d %d\n",
	   val1, val2, val3, val4));
  if((val1 > val2) && (val1 > val3) && (val1 > val4))
  {
    maxValue = val1;
  }
  else if((val2 > val3) && (val2 > val4))
  {
    maxValue = val2;
  }
  else if(val3 > val4)
  {
    maxValue = val3;
  }
  else
  {
    maxValue = val4;
  }
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecMaxOf4I FX %d\n",
	  maxValue));
  return(maxValue);
}

/************************************************************************
* Function:	RecUCharToShort					
* Returns:	void						
* Purpose:	Copies the given unsigned chars to shorts using the
*		given buffers.					
* Global refs:	-						
* Parameters:	short *sBuf:		Short buffer.		
*		unsigned char *ucBuf:	Unsigned char buffer.	
*		int count:		Number of elements in buffer.
************************************************************************/
void		RecUCharToShort(short *sBuf, unsigned char *ucBuf, int count)
{
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecUCharToShort FE 0x%lx 0x%lx %d\n",
	   (unsigned long )sBuf, (unsigned long )ucBuf, count));
  while(count-- > 0)
  {
    *sBuf++ = (short )(*ucBuf++);
  }
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecUCharToShort FX\n"));
}

/************************************************************************
* Function:	RecUCharToInt					
* Returns:	void						
* Purpose:	Copies the given unsigned chars to ints using the
*		given buffers.					
* Global refs:	-						
* Parameters:	int *iBuf:		Int buffer.		
*		unsigned char *ucBuf:	Unsigned char buffer.	
*		int count:		Number of elements in buffer.
************************************************************************/
void		RecUCharToInt(int *iBuf, unsigned char *ucBuf, int count)
{
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecUCharToInt FE 0x%lx 0x%lx %d\n",
	   (unsigned long )iBuf, (unsigned long )ucBuf, count));
  while(count-- > 0)
  {
    *iBuf++ = (int )(*ucBuf++);
  }
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecUCharToInt FX\n"));
}

/************************************************************************
* Function:	RecErrorToStr					
* Returns:	const char *:		Error string.		
* Purpose:	Returns a pointer to a const error string appropriate
*		to the given error.				
* Global refs:	-						
* Parameters:	RecError errFlag:	Given error.		
************************************************************************/
const char	*RecErrorToStr(RecError errFlag)
{
  const char	*errStr;
  const char	*errNoneStr =	"None",
  		*errUsageStr =	"Command line usage",
		*errArgsStr =	"Command line arguments",
		*errReadStr =	"Read failure",
		*errWriteStr =	"Write failure",
		*errMallocStr =	"Memory allocation failure",
		*errUnimplStr =	"Unimplemented feature",
		*errWlzStr =	"Woolz image processing error",
		*errSyntaxStr =	"File syntax error",
		*errFuncStr =	"Function parameters invalid",
		*errCancelStr =	"Operation canceled",
		*errListStr =	"Section list invalid",
		*errStrDef =	"Unknown error";

  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecErrorToStr FE %d\n",
	   (int )errFlag));
  errStr = errStrDef;
  switch(errFlag)
  {
    case REC_ERR_NONE:
      errStr = errNoneStr;
      break;
    case REC_ERR_USAGE:
      errStr = errUsageStr;
      break;
    case REC_ERR_ARGS:
      errStr = errArgsStr;
      break;
    case REC_ERR_READ:
      errStr = errReadStr;
      break;
    case REC_ERR_WRITE:
      errStr = errWriteStr;
      break;
    case REC_ERR_MALLOC:
      errStr = errMallocStr;
      break;
    case REC_ERR_UNIMPL:
      errStr = errUnimplStr;
      break;
    case REC_ERR_WLZ:
      errStr = errWlzStr;
      break;
    case REC_ERR_SYNTAX:
      errStr = errSyntaxStr;
      break;
    case REC_ERR_FUNC:
      errStr = errFuncStr;
      break;
    case REC_ERR_CANCEL:
      errStr = errCancelStr;
      break;
    case REC_ERR_LIST:
      errStr = errListStr;
      break;
  }
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecErrorToStr FX %s\n", errStr));
    return(errStr);
}

/************************************************************************
* Function:	RecMethodToStr					
* Returns:	const char *:		Method string.		
* Purpose:	Returns a pointer to a const error string appropriate
*		to the given error.				
* Global refs:	-						
* Parameters:	RecMethod method:	Given method.		
************************************************************************/
const char	*RecMethodToStr(RecMethod method)
{
  const char	*mthdStr;
  const char	*mthdNoneStr =	"Not registering sections",
  		*mthdPrinc =	"Principle axes (with centre of mass)",
		*mthdTrans =	"Translation match",
		*mthdRotate =	"Rotation match",
		*mthdIdentity =	"Identity",
		*mthdStrDef =	"Unknown method";

  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_1),
	  ("RecMethodToStr FE %d\n",
	   (int )method));
  mthdStr = mthdStrDef;
  switch(method)
  {
    case REC_MTHD_NONE:
      mthdStr = mthdNoneStr;
      break;
    case REC_MTHD_PRINC:
      mthdStr = mthdPrinc;
      break;
    case REC_MTHD_TRANS:
      mthdStr = mthdTrans;
      break;
    case REC_MTHD_ROTATE:
      mthdStr = mthdRotate;
      break;
    case REC_MTHD_IDENTITY:
      mthdStr = mthdIdentity;
      break;
  }
  REC_DBG((REC_DBG_MISC|REC_DBG_LVL_FN|REC_DBG_LVL_2),
	  ("RecMethodToStr FX %s\n", mthdStr));
    return(mthdStr);
}

/************************************************************************
* Function:	RecErrorFromWlz					
* Returns:	RecError:		Reconstruct error code.	
* Purpose:	Converts the given Woolz error code into the 	
*		appropriate reconstruct error code.		
* Global refs:	-						
* Parameters:	WlzErrorNum wlzErr:	Given woolz error code.	
************************************************************************/
RecError	RecErrorFromWlz(WlzErrorNum wlzErr)
{
  RecError	recErrNum  = REC_ERR_WLZ;

  switch(wlzErr)
  {
    case WLZ_ERR_NONE:
      recErrNum = REC_ERR_NONE;
      break;
    case WLZ_ERR_READ_EOF:
    case WLZ_ERR_READ_INCOMPLETE:
      recErrNum = REC_ERR_READ;
      break;
    case WLZ_ERR_WRITE_EOF:
    case WLZ_ERR_WRITE_INCOMPLETE:
      recErrNum = REC_ERR_WRITE;
      break;
    case WLZ_ERR_MEM_ALLOC:
    case WLZ_ERR_MEM_FREE:
      recErrNum = REC_ERR_MALLOC;
      break;
    default: 		              /* Call everything else a Woolz error! */
      recErrNum  = REC_ERR_WLZ;
      break;
  }
  return(recErrNum);
}
