#pragma ident "MRC HGU $Id$"
/*!
* \file         WlzMakeProperties.c
* \author       Richard Baldock <Richard.Baldock@hgu.mrc.ac.uk>
* \date         March 1999
* \version      MRC HGU $Id$
*               $Revision$
*               $Name$
* \par Copyright:
*               1994-2002 Medical Research Council, UK.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \ingroup      WlzProperty
* \brief        Allocation and freeing routines for woolz properties.
*               
* \todo         -
* \bug          None known
*
* Maintenance log with most recent changes at top of list.
*/

#include <pwd.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>

#include <Wlz.h>

/* function:     WlzMakeSimpleProperty    */
/*! 
* \ingroup      WlzProperty
* \brief        Allocate space for a WlzSimpleProperty which is a
structure with size bytes allocated for data.
*
* \return       WlzSimpleProperty *
* \param    size	Size of space allocated for the property.
* \param    dstErr	error return values: WLZ_ERR_NONE, WLZ_ERR_MEM_ALLOC.
* \par      Source:
*                WlzMakeProperties.c
*/
WlzSimpleProperty *
WlzMakeSimpleProperty(int size, WlzErrorNum *dstErr)
{
  WlzSimpleProperty	*p=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;

  if( (p = (WlzSimpleProperty *)
       AlcCalloc(1,sizeof(WlzSimpleProperty))) == NULL ){
    errNum = WLZ_ERR_MEM_ALLOC;
  }
  else {
    p->type = WLZ_PROPERTY_SIMPLE;
    p->linkcount = 0;
    p->size = size;

    if( (p->prop = AlcMalloc(size)) == NULL ){
      AlcFree( p );
      errNum = WLZ_ERR_MEM_ALLOC;
    }
    else {
      p->freeptr = AlcFreeStackPush(p->freeptr, p->prop, NULL);
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return( p );
}


/* function:     WlzFreeSimpleProperty    */
/*! 
* \ingroup      WlzProperty
* \brief        Free space allocated for a WlzSimpleProperty.
*
* \return       WkzErrorNum
* \param    prop	property to be freed
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzFreeSimpleProperty(WlzSimpleProperty *prop)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (prop == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(prop->linkcount), &errNum) ){
    if( prop->freeptr != NULL ){
      AlcFreeStackFree( prop->freeptr );
    }
    AlcFree( (void *) prop );
  }

  return errNum;
}

/* function:     WlzMakeEMAPProperty    */
/*! 
* \ingroup      WlzProperty
* \brief        Make an EMAP property structure.
*
* \return       WlzEMAPProperty *
* \param    type	required EMAP property type:
 WLZ_EMAP_PROPERTY_GREY_MODEL, WLZ_EMAP_PROPERTY_ANATOMY_DOMAIN or
 WLZ_EMAP_PROPERTY_OTHER_DOMAIN.
* \param    theilerStage	Theiler stage of the reconstruction
 or of the reconstruction for which this is a domain
* \param    modelName	Standard EMAP reconstruction name
* \param    version	version of the model
* \param    fileName	original filename when written using MAPaint
* \param    comment	text comment
* \param    dstErr	error return values: WLZ_ERR_NONE,
 WLZ_ERR_PARAM_DATA, WLZ_ERR_MEM_ALLOC
* \par      Source:
*                WlzMakeProperties.c
*/
WlzEMAPProperty *WlzMakeEMAPProperty(
  WlzEMAPPropertyType	type,
  int			theilerStage,
  char			*modelName,
  char			*version,
  char			*fileName,
  char			*comment,
  WlzErrorNum	*dstErr)
{
  WlzEMAPProperty	*rtnProp=NULL;
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  struct passwd		*pwdStruct;

  /* check the calling parameters */
  switch( type ){
  case WLZ_EMAP_PROPERTY_GREY_MODEL:
  case WLZ_EMAP_PROPERTY_ANATOMY_DOMAIN:
  case WLZ_EMAP_PROPERTY_OTHER_DOMAIN:
    break;

  default:
    errNum = WLZ_ERR_PARAM_DATA;
    break;
  }
  if((theilerStage > 28) || (theilerStage < 1)){
    errNum = WLZ_ERR_PARAM_DATA;
  }

  /* allocate space and copy in properties */
  if( errNum == WLZ_ERR_NONE ){
    if( rtnProp = (WlzEMAPProperty *) 
       AlcCalloc(1, sizeof(WlzEMAPProperty)) ){
      rtnProp->type = WLZ_PROPERTY_EMAP;
      rtnProp->emapType = type;
      rtnProp->creationTime = time(NULL);
      rtnProp->modificationTime = rtnProp->creationTime;
      pwdStruct = getpwuid(getuid());
      if((strlen(pwdStruct->pw_name) + strlen(pwdStruct->pw_gecos)
	  + 4) < EMAP_PROPERTY_AUTHORNAME_LENGTH){
	sprintf(rtnProp->creationAuthor, "%s (%s)",
		pwdStruct->pw_name,
		pwdStruct->pw_gecos?pwdStruct->pw_gecos:"");
      }
      else {
	strncpy(rtnProp->creationAuthor, pwdStruct->pw_name,
		EMAP_PROPERTY_AUTHORNAME_LENGTH - 1);
      }
      (void) gethostname(&(rtnProp->creationMachineName[0]),
			 EMAP_PROPERTY_AUTHORNAME_LENGTH);
      rtnProp->creationMachineName[EMAP_PROPERTY_AUTHORNAME_LENGTH-1]
	= '\0';
      strcpy(rtnProp->modificationAuthor, rtnProp->creationAuthor);

      rtnProp->theilerStage = theilerStage;
      if( modelName ){
	strncpy(rtnProp->modelName, modelName,
		EMAP_PROPERTY_MODELNAME_LENGTH - 1);
      }
      if( version ){
	strncpy(rtnProp->version, version,
		EMAP_PROPERTY_VERSION_LENGTH - 1);
      }
      if( fileName ){
	rtnProp->fileName = AlcStrDup(fileName);
	rtnProp->freeptr = AlcFreeStackPush(rtnProp->freeptr,
					    rtnProp->fileName,
					    NULL);
      }
      if( comment ){
	rtnProp->comment = AlcStrDup(comment);
	rtnProp->freeptr = AlcFreeStackPush(rtnProp->freeptr,
					    rtnProp->comment,
					    NULL);
      }
    }
    else {
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnProp;
}

/* function:     WlzChangeEMAPProperty    */
/*! 
* \ingroup      WlzProperty
* \brief        Change the values of an EMAP property. Each given
value will be checked against the current value to see if a change is
necessary. If the property is changed the modification time and author
will be updated. NULL for any of modelName, version, fileName, comment
implies no change. To change these please modify directly, updating the
modification time appropriately.
*
* \return       woolz error
* \param    prop	property to be changed
* \param    type	new type
* \param    theilerStage	new Theiler stage
* \param    modelName	new model name
* \param    version	new model version
* \param    fileName	new filename
* \param    comment	new comment
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzChangeEMAPProperty(
  WlzEMAPProperty	*prop,
  WlzEMAPPropertyType	type,
  int			theilerStage,
  char			*modelName,
  char			*version,
  char			*fileName,
  char			*comment)
{
  WlzErrorNum		errNum=WLZ_ERR_NONE;
  struct passwd		*pwdStruct;
  int			modifiedFlg=0;

  /* check the calling parameters */
  if( prop == NULL ){
    errNum = WLZ_ERR_PROPERTY_NULL;
  }
  else if( prop->type != WLZ_PROPERTY_EMAP ){
    errNum = WLZ_ERR_PROPERTY_TYPE;
  }
  else {
    switch( type ){
    case WLZ_EMAP_PROPERTY_GREY_MODEL:
    case WLZ_EMAP_PROPERTY_ANATOMY_DOMAIN:
    case WLZ_EMAP_PROPERTY_OTHER_DOMAIN:
      if((theilerStage > 28) || (theilerStage < 1)){
	errNum = WLZ_ERR_PARAM_DATA;
      }
      break;

    default:
      errNum = WLZ_ERR_PARAM_DATA;
      break;
    }
  }

  /* copy in properties */
  if( errNum == WLZ_ERR_NONE ){

    if( prop->emapType != type ){
      prop->emapType = type;
      modifiedFlg = 1;
    }

    if( prop->theilerStage != theilerStage ){
      prop->theilerStage = theilerStage;
      modifiedFlg = 1;
    }

    if( modelName ){
      if( prop->modelName ){
	if( strcmp(modelName, prop->modelName) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->modelName, modelName,
	      EMAP_PROPERTY_MODELNAME_LENGTH - 1);
    }

    if( version ){
      if( prop->version ){
	if( strcmp(version, prop->version) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      strncpy(prop->version, version,
	      EMAP_PROPERTY_VERSION_LENGTH - 1);
    }

    if( fileName ){
      if( prop->fileName ){
	if( strcmp(fileName, fileName) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      
      prop->fileName = AlcStrDup(fileName);
      prop->freeptr = AlcFreeStackPush(prop->freeptr,
				       prop->fileName,
				       NULL);
    }
    if( comment ){
      if( prop->comment ){
	if( strcmp(comment, comment) ){
	  modifiedFlg = 1;
	}
      }
      else {
	modifiedFlg = 1;
      }
      
      prop->comment = AlcStrDup(comment);
      prop->freeptr = AlcFreeStackPush(prop->freeptr,
				       prop->comment,
				       NULL);
    }
      
    if((errNum == WLZ_ERR_NONE) && modifiedFlg){
      prop->modificationTime = time(NULL);
      pwdStruct = getpwuid(getuid());
      sprintf(prop->modificationAuthor, "%s (%s)",
	      pwdStruct->pw_name,
	      pwdStruct->pw_gecos?pwdStruct->pw_gecos:"");
    }
  }

  return errNum;
}

/* function:     WlzFreeEMAPProperty    */
/*! 
* \ingroup      WlzProperty
* \brief        Free an EMAP property.
*
* \return       woolz error number values: WLZ_ERR_NONE
* \param    prop	property to be freed
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzFreeEMAPProperty(WlzEMAPProperty *prop)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;

  /* check the object pointer and linkcount */
  if (prop == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(prop->linkcount), &errNum) ){
    if( prop->freeptr != NULL ){
      AlcFreeStackFree( prop->freeptr );
    }
    AlcFree( (void *) prop );
  }

  return errNum;
}

/* function:     WlzGetProperty    */
/*! 
* \ingroup      WlzProperty
* \brief        Get a property of a given type from a property list. It
is assumed that there is only one property of a given type in the list.
*
* \return       The required property. If the core value is NULL the
 requested property is not in the list.
* \param    plist	Property list to search
* \param    type	Property type required
* \param    dstErr	error return values: WLZ_ERR_NONE,
 WLZ_ERR_PROPERTY_NULL.
* \par      Source:
*                WlzMakeProperties.c
*/
WlzProperty WlzGetProperty(
  AlcDLPList		*plist,
  WlzObjectType		type,
  WlzErrorNum 		*dstErr)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  WlzProperty	rtnProp;
  AlcDLPItem	*item;

  rtnProp.core = NULL;

  if( plist ){
    item = plist->head;
    do {
      if( item ){
	if(item->entry &&
	   (((WlzCoreProperty *) (item->entry))->type == type) ){
	  rtnProp.core = (WlzCoreProperty *) (item->entry);
	  break;
	}
	item = item->next;
      }
    } while( item != plist->head );
  }
  else {
    errNum = WLZ_ERR_PROPERTY_NULL;
  }

  if( dstErr ){
    *dstErr = errNum;
  }
  return rtnProp;
}

/* function:     WlzRemoveProperty    */
/*! 
* \ingroup      WlzProperty
* \brief        Remove, without freeing, a property from a property list.
This is typiucally used after a call to WlzGetProperty().
*
* \return       woolz error values: WLZ_ERR_NONE
* \param    plist	Property list to search
* \param    prop	Property to be removed.
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzRemoveProperty(
  AlcDLPList		*plist,
  WlzProperty		prop)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  AlcDLPItem	*item;

  if( plist == NULL ){
    errNum = WLZ_ERR_PROPERTY_NULL;
  }
  else {
    item = plist->head;
    do {
      if( item ){
	if(item->entry &&
	   ((WlzCoreProperty *) (item->entry) == prop.core) ){
	  (void) AlcDLPItemUnlink(plist, item, 0, NULL);
	  AlcFree(item);
	  break;
	}
	item = item->next;
      }
    } while( item != plist->head );
  }

  return errNum;
}

/* function:     WlzFreeProperty    */
/*! 
* \ingroup      WlzProperty
* \brief        Free a woolz property
*
* \return       wooz error values: WLZ_ERR_NONE, WLZ_ERR_MEM_ALLOC.
* \param    prop	property to be freed
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzFreeProperty(WlzProperty prop)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  AlcErrno	alcErrNum=ALC_ER_NONE;

  if(prop.core == NULL){
    return( WLZ_ERR_NONE );
  }

  if( WlzUnlink(&(prop.core->linkcount), &errNum) ){
    if( prop.core->freeptr != NULL ){
      alcErrNum = AlcFreeStackFree(prop.core->freeptr);
    }
    AlcFree( (void *) prop.core );
    if( alcErrNum != ALC_ER_NONE ){
      errNum = WLZ_ERR_MEM_ALLOC;
    }
  }

  return errNum;
}


/* function:     WlzFreePropertyList    */
/*! 
* \ingroup      WlzProperty
* \brief        Free a complete property list (including the list
structure itself).
*
* \return       woolz error values: WLZ_ERR_NONE, WLZ_ERR_PARAM_NULL.
* \param    plist	Property list to be freed
* \par      Source:
*                WlzMakeProperties.c
*/
WlzErrorNum WlzFreePropertyList(
  AlcDLPList	*plist)
{
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  AlcDLPItem	item;

  if( plist ){
    AlcDLPListFree(plist);
  }

  return errNum;
}


/* function:     WlzFreePropertyListEntry    */
/*! 
* \ingroup      WlzProperty
* \brief        Free a property list entry. This is a procedure that
can be passed e.g. to AlcDLPListEntryAppend to enable automatic freeing of
properties held on a property list.
*
* \param    prop	property to be freed
* \par      Source:
*                WlzMakeProperties.c
*/
void WlzFreePropertyListEntry(void *prop)
{
  WlzCoreProperty *core = (WlzCoreProperty *) prop;
  WlzErrorNum	errNum=WLZ_ERR_NONE;
  AlcErrno	alcErrNum=ALC_ER_NONE;

  /* check the object pointer and linkcount */
  if (prop == NULL){
    return;
  }

  if( WlzUnlink(&(core->linkcount), &errNum) ){
    if( core->freeptr != NULL ){
      alcErrNum = AlcFreeStackFree(core->freeptr);
    }
    AlcFree( (void *) prop );
  }

  return;
}

