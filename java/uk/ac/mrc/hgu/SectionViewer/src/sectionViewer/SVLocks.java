package sectionViewer;
import sectionViewer.*;

/**
 *   Provides various locks for multi-threading synchronization.
 *   This can't be an interface (=> static final) as multiple
 *   Controller progs may want their own locks, for example
 *   tie point and warping applications.
 */
public class SVLocks {

   public Object _svpLock1 = null;
   public Object _svpLock2 = null;
   public Object _svpLock3 = null;

   public SVLocks() {
      _svpLock1 = new Object();
      _svpLock2 = new Object();
      _svpLock3 = new Object();
   }

}

