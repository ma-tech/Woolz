package sectionViewer;

import sectionViewer.*;

import java.awt.*;
import javax.swing.*;
import java.io.*;

/**
 *   Utility class for SectionViewer
 */

public class Utils {

   public static Font _font = null;

//-------------------------------------------
   public static void attachMenuFont(JComponent menu, Font font) {
      _font = font;
      attachFont(menu);
   }

//-------------------------------------------
   public static void attachFont(JComponent menu) {
      int numElements = 0;
      MenuElement elements[] = null;
      MenuElement el = null;

      /*
       * note that the elements returned by getSubElements()
       * may be JMenu, JPopupMenu or JMenuItem (but not JSeparator).
       */

      if(menu instanceof JMenuBar) {
	 ((JMenuBar)menu).setFont(_font);
	 /* get the menu elements */
	 elements = ((JMenuBar)menu).getSubElements();
	 if(elements == null) return;
	 numElements = elements.length;
	 for(int i=0; i<numElements; i++) {
	    el = elements[i];
	    if(el instanceof JPopupMenu) {
	       ((JPopupMenu)el).setFont(_font);
	       attachFont((JPopupMenu)el);
	    } else if(el instanceof JMenu) {
	       ((JMenu)el).setFont(_font);
	       attachFont((JMenu)el);
	    } else if(el instanceof SelectableMenu) {
	       ((SelectableMenu)el).setFont(_font);
	       attachFont((SelectableMenu)el);
	    } else if(el instanceof JMenuItem) {
	       ((JMenuItem)el).setFont(_font);
	    }
	 }
      } else if(menu instanceof JMenu) {
	 /* get the menu elements */
	 elements = ((JMenu)menu).getSubElements();
	 if(elements == null) return;
	 numElements = elements.length;
	 for(int i=0; i<numElements; i++) {
	    el = elements[i];
	    if(el instanceof JPopupMenu) {
	       ((JPopupMenu)el).setFont(_font);
	       attachFont((JPopupMenu)el);
	    } else if(el instanceof JMenu) {
	       ((JMenu)el).setFont(_font);
	       attachFont((JMenu)el);
	    } else if(el instanceof SelectableMenu) {
	       ((SelectableMenu)el).setFont(_font);
	       attachFont((SelectableMenu)el);
	    } else if(el instanceof JMenuItem) {
	       ((JMenuItem)el).setFont(_font);
	    }
	 }
      } else if(menu instanceof SelectableMenu) {
	 /* get the menu elements */
	 elements = ((SelectableMenu)menu).getSubElements();
	 if(elements == null) return;
	 numElements = elements.length;
	 for(int i=0; i<numElements; i++) {
	    el = elements[i];
	    if(el instanceof JPopupMenu) {
	       ((JPopupMenu)el).setFont(_font);
	       attachFont((JPopupMenu)el);
	    } else if(el instanceof JMenu) {
	       ((JMenu)el).setFont(_font);
	       attachFont((JMenu)el);
	    } else if(el instanceof SelectableMenu) {
	       ((SelectableMenu)el).setFont(_font);
	       attachFont((SelectableMenu)el);
	    } else if(el instanceof JMenuItem) {
	       ((JMenuItem)el).setFont(_font);
	    }
	 }
      } else if(menu instanceof JPopupMenu) {
	 /* get the menu elements */
	 elements = ((JPopupMenu)menu).getSubElements();
	 if(elements == null) return;
	 numElements = elements.length;
	 for(int i=0; i<numElements; i++) {
	    el = elements[i];
	    /*
	     * make sure JMenu & SelectableMenu go before JMenuItem
	     * as they are also a JMenuItem.
	     */
	    if(el instanceof JMenu) {
	       ((JMenu)el).setFont(_font);
	       attachFont((JMenu)el);
	    } else if(el instanceof SelectableMenu) {
	       ((SelectableMenu)el).setFont(_font);
	       attachFont((SelectableMenu)el);
	    } else if(el instanceof JMenuItem) {
	       ((JMenuItem)el).setFont(_font);
	    }
	 }
      }
   } // attachMenuFont()
//-------------------------------------------
} // class Utils
