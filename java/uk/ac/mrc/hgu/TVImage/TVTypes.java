package uk.ac.mrc.hgu.TVImage;

import java.awt.Color;

/**
 *   Provides various constants for the TreeView application
 */
public interface TVTypes {

   public static final Color BGCOLOR = Color.white;
   public static final boolean MULTITHREADED = false;
   public static final boolean PRELOAD = false;

   //-----------------------------------
   /*
    * These ensure that only one scroll pane of
    * a particular type can be included in a scroll group.
    */
   public static final int HEATMAP = 0;
   public static final int RAW = 1;
   public static final int ANNOTATED = 2;
   public static final int SUMMARY = 3;
   public static final int NODE = 10;
   public static final int NULL = -1;
   //-----------------------------------
   // Use Cases.
   public static final int NONE = -1;
   public static final int CORRELATION_CHANGED = 1;
   public static final int SINGLE_CLICK_NODE = 2;
   public static final int DOUBLE_CLICK_NODE = 3;
   public static final int DOUBLE_CLICK_GENE = 4;
   public static final int TREE_NODE_CLICKED = 5;
   public static final int MATRIX_CLICKED = 6;
   public static final int HEATMAP_TAB_SELECTED = 8;
   public static final int RAW_TAB_SELECTED = 9;
   public static final int ANNOTATED_TAB_SELECTED = 7;
   //-----------------------------------
   // Tab titles and description of view.
   public static final String HEATMAP_DESCRIPTION = "Heatmap";
   public static final String RAW_DESCRIPTION = "Raw";
   public static final String ANNOTATED_DESCRIPTION = "Annotated";
   public static final String SUMMARY_DESCRIPTION = "Summary";
   //-----------------------------------

}

