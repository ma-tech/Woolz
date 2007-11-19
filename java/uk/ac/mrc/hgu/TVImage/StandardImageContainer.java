package uk.ac.mrc.hgu.TVImage;

import java.awt.*;
import java.awt.image.*;
import java.awt.event.*;
import java.io.*;
import java.net.*;
import java.util.*;
import javax.imageio.*;
import javax.swing.*;
import javax.swing.event.*;

import edu.stanford.genetics.treeview.*;
import edu.stanford.genetics.treeview.dendroview.*;

public abstract class StandardImageContainer extends ImageContainer {

   protected String emageId;
   protected String geneName;
   protected XmlAppSettings appSettings;

   public StandardImageContainer(String s, BufferedImage i,
	 String c, String emageId, String geneName, int sortIndex,
	 ImageViewManager tvim) {
      super(s, i, sortIndex);

      appSettings = tvim.getAppSettingsObject();

      //need to strip the start off the emageId
      if (emageId.length() > 6) {
	 this.emageId = emageId.substring(6);
      } else {
	 this.emageId = emageId;
      }

      this.geneName = geneName;

      labels.add(new ClickableLabel("EMAGE ID: " + this.emageId,
	       appSettings.getEmagePageLocation(this.emageId),
	       tvim));
      labels.add(new ClickableLabel("Gene: " + this.geneName,
	       appSettings.getGeneNameQueryString(this.geneName),
	       tvim));

      addLabels();
   }
}
