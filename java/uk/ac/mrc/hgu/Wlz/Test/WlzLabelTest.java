import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzLabelTest
{
  public static void main (String[] args)
  {
    boolean	showUsageFlg = false;
    int		optIdx = 0,
    		exitStatus = 0;
    WlzPixelV	threshV = new WlzPixelV(128);

    while((optIdx < args.length) && (args[optIdx].startsWith("-")) &&
          (args[optIdx].length() > 1))
    {
      if(args[optIdx].equals("-h"))
      {
        showUsageFlg = true;
      }
      else if(args[optIdx].equals("-t"))
      {
        if((optIdx + 1) < args.length)
	{
	  threshV = new WlzPixelV(Integer.parseInt(args[optIdx + 1]));
	}
	++optIdx;
      }
      else
      {
        showUsageFlg = true;
	exitStatus = 1;
      }
      ++optIdx;
    }
    if(showUsageFlg)
    {
      System.err.println("Usage: WlzBoundingBox [-t threshold] [-h] " +
      			 "<input object>");
      System.err.println(
      "Reads a Woolz object, thresholds it and then labels the thresholded");
      System.err.println(
      "image. The number of labeled components is out put along with the");
      System.err.println(
      "area of each component.");
      System.exit(exitStatus);
    }
    else
    {
      try
      {
	boolean firstPass = true;

	if(firstPass)
	{
	  int idx, area;
	  String inFile = args[optIdx];
	  WlzFileStream in = new WlzFileInputStream(args[optIdx]);
	  WlzObject inObj = WlzObject.WlzReadObj(in);
	  area = WlzObject.WlzArea(inObj);
	  System.out.println("input object area = " + area);
	  System.out.println("threshold = " + threshV.getIntValue());

	  WlzObject thrObj = WlzObject.WlzThreshold(inObj, threshV,
					      WlzThresholdType.WLZ_THRESH_LOW);
	  area = WlzObject.WlzArea(thrObj);
	  System.out.println("thresholded object area = " + area);

	  System.out.println("calling the beast");
	  int [] nLObjs = {0};
	  WlzObject [][] lObjs = new WlzObject[1][];
	  WlzObject.WlzLabel(thrObj, nLObjs, lObjs, 1000, 1,
	  		     WlzConnectType.WLZ_8_CONNECTED);
	  System.out.println("number of components = " + nLObjs[0]);
	  for(idx = 0; idx < nLObjs[0]; ++idx)
	  {
	    area = WlzObject.WlzArea(lObjs[0][idx]);
	    System.out.println("idx = " + idx + ", area = " + area);
	  }
	}
      }
      catch (IOException e)
      {
	System.err.println(e);
	System.exit(1);
      }
      catch (WlzException e)
      {
	System.err.println(e);
	System.exit(1);
      }
    }
  }
}
