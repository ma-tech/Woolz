import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzBndGreyInvert
{
  public static void main (String[] args)
  {
    boolean	showUsageFlg = false;
    int		optIdx = 0,
    		exitStatus = 0;
    String	outFile = "out.wlz";

    while((optIdx < args.length) && (args[optIdx].startsWith("-")) &&
          (args[optIdx].length() > 1))
    {
      if(args[optIdx].equals("-h"))
      {
        showUsageFlg = true;
      }
      else if(args[optIdx].equals("-o"))
      {
        if((optIdx + 1) < args.length)
	{
	  outFile = args[optIdx + 1];
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
      System.err.println("Usage: WlzGreyInvert [-o output] [-h] " +
      			 "<input object>");
      System.err.println(
      "Reads a Woolz object, inverts it's grey values and then writes it");
      System.err.println(
      "out.");
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
	  WlzObject.WlzBndGreyInvert(inObj);
	  WlzFileStream out = new WlzFileOutputStream(outFile);
	  WlzObject.WlzWriteObj(out, inObj);
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
