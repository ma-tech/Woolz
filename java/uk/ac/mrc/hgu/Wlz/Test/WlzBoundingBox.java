import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzBoundingBox
{
  public static void main (String[] args)
  {
    boolean	showUsageFlg = false,
    		verboseFlg = false;
    int		optIdx = 0,
    		exitStatus = 0;
    String outFile = "-";

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
	  ++optIdx;
	  if(outFile.length() < 1)
	  {
	    showUsageFlg = true;
	    exitStatus = 1;
	  }
	}
      }
      else if(args[optIdx].equals("-v"))
      {
        verboseFlg = true;
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
      System.err.println("Usage: WlzBoundingBox [-o <out file>] [-v] [-h] " +
      			 "[<input files>]");
      System.err.println(
      "Gets the background value of the given Woolz object(s).");
      System.err.println("Example: WlzBoundingBox -o out.txt myobj.wlz");
      System.err.println(
      "The input Woolz object is read from myobj.wlz. The bounding box");
      System.err.println("is calculated and  written to out.txt.");
      System.exit(exitStatus);
    }
    else
    {
      try
      {
	boolean	firstPass = true;
	String inFile = "-";
	if(outFile.equals("-") == false)
	{
	  System.setOut(new PrintStream(
	  		new FileOutputStream(outFile)));
	}
	while(firstPass || (optIdx < args.length))
	{
          WlzFileStream in = null;
	  if((firstPass && (optIdx == args.length)) ||
	     args[optIdx].equals("-"))
          {
	    inFile = "-";
	    in = new WlzFileStdInStream();
	  }
	  else
	  {
	    inFile = args[optIdx];
	    in = new WlzFileInputStream(args[optIdx]);
	  }
          WlzObject obj0 = WlzObject.WlzReadObj(in);
	  WlzIBox3 bBox = WlzObject.WlzBoundingBox3D(obj0);
	  if(verboseFlg)
	  {
	    System.out.println("file: " + inFile);
	    System.out.println("xMin: " + bBox.xMin);
	    System.out.println("yMin: " + bBox.yMin);
	    System.out.println("zMin: " + bBox.zMin);
	    System.out.println("xMax: " + bBox.xMax);
	    System.out.println("yMax: " + bBox.yMax);
	    System.out.println("zMax: " + bBox.zMax);
	  }
	  else
	  {
	    System.out.println(bBox.xMin + " " +
	    		       bBox.yMin + " " +
	    		       bBox.zMin + " " +
	    		       bBox.xMax + " " +
	    		       bBox.yMax + " " +
	    		       bBox.zMax);
	  }
	  firstPass = false;
	  ++optIdx;
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
