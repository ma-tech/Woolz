import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzConstruct3D
{
  public static void main (String[] args)
  {
    boolean	showUsageFlg = false;
    int		optIdx = 0,
    		exitStatus = 0;
    String 	outFile = "-";

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
      else
      {
        showUsageFlg = true;
	exitStatus = 1;
      }
      ++optIdx;
    }
    if(showUsageFlg)
    {
      System.err.println("Usage: WlzConstruct3D [-o <out file>] [-v] [-h] " +
      			 "[<input files>]");
      System.err.println(
      "Constructs a 3D Woolz object from the Woolz object(s) given on the");
      System.err.println(
      "command line.");
      System.err.println(
      "Example: WlzConstruct3D -o out.wlz obj0.wlz obj1.wlz obj2.wlz");
      System.err.println(
      "The Woolz objects are read from obj0.wlz, obj1.wlz and obj2.wlz and");
      System.err.println(
      "used to create a 3D Woolz object which is then written to out.wlz");
      System.exit(exitStatus);
    }
    else
    {
      try
      {
	int		idx,
			nInObj = 0;
	boolean		firstPass = true;
	String 		inFile = "-";
// 	WlzObject 	inObj[][] = new WlzObject[1][args.length - optIdx + 2];
	WlzObject 	inObj[] = new WlzObject[args.length - optIdx + 2];
	WlzFileStream 	inFS = null,
			outFS = null;

	if(outFile.equals("-") == false)
	{
	  System.setOut(new PrintStream(
	  		new FileOutputStream(outFile)));
	}
	while(firstPass || (optIdx < args.length))
	{
	  if((firstPass && (optIdx == args.length)) ||
	     args[optIdx].equals("-"))
          {
	    inFile = "-";
	    inFS = new WlzFileStdInStream();
	  }
	  else
	  {
	    inFile = args[optIdx];
	    inFS = new WlzFileInputStream(args[optIdx]);
	  }
          inObj[nInObj] = WlzObject.WlzReadObj(inFS);
	  ++nInObj;
	  firstPass = false;
	  ++optIdx;
	}
	WlzObject obj3D = WlzObject.WlzConstruct3DObjFromObj(nInObj, inObj,
				      0, 0.0f, 0.0f, 0.0f);
	outFS = new WlzFileOutputStream(outFile);
	WlzObject.WlzWriteObj(outFS, obj3D);
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
