import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzGetBackground
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
      System.err.println("Usage: WlzGetBackground [-o <out file>] [-v] [-h] " +
      			 "[<input files>]");
      System.err.println(
      "Gets the background value of the given Woolz object(s).");
      System.err.println("Example: WlzGetBackground -o out.txt myobj.wlz");
      System.err.println(
      "The input Woolz object is read from myobj.wlz. The statistics are");
      System.err.println("calculated and  written to out.txt.");
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
	  WlzPixelV bgd = WlzObject.WlzGetBackground(obj0);
	  switch(bgd.getType())
	  {
	    case WlzGreyType.WLZ_GREY_LONG:
	      System.out.println(bgd.getLongValue() + " " +
			       WlzObject.WlzStringFromGreyType(bgd.getType()));
	      break;
	    case WlzGreyType.WLZ_GREY_INT:
	      System.out.println(bgd.getIntValue() + " " +
			       WlzObject.WlzStringFromGreyType(bgd.getType()));
	      break;
	    case WlzGreyType.WLZ_GREY_SHORT:
	      System.out.println(bgd.getShortValue() + " " +
			       WlzObject.WlzStringFromGreyType(bgd.getType()));
	      break;
	    case WlzGreyType.WLZ_GREY_UBYTE:
	      System.out.println(bgd.getByteValue() + " " +
			       WlzObject.WlzStringFromGreyType(bgd.getType()));
	      break;
	    case WlzGreyType.WLZ_GREY_FLOAT:
	      System.out.println(bgd.getFloatValue() + " " +
			       WlzObject.WlzStringFromGreyType(bgd.getType()));
	      break;
	    case WlzGreyType.WLZ_GREY_DOUBLE:
	      System.out.println(bgd.getDoubleValue() + " " +
			       WlzObject.WlzStringFromGreyType(bgd.getType()));
	      break;
	    default:
	      System.err.println("Unknown grey type " + bgd.getType());
	      break;
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
