import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzFacts
{
  public static void main (String[] args)
  {
    boolean	showUsageFlg = false;
    int		manyFacts = 0,
    		optIdx = 0,
    		exitStatus = 0;
    WlzFileStream in = null;
    WlzObject	obj0 = null;
    String	facts[] = {""};

    while((optIdx < args.length) && (args[optIdx].startsWith("-")) &&
          (args[optIdx].length() > 1))
    {
      if(args[optIdx].equals("-h"))
      {
        showUsageFlg = true;
      }
      else if(args[optIdx].equals("-f"))
      {
        manyFacts = 0;
      }
      else if(args[optIdx].equals("-m"))
      {
        manyFacts = 1;
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
      System.err.println("Usage: WlzFacts [-h] [-f] [-m] [<input files>]");
      System.err.println("Print facts about the input woolz objects.");
      System.err.println("Options are:");
      System.err.println("  -h  Help - print this usage message.");
      System.err.println("  -f  Few facts (default).");
      System.err.println("  -m  Many facts (more detailed information).");
      System.exit(exitStatus);
    }
    else
    {
      try
      {
	boolean	firstPass = true;

	while(firstPass || (optIdx < args.length))
	{
	  if((firstPass && (optIdx == args.length)) ||
	     args[optIdx].equals("-"))
          {
	    in = new WlzFileStdInStream();
	  }
	  else
	  {
	    in = new WlzFileInputStream(args[optIdx]);
	  }
	  obj0 = WlzObject.WlzReadObj(in);
	  WlzObject.WlzObjectFacts(obj0, null, facts, manyFacts);
	  System.out.println(facts[0]);
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
