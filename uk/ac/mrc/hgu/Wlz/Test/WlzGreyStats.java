import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzGreyStats
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
      System.err.println("Usage: WlzGreyStats [-o <out file>] [-v] [-h] " +
      			 "[<input files>]");
      System.err.println(
      "Calculates the area, minimum, maximum, sum, sum of squares, mean");
      System.err.println(
      "and the standard deviation of the input 2D or 3D Woolz object's");
      System.err.println("grey values.");
      System.err.println(
      "If the verbose output flag is not set then the following are written");
      System.err.println("to the output file in the following order:");
      System.err.println(
      "<area> <grey type> <min> <max> <sum> <sum of sq> <mean> <std dev>.");
      System.err.println(
      "The input object is read from stdin and values are written to stdout");
      System.err.println("unless the filenames are given.");
      System.err.println("Example: WlzGreyStats -o stats.txt myobj.wlz");
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
	  int [] gTypeI = {0};
	  double [] gMin = {0.0};
	  double [] gMax = {0.0};
	  double [] gSum = {0.0};
	  double [] gSumSq = {0.0};
	  double [] gMean = {0.0};
	  double [] gStdDev = {0.0};
          WlzObject obj0 = WlzObject.WlzReadObj(in);
	  int gArea = WlzObject.WlzGreyStats(obj0, gTypeI, gMin, gMax,
					     gSum, gSumSq, gMean, gStdDev);
	  String gTypeS = WlzObject.WlzStringFromGreyType(gTypeI[0]);
	  if(verboseFlg)
	  {
	    System.out.println("file      " + inFile);
	    System.out.println("area      " + gArea);
	    System.out.println("grey type " + gTypeS);
	    System.out.println("min       " + gMin[0]);
	    System.out.println("max       " + gMax[0]);
	    System.out.println("sum       " + gSum[0]);
	    System.out.println("sum sq    " + gSumSq[0]);
	    System.out.println("mean      " + gMean[0]);
	    System.out.println("std dev   " + gStdDev[0]);
	  }
	  else
	  {
	    System.out.println(gArea + " " +
			       gTypeS + " " +
			       gMin[0] + " " +
			       gMax[0] + " " +
			       gSum[0] + " " +
			       gSumSq[0] + " " +
			       gMean[0] + " " +
			       gStdDev[0]);
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
