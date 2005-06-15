import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzArrayStats
{
  public static void main (String[] args) throws WlzException
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
      System.err.println("Usage: WlzArrayStats [-o <out file>] [-v] [-h] " +
      			 "[<input files>]");
      System.err.println(
      "Calculates the area, minimum, maximum, sum, sum of squares, mean");
      System.err.println(
      "and the standard deviation of the values in a 2D or 3D array with");
      System.err.println(
      "the same bounding box and values as the input 2D or 3D Woolz object.");
      System.err.println(
      "If the verbose output flag is not set then the following are written");
      System.err.println("to the output file in the following order:");
      System.err.println(
      "<area> <grey type> <min> <max> <sum> <sum of sq> <mean> <std dev>.");
      System.err.println(
      "The input object is read from stdin and values are written to stdout");
      System.err.println("unless the filenames are given.");
      System.err.println("Example: WlzArrayStats -o stats.txt myobj.wlz");
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

          WlzObject obj = WlzObject.WlzReadObj(in);

	  int gType = WlzObject.WlzGetObjectType(obj);
	  int gCnt = 0;
	  double gMin = 0.0;
	  double gMax = 0.0;
	  double gSum = 0.0;
	  double gSumSq = 0.0;
	  double gMean = 0.0;
	  double gStdDev = 0.0;

	  if(gType == WlzObjectType.WLZ_2D_DOMAINOBJ)
	  {
	    WlzIBox2 bBox = WlzObject.WlzBoundingBox2I(obj);
	    WlzIVertex2 org = new WlzIVertex2(bBox.xMin, bBox.yMin);
	    WlzIVertex2 size = new WlzIVertex2(bBox.xMax - bBox.xMin + 1,
	    				       bBox.yMax - bBox.yMin + 1);
	    byte dat[][][] = new byte[1][][];
	    WlzIVertex2 datSize[] = new WlzIVertex2[1];
	    WlzObject.WlzToUArray2D(datSize, dat, obj, org, size, 0);
	    gMin = dat[0][0][0];
	    gMax = dat[0][0][0];
	    gCnt = size.vtX * size.vtY;
	    for(int idY = 0; idY < size.vtY; ++idY)
	    {
	      for(int idX = 0; idX < size.vtX; ++idX)
	      {
	        double gVal = dat[0][idY][idX];
		if(gVal < gMin)
		{
		  gMin = gVal;
		}
		else if(gVal > gMax)
		{
		  gMax = gVal;
		}
		gSum += gVal;
		gSumSq += gVal * gVal;
	      }
	    }
	  }
	  else if(gType == WlzObjectType.WLZ_3D_DOMAINOBJ)
	  {
	    WlzIBox3 bBox = WlzObject.WlzBoundingBox3I(obj);
	    WlzIVertex3 org = new WlzIVertex3(bBox.xMin, bBox.yMin, bBox.zMin);
	    WlzIVertex3 size = new WlzIVertex3(bBox.xMax - bBox.xMin + 1,
	    				       bBox.yMax - bBox.yMin + 1,
	    				       bBox.zMax - bBox.zMin + 1);
	    byte dat[][][][] = new byte[1][][][];
	    WlzIVertex3 datSize[] = new WlzIVertex3[1];
	    WlzObject.WlzToUArray3D(datSize, dat, obj, org, size, 0);
	    gMin = dat[0][0][0][0];
	    gMax = dat[0][0][0][0];
	    gCnt = size.vtX * size.vtY + size.vtZ;
	    for(int idZ = 0; idZ < size.vtZ; ++idZ)
	    {
	      for(int idY = 0; idY < size.vtY; ++idY)
	      {
		for(int idX = 0; idX < size.vtX; ++idX)
		{
		  double gVal = dat[0][idZ][idY][idX];
		  if(gVal < gMin)
		  {
		    gMin = gVal;
		  }
		  else if(gVal > gMax)
		  {
		    gMax = gVal;
		  }
		  gSum += gVal;
		  gSumSq += gVal * gVal;
		}
	      }
	    }
	  }
	  else
	  {
	    String errMsg[] = new String[1];
	    throw new WlzException(
		      WlzObject.WlzStringFromErrorNum(
	    	      WlzErrorNum.WLZ_ERR_OBJECT_TYPE, errMsg));
	  }
	  gMean = gSum / gCnt;
	  gStdDev = Math.sqrt((gSum - (gSum * gSum / gCnt)) / (gCnt - 1));
	  String gTypeS = WlzObject.WlzStringFromGreyType(gType);
	  if(verboseFlg)
	  {
	    System.out.println("file      " + inFile);
	    System.out.println("area      " + gCnt);
	    System.out.println("grey type " + gTypeS);
	    System.out.println("min       " + gMin);
	    System.out.println("max       " + gMax);
	    System.out.println("sum       " + gSum);
	    System.out.println("sum sq    " + gSumSq);
	    System.out.println("mean      " + gMean);
	    System.out.println("std dev   " + gStdDev);
	  }
	  else
	  {
	    System.out.println(gCnt + " " +
			       gTypeS + " " +
			       gMin + " " +
			       gMax + " " +
			       gSum + " " +
			       gSumSq + " " +
			       gMean + " " +
			       gStdDev);
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
