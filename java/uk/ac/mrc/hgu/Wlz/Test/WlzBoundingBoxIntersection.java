
import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzBoundingBoxIntersection 
{
  public static void main (String[] args) 
  {
    boolean	showUsageFlg = false;
    int		optIdx = 0,
    		exitStatus = 0;
    double	dist = 0.0,
    		theta = 0.0,
		phi = 0.0,
		zeta = 0.0;
    String 	outFile = "-";
    WlzThreeDViewStruct viewSt;

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
      else if(args[optIdx].equals("-d"))
      {
	dist = Double.parseDouble(args[optIdx + 1]);
	++optIdx;
      }
      else if(args[optIdx].equals("-t"))
      {
	theta = Math.toRadians(Double.parseDouble(args[optIdx + 1]));
	++optIdx;
      }
      else if(args[optIdx].equals("-p"))
      {
	phi = Math.toRadians(Double.parseDouble(args[optIdx + 1]));
	++optIdx;
      }
      else if(args[optIdx].equals("-z"))
      {
	zeta = Math.toRadians(Double.parseDouble(args[optIdx + 1]));
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
      System.err.println("Usage: WlzBoundingBoxIntersection " +
		   "[-o <out file>] [-h]\n" +
		   "[-d <dist>] [-t <theta>] [-p <phi>] [-z <zeta>]\n" +
		   "[<input files>]");
      System.err.println(
      "Gets the bounding box of the given Woolz object and it's" +
      "intersection with the plane ZZZ");
      System.err.println(
      "Example: WlzBoundingBoxIntersection -o out.txt myobj.wlz");
      System.err.println(
      "The input Woolz object is read from myobj.wlz. Both the bounding\n" +
      "box and the intersection of the bounding box with the plane of\n" +
      "section are computed and written to out.txt.");
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
          System.setOut(new PrintStream(new FileOutputStream(outFile)));
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
	  
	  System.out.println(bBox.xMin + " " +
		     bBox.yMin + " " +
		     bBox.zMin + " " +
		     bBox.xMax + " " +
		     bBox.yMax + " " +
		     bBox.zMax);

	  viewSt =  obj0.WlzMake3DViewStruct(
			      WlzObjectType.WLZ_3D_VIEW_STRUCT);
	  obj0.WlzInit3DViewStruct(viewSt, obj0);
	  obj0.Wlz3DViewSetDist(viewSt, dist);
	  obj0.Wlz3DViewSetTheta(viewSt, theta);
	  obj0.Wlz3DViewSetPhi(viewSt, phi);
	  obj0.Wlz3DViewSetZeta(viewSt, zeta);

	  int[] maxVtx = new int[1];
	  WlzDVertex3 [][] vtxArray = new WlzDVertex3[1][];

	  System.out.println("intersection:");
	  try 
	  {
	    int nVtx = obj0.Wlz3DViewGetBoundingBoxIntersectionA(
			 viewSt, maxVtx, vtxArray);
	    System.out.println("number of vertices: " + nVtx);
	    System.out.println("vertices:");
	    for (int i=0; i< maxVtx[0]; i++)
			  System.out.println("  " +
					     vtxArray[0][i].vtX + " " +
					     vtxArray[0][i].vtY + " " +
					     vtxArray[0][i].vtZ);
	  }
	  catch(Exception exp) 
	  {
	    System.out.println(exp);
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
