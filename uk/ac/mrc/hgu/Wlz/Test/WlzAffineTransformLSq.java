import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzAffineTransformLSq
{
  public static void main (String[] args)
  {
    boolean	showUsageFlg = false;
    int		optIdx = 0,
    		exitStatus = 0;
    String	inFile = "-",
    		outFile = "-";

    while((optIdx < args.length) && (args[optIdx].startsWith("-")) &&
          (args[optIdx].length() > 1))
    {
      if(args[optIdx].equals("-h"))
      {
        showUsageFlg = true;
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
      System.err.println(
      "Usage: java WlzAffineTransformLSq");
      System.err.println(
      "Computes transform from 2D tie points.");
      System.exit(exitStatus);
    }
    else
    {
      try
      {
	boolean firstPass = true;

	if(firstPass)
	{
	  WlzDVertex2 [] aryS2 = {new WlzDVertex2(0.0, 1.0),
	  			  new WlzDVertex2(4.0, 1.0),
	  			  new WlzDVertex2(4.0, 5.1),
	  			  new WlzDVertex2(0.0, 5.0)};
	  WlzDVertex2 [] aryT2 = {new WlzDVertex2(0.0, 1.0),
	  			  new WlzDVertex2(4.0, 1.0),
	  			  new WlzDVertex2(4.0, 5.0),
	  			  new WlzDVertex2(0.0, 5.0)};
	  WlzAffineTransform tr2 = WlzObject.WlzAffineTransformLSq2D(
	  				4, aryS2, 4, aryT2,
					WlzTransformType.WLZ_TRANSFORM_2D_REG);

	  WlzDVertex3 [] aryS3 = {new WlzDVertex3(0.0, 1.0, 0.0),
	  		 	  new WlzDVertex3(4.0, 1.0, 0.0),
	  			  new WlzDVertex3(4.0, 5.1, 0.1),
	  			  new WlzDVertex3(0.0, 5.0, 0.0)};
	  WlzDVertex3 [] aryT3 = {new WlzDVertex3(0.0, 1.0, 0.0),
	  			  new WlzDVertex3(4.0, 1.0, 0.0),
	  			  new WlzDVertex3(4.0, 5.0, 0.0),
	  			  new WlzDVertex3(0.0, 5.0, 0.0)};
	  WlzAffineTransform tr3 = WlzObject.WlzAffineTransformLSq3D(
	  				4, aryS3, 4, aryT3,
					WlzTransformType.WLZ_TRANSFORM_3D_REG);
	}
      }
      catch (WlzException e)
      {
	System.err.println(e);
	System.exit(1);
      }
    }
  }
}
