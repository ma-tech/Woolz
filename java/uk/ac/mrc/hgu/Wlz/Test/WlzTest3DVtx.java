import java.io.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzTest3DVtx
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
    

    System.err.println("blah blah blah");
       
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
      System.err.println("Usage: WlzTest3DVtx [-h] x y z [<input files>]");
      System.err.println("Test some transform functions for 3D vertices");
      System.err.println("Options are:");
      System.err.println("  -h  Help - print this usage message.");
      System.exit(exitStatus);
    }
    else
    {
      try
      {
        System.err.println("blah blah blah");

        Double	d;
	WlzDVertex3 vtx = new WlzDVertex3();
	System.err.println("blah blah blah");

	
	if( optIdx < args.length )
	{
	  d = Double.valueOf(args[optIdx]);
	  vtx.vtX = d.doubleValue();
	  optIdx++;
	}
	if( optIdx < args.length )
	{
	  d = Double.valueOf(args[optIdx]);
	  vtx.vtY = d.doubleValue();
	  optIdx++;
	}
	if( optIdx < args.length )
	{
	  d = Double.valueOf(args[optIdx]);
	  vtx.vtZ = d.doubleValue();
	  optIdx++;
	}
				
	if( optIdx == args.length )
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
	System.out.println(vtx.vtX + " " + vtx.vtY + " " + vtx.vtZ);
	
	/* create a viewStruct, initialize it */
	WlzIBox3 bBox3 = WlzObject.WlzBoundingBox3I(obj0);
	WlzThreeDViewStruct vs=obj0.WlzMake3DViewStruct(160);
	obj0.Wlz3DViewSetTheta(vs,Math.PI/2.0d);
	obj0.Wlz3DViewSetPhi(vs,Math.PI/2.0d);
	obj0.Wlz3DViewSetDist(vs,0.0);
	obj0.Wlz3DViewSetFixed2(vs,(bBox3.xMin+bBox3.xMax)/2.0,
				(bBox3.yMin+bBox3.yMax)/2.0,
				(bBox3.zMin+bBox3.zMax)/2.0);
	obj0.WlzInit3DViewStruct(vs,obj0);
	
	/* now set the vertex to be transformed and see what we get */
	/*vtx.vtX = 0;
	vtx.vtY = 0;
	vtx.vtZ = 0;*/
	WlzDVertex3 dstVtx[] = new WlzDVertex3[1];
	obj0.Wlz3DSectionTransformInvVtxR(vs, vtx, dstVtx);
	System.out.println("vertex = " + vtx.vtX + " " + vtx.vtY + " " +
	vtx.vtZ);
	System.out.println("new vertex = " + dstVtx[0].vtX + " " + dstVtx[0].vtY +
		" " + dstVtx[0].vtZ);
		
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
