/* file name  : TVWlzObject.java
 * authors    : Tom Perry <tperry@hgu.mrc.ac.uk>
 * created    : Fri 13 Jul 2007 14:43:23 BST
 * copyright  : MRC Human Genetics Unit
 *
 * modifications:
 *
 */

package uk.ac.mrc.hgu.TVImage;

import java.awt.image.*;
import uk.ac.mrc.hgu.Wlz.*;

/** 
 * Wrapper for a Woolz object.
 *
 * Woolz is the HGU's in-house suite of image processing tools, and
 * though it is written in C, it is accesible using automatically generated
 * JNI code in uk.ac.mrc.hgu.Wlz.
 * 
 * A Woolz wrapper class already exists (see uk.ac.mrc.hgu.Wlz.WlzObject.java)
 * but provides funcionality that is inspecific to TreeView.
 * 
 * @author Tom Perry <tperry@hgu.mrc.ac.uk>
 * @version 
 */
class TVWlzObject
{
	private WlzObject wo = null;

	public TVWlzObject(BufferedImage img) { setWlzObject(bufImg2WlzObj(img)); }

	public void setWlzObject(WlzObject wo) { this.wo = wo; }
	public WlzObject getWlzObject() { return wo; }

	/** 
	 * Converts a java.awt.image.BufferedImage to a uk.ac.mrc.hgu.Wlz.WlzObject
	 * 
	 * @param img the BufferedImage to be converted
	 * @return a WlzObject corresponding to the input BufferedImage
	 */
	private WlzObject bufImg2WlzObj(BufferedImage img)
	{
		int w = img.getWidth();
		int h = img.getHeight();
		int[][] pixels = new int[h][w];
		for (int y=0; y<h; y++)
			for (int x=0; x<w; x++)
				pixels[y][x] = img.getRGB(x,y);
		WlzObject wo = null;
		int err[] = new int[1];
		WlzIVertex2 dimensions = new WlzIVertex2(w,h);
		WlzIVertex2 org = new WlzIVertex2(0,0);

		int manyFacts = 0;
		String facts[] = {""};

		wo = WlzObject.WlzFromIArray2D(dimensions,pixels,org,err);

		return wo;
	}

	/** 
	 * TODO: does this work?
	 * 
	 * @param a 
	 * @param b 
	 * @return 
	 */
	public static double getJaccardCorrelation(TVWlzObject a, TVWlzObject b)
	{
		double corr = 0.0;
		if (a == null || b == null)
			return corr;
		try
		{
			corr = WlzObject.WlzArea(WlzObject.WlzIntersect2(a.getWlzObject(),b.getWlzObject())) / WlzObject.WlzArea(WlzObject.WlzUnion2(a.getWlzObject(),b.getWlzObject()));
		}
		catch (WlzException we) { we.printStackTrace(); }
		return corr;
	}

	/** 
	 * No longer relevant
	 *
	 * @param a 
	 * @param b 
	 * @return 
	 */
	public static double getNormalisedCorrelation(TVWlzObject a, TVWlzObject b)
	{
		double[] aStdDev = new double[1], bStdDev = new double[1];
		try
		{
			WlzObject.WlzGreyStats(a.getWlzObject(),null,null,null,null,null,aStdDev,null);
			WlzObject.WlzGreyStats(b.getWlzObject(),null,null,null,null,null,bStdDev,null);
		}
		catch (WlzException we)
		{
			we.printStackTrace();
		}
		double correlation = getWlzCorrelation(a,b);
		System.out.println("a StdDev: " + aStdDev[0]);
		System.out.println("b StdDev: " + bStdDev[0]);
		System.out.println("correlation: " + correlation);

		return correlation/(aStdDev[0]*bStdDev[0]);
	}

	/** 
	 * No longer relevant
	 * 
	 * @param a 
	 * @param b 
	 * @return 
	 */
	public static double getWlzCorrelation(TVWlzObject a, TVWlzObject b)
	{
		double corr = 0;
		try { corr = WlzObject.WlzCCorS2D(a.getWlzObject(),b.getWlzObject(),0); }
		catch (WlzException we) { we.printStackTrace(); }
		return corr;
	}
}
