import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzDitherViewer extends JFrame implements ActionListener
{
  private JPanel topPanel;
  private JMenuBar menuBar;
  private JMenu menuFile;
  private JMenuItem menuFileExit;
  private JScrollPane scroll;
  private JViewport viewport;
  private WlzDitherView wlzDthV = null;

  public WlzDitherViewer() throws IOException, WlzException
  {
    setTitle("WlzDitherViewer");
    topPanel = new JPanel();
    topPanel.setLayout(new BorderLayout());
    getContentPane().add(topPanel);
    menuBar = new JMenuBar();
    menuFile = new JMenu("File");
    menuBar.add(menuFile);
    menuFileExit = new JMenuItem("Exit", 'x');
    menuFileExit.addActionListener(this);
    menuFile.add(menuFileExit);
    topPanel.add(menuBar, BorderLayout.NORTH);
    scroll = new JScrollPane();
    wlzDthV = new WlzDitherView();
    scroll.getViewport().add(wlzDthV);
    topPanel.add(scroll, BorderLayout.CENTER);
    setVisible(true);
    addWindowListener(new WindowAdapter()
    		      {
		        public void windowClosing(WindowEvent e)
			{
			  dispose();
			  System.exit(0);
			}
		      });
  }

  public class WlzDitherView extends JPanel implements MouseListener,
    MouseMotionListener
  {
    private int	gType;
    private WlzObject obj = null;
    private WlzIBox2 bBox = null;
    private WlzGreyValueWSpace gVWSp = null;
    private BufferedImage bufImage = null;

    public WlzDitherView()
    {
      setBackground(Color.white);
      addMouseListener(this);
    }

    public void setWlzObjs(String fileName0, String fileName1)
      throws IOException, WlzException
    {
      WlzFileStream in0 = null,
      		    in1 = null;

      if(fileName0.equals("-"))
      {
	in0 = new WlzFileStdInStream();
      }
      else
      {
	in0 = new WlzFileInputStream(fileName0);
      }
      if(fileName1.equals("-"))
      {
	in1 = new WlzFileStdInStream();
      }
      else
      {
	in1 = new WlzFileInputStream(fileName1);
      }
      setWlzObjs(in0, in1);
    }

    public void setWlzObjs(WlzFileStream in0, WlzFileStream in1)
      throws WlzException
    {
      setWlzObjs(WlzObject.WlzReadObj(in0), WlzObject.WlzReadObj(in1));
    }

    public void setWlzObjs(WlzObject newObj0, WlzObject newObj1)
      throws WlzException
    {
      WlzPixelV [] lPV = {null},
      		   hPV = {null};

      WlzObject.WlzGreySetRange(newObj0,
			        new WlzPixelV((int )0),
			        new WlzPixelV((int )255),
				new WlzPixelV((int )127),
				new WlzPixelV((int )0));
      WlzObject.WlzGreySetRange(newObj1,
				new WlzPixelV((int )0),
				new WlzPixelV((int )255),
				new WlzPixelV((int )255),
				new WlzPixelV((int )0));
      obj = WlzObject.WlzImageArithmetic(newObj0,
      	      WlzObject.WlzGreyDitherObj(newObj1, 0x0080),
	      WlzBinaryOperatorType.WLZ_BO_OR, 0);
      bBox = WlzObject.WlzBoundingBox2D(obj);
      gVWSp = WlzObject.WlzGreyValueMakeWSp(obj);
      Rectangle bRec = bBox.toRectangle();
      setPreferredSize(new Dimension(bRec.width, bRec.height));
      byte[] lutR = new byte[256],
      	     lutG = new byte[256],
	     lutB = new byte[256];
      for(int idN = 0; idN < 128; ++idN)
      {
	  lutR[idN] = (byte )(254 - (2 * idN));
	  lutG[idN] = (byte )(254 - (2 * idN));
	  lutB[idN] = (byte )(254 - (2 * idN));
      }
      for(int idN = 128; idN < 256; ++idN)
      {
	  lutR[idN] = (byte )0;
	  lutG[idN] = (byte )0;
	  lutB[idN] = (byte )255;
      }
      IndexColorModel cMod = new IndexColorModel(8, 256, lutR, lutG, lutB);
      byte imageData[] = new byte[bRec.width * bRec.height];
      int valI;
      for(int idY = 0; idY < bRec.height; ++idY)
      {
	for(int idX = 0; idX < bRec.width; ++idX)
	{
	  valI = WlzObject.WlzGreyValueGetI(gVWSp,
	  				    0.0,
					    (double )(bRec.y + idY),
					    (double )(bRec.x + idX));
	  if(valI < 0)
	  {
	    valI = 0;
	  }
	  else if(valI > 255)
	  {
	    valI = 255;
	  }
	  imageData[(idY * bRec.width) + idX] = (byte )valI;
	}
      }
      gType = WlzObject.WlzGreyValueGetGreyType(gVWSp);
      DataBuffer dataBuf = new DataBufferByte(imageData,
      					      bRec.width * bRec.height);
      WritableRaster ras = Raster.createPackedRaster(dataBuf,
      						     bRec.width, bRec.height,
						     8, null);
      bufImage = new BufferedImage(cMod, ras, false, null);
      Graphics2D g2d = bufImage.createGraphics();
    }

    public void paint(Graphics g)
    {
      if(bufImage != null)
      {
	g.drawImage(bufImage, 0, 0, null);
      }
    }

    public void mouseExited(MouseEvent e) {}
    public void mouseEntered(MouseEvent e) {}
    public void mouseReleased(MouseEvent e) {}
    public void mousePressed(MouseEvent e) {}
    public void mouseMoved(MouseEvent e) {}

    public void mouseDragged(MouseEvent e)
    {
      updateStats(e.getX(), e.getY());
    }

    public void mouseClicked(MouseEvent e)
    {
      updateStats(e.getX(), e.getY());
    }

    private void updateStats(int relX, int relY)
    {
      String objStats;

      if((obj == null) || (gVWSp == null))
      {
        objStats = "";
      }
      else
      {
	Point pos = new Point(bBox.xMin + relX, bBox.yMin + relY);
	objStats = "G(" + pos.x + ", " + pos.y + ") = " +
		   WlzObject.WlzGreyValueGetI(gVWSp,
		   			      0.0,
					      (double )(pos.y),
					      (double )(pos.x));
      }
    }
  }

  public void actionPerformed(ActionEvent event)
  {
    System.exit(0);
  }

  public static void main (String[] args)
  {
    WlzDitherViewer 	wlzDV = null;

    if(args.length == 2)
    {
      try
      {
	wlzDV = new WlzDitherViewer();
	wlzDV.wlzDthV.setWlzObjs(args[0], args[1]);
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
    else
    {
      System.err.println("Usage: WlzDitherViewer file0 file1");
      System.exit(1);
    }
  }
}
