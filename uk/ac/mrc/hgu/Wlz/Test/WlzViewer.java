import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.awt.image.*;
import javax.swing.*;
import uk.ac.mrc.hgu.Wlz.*;

public class WlzViewer extends JFrame implements ActionListener
{
  private JPanel topPanel;
  private JMenuBar menuBar;
  private JMenu menuFile;
  private JMenuItem menuFileExit;
  private JScrollPane scroll;
  private JViewport viewport;
  private JPanel statusPanel;
  private JLabel statusName;
  private JLabel statusStats;
  private WlzImageView wlzImgV = null;

  public WlzViewer() throws IOException, WlzException
  {
    setTitle("WlzViewer");
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
    statusPanel = new JPanel();
    statusPanel.setLayout(new GridLayout(2, 1));
    topPanel.add(statusPanel, BorderLayout.SOUTH);
    statusName = new JLabel("", SwingConstants.LEFT);
    statusStats = new JLabel("", SwingConstants.LEFT);
    statusPanel.add(statusName);
    statusPanel.add(statusStats);
    scroll = new JScrollPane();
    wlzImgV = new WlzImageView();
    scroll.getViewport().add(wlzImgV);
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

  public class WlzImageView extends JPanel implements MouseListener,
    MouseMotionListener
  {
    private int	gType;
    private WlzObject obj = null;
    private WlzIBox2 bBox = null;
    private WlzGreyValueWSpace gVWSp = null;
    private BufferedImage bufImage = null;

    public WlzImageView()
    {
      setBackground(Color.white);
      addMouseListener(this);
    }

    public void setWlzObj(String fileName) throws IOException, WlzException
    {
      WlzFileStream in = null;

      if(fileName.equals("-"))
      {
	in = new WlzFileStdInStream();
      }
      else
      {
	in = new WlzFileInputStream(fileName);
      }
      setWlzObj(in, fileName);
    }

    public void setWlzObj(WlzFileStream in, String objName) throws WlzException
    {
      WlzObject	newObj;

      newObj = WlzObject.WlzReadObj(in);
      setWlzObj(newObj, objName);
    }

    public void setWlzObj(WlzObject newObj, String objName) throws WlzException
    {
      obj = newObj;
      bBox = WlzObject.WlzBoundingBox2I(obj);
      gVWSp = WlzObject.WlzGreyValueMakeWSp(newObj);
      Rectangle bRec = bBox.toRectangle();
      setPreferredSize(new Dimension(bRec.width, bRec.height));
      byte[] lut = new byte[256];
      for(int idN = 0; idN < 256; ++idN)
      {
        lut[idN] = (byte )idN;
      }
      IndexColorModel cMod = new IndexColorModel(8, 256, lut, lut, lut);
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
      statusName.setText(objName +
      			 " (" + WlzObject.WlzStringFromGreyType(gType) + ")");
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
      statusStats.setText(objStats);
    }
  }

  public void actionPerformed(ActionEvent event)
  {
    System.exit(0);
  }

  public static void main (String[] args)
  {
    int		optIdx = 0,
    		exitStatus = 0;
    String 	fileName = "-";
    WlzViewer 	wlzV = null;

    while(optIdx < args.length)
    {
      try
      {
	boolean	firstPass = true;

	while(firstPass || (optIdx < args.length))
	{
	  if((firstPass && (optIdx == args.length)) ||
	     args[optIdx].equals("-"))
          {
	    fileName = "-";
	  }
	  else
	  {
	    fileName = args[optIdx];
	  }
	  wlzV = new WlzViewer();
	  wlzV.wlzImgV.setWlzObj(fileName);
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
