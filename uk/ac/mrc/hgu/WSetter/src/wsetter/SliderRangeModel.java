// #pragma ident "MRC HGU $Id$"
/*!
* \file         SliderRangeModel.java
* \author       Nick Burton
* \date         May 2002
* \version      $Id$
* \note
*               Copyright
*               2002 Medical Research Council, UK.
*               All rights reserved.
*               All rights reserved.
* \par Address:
*               MRC Human Genetics Unit,
*               Western General Hospital,
*               Edinburgh, EH4 2XU, UK.
* \brief        data model for slider that holds double values
* \todo         -
* \bug          None known.
*/

package wsetter;

import wsetter.*;

import javax.swing.*;
import javax.swing.event.*;

public class SliderRangeModel implements BoundedRangeModel {
    protected ChangeEvent changeEvent = null;
    protected EventListenerList listenerList = new EventListenerList();

    protected int maximum = 255;
    protected int minimum = 0;
    protected int extent = 0;
    protected double value = 0.0;
    protected double multiplier = 1.0;
    protected boolean isAdjusting = false;
    final static boolean DEBUG = false;

    public SliderRangeModel() {
    }

    public double getMultiplier() {
        if (DEBUG) {
            System.out.println("In SliderRangeModel getMultiplier");
        }
        return multiplier;
    }

    public void setMultiplier(double multiplier) {
        if (DEBUG) {
            System.out.println("In SliderRangeModel setMultiplier");
        }
        this.multiplier = multiplier;
        fireStateChanged();
    }

    public int getMaximum() {
        if (DEBUG) {
            System.out.println("In SliderRangeModel getMaximum");
        }
        return maximum;
    }

    public void setMaximum(int newMaximum) {
        if (DEBUG) {
            System.out.println("In SliderRangeModel setMaximum");
        }
        setRangeProperties(value, extent, minimum, newMaximum, isAdjusting);
    }

    public int getMinimum() {
        return (int)minimum;
    }

    public void setMinimum(int newMinimum) {
        if (DEBUG) {
            System.out.println("In SliderRangeModel setMinimum");
        }
        setRangeProperties(value, extent, newMinimum, maximum, isAdjusting);
    }

    public int getValue() {
        if (DEBUG) {
            System.out.println("In SliderRangeModel getValue");
        }
        return (int)getDoubleValue();
    }

    public void setValue(int newValue) {
        if (DEBUG) {
            System.out.println("In SliderRangeModel setValue");
        }
        setDoubleValue((double)newValue);
    }

    public double getDoubleValue() {
        if (DEBUG) {
            System.out.println("In SliderRangeModel getDoubleValue");
        }
        return value;
    }

    public float getFloatValue() {
        if (DEBUG) {
            System.out.println("In SliderRangeModel getFloatValue");
        }
        return (float)value;
    }

    public void setDoubleValue(double newValue) {
        if (DEBUG) {
            System.out.println("In SliderRangeModel setDoubleValue");
        }
        setRangeProperties(newValue, extent, minimum, maximum, isAdjusting);
    }

    public void setFloatValue(float newValue) {
        if (DEBUG) {
            System.out.println("In SliderRangeModel setFloatValue");
        }
        setRangeProperties((double)newValue, extent, minimum, maximum, isAdjusting);
    }

    public int getExtent() {
        return (int)extent;
    }

    public void setExtent(int newExtent) {
        //Do nothing.
    }

    public boolean getValueIsAdjusting() {
        return isAdjusting;
    }

    public void setValueIsAdjusting(boolean b) {
        setRangeProperties(value, extent, minimum, maximum, b);
    }

    public void setRangeProperties(int newValue,
                                   int newExtent,
                                   int newMin,
                                   int newMax,
                                   boolean newAdjusting) {
        System.out.println("In SliderRangeModel setRangeProperties");
        setRangeProperties((double)newValue,
                           newExtent,
                           newMin,
                           newMax,
                           newAdjusting);
    }

    public void setRangeProperties(double newValue,
                                   int unusedExtent,
                                   int newMin,
                                   int newMax,
                                   boolean newAdjusting) {
        if (DEBUG) {
            System.out.println("setRangeProperties(): "
                                + "newValue = " + newValue
                                + "; newMax = " + newMax);
        }
        if (newMax <= minimum) {
            newMax = minimum + 1;
            if (DEBUG) {
                System.out.println("maximum raised by 1 to " + newMax);
            }
        }
        if (Math.round(newValue) > newMax) { //allow some rounding error
            newValue = newMax;
            if (DEBUG) {
                System.out.println("value lowered to " + newMax);
            }
        }

        boolean changeOccurred = false;
        if (newValue != value) {
            if (DEBUG) {
                System.out.println("value set to " + newValue);
            }
            value = newValue;
            changeOccurred = true; 
        }
        if (newMax != maximum) {
            if (DEBUG) {
                System.out.println("maximum set to " + newMax);
            }
            maximum = newMax;
            changeOccurred = true;
        }
        if (newMin != minimum) {
            if (DEBUG) {
                System.out.println("minimum set to " + newMin);
            }
            minimum = newMin;
            changeOccurred = true;
        }
        if (newAdjusting != isAdjusting) {
            maximum = newMax;
            isAdjusting = newAdjusting;
            changeOccurred = true;
        }

        if (changeOccurred) {
            fireStateChanged();
        }
    }

    /* 
     * The rest of this is event handling code copied from 
     * DefaultBoundedRangeModel. 
     */
    public void addChangeListener(ChangeListener l) {
        listenerList.add(ChangeListener.class, l);
    }

    public void removeChangeListener(ChangeListener l) {
        listenerList.remove(ChangeListener.class, l);
    }

    protected void fireStateChanged() {
        Object[] listeners = listenerList.getListenerList();
        for (int i = listeners.length - 2; i >= 0; i -=2 ) {
            if (listeners[i] == ChangeListener.class) {
                if (changeEvent == null) {
                    changeEvent = new ChangeEvent(this);
                }
                ((ChangeListener)listeners[i+1]).stateChanged(changeEvent);
            }
        }
    }
}
