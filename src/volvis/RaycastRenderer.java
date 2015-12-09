/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volvis;

import com.jogamp.opengl.util.texture.Texture;
import com.jogamp.opengl.util.texture.awt.AWTTextureIO;
import gui.RaycastRendererPanel;
import gui.TransferFunction2DEditor;
import gui.TransferFunctionEditor;
import java.awt.image.BufferedImage;
import java.util.ArrayList;
import javax.media.opengl.GL;
import javax.media.opengl.GL2;
import util.TFChangeListener;
import util.VectorMath;
import volume.GradientVolume;
import volume.Volume;
import volume.VoxelGradient;

/**
 *
 * @author michel
 */
public class RaycastRenderer extends Renderer implements TFChangeListener {

    private Volume volume = null;
    private GradientVolume gradients = null;
    private Visualization vis = null;
    RaycastRendererPanel panel;
    TransferFunction tFunc;
    TransferFunctionEditor tfEditor;
    TransferFunction2DEditor tfEditor2D;
    private int currentMode = 0;
    private int SLICER = 0;
    private int MIP = 1;
    private int COMPOSITE = 2;
    private int TRANSFER2D = 3;
    private double sampleDistance = 2;
    private boolean debug = false;
    //private double[] pos = new double[3];


    //private boolean lowRes = true;
    
    

    //point to interpolate
        int tempValue = 0;
        double[] p = new double[4];
        //voxels
        double[] p000 = new double[4];
        double[] p001 = new double[4];
        double[] p010 = new double[4];
        double[] p011 = new double[4];
        double[] p100 = new double[4];
        double[] p101 = new double[4];
        double[] p110 = new double[4];
        double[] p111 = new double[4];
        //C values
        double[] c0 = new double[4];
        double[] c1 = new double[4];
        double[] c2 = new double[4];
        double[] c3 = new double[4];
        double[] c4 = new double[4];
        double[] c5 = new double[4];
        double[] c6 = new double[4];
        double[] c7 = new double[4];

    private boolean shading = false;
    private TFColor[] data;


    public void setCurrentMode(int currentMode) {
        this.currentMode = currentMode;
        vis.update();
    }
    
    public void setPos(double[] posi){
        //pos = posi;
    }

    public Visualization getVis() {
        return vis;
    }

    public void switchShading() {
        if(shading){
            shading = false;
        } else {
            shading = true;
        }
    }

    private int resolution(){
        if(interactiveMode){
            return 2;
        } else {
            return 1;
        }
    }

    public void setVis(Visualization vis) {
        this.vis = vis;
    }
    public int getMiddleMaxDimension(){
            int high = 0;
            if (volume.getDimX() > high){
                high = volume.getDimX();
            }
            if (volume.getDimY() > high){
                high = volume.getDimY();
            }
            if (volume.getDimZ() > high){
                high = volume.getDimZ();
            }
           return high/2;
    }
    
    
    public RaycastRenderer() {
        panel = new RaycastRendererPanel(this);
        panel.setSpeedLabel("0");
    }

    public void setVolume(Volume vol) {
        System.out.println("Assigning volume");
        volume = vol;
        System.out.println("Computing gradients");
        gradients = new GradientVolume(vol);
        data = new TFColor[volume.getDimX() * volume.getDimY() * volume.getDimZ()];
        for (int i = 0; i < data.length; i++){
            data[i] = null;
        }

        // set up image for storing the resulting rendering
        // the image width and height are equal to the length of the volume diagonal
        int imageSize = (int) Math.floor(Math.sqrt(vol.getDimX() * vol.getDimX() + vol.getDimY() * vol.getDimY()
                + vol.getDimZ() * vol.getDimZ()));
        if (imageSize % 2 != 0) {
            imageSize = imageSize + 1;
        }
        image = new BufferedImage(imageSize, imageSize, BufferedImage.TYPE_INT_ARGB);
        // create a standard TF where lowest intensity maps to black, the highest to white, and opacity increases
        // linearly from 0.0 to 1.0 over the intensity range
        tFunc = new TransferFunction(volume.getMinimum(), volume.getMaximum());
        
        // uncomment this to initialize the TF with good starting values for the orange dataset 
        tFunc.setTestFunc();
        
        
        tFunc.addTFChangeListener(this);
        tfEditor = new TransferFunctionEditor(tFunc, volume.getHistogram());
        
        tfEditor2D = new TransferFunction2DEditor(volume, gradients);
        tfEditor2D.addTFChangeListener(this);

        System.out.println("Finished initialization of RaycastRenderer");
    }

    public RaycastRendererPanel getPanel() {
        return panel;
    }

    public TransferFunction2DEditor getTF2DPanel() {
        return tfEditor2D;
    }
    
    public TransferFunctionEditor getTFPanel() {
        return tfEditor;
    }
    TFColor triLinearInterpolation1(double[] coord){
        TFColor interpolatedColor = new TFColor();
        //Rodrigo todo      
        
        //result interpol value
        double[] interpol = new double[3];
        VectorMath.setVector(interpol, 0, 0, 0);
        //init p values
        VectorMath.setVector(p, coord[0], coord[1], coord[2]);
        VectorMath.setVector(p000, Math.floor(p[0]), Math.floor(p[1]), Math.floor(p[2]));
        VectorMath.setVector(p001, p000[0], p000[1], p000[0]+1);
        VectorMath.setVector(p010, p000[0], p000[0]+1, p000[0]);
        VectorMath.setVector(p011, p000[0], p000[0]+1, p000[0]+1);
        VectorMath.setVector(p100, p000[0]+1, p000[0], p000[0]);
        VectorMath.setVector(p101, p000[0]+1, p000[0], p000[0]+1);
        VectorMath.setVector(p110, p000[0]+1, p000[0]+1, p000[0]);
        VectorMath.setVector(p111, p000[0]+1, p000[0]+1, p000[0]+1);
        
        //Get delta increments
        double deltaX = (p[0]-p000[0])/(p111[0]-p000[0]);
        double deltaY = (p[1]-p000[1])/(p111[1]-p000[1]);
        double deltaZ = (p[2]-p000[2])/(p111[2]-p000[2]);
        
        //GET colors p(x,y,z)
        tempValue = getVoxel(p000);
        TFColor p000C = tFunc.getColor(tempValue);
        VectorMath.setVector(p000, p000C.r, p000C.g, p000C.b);
        
        tempValue = getVoxel(p001);
        TFColor p001C = tFunc.getColor(tempValue);
        VectorMath.setVector(p001, p001C.r, p001C.g, p001C.b);
        
        tempValue = getVoxel(p010);
        TFColor p010C = tFunc.getColor(tempValue);
        VectorMath.setVector(p010, p010C.r, p010C.g, p010C.b);
        
        tempValue = getVoxel(p011);
        TFColor p011C = tFunc.getColor(tempValue);
        VectorMath.setVector(p011, p011C.r, p011C.g, p011C.b);
        
        
        tempValue = getVoxel(p100);
        TFColor p100C = tFunc.getColor(tempValue);
        VectorMath.setVector(p100, p100C.r, p100C.g, p100C.b);
        
        
        tempValue = getVoxel(p101);
        TFColor p101C = tFunc.getColor(tempValue);
        VectorMath.setVector(p101, p101C.r, p101C.g, p101C.b);
        
        
        tempValue = getVoxel(p110);
        TFColor p110C = tFunc.getColor(tempValue);
        VectorMath.setVector(p110, p110C.r, p110C.g, p110C.b);
        
        
        tempValue = getVoxel(p111);
        TFColor p111C = tFunc.getColor(tempValue);
        VectorMath.setVector(p111, p111C.r, p111C.g, p111C.b);
        
        
        //Calculate C values
        c0[0] = p000[0]; c0[1] = p000[1]; c0[2] = p000[2];
        VectorMath.minusVector(p100, p000, c1);
        VectorMath.minusVector(p010, p000, c2);
        VectorMath.minusVector(p001, p000, c3);
        c4[0] = p110[0] - p010[0] - p100[0] + p000[0]; 
        c4[1] = p110[1] - p010[1] - p100[1] + p000[1]; 
        c4[2] = p110[2] - p010[2] - p100[2] + p000[2];
        c5[0] = p011[0] - p001[0] - p010[0] + p000[0];
        c5[1] = p011[1] - p001[1] - p010[1] + p000[1];
        c5[2] = p011[2] - p001[2] - p010[2] + p000[2];
        c6[0] = p101[0] - p001[0] - p100[0] + p000[0];
        c6[1] = p101[1] - p001[1] - p100[1] + p000[1];
        c6[2] = p101[2] - p001[2] - p100[2] + p000[2];
        c7[0] = p111[0] - p011[0] - p101[0] - p110[0] + p100[0] + p001[0] + p010[0] - p000[0];
        c7[1] = p111[1] - p011[1] - p101[1] - p110[1] + p100[1] + p001[1] + p010[1] - p000[1];
        c7[2] = p111[2] - p011[2] - p101[2] - p110[2] + p100[2] + p001[2] + p010[2] - p000[2];
        //get interpolated value
        interpol[0] = c0[0] + c1[0]*deltaX + c2[0]*deltaY + c3[0]*deltaZ + c4[0]*deltaX*deltaY + c5[0]*deltaY*deltaZ + c6[0]*deltaZ*deltaX + c7[0]*deltaX*deltaY*deltaZ;
        interpol[1] = c0[1] + c1[1]*deltaX + c2[1]*deltaY + c3[1]*deltaZ + c4[1]*deltaX*deltaY + c5[1]*deltaY*deltaZ + c6[1]*deltaZ*deltaX + c7[1]*deltaX*deltaY*deltaZ;
        interpol[2] = c0[2] + c1[2]*deltaX + c2[2]*deltaY + c3[2]*deltaZ + c4[2]*deltaX*deltaY + c5[2]*deltaY*deltaZ + c6[2]*deltaZ*deltaX + c7[2]*deltaX*deltaY*deltaZ;

        
        return interpolatedColor;
    }
    public short trilinearInterpolation(double[] coord){
        double[] Vxyz = new double[8]; //each vertex V000, V001, V010, V011, V100, V101, V110, V111
        double alpha = 0;
        double beta = 0;
        double gamma = 0;
        //Check if inside the box (x,y,z) else just return cero
        if(  (0<=Math.floor(coord[0]) && Math.ceil(coord[0]) < volume.getDimX())  && //check [0,volDimX]
             (0<=Math.floor(coord[1]) && Math.ceil(coord[1]) < volume.getDimY())  && //check [0,volDimY]
             (0<=Math.floor(coord[2]) && Math.ceil(coord[2]) < volume.getDimZ())     //check [0,volDimZ]
        ){ 
            Vxyz[0] = volume.getVoxel((int)Math.floor(coord[0]) , (int)Math.floor(coord[1]), (int)Math.floor(coord[2])); //V000
            Vxyz[1] = volume.getVoxel((int)Math.ceil(coord[0]) , (int)Math.floor(coord[1]), (int)Math.floor(coord[2])); //V100
            Vxyz[2] = volume.getVoxel((int)Math.floor(coord[0]) , (int)Math.ceil(coord[1]), (int)Math.floor(coord[2])); //V010
            Vxyz[3] = volume.getVoxel((int)Math.ceil(coord[0]) , (int)Math.ceil(coord[1]), (int)Math.floor(coord[2])); //V110
            Vxyz[4] = volume.getVoxel((int)Math.floor(coord[0]) , (int)Math.floor(coord[1]), (int)Math.ceil(coord[2])); //V001
            Vxyz[5] = volume.getVoxel((int)Math.ceil(coord[0]) , (int)Math.floor(coord[1]), (int)Math.ceil(coord[2])); //V101
            Vxyz[6] = volume.getVoxel((int)Math.floor(coord[0]) , (int)Math.ceil(coord[1]), (int)Math.ceil(coord[2])); //V011
            Vxyz[7] = volume.getVoxel((int)Math.ceil(coord[0]) , (int)Math.ceil(coord[1]), (int)Math.ceil(coord[2])); //V111
            //compute alpha beta and gamma
            alpha = coord[0] - (int)Math.floor(coord[0]);
            beta = coord[1] - (int)Math.floor(coord[1]);
            gamma = coord[2] - (int)Math.floor(coord[2]);
            
            return (short)((1-alpha)*(1-beta)*(1-gamma)*Vxyz[0] //Value at position (x,y,z) within the cube
                            + alpha*(1-beta)*(1-gamma)*Vxyz[1]
                            + (1-alpha)*beta*(1-gamma)*Vxyz[2] 
                            + alpha*beta*(1 - gamma)*Vxyz[3]
                            + (1-alpha)*(1-beta)*gamma*Vxyz[4]
                            + alpha*(1-beta)*gamma*Vxyz[5]
                            + (1-alpha)*beta*gamma*Vxyz[6] 
                            + alpha*beta*gamma*Vxyz[7]
                           );
        }else{//Case when outside of dimensions
            return 0;
        }
    }

    short getVoxel(double[] coord) {
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return 0;
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

       return volume.getVoxel(x, y, z);
    }
    
    TFColor getPhongColor(double[] coord, double[] viewVec, TFColor testColor){
        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        if (x + volume.getDimX()*(y + volume.getDimY() * z) < data.length || x + volume.getDimX()*(y + volume.getDimY() * z) >= 0){
            if (data[x + volume.getDimX()*(y + volume.getDimY() * z)] != null){
                return data[x + volume.getDimX()*(y + volume.getDimY() * z)];
            } else {
                setColor(coord, phongShading(testColor , coord, viewVec, 0.1, 0.7, 0.2, 10));
                return data[x + volume.getDimX()*(y + volume.getDimY() * z)];
            }
        } else {
            return new TFColor(0,0,0,0);
        }
    }
    
    public void setColor(double[] coord, TFColor value) {
        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);
        
        data[x + volume.getDimX()*(y + volume.getDimY()*z)] = value;
    }
    
    VoxelGradient getGradient(double[] coord) {
        if (coord[0] < 0 || coord[0] > volume.getDimX() || coord[1] < 0 || coord[1] > volume.getDimY()
                || coord[2] < 0 || coord[2] > volume.getDimZ()) {
            return new VoxelGradient();
        }

        int x = (int) Math.floor(coord[0]);
        int y = (int) Math.floor(coord[1]);
        int z = (int) Math.floor(coord[2]);

        return gradients.getGradient(x,y,z);
    }
    
    void MIP(double[] viewMatrix){
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0xFF00FF00);
            }
        }
        
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] originPoint = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        VectorMath.setVector(originPoint, viewMatrix[3], viewMatrix[7], viewMatrix[11]);
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        double[] volumeCenter = new double[3];
        int imageCenter = image.getWidth() / 2;
        double[] pixelCoord = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        
        for (int j = 0; j < image.getHeight(); j+=resolution()) {
            for (int i = 0; i < image.getWidth(); i+=resolution()) {
                int maxVal = 0;
                
                
                for (double k = -getMiddleMaxDimension(); k < getMiddleMaxDimension();k+= sampleDistance){
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + 1 * k * viewVec[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + 1 * k * viewVec[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + 1 * k * viewVec[2];
                    
                    int value = trilinearInterpolation(pixelCoord);
                    if (value > maxVal){
                        maxVal = value;
                    }
                }
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = maxVal/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = maxVal > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                if(interactiveMode){
                    image.setRGB(i+1, j, pixelColor);
                    image.setRGB(i, j+1, pixelColor);
                    image.setRGB(i+1, j+1, pixelColor);
                }
            }     
        }   
    }
    
    void calculateVoxelColor(){
        float dx = 0;
        float dy = 0;
        float dz = 0;
        for (double i = 0; i < volume.getDimX(); i++){     
            for (double j = 0; j < volume.getDimY(); j++){   
                for (double k = 0; k < volume.getDimZ(); k++){
                    double[] coord = new double[3];
                    coord[0] = i;
                    coord[1] = j;
                    coord[2] = k;
                    
                    int intensity = trilinearInterpolation(coord);
                    TFColor newCol = tFunc.getColor(intensity);
                    setColor(coord, newCol);
                }
            }
        }
    }
    
    void composite(double[] viewMatrix){
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0xFF00FF00);
            }
        }
        
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] originPoint = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        VectorMath.setVector(originPoint, viewMatrix[3], viewMatrix[7], viewMatrix[11]);
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        double[] volumeCenter = new double[3];
        int imageCenter = image.getWidth() / 2;
        double[] pixelCoord = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        ArrayList<TFColor> compositeColors = new ArrayList<TFColor>();
        for (int j = 0; j < image.getHeight(); j+= resolution()) {
            for (int i = 0; i < image.getWidth(); i+=resolution()) {
                compositeColors.clear();
               
                for (double k = -getMiddleMaxDimension(); k < getMiddleMaxDimension(); k+= sampleDistance){
                    //System.out.println(k * viewVec[2]);
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + 1 * k * viewVec[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + 1 * k * viewVec[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + 1 * k * viewVec[2];
                    
                    int value = trilinearInterpolation(pixelCoord);
                    TFColor newColor = new TFColor();
                    
                    newColor.a = tFunc.getColor(value).a;
                    newColor.r = tFunc.getColor(value).r;
                    newColor.g = tFunc.getColor(value).g;
                    newColor.b = tFunc.getColor(value).b;
                    //if(shading){   
                    //    newColor = phongShading(newColor , pixelCoord, viewVec, 0.1, 0.7, 0.2, 10);
                    //}

                    compositeColors.add(newColor);
                }
                // Map the intensity to a grey value by linear scaling
                
                /* double r = compositeColors.get(0).r;
                double g = compositeColors.get(0).g;
                double b = compositeColors.get(0).b;
                double a = compositeColors.get(0).a;
                for (int q = 1; q < compositeColors.size(); q++){
                    double invOpacity = 1-compositeColors.get(q).a;
                    r *= invOpacity;
                    g *= invOpacity;
                    b *= invOpacity;
                    a *= invOpacity;
                    r += compositeColors.get(q).r;
                    g += compositeColors.get(q).g;
                    b += compositeColors.get(q).b;
                    a += compositeColors.get(q).a;
                }
                voxelColor.a = a;
                voxelColor.r = r;
                voxelColor.g = g;
                voxelColor.b = b;
                */
                
                double ru = 0;
                double gu = 0;
                double bu = 0;
                double au = 0;
                for (int q = 0; q < compositeColors.size(); q++){
                    
                     double aU = compositeColors.get(q).a;
                    if (aU > 0) {
                        double rU = compositeColors.get(q).r;
                        double gU = compositeColors.get(q).g;
                        double bU = compositeColors.get(q).b;
                        ru += aU * rU * (1-au);
                        gu += aU * gU * (1-au);
                        bu += aU * bU * (1-au);
                        au += aU * (1-au);
                    }
                }
                voxelColor.a = au;
                voxelColor.r = ru;
                voxelColor.g = gu;
                voxelColor.b = bu;
                
                //voxelColor = tFunc.getColor(maxVal);
                //voxelColor.r = maxVal/max;
                //voxelColor.g = voxelColor.r;
                //voxelColor.b = voxelColor.r;
                //voxelColor.a = maxVal > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                if(interactiveMode){
                    image.setRGB(i+1, j, pixelColor);
                    image.setRGB(i, j+1, pixelColor);
                    image.setRGB(i+1, j+1, pixelColor);
                }
            }     
        }   
    }
     
    private TFColor phongShading(TFColor original, double[] coord, double[] origin, double ambient, double diff, double spec, double alpha){
        VoxelGradient gradient = getGradient(coord);
        TFColor newColor = new TFColor();
        TFColor lightSource = new TFColor(1,1,1,1);

        double[] L = new double[3]; 
        double[] H = new double[3];
        double[] N = new double[3];     
        L[0] = origin[0]-coord[0];
        L[1] = origin[1]-coord[1];
        L[2] = origin[2]-coord[2];
        double mag = Math.sqrt(L[0]*L[0]+L[1]*L[1]+L[2]*L[2]);
        if (mag == 0.0 ){
            mag = 0.0000001;
        }
        L[0] = L[0] / mag;
        L[1] = L[1] / mag;
        L[2] = L[2] / mag;
        H[0] = 2*L[0];
        H[1] = 2*L[1];
        H[2] = 2*L[2];
        mag = Math.sqrt(H[0]*H[0]+H[1]*H[1]+H[2]*H[2]);
        if (mag == 0.0){
            mag = 0.0000001;
        }
        H[0] = H[0] / mag;
        H[1] = H[1] / mag;
        H[2] = H[2] / mag;
        mag = (double) gradient.mag;
        if (mag == 0.0){
            mag = 0.0000001;
        }
        
        N[0] = -((double) gradient.x / mag);
        N[1] = -((double) gradient.y / mag);
        N[2] = -((double) gradient.z / mag);
        if(N[0] != 0 || N[1] != 0 || N[2] != 0){
            //System.out.println("[" +N[0] + "," + N[1] + "," + N[2] + "]" + " " + "[" +L[0] + "," + L[1] + "," + L[2] + "]");
        }
        double dotLN = VectorMath.dotproduct(L, N);
        if(N[0] != 0 || N[1] != 0 || N[2] != 0){
            //System.out.println(dotLN);
        }
        double dotNH = VectorMath.dotproduct(N, H);
        double dotNHa = Math.pow(dotNH, alpha);
        newColor.a = original.a;
        newColor.r = lightSource.r * ambient + original.r * diff * Math.max(0.0, dotLN) + original.r * spec * Math.max(0.0, dotNHa);
        newColor.g = lightSource.g * ambient + original.g * diff * Math.max(0.0, dotLN) + original.g * spec * Math.max(0.0, dotNHa);
        newColor.b = lightSource.b * ambient + original.b * diff * Math.max(0.0, dotLN) + original.b * spec * Math.max(0.0, dotNHa);
        if (original.r * diff * Math.max(0.0, dotLN) > 0 ||  spec * Math.max(0.0, dotNHa) > 0){
            //System.out.println(lightSource.g * ambient + " " + original.g * diff * Math.max(0.0, dotLN) + " " + original.g * spec * Math.max(0.0, dotNHa));
        }
        newColor.r *= 1.5;
        newColor.g *= 1.5;
        newColor.b *= 1.5;
        return newColor;
     }
    
    void transfer2D(double[] viewMatrix){        
        double gradient_top = tfEditor2D.triangleWidget.gradUpperBound;
        double gradient_down = tfEditor2D.triangleWidget.gradLowerBound;
        for (int i = 0; i < data.length; i++) {
            data[i] = null;
        }
        System.out.println("gradients: " + gradient_top + gradient_down);
        debug = true;

        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0xFF00FF00);
            }
        }
        
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] originPoint = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        VectorMath.setVector(originPoint, viewMatrix[3], viewMatrix[7], viewMatrix[11]);
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();
        TFColor baseColor = tfEditor2D.triangleWidget.color;
        double[] volumeCenter = new double[3];
        int imageCenter = image.getWidth() / 2;
        double[] pixelCoord = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);
        
        ArrayList<TFColor> compositeColors = new ArrayList<TFColor>();
        for (int j = 0; j < image.getHeight(); j+= resolution()) {
            for (int i = 0; i < image.getWidth(); i+=resolution()) {
                compositeColors.clear();
                
                
                for (double k = -getMiddleMaxDimension(); k < getMiddleMaxDimension(); k+= sampleDistance){
                    //System.out.println(k * viewVec[2]);
                    pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0] + 1 * k * viewVec[0];
                    pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1] + 1 * k * viewVec[1];
                    pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2] + 1 * k * viewVec[2];
                    //System.out.println(pixelCoord[0] + " " + pixelCoord[1] + " " + pixelCoord[2]);
                    int value = trilinearInterpolation(pixelCoord);
                    TFColor newColor = new TFColor();
                    newColor.a = baseColor.a;
                    newColor.r = baseColor.r;
                    newColor.g = baseColor.g;
                    newColor.b = baseColor.b;
                    VoxelGradient voxGra = getGradient(pixelCoord);
                    int fv = tfEditor2D.triangleWidget.baseIntensity;
                    double r = tfEditor2D.triangleWidget.radius;
                    float magnitude = voxGra.mag;
                    
                    if(magnitude > gradient_top || magnitude < gradient_down){//Widget modification jump next iteration of loop
                        continue;
                    }
                    
                    
                    if (magnitude == 0.0f && value == fv){
                        newColor.a = baseColor.a;
                    } else if (magnitude > 0.0f &&  fv >= value - r * magnitude  && fv <= value + r * magnitude){
                        
                        newColor.a = baseColor.a * (1 - (1/r) * Math.abs((fv - value)/magnitude));
                    } else {
                        
                        newColor.a = 0;
                    }
                    if(shading){
                        TFColor testColor = new TFColor();
                        testColor.a = newColor.a;
                        testColor.r = newColor.r;
                        testColor.g = newColor.g;
                        testColor.b = newColor.b;
                        newColor = getPhongColor(pixelCoord, viewVec, testColor);
                    }
                    compositeColors.add(newColor);
                }
                
                double ru = 0;
                double gu = 0;
                double bu = 0;
                double au = 0;
                for (int q = 0; q < compositeColors.size(); q++){
                     double aU = compositeColors.get(q).a;
                    if (aU > 0) {
                        double rU = compositeColors.get(q).r;
                        double gU = compositeColors.get(q).g;
                        double bU = compositeColors.get(q).b;
                        ru += aU * rU * (1-au);
                        gu += aU * gU * (1-au);
                        bu += aU * bU * (1-au);
                        au += aU * (1-au);
                    }
                }
                voxelColor.a = au;
                voxelColor.r = ru;
                voxelColor.g = gu;
                voxelColor.b = bu;
                
                //voxelColor = tFunc.getColor(maxVal);
                //voxelColor.r = maxVal/max;
                //voxelColor.g = voxelColor.r;
                //voxelColor.b = voxelColor.r;
                //voxelColor.a = maxVal > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
                if(interactiveMode){
                    image.setRGB(i+1, j, pixelColor);
                    image.setRGB(i, j+1, pixelColor);
                    image.setRGB(i+1, j+1, pixelColor);
                }
            }     
        }   
    }


    void slicer(double[] viewMatrix) {

        // clear image
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                image.setRGB(i, j, 0);
            }
        }

        // vector uVec and vVec define a plane through the origin, 
        // perpendicular to the view vector viewVec
        double[] viewVec = new double[3];
        double[] uVec = new double[3];
        double[] vVec = new double[3];
        double[] originPoint = new double[3];
        VectorMath.setVector(viewVec, viewMatrix[2], viewMatrix[6], viewMatrix[10]);
        VectorMath.setVector(uVec, viewMatrix[0], viewMatrix[4], viewMatrix[8]);
        VectorMath.setVector(vVec, viewMatrix[1], viewMatrix[5], viewMatrix[9]);
        VectorMath.setVector(originPoint, viewMatrix[3], viewMatrix[7], viewMatrix[11]);
        //System.out.println(viewMatrix[3] + "," +  viewMatrix[7] + "," +  viewMatrix[11]);
        // image is square
        int imageCenter = image.getWidth() / 2;

        double[] pixelCoord = new double[3];
        double[] volumeCenter = new double[3];
        VectorMath.setVector(volumeCenter, volume.getDimX() / 2, volume.getDimY() / 2, volume.getDimZ() / 2);

        // sample on a plane through the origin of the volume data
        double max = volume.getMaximum();
        TFColor voxelColor = new TFColor();

        
        for (int j = 0; j < image.getHeight(); j++) {
            for (int i = 0; i < image.getWidth(); i++) {
                pixelCoord[0] = uVec[0] * (i - imageCenter) + vVec[0] * (j - imageCenter)
                        + volumeCenter[0];
                pixelCoord[1] = uVec[1] * (i - imageCenter) + vVec[1] * (j - imageCenter)
                        + volumeCenter[1];
                pixelCoord[2] = uVec[2] * (i - imageCenter) + vVec[2] * (j - imageCenter)
                        + volumeCenter[2];

                int val = trilinearInterpolation(pixelCoord);
                
                // Map the intensity to a grey value by linear scaling
                voxelColor.r = val/max;
                voxelColor.g = voxelColor.r;
                voxelColor.b = voxelColor.r;
                voxelColor.a = val > 0 ? 1.0 : 0.0;  // this makes intensity 0 completely transparent and the rest opaque
                // Alternatively, apply the transfer function to obtain a color
                // voxelColor = tFunc.getColor(val);
                
                
                // BufferedImage expects a pixel color packed as ARGB in an int
                int c_alpha = voxelColor.a <= 1.0 ? (int) Math.floor(voxelColor.a * 255) : 255;
                int c_red = voxelColor.r <= 1.0 ? (int) Math.floor(voxelColor.r * 255) : 255;
                int c_green = voxelColor.g <= 1.0 ? (int) Math.floor(voxelColor.g * 255) : 255;
                int c_blue = voxelColor.b <= 1.0 ? (int) Math.floor(voxelColor.b * 255) : 255;
                int pixelColor = (c_alpha << 24) | (c_red << 16) | (c_green << 8) | c_blue;
                image.setRGB(i, j, pixelColor);
            }
        }

    }


    private void drawBoundingBox(GL2 gl) {
        gl.glPushAttrib(GL2.GL_CURRENT_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glColor4d(1.0, 1.0, 1.0, 1.0);
        gl.glLineWidth(1.5f);
        gl.glEnable(GL.GL_LINE_SMOOTH);
        gl.glHint(GL.GL_LINE_SMOOTH_HINT, GL.GL_NICEST);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glBegin(GL.GL_LINE_LOOP);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glVertex3d(-volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, volume.getDimZ() / 2.0);
        gl.glVertex3d(volume.getDimX() / 2.0, -volume.getDimY() / 2.0, -volume.getDimZ() / 2.0);
        gl.glEnd();

        gl.glDisable(GL.GL_LINE_SMOOTH);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL2.GL_LIGHTING);
        gl.glPopAttrib();

    }

    @Override
    public void visualize(GL2 gl) {


        if (volume == null) {
            return;
        }

        drawBoundingBox(gl);

        gl.glGetDoublev(GL2.GL_MODELVIEW_MATRIX, viewMatrix, 0);

        long startTime = System.currentTimeMillis();
        //calculateVoxelColor(); //Calculate the colors for every voxel
        switch (currentMode){
            case 0: slicer(viewMatrix);
                    break;
            case 1: MIP(viewMatrix);
                    break;
            case 2: composite(viewMatrix);
                    break;
            case 3: transfer2D(viewMatrix);
                    break;
            default:slicer(viewMatrix);
                    break;
        }
        long endTime = System.currentTimeMillis();
        double runningTime = (endTime - startTime);
        panel.setSpeedLabel(Double.toString(runningTime));

        Texture texture = AWTTextureIO.newTexture(gl.getGLProfile(), image, false);

        gl.glPushAttrib(GL2.GL_LIGHTING_BIT);
        gl.glDisable(GL2.GL_LIGHTING);
        gl.glEnable(GL.GL_BLEND);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);

        // draw rendered image as a billboard texture
        texture.enable(gl);
        texture.bind(gl);
        double halfWidth = image.getWidth() / 2.0;
        gl.glPushMatrix();
        gl.glLoadIdentity();
        gl.glBegin(GL2.GL_QUADS);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        gl.glTexCoord2d(0.0, 0.0);
        gl.glVertex3d(-halfWidth, -halfWidth, 0.0);
        gl.glTexCoord2d(0.0, 1.0);
        gl.glVertex3d(-halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 1.0);
        gl.glVertex3d(halfWidth, halfWidth, 0.0);
        gl.glTexCoord2d(1.0, 0.0);
        gl.glVertex3d(halfWidth, -halfWidth, 0.0);
        gl.glEnd();
        texture.disable(gl);
        texture.destroy(gl);
        gl.glPopMatrix();

        gl.glPopAttrib();


        if (gl.glGetError() > 0) {
            System.out.println("some OpenGL error: " + gl.glGetError());
        }

    }
    private BufferedImage image;
    private double[] viewMatrix = new double[4 * 4];

    @Override
    public void changed() {
        for (int i=0; i < listeners.size(); i++) {
            listeners.get(i).changed();
        }
    }
}
