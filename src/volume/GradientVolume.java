/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package volume;

/**
 *
 * @author michel
 */
public class GradientVolume {

    public GradientVolume(Volume vol) {
        volume = vol;
        dimX = vol.getDimX();
        dimY = vol.getDimY();
        dimZ = vol.getDimZ();
        data = new VoxelGradient[dimX * dimY * dimZ];
        compute();
        maxmag = -1.0;
    }

    public VoxelGradient getGradient(int x, int y, int z) {
        return data[x + dimX * (y + dimY * z)];
    }

    
    public void setGradient(int x, int y, int z, VoxelGradient value) {
        data[x + dimX * (y + dimY * z)] = value;
    }

    public void setVoxel(int i, VoxelGradient value) {
        data[i] = value;
    }

    public VoxelGradient getVoxel(int i) {
        return data[i];
    }

    public int getDimX() {
        return dimX;
    }

    public int getDimY() {
        return dimY;
    }

    public int getDimZ() {
        return dimZ;
    }

    private void compute() {

        // this just initializes all gradients to the vector (0,0,0)
        for (int i=0; i<data.length; i++) {
            //data[i] = zero;
        }   
        float dx = 0;
        float dy = 0;
        float dz = 0;
        for (int i = 0; i < dimX; i++){     
            for (int j = 0; j < dimY; j++){   
                for (int k = 0; k < dimZ; k++){
                    if(i == 0){
                        dx = 0.5f * volume.getVoxel(i+1, j, k);
                    } else if (i == dimX-1){
                        dx = 0.5f * - volume.getVoxel(i-1, j, k);
                    } else {
                        dx = 0.5f * (volume.getVoxel(i+1, j, k) - volume.getVoxel(i-1, j, k));
                    }
                    if(j == 0){
                        dy = 0.5f * volume.getVoxel(i, j+1, k);
                    } else if (j == dimY-1){
                        dy = 0.5f * - volume.getVoxel(i, j-1, k);
                    } else {
                        dy = 0.5f * (volume.getVoxel(i, j+1, k) - volume.getVoxel(i, j-1, k));
                    }
                    if(k == 0){
                        dz = 0.5f * volume.getVoxel(i, j, k+1);
                    } else if (k == dimZ-1){
                        dz = 0.5f * - volume.getVoxel(i, j, k-1);
                    } else {
                        dz = 0.5f * (volume.getVoxel(i, j, k+1) - volume.getVoxel(i, j, k-1));
                    }
                    VoxelGradient calcGradient = new VoxelGradient(dx,dy,dz);
                    setGradient(i,j,k,calcGradient);
                }
            }
        }  
            
        
                
    }
    
    public double getMaxGradientMagnitude() {
        if (maxmag >= 0) {
            return maxmag;
        } else {
            double magnitude = data[0].mag;
            for (int i=0; i<data.length; i++) {
                magnitude = data[i].mag > magnitude ? data[i].mag : magnitude;
            }   
            maxmag = magnitude;
            return magnitude;
        }
    }
    
    private int dimX, dimY, dimZ;
    private VoxelGradient zero = new VoxelGradient();
    VoxelGradient[] data;
    Volume volume;
    double maxmag;
}
