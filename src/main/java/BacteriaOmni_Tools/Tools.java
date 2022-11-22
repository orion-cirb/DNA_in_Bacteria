package BacteriaOmni_Tools;

import BacteriaOmni_Tools.Cellpose.CellposeTaskSettings;
import BacteriaOmni_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import edu.mines.jtk.util.AtomicFloat;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import fiji.util.gui.GenericDialogPlus;
import ij.gui.WaitForUserDialog;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import ij.process.ImageProcessor;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.ResourceBundle;
import javax.swing.ImageIcon;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.plugins.util.ImageProcessorReader;
import mcib3d.geom.Voxel3D;
import mcib3d.geom2.Object3DInt;
import mcib3d.geom2.Objects3DIntPopulation;
import mcib3d.geom2.Objects3DIntPopulationComputation;
import mcib3d.geom2.VoxelInt;
import mcib3d.geom2.measurements.MeasureCentroid;
import mcib3d.geom2.measurements.MeasureFeret;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationClosestDistance;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationColocalisation;
import mcib3d.geom2.measurementsPopulation.MeasurePopulationDistance;
import mcib3d.geom2.measurementsPopulation.PairObjects3DInt;
import mcib3d.image3d.ImageHandler;
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
import net.haesleinhuepf.clij.clearcl.ClearCLBuffer;
import net.haesleinhuepf.clij2.CLIJ2;



/**
 * @author Orion-CIRB
 */
public class Tools {
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));
    public boolean canceled = false;
    public Calibration cal = new Calibration();
    private double pixelSurf = 0;
    
     // Omnipose
    private String omniposeEnvDirPath = "/opt/miniconda3/envs/omnipose";
    private String omniposeModelsPath = System.getProperty("user.home")+"/.cellpose/models/";
    private String omniposeModel = "bact_phase_omnitorch_0";
    private int omniposeDiameter = 10;
    private int omniposeMaskThreshold = 0;
    private double omniposeFlowThreshold = 0.4;
    private boolean useGpu = true;
    
    // Bacteria
    private double minBactSurface = 0.3;
    private double maxBactSurface = 5;
    // Foci
    private double minFociSurface = 0.005;
    private double maxFociSurface = 0.5;
    private String fociTh = "Moments";
    
    // DOG parameters
    private double minFociDOG = 1;
    private double maxFociDOG = 2;
    
    private final CLIJ2 clij2 = CLIJ2.getInstance();
    
    /**
     * Display a message in the ImageJ console and status bar
     */
    public void print(String log) {
        System.out.println(log);
        IJ.showStatus(log);
    }
    
    
    /**
     * Check that needed modules are installed
     */
    public boolean checkInstalledModules() {
        ClassLoader loader = IJ.getClassLoader();
        try {
            loader.loadClass("mcib3d.geom.Object3D");
        } catch (ClassNotFoundException e) {
            IJ.showMessage("Error", "3D ImageJ Suite not installed, please install from update site");
            return false;
        }
        return true;
    }
    
    
     /**
     * Difference of Gaussians 
     * Using CLIJ2
     * @param imgCL
     * @param size1
     * @param size2
     * @return imgGauss
     */ 
    public ImagePlus DOG(ImagePlus img, double size1, double size2) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLDOG = clij2.create(imgCL);
        clij2.differenceOfGaussian3D(imgCL, imgCLDOG, size1, size1, size1, size2, size2, size2);
        clij2.release(imgCL);
        ImagePlus imgDOG = clij2.pull(imgCLDOG);
        clij2.release(imgCLDOG);
        return(imgDOG);
    }
    
    
    /**
     * Threshold 
     * USING CLIJ2
     * @param imgCL
     * @param thMed
     */
    public ImagePlus threshold(ImagePlus img, String thMed) {
        ClearCLBuffer imgCL = clij2.push(img);
        ClearCLBuffer imgCLBin = clij2.create(imgCL);
        clij2.automaticThreshold(imgCL, imgCLBin, thMed);
        clij2.release(imgCL);
        ImagePlus imgBin = clij2.pull(imgCLBin);
        clij2.release(imgCLBin);
        return(imgBin);
    }
    
    /**
     * Flush and close an image
     */
    public void flush_close(ImagePlus img) {
        img.flush();
        img.close();
    }
    
    
    /**
     * Find images extension
     */
    public String findImageType(File imagesFolder) {
        String ext = "";
        File[] files = imagesFolder.listFiles();
        for (File file: files) {
            if(file.isFile()) {
                String fileExt = FilenameUtils.getExtension(file.getName());
                switch (fileExt) {
                   case "nd" :
                       ext = fileExt;
                       break;
                    case "czi" :
                       ext = fileExt;
                       break;
                    case "lif"  :
                        ext = fileExt;
                        break;
                    case "ics2" :
                        ext = fileExt;
                        break;
                    case "tif" :
                        ext = fileExt;
                        break;
                    case "tiff" :
                        ext = fileExt;
                        break;
                }
            } else if (file.isDirectory() && !file.getName().equals("Results")) {
                ext = findImageType(file);
                if (! ext.equals(""))
                    break;
            }
        }
        return(ext);
    }
     
    
    /**
     * Find images in folder
     */
    public void findImages(String imagesFolder, String imageExt, ArrayList<String> imageFiles) {
        System.out.println(imagesFolder);
        File inDir = new File(imagesFolder);
        File[] files = inDir.listFiles();
        System.out.println(files);
        
        for (File file: files) {
            if(file.isFile()) {
                String fileExt = FilenameUtils.getExtension(file.getName());
                if (fileExt.equals(imageExt) && !file.getName().startsWith("."))
                    imageFiles.add(file.getAbsolutePath());
            } else if (file.isDirectory() && !file.getName().equals("Results")) {
                findImages(file.getAbsolutePath(), imageExt, imageFiles);
            }
        }
        Collections.sort(imageFiles);
    }
     /**
     * Find channels name
     * @param imageName
     * @param meta
     * @param reader
     * @return 
     * @throws loci.common.services.DependencyException
     * @throws loci.common.services.ServiceException
     * @throws loci.formats.FormatException
     * @throws java.io.IOException
     */
    public String[] findChannels (String imageName, IMetadata meta, ImageProcessorReader reader) throws DependencyException, ServiceException, FormatException, IOException {
        int chs = reader.getSizeC();
        String[] channels = new String[chs];
        String imageExt =  FilenameUtils.getExtension(imageName);
        switch (imageExt) {
            case "nd" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "nd2" :
                for (int n = 0; n < chs; n++) 
                {
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                }
                break;
            case "lif" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null || meta.getChannelName(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelName(0, n);
                break;
            case "czi" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelFluor(0, n);
                break;
            case "ics2" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
                break;    
            default :
                for (int n = 0; n < chs; n++)
                    channels[n] = Integer.toString(n);
        }
        return(channels);         
    }
    
    /**
     * Generate dialog box
     */
    public String[] dialog(String[] channels) {
        String[] channelsName = {"Bacteria : ", "Foci1 : ", "Foci2 : "};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String ch : channels) {
            gd.addChoice(channelsName[index], channels, ch);
            index++;
        }
        gd.addMessage("Bacteria detection", Font.getFont("Monospace"), Color.blue);
        if (IJ.isWindows()) {
            omniposeEnvDirPath = System.getProperty("user.home")+"\\miniconda3\\envs\\omnipose";
            omniposeModelsPath = System.getProperty("user.home")+"\\.cellpose\\models\\";
        }
        gd.addDirectoryField("Omnipose environment directory: ", omniposeEnvDirPath);
        gd.addDirectoryField("Omnipose models path: ", omniposeModelsPath); 
        gd.addMessage("Object size threshold ", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min bacterium surface (µm2): ", minBactSurface);
        gd.addNumericField("Max bacterium surface (µm2): ", maxBactSurface);
        gd.addNumericField("Min Foci     surface (µm2): ", minFociSurface);
        gd.addNumericField("Max Foci     surface (µm2): ", maxFociSurface);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.showDialog();
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();
        omniposeEnvDirPath = gd.getNextString();
        omniposeModelsPath = gd.getNextString();
        minBactSurface = (float) gd.getNextNumber();
        maxBactSurface = (float) gd.getNextNumber();
        minFociSurface = (float) gd.getNextNumber();
        maxFociSurface = (float) gd.getNextNumber();
        cal.pixelWidth = gd.getNextNumber();
        pixelSurf = cal.pixelWidth*cal.pixelWidth;
        if(gd.wasCanceled())
           canceled = true;
        return(ch);
    }
    
    public Objects3DIntPopulation getPopFromImage(ImagePlus img) {
        // label binary images first
        ImageLabeller labeller = new ImageLabeller();
        ImageInt labels = labeller.getLabels(ImageHandler.wrap(img));
        Objects3DIntPopulation pop = new Objects3DIntPopulation(labels);
        return pop;
    }
    
    
    
    /**
     * Do Z projection
     * @param img
     * @param param
     * @return 
     */
    public ImagePlus doZProjection(ImagePlus img, int param) {
        ZProjector zproject = new ZProjector();
        zproject.setMethod(param);
        zproject.setStartSlice(1);
        zproject.setStopSlice(img.getNSlices());
        zproject.setImage(img);
        zproject.doProjection();
       return(zproject.getProjection());
    }
    
    /**
     * Find image calibration
     * @param meta
     * @return 
     */
    public void findImageCalib(IMetadata meta) {
        cal.pixelWidth = meta.getPixelsPhysicalSizeX(0).value().doubleValue();
        cal.pixelHeight = cal.pixelWidth;
        if (meta.getPixelsPhysicalSizeZ(0) != null)
            cal.pixelDepth = meta.getPixelsPhysicalSizeZ(0).value().doubleValue();
        else
            cal.pixelDepth = 1;
        cal.setUnit("microns");
        System.out.println("XY calibration = " + cal.pixelWidth + ", Z calibration = " + cal.pixelDepth);
    }
    
    /**
     * Remove object with size < min and size > max
     * @param pop
     * @param min
     * @param max
     */
    public void popFilterSize(Objects3DIntPopulation pop, double min, double max) {
        pop.getObjects3DInt().removeIf(p -> (new MeasureVolume(p).getVolumeUnit() < min) || (new MeasureVolume(p).getVolumeUnit() > max));
        pop.resetLabels();
    }
   
    /**
    * Detect bacteria with Omnipose
    */
    public Objects3DIntPopulation omniposeDetection(ImagePlus imgBact){

        ImagePlus imgIn = new Duplicator().run(imgBact);
        imgIn.setCalibration(cal);
        // Set Omnipose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(omniposeModelsPath+omniposeModel, 1, omniposeDiameter, omniposeEnvDirPath);
        settings.setVersion("0.7");
        settings.setCluster(true);
        settings.setOmni(true);
        settings.useMxNet(false);
        settings.setCellProbTh(omniposeMaskThreshold);
        settings.setFlowTh(omniposeFlowThreshold);
        settings.useGpu(useGpu);
        
        // Run Omnipose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        ImagePlus imgOut = cellpose.run();
        imgOut.setCalibration(cal);
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        Objects3DIntPopulation popFilter = new Objects3DIntPopulationComputation(pop).getExcludeBorders(ImageHandler.wrap(imgOut), false);
        popFilterSize(popFilter, minBactSurface, maxBactSurface);
        
        // Close images
        flush_close(imgIn);
        flush_close(imgOut);
        return(popFilter);
    }
    
    
    /**
     * Find foci with DOG
     * @param img
     * @return 
     */
    
    public Objects3DIntPopulation findFoci(ImagePlus img) {
        ImagePlus imgDog = DOG(img, minFociDOG, maxFociDOG);
        ImagePlus imgBin = threshold(imgDog, fociTh);
        imgBin.setCalibration(cal);
        Objects3DIntPopulation fociPop = getPopFromImage(imgBin);
        popFilterSize(fociPop, minFociSurface, maxFociSurface);
        return(fociPop);
        
    }
    
    
     /**
     * Find coloc between foci and bacteria
     * set label of colocalized bacteria in foci object
     * @param bactPop
     * @param fociPop
     */
    public void fociBactLink(Objects3DIntPopulation bactPop, Objects3DIntPopulation fociPop) {
        if (bactPop.getNbObjects() != 0 && fociPop.getNbObjects() != 0) {
            for (Object3DInt bact : bactPop.getObjects3DInt()) {
                for (Object3DInt foci : fociPop.getObjects3DInt()) {
                    MeasureCentroid fociCenter = new MeasureCentroid(foci);
                    if (bact.contains(fociCenter.getCentroidRoundedAsVoxelInt())){
                       foci.setIdObject(bact.getLabel()); 
                    }
                }
            }
        }
        // remove foci not in bacteria
        fociPop.getObjects3DInt().removeIf(p -> p.getIdObject() == 0);
        fociPop.resetLabels();
    }
    
    /**
     * Compute min distance between foci to bacteria border
     */
    private ArrayList<Double> fociBactDistance(Objects3DIntPopulation fociPop, VoxelInt bactFeret1) {
        ArrayList<Double> minFociDist = new ArrayList<>();
        for (Object3DInt foci : fociPop.getObjects3DInt()) {
            Voxel3D fociCenter = new MeasureCentroid(foci).getCentroidAsVoxel();
            double dist = bactFeret1.distance(fociCenter)*cal.pixelWidth;
            minFociDist.add(dist);
        }
        return(minFociDist);
    }
    
    /**
     * 
    */
    private Objects3DIntPopulation findFociBact(float bactLabel, Objects3DIntPopulation fociPop) {
        Objects3DIntPopulation fociBactPop = new Objects3DIntPopulation();
        for (Object3DInt foci : fociPop.getObjects3DInt()) {
                if (foci.getIdObject() == bactLabel)
                    fociBactPop.addObject(foci);
        }
        fociBactPop.resetLabels();
        return(fociBactPop) ;           
    }
    
    /**
     * Compute min,max distances between foci1 and foci2
     */
    private double[] fociFociDistance(Objects3DIntPopulation foci1Pop, Objects3DIntPopulation foci2Pop) {
        MeasurePopulationClosestDistance allDist =  new MeasurePopulationClosestDistance​(foci1Pop, foci2Pop);
        allDist.setDistanceMax(5);
        List<PairObjects3DInt> allDistPairs = allDist.getAllPairs(true);
        double[] minMaxDist = {allDistPairs.get(0).getPairValue(), allDistPairs.get(allDistPairs.size()-1).getPairValue()};
        return(minMaxDist);
    }
    
    /**
     * Compute bacteria area
     * 
     **/
    private double bacteriaArea(Object3DInt bactObj ) {
        int pixelNb = bactObj.getObject3DPlanes().get(0).getVoxels().size();
        return(pixelNb*cal.pixelWidth*cal.pixelHeight);
    }   
    
    
    
    /**
     * Compute bacteria parameters and save them in file
     * @param bactPop
     * @param foci1Pop
     * @param foci2Pop
     * @param imgName
     * @param file
     * @throws java.io.IOException
     */
    public void saveResults(Objects3DIntPopulation bactPop, Objects3DIntPopulation foci1Pop, Objects3DIntPopulation foci2Pop, String imgName, BufferedWriter file) throws IOException {
        for (Object3DInt bact : bactPop.getObjects3DInt()) {
            float bactLabel = bact.getLabel();
            //double bactSurf = new MeasureVolume(bact).getValueMeasurement(MeasureVolume.VOLUME_UNIT);
            double bactSurf = bacteriaArea(bact);
            VoxelInt feret1Unit = new MeasureFeret(bact).getFeret1Unit();
            VoxelInt feret2Unit = new MeasureFeret(bact).getFeret2Unit();
            double bactLength = feret1Unit.distance(feret2Unit)*cal.pixelWidth;
            
            Objects3DIntPopulation foci1BactPop = findFociBact(bactLabel, foci1Pop);
            int foci1Nb = foci1BactPop.getNbObjects();
            ArrayList<Double> minFoci1Dist = (foci1Nb == 0) ? null : fociBactDistance(foci1BactPop, feret1Unit);
            Objects3DIntPopulation foci2BactPop = findFociBact(bactLabel, foci2Pop);  
            int foci2Nb = foci2BactPop.getNbObjects();
            ArrayList<Double> minFoci2Dist = (foci2Nb == 0) ? null : fociBactDistance(foci2BactPop, feret1Unit);
            int fociMax = Math.max(foci1Nb, foci2Nb);
            double[] foci1Foci2Distance = (foci1Nb == 0 || foci2Nb == 0) ? null : fociFociDistance(foci1BactPop, foci2BactPop);
            file.write(imgName+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+foci1Nb+"\t"+foci2Nb+"\t");
            if (fociMax != 0) {
                for (int i = 0; i < fociMax; i++) {
                    if (i != 0)
                        file.write("\t\t\t\t\t\t");
                    if (minFoci1Dist != null && i < minFoci1Dist.size())
                        file.write(minFoci1Dist.get(i)+"\t");
                    else
                        file.write("\t");
                    if (minFoci2Dist != null && i < minFoci2Dist.size())
                        file.write(minFoci2Dist.get(i)+"\t");
                    else
                        file.write("\t");
                    if (foci1Foci2Distance != null && i == 0)
                        file.write(foci1Foci2Distance[0]+"\t"+foci1Foci2Distance[1]+"\n");
                    else
                        file.write("\t\t\n");
                }
            }
            else
                file.write("\t\t\t\t\n");
            file.flush();
        }
    }
   
    
    // Save objects image
    public void drawResults(ImagePlus img, Objects3DIntPopulation bactPop, Objects3DIntPopulation foci1Pop, Objects3DIntPopulation foci2Pop,String imgName, String outDir) {
        ImageHandler imgBact = ImageHandler.wrap(img).createSameDimensions();
        bactPop.drawInImage(imgBact);
        IJ.run(imgBact.getImagePlus(), "glasbey on dark", "");
        imgBact.getImagePlus().setCalibration(cal);
        FileSaver ImgObjectsFile = new FileSaver(imgBact.getImagePlus());
        ImgObjectsFile.saveAsTiff(outDir+imgName+"_bacteria.tif");
        
        ImageHandler imgFoci1 = ImageHandler.wrap(img).createSameDimensions();
        foci1Pop.drawInImage(imgFoci1);
        ImageHandler imgFoci2 = ImageHandler.wrap(img).createSameDimensions();
        foci2Pop.drawInImage(imgFoci2);
        ImagePlus[] imgColors = {imgFoci1.getImagePlus(), imgFoci2.getImagePlus(), null, img};
        ImagePlus imgOut = new RGBStackMerge().mergeHyperstacks(imgColors, true);
        imgOut.setCalibration(cal);
        FileSaver ImgObjectsFile2 = new FileSaver(imgOut);
        ImgObjectsFile2.saveAsTiff(outDir + imgName + "_overlay.tif");
        flush_close(imgBact.getImagePlus());
        flush_close(imgOut);
    }
    
}
