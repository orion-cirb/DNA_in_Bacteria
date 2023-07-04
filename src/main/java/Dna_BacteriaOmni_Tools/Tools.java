package Dna_BacteriaOmni_Tools;

import DNA_BacteriaOmni_Tools.Cellpose.CellposeTaskSettings;
import DNA_BacteriaOmni_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import fiji.util.gui.GenericDialogPlus;
import ij.plugin.RGBStackMerge;
import ij.plugin.ZProjector;
import java.awt.Color;
import java.awt.Font;
import java.io.BufferedWriter;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
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
import mcib3d.geom2.measurements.MeasureIntensity;
import mcib3d.geom2.measurements.MeasureVolume;
import mcib3d.image3d.ImageHandler;
import org.apache.commons.io.FilenameUtils;


/**
 * @author Orion-CIRB
 */
public class Tools {
    private final ImageIcon icon = new ImageIcon(this.getClass().getResource("/Orion_icon.png"));

    public Calibration cal = new Calibration();
    private double pixelSurf = 0;
    String[] channelsName = {"Bacteria: ", "DNA: "};
    
     // Omnipose
    private String omniposeEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\omnipose\\" : 
            "/opt/miniconda3/envs/omnipose/";
    private String omniposeModelsPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\.cellpose\\models\\" :
            System.getProperty("user.home")+"/.cellpose/models/";
    public String omniposeBactModel = "bact_phase_omnitorch_0";
    public String omniposeDnaModel = "bact_fluor_omnitorch_0";
     
    private int omniposeDiameter = 18;
    private int omniposeMaskThreshold = 0;
    private double omniposeFlowThreshold = 0;
    private boolean useGpu = true;
    
    // Bacteria
    public double minBactSurface = 0.4;
    public double maxBactSurface = 20;
    
    // DNA
    public double minDnaSurface = 0.4;
    public double maxDnaSurface = 20;
    
    
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
        String[] files = imagesFolder.list();
        for (String name : files) {
            String fileExt = FilenameUtils.getExtension(name);
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
                case "ics" :
                    ext = fileExt;
                    break;
                case "ics2" :
                    ext = fileExt;
                    break;
                case "lsm" :
                    ext = fileExt;
                    break;
                case "tif" :
                    ext = fileExt;
                    break;
                case "tiff" :
                    ext = fileExt;
                    break;
            }
        }
        return(ext);
    }
     
    
    /**
     * Find images in folder
     */
    public ArrayList<String> findImages(String imagesFolder, String imageExt) {
        File inDir = new File(imagesFolder);
        String[] files = inDir.list();
        if (files == null) {
            System.out.println("No image found in " + imagesFolder);
            return null;
        }
        ArrayList<String> images = new ArrayList();
        for (String f : files) {
            // Find images with extension
            String fileExt = FilenameUtils.getExtension(f);
            if (fileExt.equals(imageExt) && !f.startsWith("."))
                images.add(imagesFolder + File.separator + f);
        }
        Collections.sort(images);
        return(images);
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
     * Find channels name
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
            case "ics" :
                for (int n = 0; n < chs; n++) 
                    if (meta.getChannelID(0, n) == null)
                        channels[n] = Integer.toString(n);
                    else 
                        channels[n] = meta.getChannelExcitationWavelength(0, n).value().toString();
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
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 160, 0);
        gd.addImage(icon);
        
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String ch : channelsName) {
            gd.addChoice(ch, channels, channels[index]);
            index++;
        }
        
        gd.addMessage("Bacteria and DNA detection", Font.getFont("Monospace"), Color.blue);
        gd.addDirectoryField("Omnipose environment directory: ", omniposeEnvDirPath);
        gd.addDirectoryField("Omnipose models path: ", omniposeModelsPath); 
        gd.addMessage("Object size threshold ", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min bacterium surface (µm2): ", minBactSurface);
        gd.addNumericField("Max bacterium surface (µm2): ", maxBactSurface);
        gd.addNumericField("Min DNA surface (µm2): ", minDnaSurface);
        gd.addNumericField("Max DNA surface (µm2): ", maxDnaSurface);
        
        gd.addMessage("Image calibration", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("XY calibration (µm):", cal.pixelWidth);
        gd.showDialog();
        
        String[] ch = new String[channelsName.length];
        for (int i = 0; i < channelsName.length; i++)
            ch[i] = gd.getNextChoice();
        if(gd.wasCanceled())
           ch = null;
                
        omniposeEnvDirPath = gd.getNextString();
        omniposeModelsPath = gd.getNextString();
        minBactSurface = (float) gd.getNextNumber();
        maxBactSurface = (float) gd.getNextNumber();
        minDnaSurface = (float) gd.getNextNumber();
        maxDnaSurface = (float) gd.getNextNumber();
        
        cal.pixelWidth = cal.pixelHeight = gd.getNextNumber();
        cal.pixelDepth = 1;
        pixelSurf = cal.pixelWidth*cal.pixelWidth;

        return(ch);
    }
    
    
    /**
     * Do Z projection
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
    * Detect bacteria with Omnipose
    */
    public Objects3DIntPopulation omniposeDetection(ImagePlus imgBact, String model, double min, double max, boolean excludeBorders){
        ImagePlus imgIn = new Duplicator().run(imgBact);
        imgIn.setCalibration(cal);
        
        // Set Omnipose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(omniposeModelsPath+model, 1, omniposeDiameter, omniposeEnvDirPath);
        settings.setVersion("0.7");
        settings.setOmni(true);
        settings.useMxNet(false);
        settings.setCluster(true);
        settings.setCellProbTh(omniposeMaskThreshold);
        settings.setFlowTh(omniposeFlowThreshold);
        settings.useGpu(useGpu);
        
        // Run Omnipose
        CellposeSegmentImgPlusAdvanced cellpose = new CellposeSegmentImgPlusAdvanced(settings, imgIn);
        ImagePlus imgOut = cellpose.run();
        imgOut.setCalibration(cal);
        
        Objects3DIntPopulation pop = new Objects3DIntPopulation(ImageHandler.wrap(imgOut));
        if (excludeBorders)
            pop = new Objects3DIntPopulationComputation(pop).getExcludeBorders(ImageHandler.wrap(imgOut), false);
        pop = new Objects3DIntPopulationComputation(pop).getFilterSize(min/pixelSurf, max/pixelSurf);
        pop.resetLabels();
        
        // Close images
        flush_close(imgIn);
        flush_close(imgOut);
        
        return(pop);
    }
    

    /**
     * Find DNA in bacteria
     * Set label of mother bacteria in DNA object
     */
    public void dnaBactLink(Objects3DIntPopulation bactPop, Objects3DIntPopulation dnaPop) {
        if (bactPop.getNbObjects() != 0 && dnaPop.getNbObjects() != 0) {
            for (Object3DInt bact : bactPop.getObjects3DInt()) {
                for (Object3DInt dna : dnaPop.getObjects3DInt()) {
                    MeasureCentroid dnaCenter = new MeasureCentroid(dna);
                    if (bact.contains(dnaCenter.getCentroidRoundedAsVoxelInt())){
                       dna.setIdObject(bact.getLabel()); 
                    }
                }
            }
        }
        // Remove DNA not in bacteria
        dnaPop.getObjects3DInt().removeIf(p -> p.getIdObject() == 0);
        dnaPop.resetLabels();
    }
    
   
    /**
     * Compute bacteria parameters and save them in file
     * @throws java.io.IOException
     */
    public void saveResults(Objects3DIntPopulation bactPop, Objects3DIntPopulation dnaPop, ImagePlus dnaImg, String imgName, int time, BufferedWriter file) throws IOException {
        for (Object3DInt bact : bactPop.getObjects3DInt()) {
            float bactLabel = bact.getLabel();
            double bactSurf = new MeasureVolume(bact).getValueMeasurement(MeasureVolume.VOLUME_UNIT);
            VoxelInt feret1Unit = new MeasureFeret(bact).getFeret1Unit();
            VoxelInt feret2Unit = new MeasureFeret(bact).getFeret2Unit();
            double bactLength = feret1Unit.distance(feret2Unit)*cal.pixelWidth;
            
            Objects3DIntPopulation dnaBactPop = findDnaInBact(bactLabel, dnaPop);
            int dnaNb = dnaBactPop.getNbObjects();
            if (dnaNb == 0) {
                file.write(imgName+"\t"+time+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+dnaNb+"\n");
            } else {
                for (Object3DInt dna: dnaBactPop.getObjects3DInt()) {
                    double dnaSurf = new MeasureVolume(dna).getValueMeasurement(MeasureVolume.VOLUME_UNIT);
                    double dnaInt = new MeasureIntensity(dna, ImageHandler.wrap(dnaImg)).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                    double dnaDist = dnaBactDistance(dna, bact);
                    file.write(imgName+"\t"+time+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+dnaNb+"\t"+dna.getLabel()+"\t"+dnaSurf+"\t"+dnaInt+"\t"+dnaDist+"\n");
                }
            }
            file.flush();
        }
    }
    
    
    /**
     * Get DNA in bacterium
     */
    private Objects3DIntPopulation findDnaInBact(float bactLabel, Objects3DIntPopulation dnaPop) {
        Objects3DIntPopulation dnaBactPop = new Objects3DIntPopulation();
        for (Object3DInt dna : dnaPop.getObjects3DInt()) {
                if (dna.getIdObject() == bactLabel)
                    dnaBactPop.addObject(dna);
        }
        dnaBactPop.resetLabels();
        return(dnaBactPop) ;           
    }
    
    
    /**
     * Compute distance between DNA center and bacterium center
     */
    private double dnaBactDistance(Object3DInt dna, Object3DInt bact) {
        Voxel3D bactCenter = new MeasureCentroid(bact).getCentroidAsVoxel();
        Voxel3D dnaCenter = new MeasureCentroid(dna).getCentroidAsVoxel();
        double dist = bactCenter.distance(dnaCenter)*cal.pixelWidth;
        return(dist);
    }
    
    
    /**
     * Save results in images
     */
    public void drawResults(ImagePlus img1, ImagePlus img2, Objects3DIntPopulation bactPop, Objects3DIntPopulation dnaPop, String imgName, String outDir) {
        ImageHandler imgBact = ImageHandler.wrap(img1).createSameDimensions();
        bactPop.drawInImage(imgBact);
        IJ.run(imgBact.getImagePlus(), "glasbey on dark", "");
        ImagePlus[] imgColors1 = {imgBact.getImagePlus(), null, null, img1};
        ImagePlus imgOut1 = new RGBStackMerge().mergeHyperstacks(imgColors1, false);
        imgOut1.setCalibration(cal);
        FileSaver ImgObjectsFile1 = new FileSaver(imgOut1);
        ImgObjectsFile1.saveAsTiff(outDir+imgName+"_bacteria.tif");      
        
        ImageHandler imgDna = ImageHandler.wrap(img2).createSameDimensions();
        dnaPop.drawInImage(imgDna);
        ImagePlus[] imgColors2 = {imgDna.getImagePlus(), null, null, img2};
        ImagePlus imgOut2 = new RGBStackMerge().mergeHyperstacks(imgColors2, false);
        imgOut2.setCalibration(cal);
        FileSaver ImgObjectsFile2 = new FileSaver(imgOut2);
        ImgObjectsFile2.saveAsTiff(outDir+imgName+"_DNA.tif");
        
        flush_close(imgBact.getImagePlus());
        flush_close(imgDna.getImagePlus());
        flush_close(imgOut1);
        flush_close(imgOut2);
    }
    
}
