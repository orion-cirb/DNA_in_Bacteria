package Adn_BacteriaOmni_Tools;

import Adn_BacteriaOmni_Tools.Cellpose.CellposeTaskSettings;
import Adn_BacteriaOmni_Tools.Cellpose.CellposeSegmentImgPlusAdvanced;
import ij.IJ;
import ij.ImagePlus;
import ij.io.FileSaver;
import ij.measure.Calibration;
import ij.plugin.Duplicator;
import fiji.util.gui.GenericDialogPlus;
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
import mcib3d.image3d.ImageInt;
import mcib3d.image3d.ImageLabeller;
import org.apache.commons.io.FilenameUtils;
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
    private String omniposeEnvDirPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\miniconda3\\envs\\omnipose\\" : 
            "/opt/miniconda3/envs/omnipose/";
    private String omniposeModelsPath = (IJ.isWindows()) ? System.getProperty("user.home")+"\\.cellpose\\models\\" :
            System.getProperty("user.home")+"/.cellpose/models/";
    public String omniposeBactModel = "bact_phase_omnitorch_0";
    public String omniposeAdnModel = "bact_fluor_omnitorch_0";
     
    private int omniposeDiameter = 10;
    private int omniposeMaskThreshold = 0;
    private double omniposeFlowThreshold = 0.4;
    private boolean useGpu = true;
    
    // Bacteria
    public double minBactSurface = 0.3;
    public double maxBactSurface = 5;
    
    // Adn
    public double minAdnSurface = 0.10;
    public double maxAdnSurface = 5;
    
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
        String[] channelsName = {"Bacteria : ", "Adn : "};
        GenericDialogPlus gd = new GenericDialogPlus("Parameters");
        gd.setInsets​(0, 80, 0);
        gd.addImage(icon);
        gd.addMessage("Channels", Font.getFont("Monospace"), Color.blue);
        int index = 0;
        for (String ch : channelsName) {
            gd.addChoice(ch, channels, channels[index]);
            index++;
        }
        gd.addMessage("Bacteria / ADN detection", Font.getFont("Monospace"), Color.blue);
        gd.addDirectoryField("Omnipose environment directory: ", omniposeEnvDirPath);
        gd.addDirectoryField("Omnipose models path: ", omniposeModelsPath); 
        gd.addMessage("Object size threshold ", Font.getFont("Monospace"), Color.blue);
        gd.addNumericField("Min bacterium surface (µm2): ", minBactSurface);
        gd.addNumericField("Max bacterium surface (µm2): ", maxBactSurface);
        gd.addNumericField("Min Adn surface (µm2)      : ", minAdnSurface);
        gd.addNumericField("Max Adn surface (µm2)      : ", maxAdnSurface);
        
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
        minAdnSurface = (float) gd.getNextNumber();
        maxAdnSurface = (float) gd.getNextNumber();
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
    public Objects3DIntPopulation omniposeDetection(ImagePlus imgBact, String model, double min, double max){

        ImagePlus imgIn = new Duplicator().run(imgBact);
        imgIn.setCalibration(cal);
        // Set Omnipose settings
        CellposeTaskSettings settings = new CellposeTaskSettings(omniposeModelsPath+model, 1, omniposeDiameter, omniposeEnvDirPath);
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
        popFilterSize(popFilter, min, max);
        
        // Close images
        flush_close(imgIn);
        flush_close(imgOut);
        return(popFilter);
    }
    
    
     /**
     * Find coloc between Adn and bacteria
     * set label of colocalized bacteria in adn ID object
     * @param bactPop
     * @param adnPop
     */
    public void adnBactLink(Objects3DIntPopulation bactPop, Objects3DIntPopulation adnPop) {
        if (bactPop.getNbObjects() != 0 && adnPop.getNbObjects() != 0) {
            for (Object3DInt bact : bactPop.getObjects3DInt()) {
                for (Object3DInt adn : adnPop.getObjects3DInt()) {
                    MeasureCentroid adnCenter = new MeasureCentroid(adn);
                    if (bact.contains(adnCenter.getCentroidRoundedAsVoxelInt())){
                       adn.setIdObject(bact.getLabel()); 
                    }
                }
            }
        }
        // remove foci not in bacteria
        adnPop.getObjects3DInt().removeIf(p -> p.getIdObject() == 0);
        adnPop.resetLabels();
    }
    
    /**
     * Compute distance between adn center to bacteria center
     */
    private double adnBactDistance(Object3DInt adn, Object3DInt bact) {
        Voxel3D bactCenter = new MeasureCentroid(bact).getCentroidAsVoxel();
        Voxel3D adnCenter = new MeasureCentroid(adn).getCentroidAsVoxel();
        double dist = bactCenter.distance(adnCenter)*cal.pixelWidth;
        return(dist);
    }
    
    /**
     * 
    */
    private Objects3DIntPopulation findAdnBact(float bactLabel, Objects3DIntPopulation adnPop) {
        Objects3DIntPopulation adnBactPop = new Objects3DIntPopulation();
        for (Object3DInt adn : adnPop.getObjects3DInt()) {
                if (adn.getIdObject() == bactLabel)
                    adnBactPop.addObject(adn);
        }
        adnBactPop.resetLabels();
        return(adnBactPop) ;           
    }
    
    
    /**
     * Compute bacteria area
     * 
     **/
    private double bacteriaSurface(Object3DInt bactObj ) {
        int pixelNb = bactObj.getObject3DPlanes().get(0).getVoxels().size();
        return(pixelNb*pixelSurf);
    }   
    
   
    
    /**
     * Compute bacteria parameters and save them in file
     * @param bactPop
     * @param adnPop
     * @param imgName
     * @param file
     * @throws java.io.IOException
     */
    public void saveResults(Objects3DIntPopulation bactPop, Objects3DIntPopulation adnPop, ImagePlus adnImg, String imgName, BufferedWriter file) throws IOException {
        for (Object3DInt bact : bactPop.getObjects3DInt()) {
            float bactLabel = bact.getLabel();
            double bactSurf = bacteriaSurface(bact);
            VoxelInt feret1Unit = new MeasureFeret(bact).getFeret1Unit();
            VoxelInt feret2Unit = new MeasureFeret(bact).getFeret2Unit();
            double bactLength = feret1Unit.distance(feret2Unit)*cal.pixelWidth;
            
            Objects3DIntPopulation adnBactPop = findAdnBact(bactLabel, adnPop);
            int adnNb = adnBactPop.getNbObjects();
            if (adnNb == 0) {
                file.write(imgName+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+adnNb+"\n");
                file.flush();
            }
            else {
                file.write(imgName+"\t"+bactLabel+"\t"+bactSurf+"\t"+bactLength+"\t"+adnNb+"\t");
                for (Object3DInt adn : adnBactPop.getObjects3DInt()) {
                    int i = (int)adn.getLabel();
                    double adnSurf = bacteriaSurface(adn);
                    double adnInt = new MeasureIntensity(adn, ImageHandler.wrap(adnImg)).getValueMeasurement(MeasureIntensity.INTENSITY_SUM);
                    double adnDist = adnBactDistance(adn, bact);
                    if (i == 1)
                        file.write(i+"\t"+adnSurf+"\t"+adnInt+"\t"+adnDist+"\n");
                    else
                        file.write("\t\t\t\t\t"+i+"\t"+adnSurf+"\t"+adnInt+"\t"+adnDist+"\n");
                }
                file.flush();
            }
        }
    }
   
    
    // Save objects image
    public void drawResults(ImagePlus img, Objects3DIntPopulation bactPop, Objects3DIntPopulation adnPop, String imgName, String outDir) {
        ImageHandler imgBact = ImageHandler.wrap(img).createSameDimensions();
        bactPop.drawInImage(imgBact);
        IJ.run(imgBact.getImagePlus(), "glasbey on dark", "");
        imgBact.getImagePlus().setCalibration(cal);
        FileSaver ImgBactObjectsFile = new FileSaver(imgBact.getImagePlus());
        ImgBactObjectsFile.saveAsTiff(outDir+imgName+"_bacteria.tif");
        
        ImageHandler imgAdn = ImageHandler.wrap(img).createSameDimensions();
        adnPop.resetLabels();
        adnPop.drawInImage(imgAdn);
        imgAdn.getImagePlus().setCalibration(cal);
        IJ.run(imgAdn.getImagePlus(), "glasbey on dark", "");
        FileSaver ImgAdnObjectsFile = new FileSaver(imgAdn.getImagePlus());
        ImgAdnObjectsFile.saveAsTiff(outDir+imgName+"_Adn.tif");
        flush_close(imgBact.getImagePlus());
        flush_close(imgAdn.getImagePlus());
    }
    
}
