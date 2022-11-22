import BacteriaOmni_Tools.Tools;
import ij.*;
import ij.plugin.PlugIn;
import ij.plugin.ZProjector;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import loci.common.services.DependencyException;
import loci.common.services.ServiceException;
import loci.common.services.ServiceFactory;
import loci.formats.FormatException;
import loci.formats.meta.IMetadata;
import loci.formats.services.OMEXMLService;
import loci.plugins.BF;
import loci.plugins.util.ImageProcessorReader;
import loci.plugins.in.ImporterOptions;
import mcib3d.geom2.Objects3DIntPopulation;
import org.apache.commons.io.FilenameUtils;
import org.scijava.util.ArrayUtils;


/**
 * Detect foci bacteria with OmniPose
 * @author Orion-CIRB
 */
public class Foci_Bacteria implements PlugIn {
    
    Tools tools = new Tools();
    private String imageDir = "";
    public String outDirResults = "";
    public BufferedWriter results;
   
    
    public void run(String arg) {
        try {
            if (!tools.checkInstalledModules()) {
                return;
            } 
            
            imageDir = IJ.getDirectory("Choose directory containing image files...");
            if (imageDir == null) {
                return;
            }   
            // Find images with extension
            String file_ext = tools.findImageType(new File(imageDir));
            
            ArrayList<String> imageFiles = new ArrayList();
            tools.findImages(imageDir, file_ext, imageFiles);
            if (imageFiles.isEmpty()) {
                IJ.showMessage("Error", "No images found with " + file_ext + " extension");
                return;
            }
            
            // Create output folder
            outDirResults = imageDir + File.separator + "Results" + File.separator;
            File outDir = new File(outDirResults);
            if (!Files.exists(Paths.get(outDirResults))) {
                outDir.mkdir();
            }
            // Write header in results file
             String header = "Image name\t# bacterium\tBacterium surface (µm2)\tBacterium length (µm)\tFoci1 number\tFoci2 number\tFoci1 min. distance to bacterium\t"
            + "Foci2 min. distance to bacterium\tFoci1-foci2 min distance\tFoci1-foci2 max distance\n";
   
            FileWriter fwResults = new FileWriter(outDirResults + "results.xls", false);
            results = new BufferedWriter(fwResults);
            results.write(header);
            results.flush();
            
            // Create OME-XML metadata store of the latest schema version
            ServiceFactory factory;
            factory = new ServiceFactory();
            OMEXMLService service = factory.getInstance(OMEXMLService.class);
            IMetadata meta = service.createOMEXMLMetadata();
            ImageProcessorReader reader = new ImageProcessorReader();
            reader.setMetadataStore(meta);
            reader.setId(imageFiles.get(0));
            String[] channels = tools.findChannels(imageFiles.get(0), meta, reader);
            // Image calibration
            tools.findImageCalib(meta);
            
            // Dialog box
            String[] chs = tools.dialog(channels);
            if (tools.canceled || chs == null) {
                IJ.showMessage("Error", "Plugin canceled");
                return;
            }
            
            
            for (String f : imageFiles) {
                reader.setId(f);
                String rootName = FilenameUtils.getBaseName(f);
                tools.print("--- ANALYZING IMAGE " + rootName + " ------");
                
                ImporterOptions options = new ImporterOptions();
                options.setId(f);
                options.setQuiet(true);
                options.setColorMode(ImporterOptions.COLOR_MODE_GRAYSCALE);
                options.setSplitChannels(true);
                
                // Open phase channel
                int indexCh = ArrayUtils.indexOf(channels, chs[0]);
                System.out.println("Opening phase channel "+chs[0] );
                ImagePlus bactStack = BF.openImagePlus(options)[indexCh];
                
                // Detect bacteria with Omnipose
                tools.print("- Detecting bacteria -");
                ImagePlus imgBact = tools.doZProjection(bactStack, ZProjector.AVG_METHOD);
                tools.flush_close(bactStack);
                Objects3DIntPopulation bactPop = tools.omniposeDetection(imgBact);
                System.out.println(bactPop.getNbObjects() + " bacteria found");
                
                
                // Open foci1 channel
                indexCh = ArrayUtils.indexOf(channels, chs[1]);
                System.out.println("Opening foci1 channel "+chs[1]);
                ImagePlus foci1Stack = BF.openImagePlus(options)[indexCh];
                ImagePlus imgFoci1 = tools.doZProjection(foci1Stack, ZProjector.MAX_METHOD);
                tools.flush_close(foci1Stack);
                tools.print("- Detecting foci1 -");
                Objects3DIntPopulation foci1Pop = tools.findFoci(imgFoci1);
                System.out.println(foci1Pop.getNbObjects() + " foci1 found");
                tools.fociBactLink(bactPop, foci1Pop);
                System.out.println(foci1Pop.getNbObjects() + " foci1 found in bacteria");
                tools.flush_close(imgFoci1);
                
                
                // Open foci2 channel
                indexCh = ArrayUtils.indexOf(channels, chs[2]);
                System.out.println("Opening foci2 channel "+chs[2]);
                ImagePlus foci2Stack = BF.openImagePlus(options)[indexCh];
                ImagePlus imgFoci2 = tools.doZProjection(foci2Stack, ZProjector.MAX_METHOD);
                tools.flush_close(foci2Stack);
                tools.print("- Detecting foci2 -");
                Objects3DIntPopulation foci2Pop = tools.findFoci(imgFoci2);
                System.out.println(foci2Pop.getNbObjects() + " foci2 found");
                tools.fociBactLink(bactPop, foci2Pop);
                System.out.println(foci2Pop.getNbObjects() + " foci2 found in bacteria");
                tools.flush_close(imgFoci2);
                                
                // Save results
                tools.print("- Saving results -");
                tools.saveResults(bactPop, foci1Pop, foci2Pop, rootName, results);
                
                // Save images
                tools.drawResults(imgBact, bactPop, foci1Pop, foci2Pop, rootName, outDirResults);
                tools.flush_close(imgBact);
            }
        
            tools.print("--- All done! ---");
            
        }   catch (IOException | FormatException | DependencyException | ServiceException ex) {
            Logger.getLogger(Foci_Bacteria.class.getName()).log(Level.SEVERE, null, ex);
        }  
    }
}    
