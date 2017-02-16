package org.lemming.tests;

import static org.junit.Assert.*;

import java.io.File;
import java.util.Map;
import java.util.concurrent.Executors;

import org.gstorm.interfaces.Store;
import org.gstorm.modules.ImageLoader;
import org.gstorm.modules.SaveLocalizations;
import org.gstorm.pipeline.AbstractModule;
import org.gstorm.pipeline.Manager;
import org.gstorm.plugins.CUDA_Gaussian_Fitter;
import org.gstorm.plugins.NMSDetector;
import org.gstorm.tools.Utils;
import org.junit.Before;
import org.junit.Test;

import ij.ImagePlus;
import ij.plugin.FileInfoVirtualStack;
import ij.plugin.FolderOpener;

@SuppressWarnings("rawtypes")
public class GPUFitterTest {

	private Manager pipe;
	private Map<Integer, Store> storeMap;
	private ImagePlus loc_im;
	
	
	@Before
	public void setUp() throws Exception {
		
        //File file = new File("D:/Images/DRG_KO_5_1.tif");
		//File file = new File("D:/Images/DRG_WT_MT_A647_1.tif");
        File file = new File("H:/Images/test9000.tif");
		//File file = new File("D:/ownCloud/Tubulin1.tif");
        //File file = new File(System.getProperty("user.home")+"/ownCloud/exp-images.tif");

        
		if (file.isDirectory()){
        	FolderOpener fo = new FolderOpener();
        	fo.openAsVirtualStack(true);
        	loc_im = fo.openFolder(file.getAbsolutePath());
        }
        
        if (file.isFile()){
        	loc_im = FileInfoVirtualStack.openVirtual(file.getAbsolutePath());
        }
	
	    if (loc_im ==null)
		    throw new Exception("File not found");
	    
		AbstractModule tif = new ImageLoader(loc_im, Utils.readCameraSettings("camera.props"));
		//AbstractModule peak = new NMSDetector(70,7);
		AbstractModule peak = new NMSDetector(70,7,0); //
		//AbstractModule peak = new DoGFinder(4.5f,13); //DRG_KO_5_1.tif
		//AbstractModule peak = new NMSDetector(2000,5); //DRG_WT_MT_A647_1.tif
		//AbstractModule fitter = new MLE_Fitter<>(6,1152*8);
		AbstractModule fitter = new CUDA_Gaussian_Fitter<>(6,1152*8);
		//AbstractModule saver = new SaveLocalizations(new File("H:/Images/test9000-mle.tif"));
		AbstractModule saver = new SaveLocalizations(new File("H:/Images/test9000-cudaG.csv"));
		

		pipe = new Manager(Executors.newCachedThreadPool());
		pipe.add(tif);
		pipe.add(peak);
		pipe.add(fitter);
		pipe.add(saver);

		pipe.linkModules(tif, peak, true, loc_im.getStackSize());
		pipe.linkModules(peak,fitter);
		pipe.linkModules(fitter,saver,false, 128);
		storeMap = pipe.getMap();
	}

	@Test
	public void test() {
		long start = System.currentTimeMillis();
		pipe.run();
		long end = System.currentTimeMillis();
		System.out.println("Overall: " + (end-start) +"ms");
		assertEquals(true,storeMap.values().iterator().next().isEmpty());
		assertEquals(true,storeMap.values().iterator().next().isEmpty());
	}

}
