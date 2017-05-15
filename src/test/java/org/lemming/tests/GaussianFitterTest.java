package org.lemming.tests;

import static org.junit.Assert.*;

import java.io.File;
import java.util.Map;

import org.gstorm.interfaces.Store;
import org.gstorm.modules.ImageLoader;
import org.gstorm.modules.SaveLocalizations;
import org.gstorm.pipeline.AbstractModule;
import org.gstorm.pipeline.Localization;
import org.gstorm.pipeline.Manager;
import org.gstorm.plugins.GaussianFitter;
import org.gstorm.plugins.NMSDetector;
import org.gstorm.plugins.SetPeak;
import org.gstorm.tools.Utils;
import org.junit.Before;
import org.junit.Test;

import ij.ImagePlus;
import ij.plugin.FileInfoVirtualStack;
import ij.plugin.FolderOpener;

@SuppressWarnings("rawtypes")
public class GaussianFitterTest {

	private Manager pipe;
	private Map<Integer, Store> storeMap;
	private ImagePlus loc_im;
	
	@Before
	public void setUp() throws Exception {
		
        //File file = new File("H:\\Images\\test9000.tif");
        File file = new File("/home/ronny/Bilder/sequence-as-stack-MT1.0-2D-Exp.tif");
        
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
		AbstractModule peak = new NMSDetector(25,3,1);
		//AbstractModule peak = new SetPeak(new Localization(17.0,20.0,10.0,1L));
		AbstractModule fitter = new GaussianFitter<>(4, Utils.readCSV("/home/ronny/Bilder/set1-calt.csv"));
		//AbstractModule saver = new SaveLocalizations(new File("H:\\Images\\test9000-g.csv"));
		AbstractModule saver = new SaveLocalizations(new File("/home/ronny/Bilder/sequence-as-stack-MT1.0-2D-Exp-G.csv"));
		
		pipe = new Manager();
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
