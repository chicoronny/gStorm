package org.lemming.tests;

import static org.junit.Assert.*;

import java.io.File;
import java.util.Map;

import org.junit.Before;
import org.junit.Test;
import org.lemming.interfaces.Store;
import org.lemming.modules.ImageLoader;
import org.lemming.modules.SaveLocalizations;
import org.lemming.pipeline.AbstractModule;
import org.lemming.pipeline.Localization;
import org.lemming.pipeline.Manager;
import org.lemming.plugins.GaussianFitter;
import org.lemming.plugins.SetPeak;
//import org.lemming.plugins.SymmetricGaussianFitter;
import org.lemming.tools.LemmingUtils;

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
		
        File file = new File("H:\\Images\\test9000small.tif");
        
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
		
		AbstractModule tif = new ImageLoader(loc_im, LemmingUtils.readCameraSettings("camera.props"));
		AbstractModule peak = new SetPeak(new Localization(7.0,7.0,10.0,1L));
		AbstractModule fitter = new GaussianFitter<>(7, LemmingUtils.readCSV("H:\\Images\\set1-calt.csv"));
		AbstractModule saver = new SaveLocalizations(new File("H:\\Images\\test9000-g.csv"));
		
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