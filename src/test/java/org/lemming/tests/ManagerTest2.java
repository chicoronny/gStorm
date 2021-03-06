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
import org.gstorm.plugins.CentroidFitter;
import org.gstorm.plugins.NMSDetector;
import org.gstorm.tools.Utils;
import org.junit.Before;
import org.junit.Test;

import ij.ImagePlus;
import ij.plugin.FileInfoVirtualStack;
import ij.plugin.FolderOpener;

@SuppressWarnings("rawtypes")
public class ManagerTest2 {

	private Manager pipe;
	private Map<Integer, Store> storeMap;
	private ImagePlus loc_im;
	
	@Before
	public void setUp() throws Exception {
		
        File file = new File(System.getProperty("user.home")+"/ownCloud/storm/TubulinAF647.tif");
        
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

		AbstractModule peak = new NMSDetector(700,7,0);
		//Fitter fitter = new QuadraticFitter(10);
		//@SuppressWarnings("unchecked")
		//Fitter fitter = new Fitter(7, Settings.readCSV(System.getProperty("user.home")+"/ownCloud/storm/calTest.csv").get("param"));
		AbstractModule fitter = new CentroidFitter(10, 700);

		AbstractModule saver = new SaveLocalizations(new File(System.getProperty("user.home")+"/ownCloud/storm/fitted2.csv"));
		
		pipe = new Manager(Executors.newCachedThreadPool());
		pipe.add(tif);
		pipe.add(peak);
		pipe.add(fitter);
		pipe.add(saver);
		
		pipe.linkModules(tif, peak, true, loc_im.getStackSize());
		pipe.linkModules(peak,fitter);
		pipe.linkModules(fitter,saver,false, 100);
		storeMap = pipe.getMap();
	}

	@Test
	public void test() {
		pipe.run();
		assertEquals(true,storeMap.values().iterator().next().isEmpty());
		assertEquals(true,storeMap.values().iterator().next().isEmpty());
	}

}
