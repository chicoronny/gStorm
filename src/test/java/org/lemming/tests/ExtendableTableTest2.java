package org.lemming.tests;

import java.io.File;
import java.util.concurrent.Executors;

import org.gstorm.modules.DataTable;
import org.gstorm.modules.ImageLoader;
import org.gstorm.pipeline.Manager;
import org.gstorm.plugins.PeakFinder;
import org.gstorm.plugins.QuadraticFitter;
import org.gstorm.tools.Utils;
import org.junit.Before;
import org.junit.Test;

import ij.ImagePlus;
import ij.plugin.FileInfoVirtualStack;
import ij.plugin.FolderOpener;

public class ExtendableTableTest2 {
	
	private Manager pipe;
	private DataTable dt;
	private ImagePlus loc_im;

	@SuppressWarnings("rawtypes")
	@Before
	public void setUp() throws Exception {
		pipe = new Manager(Executors.newCachedThreadPool());
		
		File file = new File("/Users/ronny/ownCloud/storm/p500ast.tif");
        
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
		
		ImageLoader tif = new ImageLoader(loc_im, Utils.readCameraSettings("camera.props"));	
		
		PeakFinder peak = new PeakFinder(700,6,0);
		QuadraticFitter fitter = new QuadraticFitter(10);
		dt = new DataTable(); 
		
		pipe.add(tif);
		pipe.add(peak);
		pipe.add(fitter);
		pipe.add(dt);
		
		pipe.linkModules(tif, peak, true, loc_im.getStackSize());
		pipe.linkModules(peak,fitter);
		pipe.linkModules(fitter,dt);
		
	}

	@Test
	public void test() {
		pipe.run();
		System.out.println(dt.getTable().columnNames().toString());
		System.out.println(dt.getTable().getNumberOfRows());
	}

}
