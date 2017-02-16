package org.lemming.tests;

import static org.junit.Assert.*;

import java.io.File;
import java.util.Map;
import java.util.concurrent.Executors;

import org.gstorm.interfaces.Store;
import org.gstorm.modules.ImageLoader;
import org.gstorm.modules.SaveLocalizations;
import org.gstorm.modules.UnpackElements;
import org.gstorm.pipeline.Manager;
import org.gstorm.plugins.NMSDetector;
import org.gstorm.tools.Utils;
import org.junit.Before;
import org.junit.Test;

import ij.ImagePlus;

@SuppressWarnings("rawtypes")
public class NMSFinderTest {

	private Manager pipe;
	private Map<Integer, Store> map;
	
	@Before
	public void setUp() {
		pipe = new Manager(Executors.newCachedThreadPool());
		final ImagePlus image = new ImagePlus(System.getProperty("user.home")+"/ownCloud/p500ast_.tif");
		ImageLoader tif = new ImageLoader(image, Utils.readCameraSettings("camera.props"));
		pipe.add(tif);

		NMSDetector peak = new NMSDetector(700, 9, 0);
		pipe.add(peak);

		UnpackElements unpacker = new UnpackElements();
		pipe.add(unpacker);

		SaveLocalizations saver = new SaveLocalizations(new File(System.getProperty("user.home") + "/ownCloud/nmsfinder.csv"));
		pipe.add(saver);
		
		pipe.linkModules(tif, peak, true, image.getStackSize());
		pipe.linkModules(peak, unpacker);
		pipe.linkModules(unpacker, saver);
		map = pipe.getMap();
	}

	@Test
	public void test() {
		pipe.run();
		System.out.println("");		
		assertEquals(true,map.values().iterator().next().isEmpty());
	}

}
