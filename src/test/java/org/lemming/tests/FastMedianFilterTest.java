package org.lemming.tests;

import static org.junit.Assert.*;

import java.util.Map;
import java.util.concurrent.Executors;

import org.gstorm.interfaces.Store;
import org.gstorm.modules.ImageLoader;
import org.gstorm.modules.SaveImages;
import org.gstorm.pipeline.Manager;
import org.gstorm.plugins.NMSFastMedian;
import org.gstorm.tools.Utils;
import org.junit.Before;
import org.junit.Test;

import ij.ImagePlus;

@SuppressWarnings("rawtypes")
public class FastMedianFilterTest {

	private Manager pipe;
	private Map<Integer, Store> map;
	
	@Before
	public void setUp() {
		pipe = new Manager(Executors.newCachedThreadPool());

		ImageLoader tif = new ImageLoader(new ImagePlus("/Users/ronny/Documents/TubulinAF647.tif"), Utils.readCameraSettings("camera.props"));
		pipe.add(tif);

		NMSFastMedian fmf = new NMSFastMedian(50, true, 1, 15);

		pipe.add(fmf);

		SaveImages saver = new SaveImages("/home/ronny/Bilder/out.tif");
		pipe.add(saver);
		
		pipe.linkModules(tif, fmf);
		pipe.linkModules(fmf, saver);
		
		map = pipe.getMap();
	}

	@Test
	public void test() {
		pipe.run();
		assertEquals(true, map.values().iterator().next().isEmpty());
	}

}
