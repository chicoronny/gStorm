package org.gstorm.modules;

import java.util.List;

import org.gstorm.interfaces.Element;
import org.gstorm.pipeline.MultiRunModule;
import org.gstorm.tools.Utils;

import ij.ImagePlus;
import net.imglib2.img.Img;
import net.imglib2.img.ImgFactory;
import net.imglib2.img.array.ArrayImgFactory;
import net.imglib2.img.display.imagej.ImageJFunctions;
import net.imglib2.type.numeric.real.FloatType;

/**
 * base class for all renderer plug-ins
 * 
 * @author Ronny Sczech
 *
 */
public abstract class Renderer extends MultiRunModule {
	
	protected final Img< FloatType > img;
	protected final int yBins;
	protected final int xBins;

	protected Renderer(final int width, final int height) {
		xBins=width;
		yBins=height;
		final ImgFactory< FloatType > imgFactory = new ArrayImgFactory< FloatType >(new FloatType());
		img = imgFactory.create( new int[]{ width, height});
	}
	
	public ImagePlus getImage(){
		final ImagePlus ip = ImageJFunctions.wrap(img, "Renderer Window");
		final double max = ip.getStatistics().histMax;
		ip.getProcessor().setMinAndMax(0, max);
		return ip;
	}
	
	@Override
	public boolean check() {
		return inputs.size()==1;
	}
	
	@Override
	public void afterRun(){
		System.out.println("Rendering done in "
				+ (System.currentTimeMillis() - start) + "ms.");
	}
	
	public void show(){
		if (img!=null){
			final ImagePlus ip = ImageJFunctions.show( img );
			ip.getProcessor().setColorModel(Utils.Ice());
			double max = ip.getStatistics().histMax;
			ip.getProcessor().setMinAndMax(0, max);
			ip.updateAndDraw();
			while (ip.isVisible())
				pause(10);
		}
	}
	
	public void preview(List<Element> previewList) {
		for(Element l:previewList) newOutput(l);
		//double max = ip.getStatistics().histMax;
		//ip.getProcessor().setMinAndMax(0, max);
		//ip.updateAndRepaintWindow();
	}

}
