package org.gstorm.plugins;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;

import ij.IJ;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.gstorm.factories.FitterFactory;
import org.gstorm.gui.ConfigurationPanel;
import org.gstorm.gui.FitterPanel;
import org.gstorm.interfaces.Element;
import org.gstorm.interfaces.Frame;
import org.gstorm.math.Gaussian2DFitter;
import org.gstorm.modules.CPU_Fitter;
import org.gstorm.modules.Fitter;
import org.gstorm.pipeline.Localization;
import org.gstorm.pipeline.LocalizationPrecision3D;
import org.gstorm.tools.Utils;
import org.scijava.plugin.Plugin;

import net.imglib2.Interval;
import net.imglib2.RandomAccessibleInterval;
import net.imglib2.type.numeric.RealType;
import net.imglib2.view.Views;

public class GaussianFitter<T extends RealType<T>> extends CPU_Fitter<T> {

	private static final String NAME = "Gaussian";
	private static final String KEY = "GAUSSIANFITTER";
	private static final String INFO_TEXT = "<html>" + "Gaussian Fitter Plugin (with sx and sy)" + "</html>";

	private final Map<String, Object> params;

	public GaussianFitter(int windowSize, final Map<String,Object> params) {
		super(windowSize);
		this.params=params;
	}

	@Override
	public List<Element> fit(final List<Element> sliceLocs, Frame<T> frame, final long halfKernel) {
		final double pixelDepth = frame.getPixelDepth();
		final RandomAccessibleInterval<T> pixels = frame.getPixels();
		final List<Element> found = new ArrayList<>();
		long[] imageMin = new long[2];
		long[] imageMax = new long[2];
		for (Element el : sliceLocs) {
			final Localization loc = (Localization) el;
			final long x = Math.round(loc.getX().doubleValue()/pixelDepth);
			final long y = Math.round(loc.getY().doubleValue()/pixelDepth);
			pixels.min(imageMin);
			pixels.max(imageMax);
			final Interval roi = cropInterval(imageMin,imageMax,new long[]{x - halfKernel,y - halfKernel},new long[]{x + halfKernel,y + halfKernel});
			final Gaussian2DFitter<T> gf = new Gaussian2DFitter<>(Views.interval(pixels, roi), 200, 200);
			double[] result;
			result = gf.fit();
			if (result != null) {
				double SxSy = result[2] * result[2] - result[3] * result[3];
				result[0] *= pixelDepth;
				result[1] *= pixelDepth;
				result[6] *= pixelDepth;
				result[7] *= pixelDepth;
				found.add(new LocalizationPrecision3D( result[0], result[1], calculateZ(SxSy)
					,result[6], result[7], result[8], result[4], loc.getFrame()));
			}
		}
		return found;
	}

	private double calculateZ(final double SxSy) {
		final double[] zgrid = (double[]) params.get("zgrid");
		return calcIterZ(SxSy, zgrid[0], zgrid[zgrid.length-1], 1e-4);
	}

	private double calcIterZ(double SxSy, double start_, double end, double precision) {
		if(start_ <= end) return start_;
		final PolynomialSplineFunction psx = (PolynomialSplineFunction) params.get("psx");
		final PolynomialSplineFunction psy = (PolynomialSplineFunction) params.get("psy");
		final double zStep = (end-start_)/10;
		double curveWx = psx.value(start_);
		double curveWy = psy.value(start_);
		double calib = curveWx * curveWx - curveWy * curveWy;
		double distance = Math.abs(calib - SxSy);
		double idx = start_;
		for (double c = start_ + zStep; c <= end; c += zStep) {
			if(!psx.isValidPoint(c) || !psy.isValidPoint(c)){
				idx = c-zStep;
				break;
			}
			curveWx = psx.value(c);
			curveWy = psy.value(c);
			calib = curveWx * curveWx - curveWy * curveWy;
			double cdistance = Math.abs(calib - SxSy);
			if (cdistance < distance) {
				idx = c;
				distance = cdistance;
			}
		}
		if (zStep <= precision) {
			return idx;
		}
		return calcIterZ(SxSy, idx - zStep, idx + zStep, precision);
	}


	@Plugin(type = FitterFactory.class)
	public static class Factory implements FitterFactory {

		private Map<String, Object> settings;
		private final ConfigurationPanel configPanel = new FitterPanel();

		@Override
		public String getInfoText() {
			return INFO_TEXT;
		}

		@Override
		public String getKey() {
			return KEY;
		}

		@Override
		public String getName() {
			return NAME;
		}

		@Override
		public boolean setAndCheckSettings(Map<String, Object> settings) {
			this.settings = settings;
			return settings.get(FitterPanel.KEY_CALIBRATION_FILENAME) != null;
		}

		@Override
		public <T extends RealType<T>> Fitter<T> getFitter() {
			final int windowSize = (int) settings.get(FitterPanel.KEY_WINDOW_SIZE);
			final String calibFileName = (String) settings.get(FitterPanel.KEY_CALIBRATION_FILENAME);
			if (calibFileName == null) {
				IJ.error("No Calibration File!");
				return null;
			}
			Map<String, Object> cal = Utils.readCSV(calibFileName);
			return new GaussianFitter<>(windowSize, cal);
		}

		@Override
		public ConfigurationPanel getConfigurationPanel() {
			configPanel.setName(KEY);
			return configPanel;
		}

		@Override
		public int getHalfKernel() {
			return size;
		}
		
		@Override
		public boolean hasGPU() {
			return false;
		}
	}
}
