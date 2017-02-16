package org.gstorm.factories;

import java.util.Map;

import org.gstorm.gui.ConfigurationPanel;
import org.gstorm.interfaces.Detector;
import org.gstorm.interfaces.PluginInterface;

import net.imglib2.type.NativeType;
import net.imglib2.type.numeric.RealType;

/**
 * Factory for detectors
 * 
 * @author Ronny Sczech
 *
 */
public interface DetectorFactory extends PluginInterface{

	/**
	 * Check that the given settings map is suitable for target detector.
	 *
	 * @param settings 
	 * the map to test.
	 * @return <code>true</code> if the settings map is valid.
	 */
	boolean setAndCheckSettings(final Map<String, Object> settings);
	
	/**
	 * @param <T> data type
	 * @return  Module to process
	 */
	<T extends RealType<T> & NativeType<T>> Detector<T> getDetector();
	
	/**
	 * @return  getConfigurationPanel Returns a new GUI panel able to configure the settings suitable for this
	 * specific detector factory.
	 */
	ConfigurationPanel getConfigurationPanel();
	
	boolean hasPreProcessing();
	
}
