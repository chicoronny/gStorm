package org.gstorm.factories;

import java.util.Map;

import org.gstorm.gui.ConfigurationPanel;
import org.gstorm.interfaces.PluginInterface;
import org.gstorm.modules.Renderer;

/**
 * Factory for rendering localization data
 * 
 * @author Ronny Sczech
 *
 */
public interface RendererFactory extends PluginInterface {
	
	/**
	 * Check that the given settings map is suitable for target detector.
	 *
	 * @param settings 
	 * the map to test.
	 * @return <code>true</code> if the settings map is valid.
	 */
	boolean setAndCheckSettings(final Map<String, Object> settings);
	
	/**
	 *  @return  Renderer to process
	 */
	Renderer getRenderer();
	
	/**
	 * @return  getConfigurationPanel Returns a new GUI panel able to configure the settings suitable for this
	 * specific factory.
	 */
	ConfigurationPanel getConfigurationPanel();
	
	Map<String, Object> getInitialSettings();
	
	String KEY_xmin = "xmin";
	String KEY_xmax = "xmax";
	String KEY_ymin = "ymin";
	String KEY_ymax = "ymax";
	String KEY_xBins = "xbins";
	String KEY_yBins = "ybins";
	String KEY_zmin = "zmin";
	String KEY_zmax = "zmax";
}
