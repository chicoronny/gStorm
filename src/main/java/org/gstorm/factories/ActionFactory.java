package org.gstorm.factories;

import org.gstorm.interfaces.PluginInterface;

/**
 * Factory for manipulating localization data
 * 
 * @author Ronny Sczech
 *
 */
public interface ActionFactory extends PluginInterface {
	
	Runnable create();
}
