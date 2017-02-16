package org.gstorm.gui;

import java.util.Map;
import java.util.concurrent.Executors;
import java.util.concurrent.ScheduledExecutorService;

import javax.swing.JPanel;

public abstract class ConfigurationPanel extends JPanel {
	
	private static final long serialVersionUID = 3160662804934210143L;
	
	public static final String propertyName = "CONFIG_PANEL";
	
	protected final ScheduledExecutorService schedule = Executors.newScheduledThreadPool(4);
	
	/* Echo the parameters of the given settings on this panel. */
	public abstract void setSettings(final Map<String, Object> settings);
	
	/**
	 * @return  a new settings map object with its values set by this panel.
	 */
	public abstract Map<String, Object> getSettings();
		
	void fireChanged() {
		firePropertyChange(propertyName, null, getSettings());
	}

	public void close(){
		schedule.shutdown();
	}

}
