package org.gstorm;

import javax.swing.SwingUtilities;

import org.gstorm.gui.Controller;

import ij.plugin.PlugIn;

/**
 * Main entry point of the program
 * 
 * @author Ronny Sczech
 *
 */
@SuppressWarnings("rawtypes")
class gStorm implements PlugIn {

	@Override
	public void run(String arg) {
		Controller frame = new Controller();
		frame.setVisible(true);
	}

	public static void main(String[] args) {
		SwingUtilities.invokeLater(new Runnable() {
            
			@Override
			public void run() {
			try {
                Controller frame = new Controller();
                frame.setVisible(true);
            } catch (Exception e) {
                e.printStackTrace();
            }
			}
        });
	}

}
