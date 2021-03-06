package org.gstorm.tools;

import java.awt.event.KeyEvent;
import java.awt.event.KeyListener;
import java.util.concurrent.ScheduledExecutorService;
import java.util.concurrent.ScheduledFuture;
import java.util.concurrent.TimeUnit;

/**
 * A  that waits a certain before reading the input from a widget.
 * 
 * @author Ronny Sczech
 *
 */
public class WaitForKeyListener implements KeyListener {
	
	private long delay = 1000;
	private final ScheduledExecutorService executor;
	private ScheduledFuture< ? > future;
	private final Runnable command;


	public WaitForKeyListener(ScheduledExecutorService executor, long delay, Runnable command) {
		this.executor = executor;
		this.delay = delay;
		this.command = command;
	}

	@Override
	public void keyTyped(KeyEvent event) {
		if (future != null && !future.isDone()) {
			future.cancel(false);
		}
		future = executor.schedule(command, delay, TimeUnit.MILLISECONDS);
	}

	@Override
	public void keyPressed(KeyEvent e) {
	}

	@Override
	public void keyReleased(KeyEvent e) {

	}

}
