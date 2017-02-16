package org.gstorm.pipeline;

import java.util.concurrent.LinkedBlockingQueue;

import org.gstorm.interfaces.Element;
import org.gstorm.interfaces.Store;

/**
 * a queue using a linked list as structure
 * 
 * @author Ronny Sczech
 *
 */
public class LinkedStore extends LinkedBlockingQueue<Element> implements Store {

	/**
	 * 
	 */
	private static final long serialVersionUID = -5362955780141552726L;

	public LinkedStore(int capacity) {
		super(capacity);
	}

}
