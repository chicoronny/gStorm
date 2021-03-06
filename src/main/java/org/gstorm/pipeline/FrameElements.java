package org.gstorm.pipeline;

import java.util.List;

import org.gstorm.interfaces.Element;
import org.gstorm.interfaces.Frame;

/**
 * all objects connected to a frame
 * 
 * @author Ronny Sczech
 *
 * @param <T> - data type
 */
public class FrameElements<T> implements Element {
	
	private boolean isLast;
	private final List<Element> list;
	private final Frame<T> frame;

	public FrameElements(List<Element> list_, Frame<T> f) {
		list = list_;
		frame = f;
	}

	@Override
	public boolean isLast() {
		return isLast;
	}

	@Override
	public void setLast(boolean b) {
		isLast = b;
	}

	public Frame<T> getFrame(){
		return frame;
	}
	
	public List<Element> getList(){
		return list;
	}
	
	public String toString(){
		return ""+frame.toString() + " " + list.toString();
	}
}
