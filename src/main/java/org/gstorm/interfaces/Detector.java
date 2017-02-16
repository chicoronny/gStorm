package org.gstorm.interfaces;

import org.gstorm.pipeline.FrameElements;

import net.imglib2.type.numeric.RealType;

/**
 * interface for all detector plug-ins
 * 
 * @author Ronny Sczech
 *
 * @param <T> - data type
 */
public interface Detector<T extends RealType<T>> extends ModuleInterface {

	FrameElements<T> detect(Frame<T> frame);

}
