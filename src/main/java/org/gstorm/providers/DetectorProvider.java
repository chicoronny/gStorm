package org.gstorm.providers;

import org.gstorm.factories.DetectorFactory;

public class DetectorProvider extends AbstractProvider<DetectorFactory> {

	public DetectorProvider() {
		super(DetectorFactory.class);
	}
	
	public static void main( final String[] args ){
		final DetectorProvider provider = new DetectorProvider();
		System.out.println( provider.echo() );
	}


}
