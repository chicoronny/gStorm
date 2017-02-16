package org.gstorm.providers;

import org.gstorm.factories.RendererFactory;

public class RendererProvider extends AbstractProvider<RendererFactory> {

	public RendererProvider() {
		super(RendererFactory.class);
	}

	public static void main(String[] args) {
		final RendererProvider provider = new RendererProvider();
		System.out.println( provider.echo() );
	}

}
