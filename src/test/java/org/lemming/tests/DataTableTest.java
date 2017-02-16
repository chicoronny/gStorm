package org.lemming.tests;

import static org.junit.Assert.*;

import java.io.File;
import java.util.Map;
import java.util.concurrent.Executors;

import org.gstorm.interfaces.Store;
import org.gstorm.modules.DataTable;
import org.gstorm.modules.ReadLocalizationPrecision3D;
import org.gstorm.pipeline.ExtendableTable;
import org.gstorm.pipeline.Manager;
import org.junit.Before;
import org.junit.Test;

public class DataTableTest {
	
	private Manager pipe;
	private ExtendableTable table;
	private Map<Integer, Store> map;

	@Before
	public void setUp() {
		pipe = new Manager(Executors.newCachedThreadPool());
		
		ReadLocalizationPrecision3D reader = new ReadLocalizationPrecision3D(new File(System.getProperty("user.home")+"/ownCloud/Tubulin1.csv"),",");
		pipe.add(reader);
		
		DataTable workspace = new DataTable();
		table = workspace.getTable();
		pipe.add(workspace);
		
		pipe.linkModules(reader, workspace, true, 128);
		map = pipe.getMap();
	}

	@Test
	public void test() {
		pipe.run();
		System.out.println("Table rows: " + table.getNumberOfRows());
		assertEquals(true,map.values().iterator().next().isEmpty());
	}

}
