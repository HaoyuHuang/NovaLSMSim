package edu.usc.nova;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import edu.usc.nova.TinyRangeSim.LoadImbalanceMetric;
import edu.usc.nova.TinyRangeSim.Range;
import edu.usc.nova.TinyRangeSim.SubRange;

public class WriteReduction {

	public static final String ZIPFIAN_FILE = "/tmp/zipfian";
	public static final int numberOfKeys = 10000000;

	private static void ZipfianReduction(int keys) {
		BufferedReader br;
		try {
			br = new BufferedReader(new FileReader(new File(ZIPFIAN_FILE)));
			String line;
			double[] ref = new double[numberOfKeys];
			int key = 0;
			double sum = 0;
			while ((line = br.readLine()) != null) {
				ref[key] = Double.parseDouble(line);
				sum += ref[key];
				key++;
			}
			double reduces = 0;
			for (int i = 0; i < keys; i++) {
				reduces += ref[i];
			}
			System.out.println(String.format("%.2f", reduces / sum));
			br.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public static void main(String[] args) {
		ZipfianReduction(313);
	}

}
