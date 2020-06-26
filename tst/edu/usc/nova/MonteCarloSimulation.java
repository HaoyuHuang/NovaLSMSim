package edu.usc.nova;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Collections;
import java.util.Comparator;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import edu.usc.distributions.ZipfianGenerator;
import edu.usc.workload.KeySizeGenerator;
import edu.usc.workload.ValueSizeGenerator;

public class MonteCarloSimulation {

	public static class Record {
		int recordId = 0;
		int tenantId = 0;
		int refCount = 0;
		double accessFreq = 0;

		public Record() {

		}

		public Record(int recordId, int tenantId, int refCount,
				double accessFreq) {
			super();
			this.recordId = recordId;
			this.tenantId = tenantId;
			this.refCount = refCount;
			this.accessFreq = accessFreq;
		}

	}

	public static void main(String[] args) throws Exception {
		computeRefCounts(0.99);
	}

	public static void computeRefCounts(double zipfConstant) throws Exception {
		System.out
				.println("Monte Carlo simulation on constant " + zipfConstant);
		int maxLoop = 100000000;
		int nrecords = 10000000;
		Random rand = new Random();
		int first_range = nrecords / 10 / 64;
		int first_server = nrecords / 10;
		int refs_first_range = 0;
		int refs_first_server = 0;

		ZipfianGenerator zipf = new ZipfianGenerator(nrecords, zipfConstant,
				rand);
		for (int i = 0; i < maxLoop; i++) {
			int key = zipf.nextValue().intValue();
			if (key < first_range) {
				refs_first_range++;
			}
			if (key < first_server) {
				refs_first_server++;
			}
		}
		System.out.println(String.format("%f:%f",
				(double) refs_first_range / (double) maxLoop,
				(double) refs_first_server / (double) maxLoop));

	}

	public static void computeRefCount(double zipfConstant) throws Exception {
		System.out
				.println("Monte Carlo simulation on constant " + zipfConstant);
		int maxLoop = 100000000;
		int nrecords = 1000000000;
		Random rand = new Random();

		int[] refCount = new int[nrecords];
		ZipfianGenerator zipf = new ZipfianGenerator(nrecords, zipfConstant,
				rand);
		for (int i = 0; i < maxLoop; i++) {
			int key = zipf.nextValue().intValue();
			refCount[key] += 1;
		}
		BufferedWriter bw = new BufferedWriter(
				new FileWriter(new File("/tmp/zipfian")));
		for (int i = 0; i < refCount.length; i++) {
			bw.write(String.format("%f\n", (double) refCount[i]));
		}
		bw.close();
	}

	private static void writeToFile(String filename, int[] refCount)
			throws IOException {

		double sum = 0;
		for (int i = 0; i < refCount.length; i++) {
			sum += refCount[i];
		}
		double percentile = 0.0;

	}

}
