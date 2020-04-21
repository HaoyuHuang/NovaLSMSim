package edu.usc.nova;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.List;
import java.util.Random;

import com.google.common.collect.Lists;

import edu.usc.distributions.ZipfianGenerator;

public class LRUCacheMain {

	public static void cdf() throws Exception {
		String file = "/home/haoyu/Documents/cache/NovaSimulator/config/responsetime";
		BufferedReader br = new BufferedReader(new FileReader(new File(file)));
		String line = null;
		List<List<Integer>> lists = Lists.newArrayList();
		lists.add(Lists.newArrayList());
		lists.add(Lists.newArrayList());
		lists.add(Lists.newArrayList());
		int min = 0;
		int max = 0;
		while ((line = br.readLine()) != null) {
			String[] ems = line.split(" ");
			for (int i = 0; i < 3; i++) {
				int value = Integer.parseInt(ems[i + 2]);
				lists.get(i).add(value);
				max = Math.max(max, value);
				min = Math.min(min, value);
			}
		}

		List<List<Integer>> cdfs = Lists.newArrayList();
		cdfs.add(Lists.newArrayList());
		cdfs.add(Lists.newArrayList());
		cdfs.add(Lists.newArrayList());

		for (int i = 0; i < lists.size(); i++) {
			for (int j = 0; j < lists.get(i).size(); j++) {
				System.out.println(
						String.format("%d,%d", lists.get(i).get(j), j));
			}
		}

		br.close();
	}

	public static double[] optimalHitRate(int capacities[]) throws Exception {
		String file = "/home/haoyu/Documents/cache/NovaSimulator/config/record_access_freq_records_100000_zipf_0.990000.txt";
		BufferedReader br = new BufferedReader(new FileReader(new File(file)));
		String line = null;
		long[] refs = new long[100000];
		int index = 0;
		double sum = 0;
		double[] hitRate = new double[capacities.length];
		while ((line = br.readLine()) != null) {
			refs[index] = Integer.parseInt(line.split(",")[2]);
			sum += refs[index];
			for (int i = 0; i < capacities.length; i++) {
				if (index < capacities[i]) {
					hitRate[i] += refs[index];
				}
			}
			index++;
		}

		for (int i = 0; i < capacities.length; i++) {
			hitRate[i] = hitRate[i] / (double) (sum);
		}
		br.close();
		return hitRate;
	}

	public static void main(String[] args) throws Exception {
		cdf();
//		
//		long records = 100000;
//		Random rand = new Random(0);
//		ZipfianGenerator zip = new ZipfianGenerator(records, rand);
//		int capacities[] = new int[] { 1, 10, 100, 1000, 10000, 100000 };
//		double[] optimalHitRate = optimalHitRate(capacities);
//		int maxLoop = 1000000;
//		int partitions = 10;
//
//		for (int i = 0; i < capacities.length; i++) {
//			double hits = 0;
//			PartitionLRUCache cache = new PartitionLRUCache(partitions,
//					capacities[i] / partitions, (int) records);
//			// warm up
//			for (int j = 0; j < capacities[i]; j++) {
//				int key = capacities[i] - j - 1;
//				cache.caches.get(cache.hash(key)).set(key, key, 1);
//			}
//			for (int j = 0; j < maxLoop; j++) {
//				int key = zip.nextValue().intValue();
//				int partition = cache.hash(key);
//
//				if (cache.caches.get(partition).get(key) == -1) {
//					boolean miss = true;
//					if (partition == 0) {
//						partition = (key % (cache.caches.size() - 1)) + 1;
//						if (cache.caches.get(partition).get(key) != -1) {
//							hits += 1;
//							miss = false;
//						}
//					}
//
//					// miss in both caches.
//					if (miss) {
//						cache.caches.get(partition).set(key, key, 1);
//					}
//				} else {
//					hits += 1;
//				}
//
//			}
//			System.out.println(String.format("%d:%f:%f", capacities[i],
//					hits / (double) (maxLoop), optimalHitRate[i]));
//		}

	}

}
