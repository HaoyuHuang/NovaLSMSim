package edu.usc.nova;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Comparator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import com.google.common.collect.Lists;

import edu.usc.distributions.NumberGenerator;
import edu.usc.distributions.UniformGenerator;
import edu.usc.distributions.ZipfianGenerator;

public class SubrangeSamplingSim {
	public static int numberOfMemTables;
	public static int numberOfKeys;
	public static int numberOfSubRanges;
	public static int maxSampleSize = 10000;

	public static Random r;
	public static NumberGenerator gen = null;
	public static List<SubRange> subranges = Lists.newArrayList();
	public static double totalNumberOfInserts = 0;
	public static double fairShare = 0.0;

	public static int numberOfTotalMajor = 0;
	public static int numberOfPerformedMajor = 0;

	public static class InternalKey {
		int key;
		long sequenceNumber;
	}

	static Comparator<InternalKey> comp = new Comparator<SubrangeSamplingSim.InternalKey>() {
		@Override
		public int compare(InternalKey o1, InternalKey o2) {
			if (o1.key == o2.key) {
				if (o1.sequenceNumber > o2.sequenceNumber) {
					return -1;
				} else {
					return 1;
				}
			} else if (o1.key < o2.key) {
				return -1;
			} else {
				return 1;
			}
		}
	};

	static int samplingRatio = 1;
	static LinkedList<Integer> sampleList = new LinkedList<>();
	static TreeMap<Integer, Integer> sampleKeyFreq = new TreeMap<>();
	static int maximumSamples = 0;

	public static class SubRange {
		int lower;
		int upper;
		double numberOfInserts;
		double currentShare;
	}

	static int binarySearch(int key) {
		int l = 0, r = subranges.size() - 1;
		while (l <= r) {
			int m = l + (r - l) / 2;

			// Check if x is present at mid
			if (key < subranges.get(m).lower) {
				r = m - 1;
			} else if (key > subranges.get(m).upper) {
				l = m + 1;
			} else {
				return m;
			}
		}
		// if we reach here, then element was
		// not present
		return -1;
	}

	public static boolean majorRebalance() {
		numberOfTotalMajor++;

		if (sampleKeyFreq.size() <= subranges.size() * 2) {
			return false;
		}
		assert sampleKeyFreq.size() > subranges.size() * 2;
		numberOfPerformedMajor++;

		int sharePerSubRange = sampleList.size() / subranges.size();
		int index = 0;
		int sum = 0;
		subranges.get(0).lower = 0;
		subranges.get(subranges.size() - 1).upper = numberOfKeys;
		for (Entry<Integer, Integer> entry : sampleKeyFreq.entrySet()) {
			sum += entry.getValue();
			if (sum >= sharePerSubRange) {
				subranges.get(index).upper = entry.getKey() + 1;
				subranges.get(index + 1).lower = entry.getKey() + 1;
				index++;
				if (index == subranges.size() - 1) {
					break;
				}
				sum = 0;
			}
		}

		int priorUpper = -1;
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			if (priorUpper != -1) {
				if (sr.lower < priorUpper) {
					sr.lower = priorUpper;
				}
			}

			if (sr.upper <= sr.lower) {
				sr.upper = sr.lower + 1;
				priorUpper = sr.upper;
			}
			sr.numberOfInserts = totalNumberOfInserts / subranges.size();
			sr.currentShare = sr.numberOfInserts / totalNumberOfInserts;
		}
		subranges.get(0).lower = 0;
		subranges.get(subranges.size() - 1).upper = numberOfKeys;
		return true;
	}

	public static void addKey(InternalKey ik) {
		int index = binarySearch(ik.key);
		SubRange subrange = null;
		if (index == -1) {
			if (subranges.size() < numberOfSubRanges) {
				if (subranges.isEmpty()) {
					subrange = new SubRange();
					subrange.lower = ik.key;
					subrange.upper = ik.key + 1;
					subranges.add(subrange);
				} else if (ik.key < subranges.get(0).lower) {
					subrange = new SubRange();
					subrange.lower = ik.key;
					subrange.upper = subranges.get(0).lower;
					subranges.add(0, subrange);
				} else {
					assert ik.key >= subranges.get(subranges.size() - 1).upper;
					subrange = new SubRange();
					subrange.lower = subranges.get(subranges.size() - 1).upper;
					subrange.upper = ik.key + 1;
					subranges.add(subrange);
				}
			} else {
				assert subranges.size() == numberOfSubRanges;
				if (ik.key < subranges.get(0).lower) {
					subranges.get(0).lower = ik.key;
					subrange = subranges.get(0);
				} else {
					assert ik.key >= subranges.get(subranges.size() - 1).upper;
					subranges.get(subranges.size() - 1).upper = ik.key + 1;
					subrange = subranges.get(subranges.size() - 1);
				}
			}
		} else {
			// split a range into two.
			subrange = subranges.get(index);
			if (ik.key != subrange.lower
					&& subranges.size() < numberOfSubRanges) {
				SubRange newSubrange = new SubRange();
				newSubrange.lower = ik.key;
				newSubrange.upper = subrange.upper;
				subrange.upper = newSubrange.lower;
				subranges.add(index + 1, newSubrange);
				subrange.numberOfInserts /= 2;
				newSubrange.numberOfInserts = subrange.numberOfInserts;
				newSubrange.currentShare = newSubrange.numberOfInserts
						/ totalNumberOfInserts;
				subrange = newSubrange;
			}
		}
		assert subrange != null;

		subrange.numberOfInserts += 1;
		totalNumberOfInserts += 1;
		subrange.currentShare = subrange.numberOfInserts / totalNumberOfInserts;

		double sum = 0;
		for (int i = 0; i < subranges.size(); i++) {
			sum += subranges.get(i).numberOfInserts;
		}
		if (Math.abs(sum - totalNumberOfInserts) >= 1) {
			assert sum == totalNumberOfInserts;
		}

		if (samplingRatio > 1 && r.nextInt(samplingRatio) == 0) {
			sampleList.add(ik.key);
			if (sampleList.size() > maximumSamples) {
				int first = sampleList.pollFirst();
				sampleKeyFreq.compute(first, (k, v) -> {
					assert v != null;
					if (v - 1 == 0) {
						return null;
					}
					return v - 1;
				});
			}
			sampleKeyFreq.compute(ik.key, (k, v) -> {
				if (v == null) {
					return 1;
				}
				return v + 1;
			});

		}

		if (subranges.size() != numberOfSubRanges) {
			return;
		}

		if (totalNumberOfInserts < 100000) {
			return;
		}

		// check.
		double unfairRanges = 0;
		for (int i = 0; i < numberOfSubRanges; i++) {
			SubRange sr = subranges.get(i);
			double diff = (sr.currentShare - fairShare) * 100.0 / fairShare;
			if (diff > 10 && sr.upper - sr.lower > 1) {
				unfairRanges += 1;
			}
		}
		String prefix = "";
		if (unfairRanges > 0) {
			prefix = String.format("major at %.2f%% of iterations:",
					(totalNumberOfInserts / iterations) * 100.0);
			if (!majorRebalance()) {
				return;
			}
			printRanges(prefix);
		} else {
			return;
		}
	}

	private static void printRanges(String prefix) {
		StringBuilder builder = new StringBuilder();
		builder.append(prefix);
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			builder.append(String.format("[%d:%d):%.2f,%f", sr.lower, sr.upper,
					sr.currentShare, sr.numberOfInserts));
			builder.append(",");
		}
		builder = builder.deleteCharAt(builder.length() - 1);
		System.out.println(builder.toString());
	}

	private static void printFinalRanges() {
		StringBuilder builder = new StringBuilder();
		builder.append("final:");
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			builder.append(String.format("[%d:%d):%.2f%%", sr.lower, sr.upper,
					(sr.numberOfInserts / totalNumberOfInserts) * 100.0));
			builder.append(",");
		}
		builder = builder.deleteCharAt(builder.length() - 1);
		System.out.println(builder.toString());
	}

	public static void printStats() {
		// Uniform.
		double diff = 0;
		if ("uniform".equals(dist)) {
			diff = UniformLoadImbalance();
		} else {
			diff = ZipfianLoadImbalance();
		}
		System.out.println(
				String.format("reorg,%d,%d,%.2f", numberOfPerformedMajor,
						numberOfTotalMajor - numberOfPerformedMajor, diff));
	}

	public static final String ZIPFIAN_FILE = "/tmp/zipfian";

	private static double ZipfianLoadImbalance() {
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
			assert key == numberOfKeys;
			List<SubRange> idealRanges = Lists.newArrayList();
			double remainingShare = sum;
			int remainingRanges = numberOfSubRanges;
			key = 0;
			while (remainingShare > 0 && remainingRanges > 0) {
				double fairShare = remainingShare / remainingRanges;
				double share = 0;
				int lower = key;
				int upper = key;
				while (share < fairShare && key < numberOfKeys) {
					share += ref[key];
					key++;
				}
				upper = key;
				SubRange range = new SubRange();
				range.lower = lower;
				range.upper = upper;
				range.currentShare = share / sum;
				remainingShare -= share;
				remainingRanges -= 1;
				idealRanges.add(range);

				if (key == numberOfKeys) {
					break;
				}
			}

			StringBuilder builder = new StringBuilder();
			builder.append("ideal:");
			for (int i = 0; i < idealRanges.size(); i++) {
				SubRange sr = idealRanges.get(i);
				builder.append(String.format("[%d:%d):%.2f%%", sr.lower,
						sr.upper, sr.currentShare * 100.0));
				builder.append(",");
			}
			builder = builder.deleteCharAt(builder.length() - 1);
			System.out.println(builder.toString());

			double diff = 0.0;
			for (int i = 0; i < idealRanges.size(); i++) {
				double actualShare = 0;
				SubRange as = idealRanges.get(i);
				for (int k = as.lower; k < as.upper; k++) {
					actualShare += ref[k];
				}
				actualShare /= sum;
				if (actualShare > fairShare) {
					diff += Math.abs(actualShare - fairShare);
				}
			}
			System.out.println("ideal-load:" + diff * 100.0);

			diff = 0.0;
			for (int i = 0; i < subranges.size(); i++) {
				double actualShare = 0;
				SubRange as = subranges.get(i);
				for (int k = as.lower; k < as.upper; k++) {
					actualShare += ref[k];
				}
				actualShare /= sum;
				if (actualShare > fairShare) {
					diff += Math.abs(actualShare - fairShare);
				}
			}
			diff *= 100;
			br.close();
			return diff;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return -1;
	}

	private static double UniformLoadImbalance() {
		double diff = 0.0;
		for (int i = 0; i < subranges.size(); i++) {
			SubRange actual = subranges.get(i);
			double actualShare = (double) (actual.upper - actual.lower)
					/ (double) numberOfKeys;
			if (actualShare > fairShare) {
				diff += Math.abs(actualShare - fairShare);
			}
		}
		diff *= 100;
		return diff;
	}

	public static long iterations = 100000000;
	public static String dist = "uniform";

	public static void run() {
		dist = "zipfian";
		numberOfMemTables = 128;
		numberOfKeys = 10000000;
		numberOfSubRanges = 32;
		r = new Random(0);
		fairShare = 1.0 / (double) numberOfSubRanges;
		samplingRatio = 1;
		maximumSamples = 100000;
		iterations = 100000000;
		System.out.println(fairShare);
		ZipfianLoadImbalance();

		NumberGenerator gen = null;
		if ("zipfian".equals(dist)) {
			gen = new ZipfianGenerator(numberOfKeys, r);
		} else {
			gen = new UniformGenerator(numberOfKeys, r);
		}
		for (long seq = 0; seq < iterations; seq++) {
			int key = gen.nextValue().intValue();
			InternalKey ik = new InternalKey();
			ik.key = key;
			ik.sequenceNumber = seq;
			addKey(ik);
		}
		printFinalRanges();
		printStats();
	}

	public static void main(String[] args) {
		run();

		dist = args[0];
		numberOfMemTables = Integer.parseInt(args[1]);
		numberOfKeys = Integer.parseInt(args[2]);
		numberOfSubRanges = Integer.parseInt(args[3]);
		iterations = Long.parseLong(args[4]);
		samplingRatio = Integer.parseInt(args[5]);
		maximumSamples = Integer.parseInt(args[6]);
		r = new Random(0);
		fairShare = 1.0 / (double) numberOfSubRanges;

		NumberGenerator gen = null;
		if ("zipfian".equals(dist)) {
			gen = new ZipfianGenerator(numberOfKeys, r);
		} else {
			gen = new UniformGenerator(numberOfKeys, r);
		}
		for (long seq = 0; seq < iterations; seq++) {
			int key = gen.nextValue().intValue();
			InternalKey ik = new InternalKey();
			ik.key = key;
			ik.sequenceNumber = seq;
			addKey(ik);
		}
		printFinalRanges();
		printStats();
	}

}
