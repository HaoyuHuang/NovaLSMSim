package edu.usc.nova;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;

import edu.usc.distributions.NumberGenerator;
import edu.usc.distributions.UniformGenerator;
import edu.usc.distributions.ZipfianGenerator;

public class SubrangeSim {
	public static int numberOfMemTables;
	public static int numberOfKeys;
	public static int numberOfSubRanges;
	public static int maxSampleSize = 10000;

	public static Random r;
	public static NumberGenerator gen = null;
	public static List<SubRange> subranges = Lists.newArrayList();
	public static double totalNumberOfInserts = 0;
	public static double totalNumberOfInsertsSinceLastMajor = 0;
	public static double fairShare = 0.0;

	public static int numberOfTotalMajor = 0;
	public static int numberOfTotalMinor = 0;
	public static int numberOfPerformedMajor = 0;
	public static int numberOfPerformedMinor = 0;

	public static class InternalKey {
		int key;
		long sequenceNumber;
	}

	static Comparator<InternalKey> comp = new Comparator<SubrangeSim.InternalKey>() {
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

	public static class SubRange {
		int lower;
		int upper;
		double numberOfInserts;
		double currentShare;
		double cumulativeNumberOfInserts = 0;

		public SubRange() {
			for (int i = 0; i < numberOfMemTables / numberOfSubRanges; i++) {
				memtables.add(new TreeSet<>(comp));
			}
		}

		List<TreeSet<InternalKey>> memtables = Lists.newArrayList();

		public void addKey(InternalKey ik) {
			int index = 0; // r.nextInt(memtables.size());
			if (memtables.get(index).size() == 16384) {
				memtables.get(index).clear();
			}
			memtables.get(index).add(ik);
		}

		int ninternalKeys = 0;

		public TreeMap<Integer, Integer> sampleKeys(int sampleInterval) {
			// sample my keys.
			ninternalKeys = 0;
			TreeMap<Integer, Integer> sortedMap = new TreeMap<>();
			Iterator<InternalKey> it = memtables.stream()
					.flatMap(TreeSet::stream).iterator();
			int i = -1;
			while (it.hasNext()) {
				InternalKey ik = it.next();
				i++;
				if (i % sampleInterval != 0) {
					continue;
				}

				if (samplingRatio > 1 && r.nextInt(samplingRatio) != 0) {
					continue;
				}

				ninternalKeys++;
				sortedMap.compute(ik.key, (k, v) -> {
					if (v == null) {
						return 1;
					}
					return v + 1;
				});
			}
			return sortedMap;
		}

		public int keys() {
			return memtables.stream().mapToInt(TreeSet::size).sum();
		}

		public TreeMap<Integer, Integer> sampleKeys() {
			// sample my keys.
			ninternalKeys = 0;
			TreeMap<Integer, Integer> sortedMap = new TreeMap<>();
			for (int i = 0; i < memtables.size(); i++) {
				for (InternalKey ik : memtables.get(i)) {
					if (samplingRatio > 1 && r.nextInt(samplingRatio) != 0) {
						continue;
					}

					ninternalKeys++;
					sortedMap.compute(ik.key, (k, v) -> {
						if (v == null) {
							return 1;
						}
						return v + 1;
					});
				}
			}
			return sortedMap;
		}
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

	public static double last_major_reorg_seq = 0;
	public static double last_minor_reorg_seq = 0;

	public static int majorRebalance() {
		TreeMap<Integer, Double> sortedMap = new TreeMap<>();
		double totalShares = 0.0;
		int minKeys = Integer.MAX_VALUE;
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			if (sr.keys() > 100) {
				minKeys = Math.min(subranges.get(i).keys(), minKeys);
			}
		}

		numberOfTotalMajor++;
		for (int i = 0; i < subranges.size(); i++) {
			// sample all keys
			int puts = subranges.get(i).keys();
			TreeMap<Integer, Integer> map;
			if (puts <= 100) {
				map = subranges.get(i).sampleKeys(1);
			} else {
				map = subranges.get(i).sampleKeys(puts / minKeys);
			}

			double insertionRatio = subranges.get(i).currentShare;
			for (Entry<Integer, Integer> entry : map.entrySet()) {
				sortedMap.compute(entry.getKey(), (k, v) -> {
					return entry.getValue().doubleValue() * insertionRatio;
				});
			}
		}
		for (Double value : sortedMap.values()) {
			totalShares += value;
		}

		System.out.println(sortedMap.size());

		if (sortedMap.size() <= subranges.size() * 2) {
			return 0;
		}

		printRanges("");

		assert sortedMap.size() > subranges.size() * 2;
		numberOfPerformedMajor++;
		last_major_reorg_seq = totalNumberOfInserts;

		double sharePerSubRange = totalShares / subranges.size();
		int index = 0;
		double sum = 0;
		subranges.get(0).lower = 0;
		subranges.get(subranges.size() - 1).upper = numberOfKeys;
		for (Entry<Integer, Double> entry : sortedMap.entrySet()) {
			sum += entry.getValue();
			totalShares -= entry.getValue();
			if (sum > sharePerSubRange) {
				subranges.get(index).upper = entry.getKey() + 1;
				subranges.get(index + 1).lower = entry.getKey() + 1;
				index++;
				if (index == subranges.size() - 1) {
					break;
				}
				sum = 0;
				sharePerSubRange = totalShares / (subranges.size() - index);
			}
		}

		int priorUpper = -1;
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			if (priorUpper != -1) {
				assert sr.lower > priorUpper;
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
			sr.ninternalKeys = 0;
		}
		subranges.get(0).lower = 0;
		subranges.get(subranges.size() - 1).upper = numberOfKeys;
		fixedMinor.clear();
		return sortedMap.size();
	}

	public static enum UpdateBoundary {
		LOWER, UPPER, NONE
	}

	public static Set<Integer> fixedMinor = Sets.newHashSet();

	public static boolean minorRebalanceHigherShare(int index) {
		fixedMinor.add(index);
		SubRange sr = subranges.get(index);
		assert sr.currentShare > fairShare;
		sr.ninternalKeys = 0;
		SubRange other = null;

		UpdateBoundary update_boundary = UpdateBoundary.NONE;
		for (int i = -1; i <= 1; i += 2) {
			if (index + i >= 0 && index + i < subranges.size()) {
				other = subranges.get(index + i);
				if (i == -1) {
					update_boundary = UpdateBoundary.LOWER;
				} else {
					update_boundary = UpdateBoundary.UPPER;
				}
				break;
			}
		}

		if (UpdateBoundary.NONE.equals(update_boundary)) {
			return false;
		}

		// Distribute the load across adjacent subranges.

		// higher share.
		TreeMap<Integer, Integer> samples = sr.sampleKeys();
		if (samples.isEmpty()) {
			return false;
		}
		System.out.println("!!!!!!!minor-" + index + "," + samples.size());
		last_minor_reorg_seq = totalNumberOfInserts;
		numberOfTotalMinor++;
		double inserts = (sr.currentShare - fairShare) * totalNumberOfInserts;
		assert inserts < sr.numberOfInserts;
		double removeShare = (double) ((sr.currentShare - fairShare)
				* sr.ninternalKeys);

		if (UpdateBoundary.LOWER.equals(update_boundary)) {
			int newLower = sr.lower;
			for (Entry<Integer, Integer> entry : samples.entrySet()) {
				if (removeShare <= 0) {
					break;
				}
				removeShare -= entry.getValue();
				newLower = entry.getKey();
			}
			if (newLower > sr.lower) {
				sr.lower = newLower;
				sr.numberOfInserts -= inserts;
				sr.currentShare = sr.numberOfInserts / totalNumberOfInserts;

				other.numberOfInserts += inserts;
				other.upper = sr.lower;
				other.currentShare = other.numberOfInserts
						/ totalNumberOfInserts;
				numberOfPerformedMinor++;
				return true;
			}
		} else {
			int newUpper = sr.upper;
			for (Entry<Integer, Integer> entry : samples.descendingMap()
					.entrySet()) {
				if (removeShare <= 0) {
					break;
				}
				removeShare -= entry.getValue();
				newUpper = entry.getKey();
			}
			if (newUpper < sr.upper) {
				sr.upper = newUpper + 1;
				sr.numberOfInserts -= inserts;
				sr.currentShare = sr.numberOfInserts / totalNumberOfInserts;

				other.numberOfInserts += inserts;
				other.lower = sr.upper;
				other.currentShare = other.numberOfInserts
						/ totalNumberOfInserts;
				numberOfPerformedMinor++;
				return true;
			}
		}
		return false;
	}

	public static boolean minorRebalanceHigherShareDistribute(int index) {
		fixedMinor.add(index);
		SubRange sr = subranges.get(index);
		assert sr.currentShare > fairShare;
		sr.ninternalKeys = 0;

		// Distribute the load across adjacent subranges.
		numberOfTotalMinor++;
		// higher share.
		TreeMap<Integer, Integer> samples = sr.sampleKeys();
		if (samples.size() < 100) {
			return false;
		}
		System.out.println(
				String.format("!!!!!!minor-%d-%d-%.2f", index, samples.size(),
						(sr.currentShare - fairShare) * 100.0 / fairShare));
		last_minor_reorg_seq = totalNumberOfInserts;
		double inserts = (sr.currentShare - fairShare) * totalNumberOfInserts;
		assert inserts < sr.numberOfInserts;
		double removeShare = (double) ((sr.currentShare - fairShare)
				* sr.ninternalKeys);
		if (index != 0 && index != subranges.size() - 1) {
			inserts /= 2;
			removeShare /= 2;
		}
		numberOfPerformedMinor++;
		// update lower
		if (index > 0) {
			int newLower = sr.lower;
			double removed_share = removeShare;
			Set<Integer> removed = Sets.newHashSet();
			for (Entry<Integer, Integer> entry : samples.entrySet()) {
				if (removed_share <= 0) {
					break;
				}
				removed_share -= entry.getValue();
				newLower = entry.getKey();
				removed.add(entry.getKey());
			}

			removed.forEach(k -> {
				samples.remove(k);
			});

			if (newLower > sr.lower) {
				sr.lower = newLower;
				sr.numberOfInserts -= inserts;
				sr.currentShare = sr.numberOfInserts / totalNumberOfInserts;
				SubRange other = subranges.get(index - 1);
				other.numberOfInserts += inserts;
				other.upper = sr.lower;
				other.currentShare = other.numberOfInserts
						/ totalNumberOfInserts;

			}
		}
		// update upper.
		if (index < subranges.size() - 1) {
			double removed_share = removeShare;
			int newUpper = sr.upper;
			for (Entry<Integer, Integer> entry : samples.descendingMap()
					.entrySet()) {
				if (removed_share <= 0) {
					break;
				}
				removed_share -= entry.getValue();
				newUpper = entry.getKey();
			}
			if (newUpper < sr.upper) {
				sr.upper = newUpper + 1;
				sr.numberOfInserts -= inserts;
				sr.currentShare = sr.numberOfInserts / totalNumberOfInserts;
				SubRange other = subranges.get(index + 1);
				other.numberOfInserts += inserts;
				other.lower = sr.upper;
				other.currentShare = other.numberOfInserts
						/ totalNumberOfInserts;
			}
		}
		return false;
	}

	public static double num_iteration_grow = 0;

	public static void addKey(InternalKey ik) {
		int index = binarySearch(ik.key);
		SubRange subrange = null;
		if (index == -1) {
			if (subranges.size() < numberOfSubRanges) {
				num_iteration_grow = totalNumberOfInserts;
				if (subranges.isEmpty()) {
					subrange = new SubRange();
					subrange.lower = ik.key;
					subrange.upper = ik.key + 1;
					subranges.add(subrange);
				} else {
					if (subranges.size() == 1 && subranges.get(0).upper
							- subranges.get(0).lower == 1) {
						// one subrange and it is a point.
						if (ik.key < subranges.get(0).lower) {
							subranges.get(0).lower = ik.key;
							subrange = subranges.get(0);
						} else {
							subranges.get(0).upper = ik.key + 1;
							subrange = subranges.get(0);
						}
					} else if (ik.key < subranges.get(0).lower) {
						subrange = new SubRange();
						subrange.lower = ik.key;
						subrange.upper = subranges.get(0).lower;
						subranges.add(0, subrange);
					} else {
						assert ik.key >= subranges
								.get(subranges.size() - 1).upper;
						subrange = new SubRange();
						subrange.lower = subranges
								.get(subranges.size() - 1).upper;
						subrange.upper = ik.key + 1;
						subranges.add(subrange);
					}

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
				num_iteration_grow = totalNumberOfInserts;
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
		subrange.addKey(ik);
		subrange.numberOfInserts += 1;
		subrange.cumulativeNumberOfInserts += 1;
		totalNumberOfInserts += 1;
		subrange.currentShare = subrange.numberOfInserts / totalNumberOfInserts;

		double sum = 0;
		for (int i = 0; i < subranges.size(); i++) {
			sum += subranges.get(i).numberOfInserts;
		}
		if (Math.abs(sum - totalNumberOfInserts) >= 1) {
			assert sum == totalNumberOfInserts;
		}

		if (subranges.size() != numberOfSubRanges) {
			return;
		}

		if (totalNumberOfInserts < 1000000) {
			return;
		}

		// check.
		double unfairRanges = 0;
		int mostUnfairRange = -1;
		double mostUnfair = 0;
		for (int i = 0; i < numberOfSubRanges; i++) {
			SubRange sr = subranges.get(i);
			double diff = (sr.currentShare - fairShare) * 100.0 / fairShare;
			if (Math.abs(diff) > 20 && sr.upper - sr.lower > 1) {
				unfairRanges += 1;
			}

			if (diff > 20 && sr.upper - sr.lower > 1) {
				if (diff > mostUnfair) {
					mostUnfairRange = i;
					mostUnfair = diff;
				}
			}
		}
		String prefix = "";
		if (unfairRanges / (double) (numberOfSubRanges) > 0.1) {
			int sampledKeys = 0;
			if (totalNumberOfInserts - last_major_reorg_seq > 1000000) {
				sampledKeys = majorRebalance();
			}
			prefix = String.format("major at %.2f%% of iterations samples %d:",
					(totalNumberOfInserts / iterations) * 100.0, sampledKeys);
			if (sampledKeys == 0) {
				if (enableMinor && mostUnfairRange != -1 && totalNumberOfInserts
						- last_minor_reorg_seq > 100000) {
					prefix = "minor-" + mostUnfairRange + ":";
					minorRebalanceHigherShareDistribute(mostUnfairRange);
					return;
				}
				return;
			}
		} else if (enableMinor && mostUnfairRange != -1) {
//			if (totalNumberOfInserts - last_minor_reorg_seq <= 100000) {
//				return;
//			}

			prefix = "minor-" + mostUnfairRange + ":";
			if (!minorRebalanceHigherShareDistribute(mostUnfairRange)) {
				return;
			}
		} else {
			return;
		}

		for (int i = 0; i < subranges.size(); i++) {
			subranges.get(i).cumulativeNumberOfInserts = 0;
			history.get(i).cumulativeNumberOfInserts = 0;
		}
		printRanges(prefix);
	}

	public static boolean enableMinor = true;

	private static void printRanges(String prefix) {
		StringBuilder builder = new StringBuilder();
		builder.append(prefix);
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			builder.append(String.format("[%d:%d):%.2f,%f", sr.lower, sr.upper,
					sr.currentShare, sr.numberOfInserts));
			builder.append(",");
			sr.ninternalKeys = 0;
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
			subranges.get(i).ninternalKeys = 0;
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
		double last_major_percent = ((double) last_major_reorg_seq
				/ (double) iterations) * 100.0;
		double last_minor_percent = ((double) last_minor_reorg_seq
				/ (double) iterations) * 100.0;

		System.out.println(String.format("reorg,%d,%d,%d,%d,%.2f,%.2f,%d,%.2f",
				numberOfPerformedMajor,
				numberOfTotalMajor - numberOfPerformedMajor,
				numberOfPerformedMinor,
				numberOfTotalMinor - numberOfPerformedMinor, last_major_percent,
				last_minor_percent, (int) num_iteration_grow, diff));
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
			double fs = sum / numberOfSubRanges;
			{
				double share = 0;
				int lower = 0;
				for (int i = 0; i < numberOfKeys; i++) {
					double percent = ref[i];
					if (share + percent > fs) {
						SubRange subrange = new SubRange();
						subrange.lower = lower;
						subrange.upper = i;
						subrange.currentShare = share / sum;
						lower = i;
						idealRanges.add(subrange);
						share = 0;
					}
					share += percent;
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

			double max = 0.0;
			double diff = 0.0;
			for (int i = 0; i < idealRanges.size(); i++) {
				double actualShare = 0;
				SubRange as = idealRanges.get(i);
				for (int k = as.lower; k < as.upper; k++) {
					actualShare += ref[k];
				}
				actualShare /= sum;
				max = Math.max(max, actualShare);
			}
			diff = max - fairShare;
			diff /= fairShare;
			System.out.println("ideal-load:" + (diff) * 100.0);

			diff = 0.0;
			max = 0.0;
			for (int i = 0; i < subranges.size(); i++) {
				double actualShare = 0;
				SubRange as = subranges.get(i);
				for (int k = as.lower; k < as.upper; k++) {
					actualShare += ref[k];
				}
				actualShare /= sum;
				max = Math.max(max, actualShare);
			}
			diff = max - fairShare;
			diff /= fairShare;
			diff *= 100;
			br.close();
			return diff;
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return -1;
	}

	private static void evenDist(double[] ref, double sum,
			List<SubRange> idealRanges) {
		int key;
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
	}

	private static double UniformLoadImbalance() {
		double diff = 0.0;
		double max = 0.0;
		for (int i = 0; i < subranges.size(); i++) {
			SubRange actual = subranges.get(i);
			double actualShare = (double) (actual.upper - actual.lower)
					/ (double) numberOfKeys;
			max = Math.max(actualShare, max);
		}
		diff = max - fairShare;
		diff /= fairShare;
		diff *= 100;
		return diff;
	}

	public static long iterations = 100000000;
	public static String dist = "uniform";

	public static List<SubRange> history = Lists.newArrayList();

	public static void run() {
		dist = "uniform";
		numberOfMemTables = 128;
		numberOfKeys = 10000000;
		numberOfSubRanges = 64;
		r = new Random(0);
		enableMinor = true;
		fairShare = 1.0 / numberOfSubRanges;
		samplingRatio = 10;
		iterations = 100000000;

		NumberGenerator gen = null;
		if ("zipfian".equals(dist)) {
			gen = new ZipfianGenerator(numberOfKeys, r);
		} else {
			gen = new UniformGenerator(numberOfKeys, r);
		}
		runSim(gen);
	}

	public static void main(String[] args) {
		run();

		dist = args[0];
		numberOfMemTables = Integer.parseInt(args[1]);
		numberOfKeys = Integer.parseInt(args[2]);
		numberOfSubRanges = Integer.parseInt(args[3]);
		iterations = Long.parseLong(args[4]);
		samplingRatio = Integer.parseInt(args[5]);
		enableMinor = Boolean.parseBoolean(args[6]);

		r = new Random(0);
		fairShare = 1.0 / numberOfSubRanges;

		NumberGenerator gen = null;
		if ("zipfian".equals(dist)) {
			gen = new ZipfianGenerator(numberOfKeys, r);
		} else {
			gen = new UniformGenerator(numberOfKeys, r);
		}
		runSim(gen);
	}

	private static void runSim(NumberGenerator gen) {
		for (int i = 0; i < numberOfSubRanges; i++) {
			SubRange newS = new SubRange();
			history.add(newS);
		}
		long last_seq = 0;
		int interval = 100000;
		for (long seq = 0; seq < iterations; seq++) {
			int key = gen.nextValue().intValue();
			InternalKey ik = new InternalKey();
			ik.key = key;
			ik.sequenceNumber = seq;
			addKey(ik);
			if (subranges.size() == numberOfSubRanges && seq % interval == 0) {
				double max_inbalance = 0;
				double sum = 0;
				for (int i = 0; i < numberOfSubRanges; i++) {
					double diff = subranges.get(i).cumulativeNumberOfInserts
							- history.get(i).cumulativeNumberOfInserts;
					sum += diff;
					history.get(i).cumulativeNumberOfInserts = subranges
							.get(i).cumulativeNumberOfInserts;
					max_inbalance = Math.max(max_inbalance, diff);
				}
				if (sum > 10000) {
					double fair = 1.0 / numberOfSubRanges;
					double perentage = max_inbalance / sum;
					double imb = (perentage - fair) / fair;
					System.out.println(
							String.format("seq,%d,%.2f", seq, imb * 100.0));
				}
			}
		}
		printFinalRanges();
		printStats();
	}

}
