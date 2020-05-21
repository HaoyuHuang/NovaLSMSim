package edu.usc.nova;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Comparator;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.TreeMap;
import java.util.TreeSet;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;

import edu.usc.distributions.NumberGenerator;
import edu.usc.distributions.UniformGenerator;
import edu.usc.distributions.ZipfianGenerator;

public class TinyRangeSim {

	public static final class LoadImbalanceMetric {
		double maximum_diff_percent = 0;
		double standard_deviation = 0;

		@Override
		public String toString() {
			return String.format("%.2f,%.8f", maximum_diff_percent,
					standard_deviation);
		}

		public static LoadImbalanceMetric compute(List<Double> loads) {
			double fair = 1.0 / loads.size();

			LoadImbalanceMetric metric = new LoadImbalanceMetric();
			double sum = 0.0;
			double highest_load = 0.0;
			double stdev = 0.0;

			for (int i = 0; i < loads.size(); i++) {
				double load = loads.get(i);
				sum += load;
				highest_load = Math.max(highest_load, load);

				double diff = (load - fair) * 100.0;
				stdev += Math.pow(diff, 2);
			}
			assert Math.abs(sum - 1.0) <= 0.01;
			metric.maximum_diff_percent = (highest_load - fair) * 100.0 / fair;
			metric.standard_deviation = Math.sqrt(stdev / (double) loads.size())
					/ 100.0;
			return metric;
		}
	}

	public static int numberOfMemTables;
	public static int numberOfKeys;
	public static int numberOfSubRanges;
	public static int numberOfTinyRangesPerSubRange;

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
	public static int numberOfPerformedMinorSampling = 0;
	public static int numberOfPerformedMinorDup = 0;

	public static class InternalKey {
		int key;
		long sequenceNumber;

		public int byteSize() {
			return String.valueOf(key).length() + 8 + 1 + 1024;
		}
	}

	static Comparator<InternalKey> comp = new Comparator<TinyRangeSim.InternalKey>() {
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

	public static class Range {
		int lower;
		int upper;
		double numberOfInserts;
		double currentShare;
		boolean duplicatedRange = false;
		int priorSubrangeId = -1;

		public int byteSize() {
			return String.valueOf(lower).length()
					+ String.valueOf(upper).length() + 8 + 8 + 8;
		}

		public void reset() {
			numberOfInserts = 0;
			currentShare = 0;
		}

		@Override
		public int hashCode() {
			final int prime = 31;
			int result = 1;
			result = prime * result + lower;
			result = prime * result + upper;
			return result;
		}

		@Override
		public boolean equals(Object obj) {
			if (this == obj)
				return true;
			if (obj == null)
				return false;
			if (getClass() != obj.getClass())
				return false;
			Range other = (Range) obj;
			if (lower != other.lower)
				return false;
			if (upper != other.upper)
				return false;
			return true;
		}
	}

	public static class SubRange {
		LinkedList<Range> tinyRanges = Lists.newLinkedList();
		List<TreeSet<InternalKey>> memtables = Lists.newArrayList();
		double cumulativeNumberOfInserts = 0;

		public SubRange() {
			for (int i = 0; i < numberOfMemTables / numberOfSubRanges; i++) {
				memtables.add(new TreeSet<>(comp));
			}
		}

		public int byteSize() {
			return tinyRanges.stream().mapToInt(Range::byteSize).sum();
		}

		@Override
		public String toString() {
			SubRange sr = this;
			StringBuilder builder = new StringBuilder();
			builder.append(String.format("[%d:%d):%.2f,%d,%d", sr.lower(),
					sr.upper(), sr.currentShare(), (int) sr.numberOfInserts(),
					sr.upper() - sr.lower()));
			builder.append(" tiny:");
			for (int j = 0; j < sr.tinyRanges.size(); j++) {
				builder.append(String.format("%d:[%d,%d),",
						sr.tinyRanges.get(j).upper - sr.tinyRanges.get(j).lower,
						sr.tinyRanges.get(j).lower,
						sr.tinyRanges.get(j).upper));
			}
			builder = builder.deleteCharAt(builder.length() - 1);
			return builder.toString();
		}

		public void clearPriorOnPushLeft() {
			for (int i = 0; i < tinyRanges.size() - 1; i++) {
				tinyRanges.get(i).priorSubrangeId = -1;
			}
		}

		public void clearPriorOnPushRight() {
			for (int i = 1; i < tinyRanges.size(); i++) {
				tinyRanges.get(i).priorSubrangeId = -1;
			}
		}

		public boolean freezed(int index) {
			if (tinyRanges.size() == 1) {
				return true;
			}
			int l = index - 1;
			int r = index + 1;
			boolean leftF = (l < 0)
					|| (tinyRanges.getFirst().priorSubrangeId == l)
					|| (subranges.get(l).duplicatedRange());
			boolean rightF = (r == subranges.size())
					|| (tinyRanges.getLast().priorSubrangeId == r)
					|| (subranges.get(r).duplicatedRange());
			return leftF && rightF;

		}

		public Range tinyRange(int key) {
			int index = -1;
			for (int i = 0; i < tinyRanges.size(); i++) {
				Range r = tinyRanges.get(i);
				if (r.lower <= key && key < r.upper) {
					index = i;
					break;
				}
			}
			assert index != -1;
			return tinyRanges.get(index);
		}

		public int lower() {
			return tinyRanges.getFirst().lower;
		}

		public int upper() {
			return tinyRanges.getLast().upper;
		}

		public boolean IsAPoint() {
			if (tinyRanges.size() > 1) {
				return false;
			}
			return tinyRanges.getFirst().upper
					- tinyRanges.getFirst().lower == 1;
		}

		public double currentShare() {
			return tinyRanges.stream().mapToDouble(i -> {
				return i.currentShare;
			}).sum();
		}

		public boolean duplicatedRange() {
			if (tinyRanges.getFirst().duplicatedRange) {
				assert tinyRanges.size() == 1;
				return true;
			}
			return false;
		}

		public double numberOfInserts() {
			return tinyRanges.stream().mapToDouble(i -> {
				return i.numberOfInserts;
			}).sum();
		}

		public void addKey(InternalKey ik) {
			boolean success = false;
			for (int i = 0; i < memtables.size(); i++) {
				TreeSet<InternalKey> memtable = memtables.get(i);
				if (memtable.size() < 16384) {
					memtable.add(ik);
					success = true;
					break;
				}
			}
			if (!success) {
				// All memtables are full.
				memtables.get(0).clear();
				memtables.get(0).add(ik);
			}
		}

		int ninternalKeys = 0;

		public TreeMap<Integer, Integer> sampleKeys(int sampleInterval) {
			// sample my keys.
			ninternalKeys = 0;
			TreeMap<Integer, Integer> sortedMap = new TreeMap<>();
			Iterator<InternalKey> it = memtables.stream()
					.flatMap(TreeSet::stream).iterator();
			int i = -1;
			int lower = lower();
			int upper = upper();
			while (it.hasNext()) {
				InternalKey ik = it.next();
				i++;
				if (ik.key < lower) {
					continue;
				}
				if (ik.key >= upper) {
					continue;
				}
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

	static void assertSubRanges() {
		List<Range> ranges = Lists.newArrayList();
		subranges.forEach(subrange -> {
			assert !subrange.tinyRanges.isEmpty();
			if (subrange.duplicatedRange()) {
				assert subrange.tinyRanges.size() == 1;
				assert subrange.IsAPoint();
			}

			ranges.addAll(subrange.tinyRanges);
		});

		if (ranges.isEmpty()) {
			return;
		}
		int priorLower = ranges.get(0).lower;
		int priorUpper = ranges.get(0).upper;
		boolean isPriorDup = ranges.get(0).duplicatedRange;
		for (int i = 1; i < ranges.size(); i++) {
			Range r = ranges.get(i);
			if (r.lower == priorLower) {
				assert r.upper == priorUpper;
				assert r.duplicatedRange = true;
				assert isPriorDup;
			} else {
				assert r.lower >= priorUpper;
				assert r.upper > r.lower;
				priorLower = r.lower;
				priorUpper = r.upper;
				isPriorDup = r.duplicatedRange;
			}
		}
	}

	static boolean rangeContains(int key, SubRange r) {
		if (key < r.lower()) {
			return false;
		}
		if (key >= r.upper()) {
			return false;
		}
		return true;
	}

	static int binarySearchWithDuplicate(int key) {
		int index = binarySearch(key);
		if (index == -1) {
			return index;
		}

		List<Integer> possibleRanges = Lists.newArrayList();
		possibleRanges.add(index);
		for (int i = 1; i <= subranges.size(); i++) {
			int low = index - i;
			int high = index + i;
			if (low < 0 || high >= subranges.size()) {
				break;
			}
			boolean contains = false;
			if (low >= 0 && rangeContains(key, subranges.get(low))) {
				possibleRanges.add(low);
				contains = true;
			}
			if (high < subranges.size()
					&& rangeContains(key, subranges.get(high))) {
				possibleRanges.add(high);
				contains = true;
			}
			if (!contains) {
				break;
			}
		}

		int i = r.nextInt(possibleRanges.size());
		index = possibleRanges.get(i);
		assert rangeContains(key, subranges.get(index));
		return index;
	}

	static int binarySearch(int key) {
		int l = 0, r = subranges.size() - 1;
		while (l <= r) {
			int m = l + (r - l) / 2;

			// Check if x is present at mid
			if (key < subranges.get(m).lower()) {
				r = m - 1;
			} else if (key >= subranges.get(m).upper()) {
				l = m + 1;
			} else {
				return m;
			}
		}
		return -1;
	}

	public static double last_major_reorg_seq = 0;
	public static double last_minor_reorg_seq = 0;

	public static int majorRebalance() {
		TreeMap<Integer, Double> userkeyFreqMap = new TreeMap<>();
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

			double insertionRatio = subranges.get(i).currentShare();
			for (Entry<Integer, Integer> entry : map.entrySet()) {
				userkeyFreqMap.compute(entry.getKey(), (k, v) -> {
					return entry.getValue().doubleValue() * insertionRatio;
				});
			}
		}
		if (userkeyFreqMap.size() <= numberOfSubRanges * 2) {
			return 0;
		}

		assert userkeyFreqMap.size() > numberOfSubRanges * 2;

		numberOfPerformedMajor++;
		last_major_reorg_seq = totalNumberOfInserts;
		last_minor_reorg_seq = totalNumberOfInserts;
		for (Double value : userkeyFreqMap.values()) {
			totalShares += value;
		}

		subranges.clear();
		// First, construct subranges with each subrange containing one tiny
		// range.
		List<Range> tmp_subranges = constructRanges(userkeyFreqMap, totalShares,
				0, numberOfKeys, numberOfSubRanges, true);
		tmp_subranges.forEach(r -> {
			r.reset();
		});
		// Second, break each subrange that contains more than one value into
		// alpha tiny ranges.
		for (int i = 0; i < tmp_subranges.size(); i++) {
			TreeMap<Integer, Double> subUserkeyFreqMap = Maps.newTreeMap();
			int lower = tmp_subranges.get(i).lower;
			int upper = tmp_subranges.get(i).upper;
			List<Range> tinyRanges = Lists.newArrayList(tmp_subranges.get(i));
			if (upper - lower > 1) {
				double subTotalShare = 0;
				for (Entry<Integer, Double> it : userkeyFreqMap.entrySet()) {
					if (it.getKey() < lower) {
						continue;
					}
					if (it.getKey() >= upper) {
						continue;
					}
					subTotalShare += it.getValue();
					subUserkeyFreqMap.put(it.getKey(), it.getValue());
				}
				tinyRanges = constructRanges(subUserkeyFreqMap, subTotalShare,
						lower, upper, numberOfTinyRangesPerSubRange, false);
				tinyRanges.forEach(r -> {
					r.reset();
				});
			}

			SubRange sr = new SubRange();
			sr.tinyRanges.addAll(tinyRanges);
			subranges.add(sr);
		}

		totalNumberOfInsertsSinceLastMajor = 0;
		return userkeyFreqMap.size();
	}

	private static List<Range> constructRanges(
			TreeMap<Integer, Double> sortedMap, double totalShares, int lower,
			int upper, int numberOfRanges, boolean isConstructingSubRanges) {
		assert upper - lower > 1;
		assert numberOfRanges > 1;
		List<Range> tmp_subranges = Lists.newArrayList();
		double sharePerSubRange = totalShares / numberOfRanges;
		double fairShare = totalShares / numberOfRanges;
		double total = totalShares;
		double currentRate = 0;

		Range currentRange = new Range();
		currentRange.lower = lower;
		for (Entry<Integer, Double> entry : sortedMap.entrySet()) {
			assert entry.getKey() >= lower;
			assert entry.getKey() < upper;

			double rate = entry.getValue();
			if (rate >= fairShare && isConstructingSubRanges) {
				System.out.println(String.format("hot key %d:%.2f:%.2f",
						entry.getKey(), rate / total, fairShare / total));
				// close the current subrange.
				if (currentRange.lower < entry.getKey()) {
					currentRange.upper = entry.getKey();
					currentRange.currentShare = currentRate / total;
					tmp_subranges.add(currentRange);
					currentRange = new Range();
				}

				int nDuplicates = (int) Math.ceil(rate / fairShare);
				for (int i = 0; i < nDuplicates; i++) {
					currentRange.duplicatedRange = true;
					currentRange.lower = entry.getKey();
					currentRange.upper = entry.getKey() + 1;
					currentRange.currentShare = (rate / nDuplicates) / total;
					tmp_subranges.add(currentRange);
					currentRange = new Range();

				}
				currentRange.lower = entry.getKey() + 1;
				totalShares -= entry.getValue();
				currentRate = 0;
				sharePerSubRange = totalShares
						/ (numberOfRanges - tmp_subranges.size());
				continue;
			}

			if (currentRate + rate > sharePerSubRange) {
				if (currentRange.lower == entry.getKey()) {
					currentRate += rate;
					currentRange.currentShare = currentRate / total;
					currentRange.upper = entry.getKey() + 1;
					currentRange.currentShare = currentRate / total;
					tmp_subranges.add(currentRange);
					currentRange = new Range();

					currentRange.lower = entry.getKey() + 1;
					if (tmp_subranges.size() + 1 == numberOfRanges) {
						break;
					}
					currentRate = 0;
					totalShares -= rate;
					sharePerSubRange = totalShares
							/ (numberOfRanges - tmp_subranges.size());
					continue;
				} else {
					currentRange.currentShare = currentRate / total;
					currentRange.upper = entry.getKey();

					tmp_subranges.add(currentRange);
					currentRange = new Range();

					currentRange.lower = entry.getKey();
					if (tmp_subranges.size() + 1 == numberOfRanges) {
						break;
					}
					currentRate = 0;
					sharePerSubRange = totalShares
							/ (numberOfRanges - tmp_subranges.size());
				}
			}
			currentRate += rate;
			totalShares -= rate;
		}

		if (isConstructingSubRanges) {
			currentRange.currentShare = currentRate / total;
			tmp_subranges.add(currentRange);
			assert tmp_subranges.size() == numberOfRanges;
		} else {
			if (currentRange.lower < upper) {
				currentRange.currentShare = currentRate / total;
				tmp_subranges.add(currentRange);
			}
			assert tmp_subranges.size() <= numberOfRanges;
		}

		tmp_subranges.get(0).lower = lower;
		tmp_subranges.get(tmp_subranges.size() - 1).upper = upper;
		return tmp_subranges;
	}

	public static enum UpdateBoundary {
		LOWER, UPPER, NONE
	}

	public static void moveShare(int index) {
		SubRange sr = subranges.get(index);
		int lower = sr.lower();
		int nDuplicates = 0;
		int start = -1;
		int end = -1;
		for (int i = 0; i < subranges.size(); i++) {
			SubRange r = subranges.get(i);
			if (!r.duplicatedRange()) {
				continue;
			}
			assert r.tinyRanges.size() == 1;
			if (r.lower() != lower) {
				continue;
			}
			end = i;
			if (start == -1) {
				start = i;
			}
		}

		nDuplicates = end - start;
		double share = sr.numberOfInserts() / nDuplicates;
		for (int i = start; i <= end; i++) {
			SubRange r = subranges.get(i);
			if (nDuplicates == 1) {
				r.tinyRanges.getFirst().duplicatedRange = false;
			}
			r.tinyRanges.getFirst().numberOfInserts += share;
			r.tinyRanges.getFirst().currentShare = r.tinyRanges
					.get(0).numberOfInserts
					/ totalNumberOfInsertsSinceLastMajor;
		}
	}

	public static boolean destroyDuplicates(int index, boolean force) {
		SubRange sr = subranges.get(index);
		if (!sr.duplicatedRange()) {
			return false;
		}

		if (force) {
			moveShare(index);
			subranges.remove(index);
			return true;
		}

		double percent = sr.currentShare() / fairShare;
		if (percent >= 0.5) {
			return false;
		}
		numberOfPerformedMinorDup++;
		System.out.println("Destroy subrange " + index);
		// destroy this duplicate.
		moveShare(index);
		subranges.remove(index);
		return true;
	}

	public static boolean minorRebalanceDuplicate(int index) {
		SubRange sr = subranges.get(index);
		sr.ninternalKeys = 0;
		if (sr.tinyRanges.size() != 1) {
			return false;
		}
		int lower = sr.tinyRanges.getFirst().lower;
		int upper = sr.tinyRanges.getFirst().upper;

		if (upper - lower != 1) {
			return false;
		}
		if (sr.currentShare() <= 1.5 * fairShare) {
			return false;
		}
		int nDuplicates = (int) Math.floor(sr.currentShare() / fairShare);
		if (nDuplicates == 0) {
			return false;
		}

		last_minor_reorg_seq = totalNumberOfInserts;
		numberOfPerformedMinorDup++;

		// Create new duplicate subranges.
		sr.tinyRanges.getFirst().duplicatedRange = true;
		double remainingSum = sr.numberOfInserts();
		for (int i = 0; i < nDuplicates; i++) {
			SubRange newSR = new SubRange();
			Range tinyRange = new Range();

			tinyRange.lower = lower;
			tinyRange.upper = upper;
			tinyRange.numberOfInserts = sr.numberOfInserts()
					/ (nDuplicates + 1);
			remainingSum -= tinyRange.numberOfInserts;
			tinyRange.currentShare = tinyRange.numberOfInserts
					/ totalNumberOfInsertsSinceLastMajor;
			tinyRange.duplicatedRange = true;
			newSR.tinyRanges.add(tinyRange);
			subranges.add(index, newSR);
		}
		sr.tinyRanges.getFirst().numberOfInserts = remainingSum;
		sr.tinyRanges.getFirst().currentShare = remainingSum
				/ totalNumberOfInsertsSinceLastMajor;

		// Remove subranges if the number of subranges exceeds max.
		// For each removed subrange, move its tiny ranges to its neighboring
		// subranges.
		while (subranges.size() > numberOfSubRanges) {
			// remove the subrange that has the lowest insertion rate.
			double minShare = Double.MAX_VALUE;
			int minRangeId = 0;
			for (int i = 0; i < subranges.size(); i++) {
				SubRange minSR = subranges.get(i);
				// Skip the new subranges.
				if (minSR.tinyRanges.size() == 1) {
					if (minSR.tinyRanges.getFirst().lower == lower) {
						continue;
					}
				}
				if (minSR.currentShare() < minShare) {
					minShare = minSR.currentShare();
					minRangeId = i;
				}
			}

			if (subranges.get(minRangeId).duplicatedRange()) {
				destroyDuplicates(minRangeId, true);
				continue;
			}

			int left = minRangeId - 1;
			int right = minRangeId + 1;
			boolean mergeLeft = true;
			if (left >= 0 && right < subranges.size()) {
				if (subranges.get(left).currentShare() < subranges.get(right)
						.currentShare()) {
					mergeLeft = true;
				} else {
					mergeLeft = false;
				}
			} else if (left >= 0) {
				mergeLeft = true;
			} else {
				mergeLeft = false;
			}

			if (mergeLeft && subranges.get(left).duplicatedRange()) {
				destroyDuplicates(left, true);
			} else if (!mergeLeft && subranges.get(right).duplicatedRange()) {
				destroyDuplicates(right, true);
			} else {
				int nranges = subranges.get(minRangeId).tinyRanges.size();
				int pushedRanges = pushTinyRanges(minRangeId, false);
				assert pushedRanges == nranges;
			}
		}
		System.out.println("minor duplicate subrange-" + index);
		return true;
	}

	private static int pushTinyRanges(int subrangeId,
			boolean stopWhenBelowFair) {
		// Push its tiny ranges to its neighbors.
		int left = subrangeId - 1;
		int right = subrangeId + 1;

		SubRange leftSR = null;
		SubRange rightSR = null;
		if (left >= 0) {
			leftSR = subranges.get(left);
		}
		if (right < subranges.size()) {
			rightSR = subranges.get(right);
		}
		int moved = 0;
		SubRange minSR = subranges.get(subrangeId);
		// move to left
		// move the remaining
		double leftShare = Double.MAX_VALUE;
		double rightShare = Double.MAX_VALUE;
		if (leftSR != null && !leftSR.duplicatedRange()) {
			leftShare = leftSR.currentShare();
		}
		if (rightSR != null && !rightSR.duplicatedRange()) {
			rightShare = rightSR.currentShare();
		}
		if (leftShare != Double.MAX_VALUE || rightShare != Double.MAX_VALUE) {
			while (!minSR.tinyRanges.isEmpty()) {
				Range first = minSR.tinyRanges.getFirst();
				Range last = minSR.tinyRanges.getLast();
				boolean pushLeft = false;

				if (leftShare + first.currentShare < rightShare
						+ last.currentShare) {
					pushLeft = true;
				}

				if (stopWhenBelowFair) {
					if (first.priorSubrangeId == left) {
						// the left pushed this tiny range to me in the
						// last reorg.
						if (rightShare != Double.MAX_VALUE) {
							pushLeft = false;
						} else {
							// I cannot push anything to the right.
						}

					}
					if (last.priorSubrangeId == right) {
						// the right pushed this tiny range to me in the
						// last reorg.
						if (leftShare != Double.MAX_VALUE) {
							pushLeft = true;
						} else {
							// I cannot push anything to the left.
						}

					}

					if (pushLeft) {
						if (minSR.currentShare()
								- first.currentShare < fairShare && moved > 0) {
							return moved;
						}

						if (first.priorSubrangeId == left) {
							first.priorSubrangeId = -1;
							return moved;

						}
						first.priorSubrangeId = subrangeId;
					} else {
						if (minSR.currentShare() - last.currentShare < fairShare
								&& moved > 0) {
							return moved;
						}
						if (last.priorSubrangeId == right) {
							last.priorSubrangeId = -1;
							return moved;
						}
						last.priorSubrangeId = subrangeId;
					}
				}

				if (pushLeft) {
					// move to left.
					leftSR.tinyRanges.addLast(minSR.tinyRanges.removeFirst());
					leftShare = leftSR.currentShare();
					leftSR.clearPriorOnPushLeft();
				} else {
					// move to right.
					rightSR.tinyRanges.addFirst(minSR.tinyRanges.removeLast());
					rightShare = rightSR.currentShare();
					rightSR.clearPriorOnPushRight();
				}
				moved++;
			}
		}
		if (minSR.tinyRanges.isEmpty()) {
			subranges.remove(subrangeId);
		}
		return moved;
	}

	private static int pullTinyRanges(int subrangeId) {
		// Pull tiny ranges from its neighbors.
		int left = subrangeId - 1;
		int right = subrangeId + 1;
		int moved = 0;
		SubRange leftSR = null;
		SubRange rightSR = null;
		if (left > 0) {
			leftSR = subranges.get(left);
		}
		if (right < subranges.size()) {
			rightSR = subranges.get(right);
		}

		SubRange minSR = subranges.get(subrangeId);
		double leftShare = Double.MAX_VALUE;
		double rightShare = Double.MAX_VALUE;
		if (leftSR != null && !leftSR.duplicatedRange()) {
			leftShare = leftSR.currentShare();
		}
		if (rightSR != null && !rightSR.duplicatedRange()) {
			rightShare = rightSR.currentShare();
		}
		if (leftShare != Double.MAX_VALUE || rightShare != Double.MAX_VALUE) {
			while (true) {
				boolean pullRight = true;
				double pullRightShare = 0;
				double pullLeftShare = 0;

				if (rightShare != Double.MAX_VALUE) {
					pullRightShare = rightSR.tinyRanges.getFirst().currentShare;
				}

				if (leftShare != Double.MAX_VALUE) {
					pullLeftShare = leftSR.tinyRanges.getLast().currentShare;
				}

				if (rightShare - pullRightShare > leftShare - pullLeftShare) {
					pullRight = true;
				} else {
					pullRight = false;
				}

				if (pullRight) {
					// pull from right.
					if (rightSR.tinyRanges.size() == 1) {
						return moved;
					}
					if (minSR.currentShare() + pullRightShare > fairShare
							&& moved > 0) {
						return moved;
					}

					minSR.tinyRanges.addLast(rightSR.tinyRanges.removeFirst());
					rightShare = rightSR.currentShare();
				} else {
					// pull from right.
					if (leftSR.tinyRanges.size() == 1) {
						return moved;
					}
					if (minSR.currentShare() + pullLeftShare > fairShare
							&& moved > 0) {
						return moved;
					}
					minSR.tinyRanges.addFirst(leftSR.tinyRanges.removeLast());
					leftShare = leftSR.currentShare();
				}

				moved++;
			}
		}
		return moved;
	}

	public static void minorRebalance(int index) {
		last_minor_reorg_seq = totalNumberOfInserts;
		numberOfTotalMinor++;

		SubRange sr = subranges.get(index);
		assert sr.currentShare() > fairShare;

		if (sr.upper() - sr.lower() > 100 && sr.tinyRanges.size() > 1) {
			double unfairTinyRanges = 0;
			double fair_share = 1.0 / sr.tinyRanges.size();
			for (int i = 0; i < sr.tinyRanges.size(); i++) {
				Range r = sr.tinyRanges.get(i);
				double diff = (r.currentShare - fair_share) * 100.0
						/ fair_share;
				if (Math.abs(diff) > 20) {
					unfairTinyRanges += 1;
				}
			}

			if (unfairTinyRanges / sr.tinyRanges.size() > 0.1) {
				double beforeShare = sr.currentShare();
				double beforeInserts = sr.numberOfInserts();
				numberOfPerformedMinorSampling++;
				TreeMap<Integer, Integer> map = sr.sampleKeys(1);
				TreeMap<Integer, Double> userkeyFreqMap = Maps.newTreeMap();
				double insertionRatio = sr.currentShare();
				double totalShare = 0.0;
				for (Entry<Integer, Integer> entry : map.entrySet()) {
					userkeyFreqMap.put(entry.getKey(),
							entry.getValue().doubleValue() * insertionRatio);
					totalShare += entry.getValue().doubleValue()
							* insertionRatio;
				}

				if (userkeyFreqMap.size() > numberOfTinyRangesPerSubRange * 2) {
					List<Range> ranges = constructRanges(userkeyFreqMap,
							totalShare, sr.lower(), sr.upper(),
							numberOfTinyRangesPerSubRange, false);
					double sum = 0;
					for (int i = 0; i < ranges.size(); i++) {
						Range r = ranges.get(i);
						r.numberOfInserts = r.currentShare * beforeInserts;
						sum += r.numberOfInserts;
					}
					ranges.get(
							ranges.size() - 1).numberOfInserts += beforeInserts
									- sum;
					for (int i = 0; i < ranges.size(); i++) {
						Range r = ranges.get(i);
						r.currentShare = r.numberOfInserts
								/ totalNumberOfInsertsSinceLastMajor;
					}

					ranges.get(0).priorSubrangeId = sr.tinyRanges
							.get(0).priorSubrangeId;
					ranges.get(ranges.size()
							- 1).priorSubrangeId = sr.tinyRanges.get(
									sr.tinyRanges.size() - 1).priorSubrangeId;

					sr.tinyRanges.clear();
					sr.tinyRanges.addAll(ranges);
					assert Math
							.abs(beforeInserts - sr.numberOfInserts()) < 0.01;
					assert Math.abs(beforeShare - sr.currentShare()) < 0.01;
				}
			}
		}

		if (minorRebalancePush(index)) {
			numberOfPerformedMinor++;
		}
	}

	public static boolean minorRebalancePush(int index) {
		SubRange sr = subranges.get(index);
		double totalRemoveInserts = (sr.currentShare() - fairShare)
				* totalNumberOfInsertsSinceLastMajor;

		if (totalRemoveInserts >= sr.numberOfInserts()) {
			return false;
		}

		if (sr.tinyRanges.size() == 1) {
			return false;
		}

		assert sr.currentShare() > fairShare;

		sr.ninternalKeys = 0;
		// Distribute the load across adjacent subranges.
		last_minor_reorg_seq = totalNumberOfInserts;

		System.out.println(index + " push minor before " + sr.toString());
		// push tiny ranges.
		if (pushTinyRanges(index, true) == 0) {
			return false;
		}
		System.out.println(index + " push minor after " + sr.toString());
		return true;
	}

	public static boolean minorRebalancePull(int index) {
		SubRange sr = subranges.get(index);
		double totalNeededInserts = (fairShare - sr.currentShare())
				* totalNumberOfInsertsSinceLastMajor;

		if (totalNeededInserts >= 0) {
			return false;
		}

		assert sr.currentShare() < fairShare;

		sr.ninternalKeys = 0;
		// Distribute the load across adjacent subranges.
		numberOfTotalMinor++;

		System.out.println(index + " pull minor before " + sr.toString());
		// push tiny ranges.
		if (pullTinyRanges(index) == 0) {
			return false;
		}
		System.out.println(index + " pull minor after " + sr.toString());
		return true;
	}

	public static double num_iteration_grow = 0;
	public static double minor_reorg_interval = 100000;
	public static double major_reorg_interval = 1000000;

	public static void addKey(InternalKey ik) {
//		assertSubRanges();

		int index = binarySearchWithDuplicate(ik.key);
		SubRange subrange = constructSubranges(ik, index);
		assert subrange != null;
		subrange.addKey(ik);
		subrange.cumulativeNumberOfInserts += 1;
		totalNumberOfInserts += 1;
		totalNumberOfInsertsSinceLastMajor += 1;

		Range tinyRange = subrange.tinyRange(ik.key);
		tinyRange.numberOfInserts += 1;
		tinyRange.currentShare = tinyRange.numberOfInserts
				/ totalNumberOfInsertsSinceLastMajor;

		double sum = 0;
		for (int i = 0; i < subranges.size(); i++) {
			sum += subranges.get(i).numberOfInserts();
		}
		if (Math.abs(sum - totalNumberOfInsertsSinceLastMajor) >= 1) {
			assert sum == totalNumberOfInsertsSinceLastMajor;
		}

		if (subranges.size() != numberOfSubRanges) {
			return;
		}

		if (totalNumberOfInserts < 1000000) {
			return;
		}

		// first try destroy/duplicate subranges with a single value.
		if (totalNumberOfInserts
				- last_minor_reorg_seq > minor_reorg_interval) {
			int pivot = 0;
			boolean success = false;
			while (pivot < subranges.size()) {
				SubRange sr = subranges.get(pivot);
				if (destroyDuplicates(pivot, false)) {
					// no need to update pivot since it already points to the
					// next range.
					success = true;
				} else if (minorRebalanceDuplicate(pivot)) {
					// go to the next subrange.
					for (int i = 0; i < subranges.size(); i++) {
						SubRange other = subranges.get(i);
						if (other.equals(sr)) {
							pivot = i;
						}
					}
					pivot++;
					success = true;
				} else {
					pivot++;
				}
			}
			if (success) {
				last_minor_reorg_seq = totalNumberOfInserts;
			}
		}

		double unfairRanges = 0;
		int mostUnfairRange = -1;
		double mostUnfair = 0;
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			double diff = (sr.currentShare() - fairShare) * 100.0 / fairShare;

			if (Math.abs(diff) > 20 && !sr.IsAPoint()) {
				unfairRanges += 1;
			}

			if (diff > 20
					&& totalNumberOfInserts
							- last_minor_reorg_seq > minor_reorg_interval
					&& !sr.IsAPoint()) {
				if (diff > mostUnfair) {
					mostUnfairRange = i;
					mostUnfair = diff;
				}
			}
		}
		String prefix = "";
		if (unfairRanges / (double) (numberOfSubRanges) > 0.1) {
			int sampledKeys = 0;
			if (totalNumberOfInserts
					- last_major_reorg_seq > major_reorg_interval) {
				sampledKeys = majorRebalance();
			}
			prefix = String.format("major at %.2f%% of iterations samples %d:",
					(totalNumberOfInserts / iterations) * 100.0, sampledKeys);
			if (sampledKeys == 0) {
				if (enableMinor && mostUnfairRange != -1) {
					prefix = "minor-" + mostUnfairRange + ":";
					minorRebalance(mostUnfairRange);
					return;
				}
				return;
			}
			printRanges(prefix);
		} else if (enableMinor && mostUnfairRange != -1) {
			prefix = "minor-" + mostUnfairRange + ":";
			minorRebalance(mostUnfairRange);
		} else {
			return;
		}

		for (int i = 0; i < subranges.size(); i++) {
			subranges.get(i).cumulativeNumberOfInserts = 0;
			history.get(i).cumulativeNumberOfInserts = 0;
		}
	}

	private static SubRange constructSubranges(InternalKey ik, int index) {
		SubRange subrange = null;
		if (index == -1) {
			if (subranges.size() < numberOfSubRanges) {
				num_iteration_grow = totalNumberOfInserts;
				if (subranges.isEmpty()) {
					subrange = new SubRange();
					Range r = new Range();
					subrange.tinyRanges.add(r);
					r.lower = ik.key;
					r.upper = ik.key + 1;
					subranges.add(subrange);
				} else {
					if (subranges.size() == 1 && subranges.get(0).IsAPoint()) {
						// one subrange and it is a point.
						if (ik.key < subranges.get(0).lower()) {
							subranges.get(0).tinyRanges.get(0).lower = ik.key;
							subrange = subranges.get(0);
						} else {
							subranges.get(0).tinyRanges.get(0).upper = ik.key
									+ 1;
							subrange = subranges.get(0);
						}
					} else if (ik.key < subranges.get(0).lower()) {
						subrange = new SubRange();
						Range r = new Range();
						r.lower = ik.key;
						r.upper = subranges.get(0).lower();
						subrange.tinyRanges.add(r);
						subranges.add(0, subrange);
					} else {
						assert ik.key >= subranges.get(subranges.size() - 1)
								.upper();
						subrange = new SubRange();
						Range r = new Range();
						subrange.tinyRanges.add(r);
						r.lower = subranges.get(subranges.size() - 1).upper();
						r.upper = ik.key + 1;
						subranges.add(subrange);
					}
				}
			} else {
				assert subranges.size() == numberOfSubRanges;
				if (ik.key < subranges.get(0).lower()) {
					subranges.get(0).tinyRanges.get(0).lower = ik.key;
					subrange = subranges.get(0);
				} else {
					SubRange last = subranges.get(subranges.size() - 1);
					assert ik.key >= last.upper();
					last.tinyRanges
							.get(last.tinyRanges.size() - 1).upper = ik.key + 1;
					subrange = subranges.get(subranges.size() - 1);
				}
			}
		} else {
			// split a range into two.
			subrange = subranges.get(index);
			if (ik.key != subrange.lower()
					&& subranges.size() < numberOfSubRanges) {
				num_iteration_grow = totalNumberOfInserts;
				SubRange newSubrange = new SubRange();
				Range r = new Range();
				newSubrange.tinyRanges.add(r);
				r.lower = ik.key;
				r.upper = subrange.upper();

				Range last = subrange.tinyRanges
						.get(subrange.tinyRanges.size() - 1);
				last.upper = newSubrange.lower();
				last.numberOfInserts /= 2;
				last.currentShare = last.numberOfInserts
						/ totalNumberOfInsertsSinceLastMajor;
				r.numberOfInserts = last.numberOfInserts;
				r.currentShare = r.numberOfInserts
						/ totalNumberOfInsertsSinceLastMajor;
				subranges.add(index + 1, newSubrange);
				subrange = newSubrange;
			}
		}
		return subrange;
	}

	public static boolean enableMinor = true;

	private static void printRanges(String prefix) {
		StringBuilder builder = new StringBuilder();
		builder.append(prefix);
		builder.append("\n");
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			builder.append(String.format("%d [%d:%d):%.2f,%d,%d", i, sr.lower(),
					sr.upper(), sr.currentShare(), (int) sr.numberOfInserts(),
					sr.upper() - sr.lower()));
			builder.append(" tiny:");
			for (int j = 0; j < sr.tinyRanges.size(); j++) {
				builder.append(
						String.format("[%d,%d),", sr.tinyRanges.get(j).lower,
								sr.tinyRanges.get(j).upper));
			}
			builder = builder.deleteCharAt(builder.length() - 1);
			builder.append("\n");
			sr.ninternalKeys = 0;
		}
		builder = builder.deleteCharAt(builder.length() - 1);
		System.out.println(builder.toString());
	}

	private static void printFinalRanges() {
		StringBuilder builder = new StringBuilder();
		builder.append("final:\n");
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			builder.append(String.format("%d [%d:%d):%.2f%%,%d", i, sr.lower(),
					sr.upper(),
					(sr.numberOfInserts() / totalNumberOfInsertsSinceLastMajor)
							* 100.0,
					sr.upper() - sr.lower()));

			builder.append(" tiny:");
			for (int j = 0; j < sr.tinyRanges.size(); j++) {
				builder.append(
						String.format("[%d,%d),", sr.tinyRanges.get(j).lower,
								sr.tinyRanges.get(j).upper));
			}
			builder = builder.deleteCharAt(builder.length() - 1);

			builder.append("\n");
			subranges.get(i).ninternalKeys = 0;
		}
		builder = builder.deleteCharAt(builder.length() - 1);
		System.out.println(builder.toString());
	}

	public static void printStats() {
		// Uniform.
		LoadImbalanceMetric diff;
		if ("uniform".equals(dist)) {
			diff = UniformLoadImbalance();
		} else {
			diff = ZipfianLoadImbalance();
		}
		double last_major_percent = ((double) last_major_reorg_seq
				/ (double) iterations) * 100.0;
		double last_minor_percent = ((double) last_minor_reorg_seq
				/ (double) iterations) * 100.0;

		double memtable_sizes = 0;
		double subrange_sizes = 0;

		for (SubRange sr : subranges) {
			for (TreeSet<InternalKey> memtable : sr.memtables) {
				for (InternalKey ik : memtable) {
					memtable_sizes += ik.byteSize();
				}
			}
			subrange_sizes += sr.byteSize();
		}

		System.out.println(String.format(
				"reorg,%d,%d,%d,%d,%d,%d,%.2f,%.2f,%d,%s,%f,%f",
				numberOfPerformedMajor,
				numberOfTotalMajor - numberOfPerformedMajor,
				numberOfPerformedMinor, numberOfPerformedMinorSampling,
				numberOfTotalMinor - numberOfPerformedMinor,
				numberOfPerformedMinorDup, last_major_percent,
				last_minor_percent, (int) num_iteration_grow, diff.toString(),
				memtable_sizes / 1024, subrange_sizes / 1024));
	}

	public static final String ZIPFIAN_FILE = "/tmp/zipfian";

	private static LoadImbalanceMetric ZipfianLoadImbalance() {
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
						Range r = new Range();
						subrange.tinyRanges.add(r);
						if (lower == i) {
							// include this key.
							r.lower = lower;
							r.upper = i + 1;
							r.currentShare = (share + percent) / sum;
							lower = i + 1;
							idealRanges.add(subrange);
							share = 0;
							continue;
						} else {
							r.lower = lower;
							r.upper = i;
							r.currentShare = share / sum;
							lower = i;
							idealRanges.add(subrange);
							share = 0;
						}
					}
					share += percent;
				}
			}

			StringBuilder builder = new StringBuilder();
			builder.append("ideal:");
			for (int i = 0; i < idealRanges.size(); i++) {
				SubRange sr = idealRanges.get(i);
				builder.append(String.format("[%d:%d):%.2f%%", sr.lower(),
						sr.upper(), sr.currentShare() * 100.0));
				builder.append(",");
			}
			builder = builder.deleteCharAt(builder.length() - 1);
			System.out.println(builder.toString());

			List<Double> loads = Lists.newArrayList();
			for (int i = 0; i < idealRanges.size(); i++) {
				double actualShare = 0;
				SubRange as = idealRanges.get(i);
				for (int k = as.lower(); k < as.upper(); k++) {
					actualShare += ref[k];
				}
				actualShare /= sum;
				loads.add(actualShare);
			}
			LoadImbalanceMetric metric = LoadImbalanceMetric.compute(loads);
			System.out.println("ideal-load," + metric.toString());

			Map<Integer, Integer> indexes = Maps.newHashMap();
			for (int i = 0; i < subranges.size(); i++) {
				SubRange as = subranges.get(i);
				if (as.duplicatedRange()) {
					for (int k = as.lower(); k < as.upper(); k++) {
						indexes.compute(k, (kk, v) -> {
							if (v == null) {
								return 1;
							}
							return v + 1;
						});
					}
				}
			}

			loads.clear();
			for (int i = 0; i < subranges.size(); i++) {
				double actualShare = 0;
				SubRange as = subranges.get(i);
				for (int k = as.lower(); k < as.upper(); k++) {
					Integer count = indexes.get(k);
					if (count != null) {
						actualShare += ref[k] / count.intValue();
					} else {
						actualShare += ref[k];
					}
				}
				actualShare /= sum;
				loads.add(actualShare);
			}
			br.close();
			return LoadImbalanceMetric.compute(loads);
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return null;
	}

	private static LoadImbalanceMetric UniformLoadImbalance() {
		List<Double> loads = Lists.newArrayList();
		for (int i = 0; i < subranges.size(); i++) {
			SubRange actual = subranges.get(i);
			double actualShare = (double) (actual.upper() - actual.lower())
					/ (double) numberOfKeys;
			loads.add(actualShare);
		}
		return LoadImbalanceMetric.compute(loads);
	}

	public static long iterations = 100000000;
	public static String dist = "uniform";

	public static List<SubRange> history = Lists.newArrayList();

	public static void run() {
		dist = "uniform";
		numberOfMemTables = 256;
		numberOfKeys = 10000000;
		numberOfSubRanges = 128;
		r = new Random(0);
		enableMinor = true;
		fairShare = 1.0 / numberOfSubRanges;
		samplingRatio = 1;
		numberOfTinyRangesPerSubRange = 10;
		iterations = 10000000;

		NumberGenerator gen = null;
		if ("zipfian".equals(dist)) {
			gen = new ZipfianGenerator(numberOfKeys, r);
		} else {
			gen = new UniformGenerator(numberOfKeys, r);
		}
		runSim(gen);
	}

	public static void main(String[] args) {
//		run();
		numberOfTinyRangesPerSubRange = 10;
		dist = args[0];
		numberOfMemTables = Integer.parseInt(args[1]);
		numberOfKeys = Integer.parseInt(args[2]);
		numberOfSubRanges = Integer.parseInt(args[3]);
		iterations = Long.parseLong(args[4]);
		samplingRatio = Integer.parseInt(args[5]);
		enableMinor = Boolean.parseBoolean(args[6]);
		long seed = Long.parseLong(args[7]);
		r = new Random(seed);
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
		int interval = 100000;
		for (long seq = 0; seq < iterations; seq++) {
			int key = gen.nextValue().intValue();
			InternalKey ik = new InternalKey();
			ik.key = key;
			ik.sequenceNumber = seq;
			addKey(ik);
			if (subranges.size() == numberOfSubRanges && seq % interval == 0) {
				double sum = 0;
				List<Double> loads = Lists.newArrayList();
				for (int i = 0; i < numberOfSubRanges; i++) {
					double diff = subranges.get(i).cumulativeNumberOfInserts
							- history.get(i).cumulativeNumberOfInserts;
					sum += diff;
					loads.add(diff);
					history.get(i).cumulativeNumberOfInserts = subranges
							.get(i).cumulativeNumberOfInserts;

				}
				for (int i = 0; i < numberOfSubRanges; i++) {
					loads.set(i, loads.get(i) / sum);
				}

				if (sum > 10000) {
					LoadImbalanceMetric m = LoadImbalanceMetric.compute(loads);
					System.out.println(
							String.format("seq,%d,%s", seq, m.toString()));
				}
			}
		}
		printFinalRanges();
		printStats();
	}

}
