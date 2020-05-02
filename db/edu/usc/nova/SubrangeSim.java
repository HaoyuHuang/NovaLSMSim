package edu.usc.nova;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.util.Comparator;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
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
	public static int numberOfPerformedMinorDup = 0;

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
		boolean duplicatedRange = false;

		public SubRange() {
			for (int i = 0; i < numberOfMemTables / numberOfSubRanges; i++) {
				memtables.add(new TreeSet<>(comp));
			}
		}

		@Override
		public String toString() {
			return "SubRange [lower=" + lower + ", upper=" + upper
					+ ", numberOfInserts=" + numberOfInserts + ", currentShare="
					+ currentShare + ", keys=" + (upper - lower) + "]";
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

	static boolean subrangeContains(int key, SubRange sr) {
		if (key < sr.lower) {
			return false;
		}
		if (key >= sr.upper) {
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
			if (low >= 0 && subrangeContains(key, subranges.get(low))) {
				possibleRanges.add(low);
				contains = true;
			}
			if (high < subranges.size()
					&& subrangeContains(key, subranges.get(high))) {
				possibleRanges.add(high);
				contains = true;
			}
			if (!contains) {
				break;
			}
		}

		int i = r.nextInt(possibleRanges.size());
		index = possibleRanges.get(i);
		assert subrangeContains(key, subranges.get(index));
		return index;
	}

	static int binarySearch(int key) {
		int l = 0, r = subranges.size() - 1;
		while (l <= r) {
			int m = l + (r - l) / 2;

			// Check if x is present at mid
			if (key < subranges.get(m).lower) {
				r = m - 1;
			} else if (key >= subranges.get(m).upper) {
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
		if (sortedMap.size() <= numberOfSubRanges * 2) {
			return 0;
		}

		assert sortedMap.size() > numberOfSubRanges * 2;
		numberOfPerformedMajor++;
		last_major_reorg_seq = totalNumberOfInserts;

		for (int i = 0; i < numberOfSubRanges; i++) {
			fixedMinor.put(i, totalNumberOfInserts);
		}

		subranges.clear();
		for (int i = 0; i < numberOfSubRanges; i++) {
			subranges.add(new SubRange());
		}

		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			sr.duplicatedRange = false;
		}

		double sharePerSubRange = totalShares / subranges.size();
		double fairShare = totalShares / subranges.size();
		double total = totalShares;
		int index = 0;
		double sum = 0;
		subranges.get(0).lower = 0;
		subranges.get(subranges.size() - 1).upper = numberOfKeys;
		for (Entry<Integer, Double> entry : sortedMap.entrySet()) {
			double rate = entry.getValue();
			if (rate >= fairShare) {
				System.out.println(String.format("hot key %d:%.2f:%.2f",
						entry.getKey(), rate / total, fairShare / total));
				// close the current subrange.
				if (subranges.get(index).lower < entry.getKey()) {
					subranges.get(index).upper = entry.getKey();
					index++;
				}

				int nDuplicates = (int) Math.ceil(rate / fairShare);
				for (int i = 0; i < nDuplicates; i++) {
					subranges.get(index).duplicatedRange = true;
					subranges.get(index).lower = entry.getKey();
					subranges.get(index).upper = entry.getKey() + 1;
					index++;
				}
				subranges.get(index).lower = entry.getKey() + 1;
				totalShares -= entry.getValue();
				sum = 0;
				sharePerSubRange = totalShares / (subranges.size() - index);
				continue;
			}

			if (sum + rate > sharePerSubRange) {
				if (subranges.get(index).lower == entry.getKey()) {
					sum += rate;
					subranges.get(index).currentShare = sum / total;
					subranges.get(index).upper = entry.getKey() + 1;
					subranges.get(index + 1).lower = entry.getKey() + 1;
					index++;
					if (index == subranges.size() - 1) {
						break;
					}
					sum = 0;
					totalShares -= rate;
					sharePerSubRange = totalShares / (subranges.size() - index);
					continue;
				} else {
					subranges.get(index).currentShare = sum / total;
					subranges.get(index).upper = entry.getKey();
					subranges.get(index + 1).lower = entry.getKey();
					index++;
					if (index == subranges.size() - 1) {
						break;
					}
					sum = 0;
					sharePerSubRange = totalShares / (subranges.size() - index);
				}
			}
			sum += rate;
			totalShares -= rate;
		}

		int priorUpper = -1;
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			if (priorUpper != -1 && priorUpper - sr.lower > 1) {
				assert sr.lower > priorUpper;
				if (sr.lower < priorUpper) {
					sr.lower = priorUpper;
				}
			}
			if (sr.upper <= sr.lower) {
				sr.upper = sr.lower + 1;
				priorUpper = sr.upper;
			}
			sr.numberOfInserts = 0;
			sr.currentShare = 0;
			sr.ninternalKeys = 0;
		}
		totalNumberOfInsertsSinceLastMajor = 0;
		subranges.get(0).lower = 0;
		subranges.get(subranges.size() - 1).upper = numberOfKeys;
		return sortedMap.size();
	}

	public static enum UpdateBoundary {
		LOWER, UPPER, NONE
	}

	public static Map<Integer, Double> fixedMinor = Maps.newHashMap();

	public static void moveShare(int index) {
		SubRange sr = subranges.get(index);
		int nDuplicates = 0;
		int start = -1;
		int end = -1;
		for (int i = 0; i < subranges.size(); i++) {
			if (subranges.get(i).lower != sr.lower) {
				continue;
			}
			end = i;
			if (start == -1) {
				start = i;
			}
		}
		nDuplicates = end - start;
		double share = sr.numberOfInserts / nDuplicates;
		for (int i = start; i <= end; i++) {
			if (nDuplicates == 1) {
				subranges.get(i).duplicatedRange = false;
			}
			subranges.get(i).numberOfInserts += share;
			subranges.get(i).currentShare = subranges.get(i).numberOfInserts
					/ totalNumberOfInsertsSinceLastMajor;
		}
	}

	public static boolean destroyDuplicates(int index) {
		SubRange sr = subranges.get(index);
		if (!sr.duplicatedRange) {
			return false;
		}

		double percent = sr.currentShare / fairShare;
		if (percent >= 0.5) {
			return false;
		}

//		numberOfTotalMinor++;
		last_minor_reorg_seq = totalNumberOfInserts;
//		numberOfPerformedMinor++;
		numberOfPerformedMinorDup++;

		for (int i = 0; i < numberOfSubRanges; i++) {
			fixedMinor.put(i, totalNumberOfInserts);
		}

		System.out.println("Destroy subrange " + index);
		// destroy this duplicate.
		moveShare(index);
		subranges.remove(index);
		return true;
	}

	public static boolean minorRebalanceDuplicate(int index) {
		SubRange sr = subranges.get(index);
		assert sr.currentShare > fairShare;
		sr.ninternalKeys = 0;
		if (sr.upper - sr.lower != 1) {
			return false;
		}

		int nDuplicates = (int) Math.floor(sr.currentShare / fairShare);
		if (nDuplicates == 0) {
			return false;
		}

		last_minor_reorg_seq = totalNumberOfInserts;
		numberOfPerformedMinorDup++;

		for (int i = 0; i < numberOfSubRanges; i++) {
			fixedMinor.put(i, totalNumberOfInserts);
		}

		sr.duplicatedRange = true;
		double remainingSum = sr.numberOfInserts;
		for (int i = 0; i < nDuplicates; i++) {
			SubRange newSR = new SubRange();
			newSR.lower = sr.lower;
			newSR.upper = sr.upper;
			newSR.numberOfInserts = sr.numberOfInserts / (nDuplicates + 1);
			remainingSum -= newSR.numberOfInserts;
			newSR.currentShare = newSR.numberOfInserts
					/ totalNumberOfInsertsSinceLastMajor;
			newSR.duplicatedRange = true;
			subranges.add(index, newSR);
		}
		sr.numberOfInserts = remainingSum;
		sr.currentShare = sr.numberOfInserts
				/ totalNumberOfInsertsSinceLastMajor;

		// if (index == 18) {
//			printRanges("debug");
//			System.out.println();
//		}
//		
		while (subranges.size() > numberOfSubRanges) {
			// remove min.
			double minShare = Double.MAX_VALUE;
			int minRangeId = 0;
			for (int i = 0; i < subranges.size(); i++) {
				SubRange minSR = subranges.get(i);
				// Skip the new subranges.
				if (minSR.lower == sr.lower) {
					continue;
				}
//				if (minSR.duplicatedRange) {
//					continue;
//				}
				if (minSR.currentShare < minShare) {
					minShare = minSR.currentShare;
					minRangeId = i;
				}
			}

			// merge min with neighboring subrange.
			int left = minRangeId - 1;
			int right = minRangeId + 1;
			boolean mergeLeft = true;
			if (left >= 0 && right < subranges.size()) {
				if (subranges.get(left).currentShare < subranges
						.get(right).currentShare) {
					mergeLeft = true;
				} else {
					mergeLeft = false;
				}
			} else if (left >= 0) {
				mergeLeft = true;
			} else {
				mergeLeft = false;
			}
			SubRange minSR = subranges.get(minRangeId);
			if (mergeLeft) {
				SubRange l = subranges.get(left);
				if (l.duplicatedRange) {
					// move its shares to other duplicates.
					moveShare(left);
					// destroy this subrange.
					subranges.remove(left);
				} else {
					l.currentShare += minSR.currentShare;
					l.numberOfInserts += minSR.numberOfInserts;
					l.upper = minSR.upper;
					subranges.remove(minRangeId);
				}
			} else {
				SubRange r = subranges.get(right);
				if (r.duplicatedRange) {
					// move its shares to other duplicates.
					moveShare(right);
					// destroy this subrange.
					subranges.remove(right);
				} else {
					r.currentShare += minSR.currentShare;
					r.numberOfInserts += minSR.numberOfInserts;
					r.lower = minSR.lower;
					subranges.remove(minRangeId);
				}
			}
		}
		printRanges("minor duplicate subrange-" + index);
		return true;
	}

	public static boolean minorRebalanceHigherShareDistribute(int index) {
		SubRange sr = subranges.get(index);
		double totalRemoveInserts = (sr.currentShare - fairShare)
				* totalNumberOfInsertsSinceLastMajor;

		if (totalRemoveInserts >= sr.numberOfInserts) {
			return false;
		}

		if (minorRebalanceDuplicate(index)) {
			return true;
		}

		assert sr.currentShare > fairShare;
		sr.ninternalKeys = 0;
		for (int i = 0; i < numberOfSubRanges; i++) {
			fixedMinor.put(i, totalNumberOfInserts);
		}
		// Distribute the load across adjacent subranges.
		// higher share.
		numberOfTotalMinor++;
		TreeMap<Integer, Integer> samples = sr.sampleKeys();
		if (samples.size() <= 1 || sr.ninternalKeys < 1000
				|| sr.ninternalKeys > sr.numberOfInserts) {
			sr.ninternalKeys = 0;
			return false;
		}

		last_minor_reorg_seq = totalNumberOfInserts;
		System.out.println(
				String.format("!!!!!!minor-%d-%d-%.2f", index, samples.size(),
						(sr.currentShare - fairShare) * 100.0 / fairShare));
		System.out.println(sr.toString());

		double totalRemovePuts = (double) (((sr.currentShare - fairShare)
				/ sr.currentShare) * sr.ninternalKeys);

		double leftPuts = totalRemovePuts;
		double rightPuts = totalRemovePuts;
		double leftInserts = totalRemoveInserts;
		double rightInsert = totalRemoveInserts;

		double totalRateToDistribute = sr.currentShare - fairShare;
		double rightRate = 0.0;
		double leftRate = 0.0;
		double remainingRate = totalRateToDistribute;

		double newTotalRightRate = 0;
		double newTotalLeftRate = 0;

		if (index != 0 && index != subranges.size() - 1
				&& !subranges.get(index - 1).duplicatedRange
				&& !subranges.get(index + 1).duplicatedRange) {
			SubRange left = subranges.get(index - 1);
			SubRange right = subranges.get(index + 1);

			if (left.currentShare > fairShare
					&& right.currentShare < fairShare) {
				double needRate = fairShare - right.currentShare;
				if (totalRateToDistribute <= needRate) {
					// give all to right.
					rightRate = totalRateToDistribute;
					remainingRate = 0;
				} else {
					rightRate = totalRateToDistribute - needRate;
					remainingRate -= rightRate;
				}
			} else if (left.currentShare < fairShare
					&& right.currentShare > fairShare) {
				double needRate = fairShare - left.currentShare;
				if (totalRateToDistribute <= needRate) {
					// give all to left.
					leftRate = totalRateToDistribute;
					remainingRate = 0;
				} else {
					leftRate = totalRateToDistribute - needRate;
					remainingRate -= leftRate;
				}
			}

			if (remainingRate > 0) {
				newTotalRightRate = right.currentShare + rightRate;
				newTotalLeftRate = left.currentShare + leftRate;

				double total = newTotalRightRate + newTotalLeftRate;
				double leftp = newTotalRightRate / total;
				double rightp = newTotalLeftRate / total;

				rightRate += remainingRate * rightp;
				leftRate += remainingRate * leftp;
			}

			assert Math
					.abs(totalRateToDistribute - rightRate - leftRate) < 0.001;

			double leftRatio = leftRate / totalRateToDistribute;
			double rightRatio = rightRate / totalRateToDistribute;

			leftPuts = totalRemovePuts * leftRatio;
			leftInserts = totalRemoveInserts * leftRatio;

			rightPuts = totalRemovePuts * rightRatio;
			rightInsert = totalRemoveInserts * rightRatio;

			if (leftPuts < samples.firstEntry().getValue()
					&& rightPuts < samples.lastEntry().getValue()) {
				// cannot move either to left or right.
				// move to right.
				rightPuts = samples.lastEntry().getValue();

				rightInsert += leftInserts;
				leftPuts = 0;
				leftInserts = 0;
			}
		}

		boolean success = false;
		// update lower
		if (index > 0 && leftPuts > 0
				&& !subranges.get(index - 1).duplicatedRange) {
			int newLower = sr.lower;
			double removed_share = leftPuts;
			Set<Integer> removed = Sets.newHashSet();
			for (Entry<Integer, Integer> entry : samples.entrySet()) {
				if (removed_share - entry.getValue() < 0) {
					break;
				}
				newLower = entry.getKey();
				removed_share -= entry.getValue();
				removed.add(entry.getKey());
			}

			removed.forEach(k -> {
				samples.remove(k);
			});

			if (newLower > sr.lower) {
				sr.lower = newLower;
				sr.numberOfInserts -= leftInserts;
				sr.currentShare = sr.numberOfInserts
						/ totalNumberOfInsertsSinceLastMajor;
				SubRange other = subranges.get(index - 1);
				other.numberOfInserts += leftInserts;
				other.upper = sr.lower;
				other.currentShare = other.numberOfInserts
						/ totalNumberOfInsertsSinceLastMajor;
				success = true;
			}
		}
		// update upper.
		if (index < subranges.size() - 1 && rightPuts > 0
				&& !subranges.get(index + 1).duplicatedRange) {
			double removed_share = rightPuts;
			int newUpper = sr.upper;
			for (Entry<Integer, Integer> entry : samples.descendingMap()
					.entrySet()) {
				if (removed_share - entry.getValue() < 0) {
					break;
				}
				newUpper = entry.getKey();
				removed_share -= entry.getValue();
			}
			if (newUpper < sr.upper && newUpper - sr.lower >= 1) {
				sr.upper = newUpper;
				sr.numberOfInserts -= rightInsert;
				sr.currentShare = sr.numberOfInserts
						/ totalNumberOfInsertsSinceLastMajor;
				SubRange other = subranges.get(index + 1);
				other.numberOfInserts += rightInsert;
				other.lower = sr.upper;
				other.currentShare = other.numberOfInserts
						/ totalNumberOfInsertsSinceLastMajor;
				success = true;
			}
		}
		if (success) {
			numberOfPerformedMinor++;
			System.out.println("after-" + sr.toString());
		}
		return success;
	}

	public static double num_iteration_grow = 0;
	public static double minor_reorg_interval = 100000;
	public static double major_reorg_interval = 1000000;

	public static void addKey(InternalKey ik) {
		int index = binarySearchWithDuplicate(ik.key);
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
						/ totalNumberOfInsertsSinceLastMajor;
				subrange = newSubrange;
			}
		}
		assert subrange != null;
		subrange.addKey(ik);
		subrange.numberOfInserts += 1;
		subrange.cumulativeNumberOfInserts += 1;
		totalNumberOfInserts += 1;
		totalNumberOfInsertsSinceLastMajor += 1;
		subrange.currentShare = subrange.numberOfInserts
				/ totalNumberOfInsertsSinceLastMajor;

		double sum = 0;
		for (int i = 0; i < subranges.size(); i++) {
			sum += subranges.get(i).numberOfInserts;
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

		// check.
		double unfairRanges = 0;
		int mostUnfairRange = -1;
		double mostUnfair = 0;
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			double diff = (sr.currentShare - fairShare) * 100.0 / fairShare;
			if (Math.abs(diff) > 20 && sr.upper - sr.lower > 1) {
				unfairRanges += 1;
			}

			if (diff > 20
					&& totalNumberOfInserts
							- fixedMinor.get(i) > minor_reorg_interval
					&& sr.numberOfInserts > minor_reorg_interval
							/ subranges.size()) {
				if (diff > mostUnfair) {
					mostUnfairRange = i;
					mostUnfair = diff;
				}
			}

			if (totalNumberOfInserts - fixedMinor.get(i) > minor_reorg_interval
					&& sr.numberOfInserts > minor_reorg_interval
							/ subranges.size()) {
				if (destroyDuplicates(i)) {
					mostUnfairRange = -1;
					mostUnfair = -1;
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
				if (enableMinor && mostUnfairRange != -1) {
					prefix = "minor-" + mostUnfairRange + ":";
					minorRebalanceHigherShareDistribute(mostUnfairRange);
					return;
				}
				return;
			}
			printRanges(prefix);
		} else if (enableMinor && mostUnfairRange != -1) {
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
	}

	public static boolean enableMinor = true;

	private static void printRanges(String prefix) {
		StringBuilder builder = new StringBuilder();
		builder.append(prefix);
		builder.append("\n");
		for (int i = 0; i < subranges.size(); i++) {
			SubRange sr = subranges.get(i);
			builder.append(String.format("%d [%d:%d):%.2f,%d,%d", i, sr.lower,
					sr.upper, sr.currentShare, (int) sr.numberOfInserts,
					sr.upper - sr.lower));
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
			builder.append(String.format("%d [%d:%d):%.2f%%,%d", i, sr.lower,
					sr.upper,
					(sr.numberOfInserts / totalNumberOfInsertsSinceLastMajor)
							* 100.0,
					sr.upper - sr.lower));
			builder.append("\n");
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

		System.out
				.println(String.format("reorg,%d,%d,%d,%d,%d,%.2f,%.2f,%d,%.2f",
						numberOfPerformedMajor,
						numberOfTotalMajor - numberOfPerformedMajor,
						numberOfPerformedMinor,
						numberOfTotalMinor - numberOfPerformedMinor,
						numberOfPerformedMinorDup, last_major_percent,
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

			Map<Integer, Integer> indexes = Maps.newHashMap();
			for (int i = 0; i < subranges.size(); i++) {
				SubRange as = subranges.get(i);
				if (as.duplicatedRange) {
					for (int k = as.lower; k < as.upper; k++) {
						indexes.compute(k, (kk, v) -> {
							if (v == null) {
								return 1;
							}
							return v + 1;
						});
					}
				}
			}

			int maxRangeId = 0;
			for (int i = 0; i < subranges.size(); i++) {
				double actualShare = 0;
				SubRange as = subranges.get(i);
				for (int k = as.lower; k < as.upper; k++) {
					Integer count = indexes.get(k);
					if (count != null) {
						actualShare += ref[k] / count.intValue();
					} else {
						actualShare += ref[k];
					}
				}
				actualShare /= sum;

				if (actualShare > max) {
					maxRangeId = i;
				}
				max = Math.max(max, actualShare);
			}

			System.out.println("max range id " + maxRangeId);

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
		dist = "zipfian";
		numberOfMemTables = 64;
		numberOfKeys = 10000000;
		numberOfSubRanges = 32;
		r = new Random(0);
		enableMinor = true;
		fairShare = 1.0 / numberOfSubRanges;
		samplingRatio = 1;
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
			fixedMinor.put(i, 0.0);
		}
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
