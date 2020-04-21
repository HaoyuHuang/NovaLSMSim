package edu.usc.distributions;

import java.util.Random;

public class UniformGenerator extends NumberGenerator {

	private final int maxNumber;
	private final Random random;

	public UniformGenerator(int maxNumber, Random r) {
		super();
		this.maxNumber = maxNumber;
		this.random = r;
	}

	@Override
	public double mean() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public Integer nextValue() {
		return this.random.nextInt(maxNumber);
	}

}
