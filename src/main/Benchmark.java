package main;

import java.util.ArrayList;

import algorithms.BruteForceJoin;
import algorithms.IRangeJoin;
import algorithms.QuickJoin;
import algorithms.SuperEGO;
import common.Table;
import common.collectors.HashCollector;


public class Benchmark {
	public static void main(String[] args) {
		
		ArrayList<Table> datasets = new ArrayList<>();
		datasets.add(Table.readCSV("datasets/iris.csv", ',', 4, true));
		
		ArrayList<IRangeJoin> algorithms = new ArrayList<>();
		algorithms.add(new BruteForceJoin(true, new HashCollector()));
		algorithms.add(new QuickJoin(10, true, new HashCollector()));
		algorithms.add(new SuperEGO(10, false, true, new HashCollector()));
		
		for (Table dataset : datasets) {
			for (IRangeJoin algorithm : algorithms) {
				// run algorithm
				System.out.println(algorithm.getClass().getName());
				algorithm.range(dataset, 25);
				
				// Print the number of neighbors for each instance 
				for (int i = 0; i < dataset.getNumRows(); i++) {
					int idx = dataset.getInstance(i).getId();
					System.out.println(idx + ": " + algorithm.getCollector().getPairsOf(idx).size());
				}
				
				// reset algorithm
				algorithm.reset();
			}
		}

	}
}
