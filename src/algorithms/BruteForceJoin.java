package algorithms;

import java.util.Collections;
import java.util.HashMap;
import java.util.Objects;
import java.util.PriorityQueue;

import common.EuclideanDistance;
import common.IDistanceFunction;
import common.Instance;
import common.Pair;
import common.Table;
import common.collectors.ICollector;

public class BruteForceJoin implements IRangeJoin {

	private boolean allowSelfSimilar;
	private IDistanceFunction dist;
	private ICollector collector;
	
	public BruteForceJoin(boolean allowSelfSimilar, ICollector collector) {
		this.allowSelfSimilar = allowSelfSimilar;
		this.dist = new EuclideanDistance();
		this.collector = collector;
	}
	
	public ICollector getCollector() {
		return collector;
	}

	public void range(Table a, double r) {
		range(a, a, r);
	}

	public void range(Table a, Table b, double r) {
		assert Objects.nonNull(a);
		assert Objects.nonNull(b);
		assert a.getNumRows() > 1;
		assert b.getNumRows() > 1;
		assert r > 0;
		
		reset();
		Instance p, q;

		for (int i = 0; i < a.getNumRows(); i++ ) {
			p =  a.getInstance(i);
			
			for (int j = 0; j < b.getNumRows(); j++) {
				if (!allowSelfSimilar && i == j) continue;				
				q =  b.getInstance(j);
				
				if (dist.compute(p.getValues(), q.getValues()) <= r) {
					collector.addPair(p.getId(), q.getId());
				}
			}
		}

	}

	public void kNN(Table t1, int k) {
		kNN(t1, t1, k);
	}

	public void kNN(Table t1, Table t2, int k) {
		assert Objects.nonNull(t1);
		assert Objects.nonNull(t2);
		assert t1.getNumRows() > 1;
		assert t2.getNumRows() > 1;
		assert k > 0;
		
		reset();
		Instance p, q;
		
		HashMap<Integer, PriorityQueue<Pair<Double, Integer>>> R = new HashMap<>();
		
		for (int i = 0; i < t1.getNumRows(); i++ ) {
			p = t1.getInstance(i);

			for (int j = 0; j < t2.getNumRows(); j++) {
				q = t2.getInstance(j);
			
				if (!allowSelfSimilar && p.getId() == q.getId())
					continue;
				
				if (!R.containsKey(p.getId())) {					
					PriorityQueue<Pair<Double, Integer>> pq;
					pq = new PriorityQueue<>(k + 1, Collections.reverseOrder()); 
					R.put(i, pq);
				}
				
				Double pqDist = dist.compute(p.getValues(), q.getValues());
				PriorityQueue<Pair<Double, Integer>> neighbors = R.get(p.getId());
				
				if (neighbors.size() < k) {
					neighbors.add(new Pair<>(pqDist, q.getId()));
					continue;
				}

				// else find max dist in p neighborhood
				Pair<Double, Integer> mostDistantIdx = neighbors.peek();
				if (pqDist < mostDistantIdx.getValue0()) {
					neighbors.poll();
					neighbors.add(new Pair<>(pqDist, q.getId()));
				}
			}
		}
		
		for (Integer idx1 : R.keySet()) {
			for (Pair<Double, Integer> distIdx : R.get(idx1)) {
				collector.addPair(idx1, distIdx.getValue1());
			}
		}
	}

	public void reset() {
		collector.clear();
	}

}
