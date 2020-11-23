package algorithms;

import java.util.Collections;
import java.util.Comparator;
import java.util.Random;

import common.EuclideanDistance;
import common.IDistanceFunction;
import common.Instance;
import common.Table;
import common.collectors.ICollector;

/*
 * Single thread implementation of SuperEGO algorithm.
 * 
 * Original Source Code:
 * https://www.ics.uci.edu/~dvk/code/SuperEGO.html
 * 
 * Publication:
 * Kalashnikov, D. V. (2013).
 * Super-EGO: Fast Multi-dimensional Similarity Join.
 * VLDB Journal, 22(4), 561–585. https://doi.org/10.1007/s00778-012-0305-7 
 */
public class SuperEGO implements IRangeJoin {

	private double[][] range1;
	private double[][] range2;
	private double[][] range3;
	
	private Table table1;
	private Table table2;
	
	private final Random rand;
	private final IDistanceFunction dist;
	private ICollector collector;	
	
	// number of dimensions
	private int numDim;

	// threshold value to execute join operation between two subsets
	private int t;
	
	// when true, perform dimension reordering for performance gains
	private final boolean reorderDim;
	
	// when true, consider single instance as self-similar pair 
	private final boolean allowSelfSimilar;
	
	// when true, means that the algorithm must consider a self-similarity join
	private boolean isSelfJoin;
	
	public SuperEGO(int joinThreshold,
					boolean reorderDim,
					boolean allowSelfSimilar,
					ICollector collector) {
		this.reorderDim = reorderDim;
		this.allowSelfSimilar = allowSelfSimilar;
		this.collector = collector;
		this.t = joinThreshold;
		this.isSelfJoin = false;
		this.rand = new Random();
		this.dist = new EuclideanDistance();
	}
	
	@Override
	public void range(Table a, double r) {
		isSelfJoin = true;
		range(a, a, r);
		isSelfJoin = false;
	}
	
	@Override
	public void range(Table a, Table b, double r) {
		assert a.getNumCols() == b.getNumCols();
		assert r > 0;
		
		reset();
		
		numDim = a.getNumCols();
		
		// ranges for SimpleJoin
		range1 = new double[numDim + 1][2];
		range2 = new double[numDim + 1][2];
		range3 = new double[numDim + 1][2];		

		if (isSelfJoin) {
			table2 = table1 = a.clone();
		} else {
			table1 = a.clone();
			table2 = b.clone();
		}
		
		// reorder dimension
		if (reorderDim)
			doDimensionReorder(r);
		
		// ego-sort
		EGOSort(table1, r);
		if (!isSelfJoin)
			EGOSort(table2, r);

		// ego-join
		int startDim = 0;
		int frA = 0;
		int toA = table1.getNumRows() - 1;
		int frB = 0;
		int toB = table2.getNumRows() - 1;
		
		EGOJoin(frA, toA, frB, toB, startDim, r);	
	}
	
	private void doDimensionReorder(double eps) {
		assert table1.getNumCols() >= 2;
		assert table2.getNumCols() >= 2;
		assert eps > 0;
		
		int sampleSize = (int) Math.ceil(table1.getNumRows() * 0.2);
		int numBuckets = (int) Math.ceil(1/eps) + 1;
		
		double[][] hA = new double[numBuckets][numDim];
		double[][] hB = new double[numBuckets][numDim];
		double[] avgDist = new double[numDim];
		
		for (int i = 0; i < sampleSize; i++) {
			// get random point from A
			int maxRowsA = (int) table1.getNumRows() - 1;
			int rndRowAId = nextInt(0, maxRowsA);
			double[] rndA = table1.getInstance(rndRowAId).getValues();
			
			// get random point from B
			int maxRowsB = (int) table2.getNumRows() - 1;
			int rndRowBId = nextInt(0, maxRowsB);
			double[] rndB = table2.getInstance(rndRowBId).getValues();
			
			// update stats
			for (int j = 0; j < numDim; j++) {
				avgDist[j] += Math.abs(rndA[j] - rndB[j]);
				int bucket_a = (int) (rndA[j]/eps);
				int bucket_b = (int) (rndB[j]/eps);
				hA[bucket_a][j] += 1;
				hB[bucket_b][j] += 1;
			}
		}
		
		// compute fail factor f in each dimension
		double[] f = new double[numDim];
		double[] g = new double[numDim];
		
		for (int i = 0; i < numDim; i++) {
			f[i] = 0;
			for (int j = 0; j < numBuckets; j++) {
				if (j == 0) {
					f[i] += hA[j][i] * (hB[j][i] + hB[j+1][i]); // assumes at least two columns (obviously)
					continue;
				} else if (j == (numBuckets - 1)) {
					f[i] += hA[j][i] * (hB[j][i] + hB[j-1][i]);
					continue;
				}
				
				f[i] += hA[j][i] * (hB[j][i] + hB[j-1][i] +hB[j+1][i]);
			}
		}
		
		// normalizing f
		for (int i = 0; i < numDim; i++) {
			avgDist[i] = avgDist[i] / sampleSize;
			f[i] = f[i]/ (sampleSize * sampleSize);
			g[i] = -avgDist[i];
		}
		
		int[] map = new int[numDim];
		
		// constracting map of remapping (inefficient)
		for (int i = 0; i < numDim; i++) {
			double min = Double.MAX_VALUE;
			int min_idx = -1;
			for (int j = 0; j < numDim; j++) {
				if (g[j] < min) {
					min = g[j];
					min_idx = j;
				}
			}
			map[i] = min_idx;
			g[min_idx] = Double.MAX_VALUE;
		}
		
		
		// reorder dimension in A (inefficient)
		for (int i = 0; i < table1.getNumRows(); i++) {
			double[] x = new double[numDim];
			for (int j = 0; j < numDim; j++)
				x[j] = table1.getInstance(i).getAt(j);
			for (int j = 0; j < numDim; j++)
				table1.setAt(i, j, x[map[j]]);
		}
		
		// reorder dimension in A (inefficient)
		if (!isSelfJoin) {
			for (int i = 0; i < table2.getNumRows(); i++) {
				double[] x = new double[numDim];
				for (int j = 0; j < numDim; j++)
					x[j] = table2.getInstance(i).getAt(j);
				for (int j = 0; j < numDim; j++)
					table2.setAt(i, j, x[map[j]]);
			}
		}
		
		// reorder stats accordingly
		double[] rs = new double[numDim];
		double[] rd = new double[numDim];
		for(int j = 0; j < numDim; j++) {
			rs[j] = 1 - f[map[j]];
			rd[j] = avgDist[map[j]];
		}
		
		
		//-- Case 1: zero inactive dimensions --
		int smallSeqSize = Math.min(1, numDim);
		for (int i = 0; i< smallSeqSize; i++) {
			range1[i][0] = 0;
			range1[i][1] = numDim - 1;
			range2[i][0] = 0;
			range2[i][1] = -1;
			range3[i][0] = 0;
			range3[i][1] = -1;
		}

		//-- Case 2: all dims are inactive --
		range1[numDim][0] = 0;
		range1[numDim][1] = numDim - 1;
		range2[numDim][0] = 0;
		range2[numDim][1] = -1;
		range3[numDim][0] = 0;
		range3[numDim][1] = -1;
	    
		//-- Case 3: remaining cases --
	    //-- find first k s.t. rd[k] < eps/2 --
	    int k = numDim;
	    
	    for (int i = 0; i < numDim; i++) {
	        if (rd[i] < eps/2) {
	            k = i;
	            break;
	        }
	    }
	    
	    for (int i = smallSeqSize; i < numDim; i++) {
	    	//-- Case I: 1-interval --
	    	if (rd[i] < eps/2) {
				range1[i][0] = 0;
				range1[i][1] = numDim - 1;
				range2[i][0] = 0;
				range2[i][1] = -1;
				range3[i][0] = 0;
				range3[i][1] = -1;
				continue;
	    	}
	    	//-- Case II: 3-intervals --
	    	if (k < numDim) {
				range1[i][0] = i;
				range1[i][1] = k - 1;
				range2[i][0] = 0;
				range2[i][1] = i - 1;
				range3[i][0] = k;
				range3[i][1] = numDim - 1;
				continue;
	    	}
	    	//-- Case III: 2-interval --
			range1[i][0] = i;
			range1[i][1] = numDim - 1;
			range2[i][0] = 0;
			range2[i][1] = i - 1;
			range3[i][0] = 0;
			range3[i][1] = -1;
	    }
	}
	
	private int nextInt(int min, int max) {
		return rand.nextInt(max - min + 1) + min;
	}
	
	private void EGOSort(Table t, double eps) {
		Collections.sort(t.asList(), new Comparator<Instance>() {
			@Override
			public int compare(Instance r1, Instance r2) {
				for (int i = 0; i < r1.getValues().length; i++) {
					int d = ((int) (r1.getAt(i)/eps)) - ((int) (r2.getAt(i)/eps));					
					if (d != 0)
						return d;
				}
				return 0;
			}
		});
	}
	
	private void EGOJoin(int frA, int toA, int frB, int toB, int startDim, double eps) {
		int szA = toA - frA + 1; 
		int szB = toB - frB + 1;

		double[] fstA = table1.getInstance(frA).getValues();
		double[] lstA = table1.getInstance(toA).getValues();
		double[] fstB = table2.getInstance(frB).getValues();
		double[] lstB = table2.getInstance(toB).getValues();
		
		// Ego-Strategy
		double loA, hiA, loB, hiB;
		for (int i = startDim; i < numDim; i++) {
			loA = (int) (fstA[i] / eps);
			hiB = (int) (lstB[i] / eps);
			if (loA > hiB + 1) return;
			loB = (int) (fstB[i] / eps);
			hiA = (int) (lstA[i] / eps);
			if (loB > hiA + 1) return;
			if ((loA < hiA) || (loB < hiB)) {
				startDim = i;
				break;
			}
		}
		
		// Ego-Join
		int midA = (int) (frA + (szA/2.0));
		int midB = (int) (frB + (szB/2.0));
		
		if ((szA < t) && (szB < t)) {
//			NaiveJoin(frA, toA, frB, toB, eps);
			SimpleJoin(frA, toA, frB, toB, eps);
			return;
		}

		if ((szA < t) && (szB >= t)) {
			EGOJoin(frA, toA, frB     , midB, startDim, eps);
			EGOJoin(frA, toA, midB + 1,  toB, startDim, eps);
			return;
		}
		
		if ((szA >= t) && (szB < t)) {
			EGOJoin(frA     , midA, frB, toB, startDim, eps);
			EGOJoin(midA + 1, toA , frB, toB, startDim, eps);
			return;
		}
		
		if ((szA >= t) && (szB >= t)) {
			EGOJoin(frA     , midA, frB     , midB, startDim, eps);
			EGOJoin(frA     , midA, midB + 1, toB , startDim, eps);
			EGOJoin(midA + 1, toA , midB + 1, toB , startDim, eps);
//			if (!isSelfJoin)
				EGOJoin(midA + 1, toA , frB     , midB, startDim, eps);
			return;
		}
	}
	
	private void NaiveJoin(int frA, int toA, int frB, int toB, double eps) {
		// TODO: This implementation must be reviewed 
		// and possibly extracted to a interface containing 
		// a join method so other versions of joins 
		// could be easily implemented
		Instance p, q;		

		for (int i = frA; i <= toA; i++) {
			p = table1.getInstance(i);
			for(int j = frB; j <= toB; j++) {
				q = table2.getInstance(j);
				
				if (!allowSelfSimilar && (p.getId() == q.getId()))
					continue;
				
				if (dist.compute(p.getValues(), q.getValues()) <= eps) {
					if (isSelfJoin) {
						collector.addPair(p.getId(), q.getId());
						collector.addPair(q.getId(), p.getId());
					} else {
						collector.addPair(p.getId(), q.getId());						
					}
				}
			}
		}
	}	
	
	private void SimpleJoin(int frA, int toA, int frB, int toB, double eps) {
		Instance p, q;
		double eps2 = eps * eps;
		double s, dx;
		boolean skip;
		
		for (int i = frA; i <= toA; i++) {
			p = table1.getInstance(i);

			for (int j = frB; j <= toB; j++) {
				q = table2.getInstance(j);
				
				if (!allowSelfSimilar && (p.getId() == q.getId()))
					continue;
				
				s = 0;
				dx = 0;
				skip = false;

				for (int k = 0; k < numDim; k++) {
					dx = p.getAt(k) - q.getAt(k);
					s += dx * dx;
					if (s >= eps2) {
						skip = true;
						break;
					}
				}
				
				if (!skip) {
					if (isSelfJoin) {
						collector.addPair(p.getId(), q.getId());
						collector.addPair(q.getId(), p.getId());
					} else {
						collector.addPair(p.getId(), q.getId());						
					}					
				}
				
				skip = false;
			}
		}
	}
	
	private void  SimpleJoinAlternative(int frA, int toA, int frB, int toB, double eps) {
		// TODO: This implementation must be reviewed 
		// and possibly extracted to a interface containing 
		// a join method so other versions of joins 
		// could be easily implemented
		Instance p, q;
		double s;
		double dx;
		boolean skip;
		final double eps2 = eps * eps;

		for (int i = frA; i <= toA; i++) {
			p = table1.getInstance(i);
			for (int j = frB; j <= toB; j++) {
				q = table2.getInstance(j);
				s = 0;
				skip = false;
				int mid = (int) Math.floor(numDim/2);
				
				for (int k = mid; !skip && k < numDim; k++) {
					dx = p.getAt(k) - q.getAt(k);
					s += dx * dx;
					if (s >= eps2) {
						skip = true;
						break;
					}
				}
				
				for (int k = 0; !skip && k < mid; k++) {
					dx = p.getAt(k) - q.getAt(k);
					s += dx * dx;
					if (s >= eps2) {						
						skip = true;
						break;
					}
				}
				
				if (!skip) {
					if (isSelfJoin) {
						collector.addPair(p.getId(), q.getId());
						collector.addPair(q.getId(), p.getId());
					} else {
						collector.addPair(p.getId(), q.getId());						
					}
				}
				
				skip = false;
			}
		}
	}

	@Override
	public void reset() {
		collector.clear();
	}

	@Override
	public ICollector getCollector() {
		return collector;
	}
		
}
