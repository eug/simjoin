package algorithms;

import java.util.ArrayList;
import java.util.Random;

import common.EuclideanDistance;
import common.IDistanceFunction;
import common.Instance;
import common.Pair;
import common.Table;
import common.collectors.ICollector;

/*
 * Single thread implementation of QuickJoin algorithm.
 * 
 * Publication:
 * Fredriksson, K., & Braithwaite, B. (2015).
 * Quicker range- and k-NN joins in metric spaces. 
 * Information Systems, 52, 189–204. https://doi.org/10.1016/j.is.2014.09.006
 */
public class QuickJoin implements IRangeJoin {
	
	private class Indices extends ArrayList<Integer> {
		private static final long serialVersionUID = 1L;
		public Table table;
		public Indices(Table t) {
			this(t, 10);
		}
		public Indices(Table t, int c){
			super(c);
			table = t;
		}
	}
	
	private class Partitions {
		public final Indices L;
		public final Indices G;
		public final Indices Lw;
		public final Indices Gw;
		
		public Partitions(Indices L, Indices G, Indices Lw, Indices Gw) {
			this.L = L;
			this.G = G;
			this.Lw = Lw;
			this.Gw = Gw;
		}
	}
	
	private Random rand = new Random();
	private double beta = 1;
	
	private int c;
	private IDistanceFunction dist;
	private ICollector collector;
	
	public QuickJoin(int joinThreshold,
					 boolean selfSimilar,
					 ICollector collector) {
		this.c = joinThreshold;
		this.dist = new EuclideanDistance();
		this.collector = collector;
	}
	
	@Override
	public void range(Table a, double r) {
		Indices s1 = getIndices(a);
		doQuickJoin(s1, r);
	}
	
	@Override
	public void range(Table a, Table b, double r) {
		Indices s1 = getIndices(a);
		Indices s2 = getIndices(b);	
		doQuickJoinWin(s1, s2, r);
	}

	private Indices getIndices(Table t) {
		Indices s = new Indices(t, t.getNumRows());
		for (int i = 0; i < t.getNumRows(); i++)
			s.add(t.getInstance(i).getId());
		return s;
	}
	
	private void doQuickJoin(Indices s, double r) {
		if (s.size() < c) {
			SelfJoinBF(s, r);
			return;
		}
		
		Pair<Instance, Instance> pivots = Pivots(s);
		Instance p1 = pivots.getValue0();
		Instance p2 = pivots.getValue1();
		
		double rho = beta * dist.compute(p1.getValues(), p2.getValues());
		
		Partitions part = Partition(s, p1, rho, r);
		
		doQuickJoinWin(part.Lw, part.Gw, r);
		doQuickJoin(part.L, r);
		doQuickJoin(part.G, r);
	}
	
	private void doQuickJoinWin(Indices s1, Indices s2, double r) {
		if (s1.size() + s2.size() <= c) {
			JoinBF(s1, s2, r);
			return;
		}
		
		Pair<Instance, Instance> pivots = Pivots(s1, s2);
		Instance p1 = pivots.getValue0();
		Instance p2 = pivots.getValue1();
		
		double rho = beta * dist.compute(p1.getValues(), p2.getValues());
		
		Partitions part1 = Partition(s1, p1, rho, r);
		Partitions part2 = Partition(s2, p1, rho, r); // in paper says p1
		
		doQuickJoinWin(part1.Lw, part2.Gw, r);
		doQuickJoinWin(part1.Gw, part2.Lw, r);
		doQuickJoinWin(part1.L, part2.L, r);
		doQuickJoinWin(part1.L, part2.G, r);
	}
	
	private Pair<Instance, Instance> Pivots(Indices s) {
		int pos1 = nextInt(0, s.size() - 1);
		int pos2 = nextInt(0, s.size() - 1);
		
		while (pos1 == pos2)
			pos1 = nextInt(0, s.size() - 1);

		Instance r1 = s.table.getInstance(s.get(pos1));
		Instance r2 = s.table.getInstance(s.get(pos2));

		return new Pair<Instance, Instance>(r1, r2);
	}
	
	private Pair<Instance, Instance> Pivots(Indices s1, Indices s2) {
		Instance r1 = null, r2 = null;
		int sz = s1.size() + s2.size() - 1;

		int pos1 = nextInt(0, sz);
		int pos2 = nextInt(0, sz);

		while (pos1 == pos2)
			pos1 = nextInt(0, sz);
		
		if ((pos1 < s1.size()) && (pos2 < s1.size())) {
			r1 = s1.table.getInstance(s1.get(pos1));
			r2 = s2.table.getInstance(s1.get(pos2));
		}
		else if ((pos1 >= s1.size()) && (pos2 >= s1.size())) {
			r1 = s2.table.getInstance(s2.get(pos1 - s1.size()));
			r2 = s2.table.getInstance(s2.get(pos2 - s1.size()));
		}
		else if ((pos1 < s1.size()) && (pos2 >= s1.size())) {
			r1 = s1.table.getInstance(s1.get(pos1));
			r2 = s2.table.getInstance(s2.get(pos2 - s1.size()));
		}		
		else if ((pos1 >= s1.size()) && (pos2 < s1.size())) {
			r1 = s2.table.getInstance(s2.get(pos1 - s1.size()));
			r2 = s1.table.getInstance(s1.get(pos2));
		}

		return new Pair<Instance, Instance>(r1, r2);
	}
	
	private Partitions Partition(Indices s, Instance p, double rho, double r) {
		Instance q;
		Indices L  = new Indices(s.table);
		Indices G  = new Indices(s.table);
		Indices Lw = new Indices(s.table);
		Indices Gw = new Indices(s.table);
		
		for (int i = 0; i < s.size(); i++) {
			q = s.table.getInstance(s.get(i));
			if (dist.compute(q.getValues(), p.getValues()) < rho) {
				L.add(s.get(i));
				if ((rho - r) <= dist.compute(q.getValues(), p.getValues()))
					Lw.add(s.get(i));
			} else {
				G.add(s.get(i));
				if (dist.compute(q.getValues(), p.getValues()) <= (rho + r))
					Gw.add(s.get(i));
			}
		}
		
		return new Partitions(L, G, Lw, Gw);
	}

	private void SelfJoinBF(Indices s, double r) {
		Instance p, q; 
		for (int i = 0; i < s.size(); i++) {
			p = s.table.getInstance(s.get(i));
			for (int j = i + 1; j < s.size(); j++) {
				q = s.table.getInstance(s.get(j));
				if (dist.compute(p.getValues(), q.getValues()) <= r) {
					collector.addPair(p.getId(), q.getId());
				}
			}
		}
	}
	
	private void JoinBF(Indices s1, Indices s2, double r) {
		Instance p, q; 
		for (int i = 0; i < s1.size(); i++) {
			p = s1.table.getInstance(s1.get(i));
			for (int j = 0; j < s2.size(); j++) {
				q = s2.table.getInstance(s2.get(j));
				if (dist.compute(p.getValues(), q.getValues()) <= r) {
					collector.addPair(p.getId(), q.getId());
				}
			}
		}
	}
	
	private void JoinPivots(Indices s1, Indices s2, double r) {
		// TODO: This implementation must be reviewed 
		// and possibly extracted to a interface containing 
		// a join method so other versions of joins 
		// could be easily implemented
		Instance p, q;
		int k = nextInt((int) 0.3 * s1.size(), (int) 0.6 * s1.size());
		double[] D = new double[k];
		double[][] P = new double[k][s2.size()];
		
		for (int i = 0; i <= k; i++) {
			for (int j = 0; j < s2.size(); j++) {
				p = s1.table.getInstance(s1.get(i));
				q = s2.table.getInstance(s2.get(j));
				P[i][j] = dist.compute(p.getValues(), q.getValues());
				if (P[i][j] <= r) {
					collector.addPair(p.getId(), q.getId());
				}
			}
		}
		
		for (int i = k + 1; i < s1.size(); i++) {
			for (int l = 0; l <= k; l++) {
				p = s1.table.getInstance(s1.get(l));
				q = s1.table.getInstance(s1.get(i)); 
				D[l] = dist.compute(p.getValues(), q.getValues());
			}
			for (int j = 0; j < s2.size(); j++) {
				boolean f = false;
				for (int l = 0; l <= k; l++) {
					if (Math.abs(P[l][j] - D[l]) > r) {
						f = true;
						break;
					}
				}
				p = s1.table.getInstance(s1.get(i));
				q = s2.table.getInstance(s2.get(j));
				double e = dist.compute(p.getValues(), q.getValues());
				if (!f && (e <= r)) {
					collector.addPair(p.getId(), q.getId());
				}
			}
		}
	}

	private void JoinDC(Indices s1, Indices s2, double r) {
		// TODO: This implementation must be reviewed 
		// and possibly extracted to a interface containing 
		// a join method so other versions of joins 
		// could be easily implemented
		double[] dl = new double[s2.table.getNumRows()];
		double[] du = new double[s2.table.getNumRows()];
		
		for (int j = 0; j < s2.size(); j++) {
			dl[j] = 0;
			du[j] = Double.MAX_VALUE;
		}
		
		Instance p, q;
		double e = Double.MAX_VALUE;
		for (int i = 0; i < s2.size(); i++) {
			
			if (i > 1) {
				p = s1.table.getInstance(s1.get(i));
				q = s1.table.getInstance(s1.get(i-1));
				e = dist.compute(p.getValues(), q.getValues());
			}
			
			for (int j = 0; j < s2.size(); j++) {
				dl[j] = Math.max(Math.max(e-du[j], dl[j]-e), 0);
				if (du[j] <= r) {
					collector.addPair(s1.get(i), s2.get(j));
				} else if (dl[j] <= r) {
					p = s1.table.getInstance(s1.get(i));
					q = s2.table.getInstance(s2.get(j));
					dl[j] = du[j] = dist.compute(p.getValues(), q.getValues());
					if (du[j] <= r) {
						collector.addPair(s1.get(i), s2.get(j));
					}
				} 
			}
		}
	}
	
	private int nextInt(int min, int max) {
		return rand.nextInt(max - min + 1) + min;
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
