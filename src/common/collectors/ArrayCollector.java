package common.collectors;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;

public class ArrayCollector implements ICollector {

	private ArrayList<Integer> R = new ArrayList<>();

	@Override
	public void addPair(Integer idx1, Integer idx2) {
		R.add(idx1);
		R.add(idx2);
	}

	@Override
	public Collection<Integer> getPairsOf(Integer idx) {
		ArrayList<Integer> pairs = new ArrayList<>();
		Iterator<Integer> iter = R.iterator();
		
		while (iter.hasNext()) {
			Integer value1 = iter.next();
			Integer value2 = iter.next();
			
			if (value1.equals(idx)) {
				pairs.add(value2);
			}
			
			if (value2.equals(idx)) {
				pairs.add(value1);
			}
		}
		return pairs;
	}

	@Override
	public void clear() {
		R.clear();
	}

}
