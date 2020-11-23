package common.collectors;

import java.util.Collection;

public interface ICollector {
	public void addPair(Integer idx1, Integer idx2);
	public Collection<Integer> getPairsOf(Integer idx);
	public void clear();
}
