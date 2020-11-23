package algorithms;

import common.Table;
import common.collectors.ICollector;

public interface IRangeJoin {
	// similarity self-join operation
	public void range(Table a, double r);
	// similarity join operation
	public void range(Table a, Table b, double r);
	// reset collector
	public void reset();
	// get collector containing similar pairs
	public ICollector getCollector();
}
