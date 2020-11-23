package common;

public class Pair<A, B> {
    private final A val0;
    private final B val1;
    
    public Pair(A val0, B val1) {
    	this.val0 = val0;
    	this.val1 = val1;
    }
    
    public A getValue0() {
        return this.val0;
    }

    public B getValue1() {
        return this.val1;
    }
}
