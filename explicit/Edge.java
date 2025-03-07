package explicit;
public class Edge {
	private int source;
	private double weight;

	public Edge(int source, double weight) {
		this.source = source;
		this.weight = weight;
	}

	public int getSource() {
		return this.source;
	}

	public double getWeight() {
		return this.weight;
	}
}
