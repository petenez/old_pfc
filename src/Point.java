package coordinator;
 
import java.util.ArrayList;

/*
Author: Petri Hirvonen, petri.hirvonen@aalto.fi
December 9, 2017
 */
 
public class Point {	// implements carbons and rings

	private int k;		// number of carbon/ring
	private double[] r;	// its position
	// for a carbon, its bros are its neighboring carbons
	// for a ring, its bros are its neighboring rings
	private ArrayList<Point> bros = new ArrayList<Point>();
	// for a carbon, its hos are its neighboring rings
	// for a ring, its hos are its neighboring carbons
	private ArrayList<Point> hos = new ArrayList<Point>();
	 
	public Point(int k, double[] r) {
		this.k = k;
		this.r = r;
	}
	
	public int k() {
		return k;
	}
	 
	public double[] r() {
		return r;
	}
	
	public void r(double[] r) {
		this.r = r;
	}
	
	public void add_bro(Point p) {
		bros.add(p);
	}
	
	public Point get_bro(int n) {
		return bros.get(n);
	}
	
	public boolean is_bro(Point p) {
		if(bros.contains(p)) return true;
		return false;
	}
	
	public int count_bros() {
		return bros.size();
	}
	
	public void add_ho(Point p) {
		hos.add(p);
	}
	
	public Point get_ho(int n) {
		return hos.get(n);
	}
	
	public boolean is_ho(Point p) {
		if(hos.contains(p)) return true;
		return false;
	}
	
	public int count_hos() {
		return hos.size();
	}
}