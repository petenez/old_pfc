package coordinator;

import java.io.*;
import java.util.*;

/*
Author: Petri Hirvonen, petri.hirvonen@aalto.fi
December 9, 2017
 */

//implements a "coordinator" for converting a PFC density field into atomic coordinates
public class Coordinator {
	
	private static int[] u = new int[] {1, 1, 0, -1, -1, -1, 0, 1};	// arrays for looping through ...
	private static int[] v = new int[] {0, 1, 1, 1, 0, -1, -1, -1};	// ... neighboring indices
	
	private int Lx;		// PFC grid dimensions
	private int Ly;
	private int Nx;		// local minima bin counts
	private int Ny;
	private double w;	// system dimensions (Å)
	private double h;
	private double lA;
	private double ux;	// length units
	private double uy;
	private ArrayList<Point> carbons = new ArrayList<Point>();	// list of carbons
	private ArrayList<Point> rings = new ArrayList<Point>();	// list of rings
	private ArrayList<ArrayList<ArrayList<Point>>> bins;		// local minima bins
	
	public Coordinator(int Lx, int Ly, double dx, double dy, double lPFC, double lA) {
		this.Lx = Lx;
		this.Ly = Ly;
		this.lA = lA;
		ux = dx*lA/lPFC;
		uy = dy*lA/lPFC;
		w = Lx*ux;
		h = Ly*uy;
		this.Nx = (int)Math.ceil(0.25*w);	// bin size ~4 Å
		this.Ny = (int)Math.ceil(0.25*h);
		bins = new ArrayList<ArrayList<ArrayList<Point>>>(Ny);	// initialize bins
		for(int m = 0; m < Ny; m++) {
			bins.add(new ArrayList<ArrayList<Point>>(Nx));
			for(int n = 0; n < Nx; n++) {
				bins.get(m).add(new ArrayList<Point>());
			}
		}
	}
	
	// wraps x-coordinates
	public int wx(int x) {
		if(x < 0) return x+Lx;
		if(x >= Lx) return x-Lx;
		return x;
	}
	
	// wraps y-coordinates
	public int wy(int y) {
		if(y < 0) return y+Ly;
		if(y >= Ly) return y-Ly;
		return y;
	}
	
	// wraps bin x-indices
	public int wnx(int nx) {
		if(nx < 0) return nx+Nx;
		if(nx >= Nx) return nx-Nx;
		return nx;
	}
	
	// wraps bin y-indices
	public int wny(int ny) {
		if(ny < 0) return ny+Ny;
		if(ny >= Ny) return ny-Ny;
		return ny;
	}
	
	// wraps a vector wrt. system dimensions
	public double[] w(double[] r) {
		double x = r[0];
		double y = r[1];
		if(x < 0.0) x += w;
		else if(x >= w) x -= w;
		if(y < 0.0) y += h;
		else if(y >= h) y -= h;
		return new double[] {x, y};
	}
	
	// returns the bin x-index corresponding to vector r
	public int nx(double[] r) {
		return (int)(r[0]/w*Nx);
	}
	
	// returns the bin y-index corresponding to vector r
	public int ny(double[] r) {
		return (int)(r[1]/h*Ny);
	}
	
	// returns the bin indexed by nx and ny
	public ArrayList<Point> bin(int nx, int ny) {
		return bins.get(ny).get(nx);
	}
	
	// returns the bin corresponding to vector r
	public ArrayList<Point> bin(double[] r) {
		return bins.get((int)ny(r)).get((int)nx(r));
	}
	
	// adds C as a member of the three rings, adds the three rings as neighbors if not
	// previously neighbors
	public void bros_n_hos(Point C, Point ith, Point jth, Point kth) {
		C.add_ho(ith);
		C.add_ho(jth);
		C.add_ho(kth);
		ith.add_ho(C);
		jth.add_ho(C);
		kth.add_ho(C);
		if(!ith.is_bro(jth)) {
			ith.add_bro(jth);
			jth.add_bro(ith);
		}
		if(!ith.is_bro(kth)) {
			ith.add_bro(kth);
			kth.add_bro(ith);
		}
		if(!jth.is_bro(kth)) {
			jth.add_bro(kth);
			kth.add_bro(jth);
		}
	}
	
	// Loads the PFC system, detects local minima, finds nearest-neighbor triplets of them,
	// places carbons in the center of each triplet, finds their neighbors and checks that
	// there's a reasonable number of neighbors.
	public void load(String input) throws IOException {
		double[][] psi = new double[Ly][Lx];
		BufferedReader reader = new BufferedReader(new FileReader(input));
		for(int j = 0; j < Ly; j++) {		// load PFC density data
			for(int i = 0; i < Lx; i++) {
				psi[j][i] = Double.parseDouble(reader.readLine());
			}
		}
		reader.close();
		Point p;
		for(int j = 0; j < Ly; j++) {
			outer:
			for(int i = 0; i < Lx; i++) {
				for(int n = 0; n < 8; n++) {
					// if neighbors density is lower, not a local minimum
					if(psi[wy(j+v[n])][wx(i+u[n])] < psi[j][i]) {
						continue outer;
					}
				}
				// density data might have limited precision whereby neighboring grid points
				// might have identical values, set current density and those of neighbors
				// very low to ensure uniqueness of current local minimum
				psi[j][i] -= 1000.0;
				for(int n = 0; n < 8; n++) {
					psi[wy(j+v[n])][wx(i+u[n])] -= 1000.0;
				}
				// use quadratic interpolation to estimate more accurate location
				double x0 = i-1;	// neighbor
				double x1 = i;	// current minimum
				double x2 = i+1;	// neighbor
				double x02 = x0*x0;	
				double x12 = x1*x1;
				double x22 = x2*x2;
				double p0 = psi[j][wx(i-1)];	// periodic boundary conditions
				double p1 = psi[j][i];
				double p2 = psi[j][wx(i+1)];
				double A = (x2*(p1-p0)+x1*(p0-p2)+x0*(p2-p1));	// y = A*x^2 + B*x + C
				double B = (x22*(p0-p1)+x12*(p2-p0)+x02*(p1-p2));
				double x = -B/(2.0*A)*ux;	// extremum at x = -B/(2A)
				if(Math.abs(A) < 1.0e-14) x = i;	// linear fit
				x0 = j-1;	// same in the other direction
				x1 = j;
				x2 = j+1;
				x02 = x0*x0;
				x12 = x1*x1;
				x22 = x2*x2;
				p0 = psi[wy(j-1)][i];
				p1 = psi[j][i];
				p2 = psi[wy(j+1)][i];
				A = (x2*(p1-p0)+x1*(p0-p2)+x0*(p2-p1));
				B = (x22*(p0-p1)+x12*(p2-p0)+x02*(p1-p2));
				double y = -B/(2.0*A)*uy;
				if(Math.abs(A) < 1.0e-14) y = j;
				p = new Point(rings.size(), w(new double[] {x, y}));
				rings.add(p);		// add local minimum to corresponding list
				bin(p.r()).add(p);	// add to corresponding bin
			}
		}
		// Find triplets of local minima that are closest neighbors to each other.
		// Practically find the triple points of three local minima (where each one
		// is equally distant) and make sure that no fourth carbon atom is closer
		// to it. Place a carbon atom at the average of the three minimas'
		// coordinates. Find the neighbors of the carbon atoms and minima (which
		// actually correspond to carbon atom rings) and check that the bonding
		// makes sense.
		Point ith, jth, kth, lth, C;	// 'Points' for minima and carbons
		ArrayList<Point> binj, bink, binl;	// bins
		int nx, ny;								// indices for bins
		double a, d2;
		double[] ij, ijp, ik, ikp, it;			// vectors
		for(int i = 0; i < rings.size(); i++) {	// go through rings
			ith = rings.get(i);
			nx = nx(ith.r());
			ny = ny(ith.r());
			for(int vj = -1; vj < 2; vj++) {	// go through its surrounding bins
				for(int uj = -1; uj < 2; uj++) {
					binj = bin(wnx(nx+uj), wny(ny+vj));
					for(int j = 0; j < binj.size(); j++) {	// go through minima in current bin
						jth = binj.get(j);
						if(jth.k() <= ith.k()) continue;	// consider each triplet only once
						// compute some vectors to find the triple point
						ij = diff(jth.r(), ith.r());
						if(V.abs(ij) > 2.0*lA) continue;	// closest neighbors can't be farther apart
						ijp = V.x(ij);
						for(int vk = -1; vk < 2; vk++) {	// again check bins
							for(int uk = -1; uk < 2; uk++) {
								bink = bin(wnx(nx+uk), wny(ny+vk));
								k_loop:
								// rings in this bin
								for(int k = 0; k < bink.size(); k++) {
									kth = bink.get(k);
									if(kth.k() <= jth.k()) continue;
									// compute vectors
									ik = diff(kth.r(), ith.r());
									if(V.abs(ik) > 2.0*lA) continue;
									if(V.abs(diff(kth.r(), jth.r())) > 2.0*lA) continue;
									ikp = V.x(ik);
									a = 0.5/ijp[0]*(ik[0]-ij[0]+ikp[0]/ikp[1]*(ij[1]-ik[1]))/(1.0-ijp[1]*ikp[0]/(ijp[0]*ikp[1]));
									it = V.sum(V.mul(0.5, ij), V.mul(a, ijp));	// position of triple point (from ith minimum)
									d2 = V.abs2(it);
									// check that no fourth minimum is closer to triple point
									for(int vl = -1; vl < 2; vl++) {
										for(int ul = -1; ul < 2; ul++) {
											binl = bin(wnx(nx+ul), wny(ny+vl));
											for(int l = 0; l < binl.size(); l++) {
												lth = binl.get(l);
												if(lth == ith || lth == jth || lth == kth) continue;
												if(V.abs2(V.sub(diff(lth.r(), ith.r()), it)) < d2) continue k_loop;
											}
										}
									}
									// add carbon atom, update neighbors
									C = new Point(carbons.size(), w(V.sum(ith.r(), V.mul(1.0/3.0, V.sum(ij, ik)))));
									carbons.add(C);
									bros_n_hos(C, ith, jth, kth);
								}
							}
						}
					}
				}
			}
		}
		// find neighboring carbons of each carbon
		Point D;
		for(int c = 0; c < carbons.size(); c++) {	// all carbons
			C = carbons.get(c);
			for(int i = 0; i < C.count_hos(); i++) {	// ith ring its a member of
				ith = C.get_ho(i);
				for(int j = i+1; j < C.count_hos(); j++) {	// jth ring its a member of
					jth = C.get_ho(j);
					for(int k = 0; k < ith.count_hos(); k++) {	// kth carbon of ith ring
						D = ith.get_ho(k);
						// if its not C, but shares two rings with it and is not yet a neighbor ...
						if(D != C && D.is_ho(jth) && !C.is_bro(D)) {
							C.add_bro(D);	// ... add as neighbors
							D.add_bro(C);
						}
					}
				}
			}
		}
		// check if every carbon has three neighbors and if each ring has 5-7 neighboring rings
		for(int c = 0; c < carbons.size(); c++) {
			C = carbons.get(c);
			if(C.count_bros() != 3) {
				System.out.println("Warning: The "+c+"th carbon has "+C.count_bros()+" bonds!");
			}
			if(C.count_hos() != 3) {
				System.out.println("Warning: The "+c+"th carbon is part of "+C.count_hos()+" rings!");
			}
		}
		for(int i = 0; i < rings.size(); i++) {
			ith = rings.get(i);
			if(ith.count_bros() < 5 || ith.count_bros() > 7) {
				System.out.println("Warning: The "+i+"th ring is a "+ith.count_bros()+"-gon!");
			}
			if(ith.count_bros() != ith.count_hos()) {
				throw new RuntimeException("Error: For "+i+"th ring, the numbers of carbons and neighboring rings don't match!");
			}
		}
	}
	
	// computes the difference vector from a to b taking periodic boundary conditions into account
	public double[] diff(double[] b, double[] a) {
		double[] ab = V.sub(b, a);
		if(ab[0] > 0.5*w) ab[0] -= w;
		else if(ab[0] < -0.5*w) ab[0] += w;
		if(ab[1] > 0.5*h) ab[1] -= h;
		else if(ab[1] < -0.5*h) ab[1] += h;
		return ab;
	}
	
	// writes out the carbon atom data
	// format: number of carbon, its coordinates and numbers of its neighbors
	public void write_carbons(String output) throws IOException {
		Point C, D;
		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		for(int c = 0; c < carbons.size(); c++) {
			C = carbons.get(c);
			writer.write(C.k()+" "+C.r()[0]+" "+C.r()[1]);
			for(int d = 0; d < C.count_bros(); d++) {
				D = C.get_bro(d);
				writer.write(" "+D.k());
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	// writes out the data of nonhexagonal rings
	// format: numbers of the member carbons
	public void write_nonhexagons(String output) throws IOException {
		Point ring, C;
		BufferedWriter writer = new BufferedWriter(new FileWriter(output));
		for(int n = 0; n < rings.size(); n++) {
			ring = rings.get(n);
			if(ring.count_bros() == 6) continue;
			for(int c = 0; c < ring.count_hos(); c++) {
				C = ring.get_ho(c);
				writer.write(C.k()+" ");
			}
			writer.write("\n");
		}
		writer.close();
	}
	
	public static void main(String[] args) throws IOException {
		String input = args[0];
		int Lx = Integer.parseInt(args[1]);
		int Ly = Integer.parseInt(args[2]);
		double dx = Double.parseDouble(args[3]);
		double dy = Double.parseDouble(args[4]);
		double lPFC = Double.parseDouble(args[5]);
		double lA = Double.parseDouble(args[6]);
		String cout = args[7];
		String nhout = args[8];
		Coordinator coordinator = new Coordinator(Lx, Ly, dx, dy, lPFC, lA);
		coordinator.load(input);
		coordinator.write_carbons(cout);
		coordinator.write_nonhexagons(nhout);
	}
}
