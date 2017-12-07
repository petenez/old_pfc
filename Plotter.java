package visualization;

/*
Author: Petri Hirvonen, petri.hirvonen@aalto.fi
December 7, 2017
*/

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import javax.imageio.ImageIO;

public class Plotter {
	
	// maps linearly x from between [min, max] to [0, 1]
	public static double map(double x, double min, double max) {
		return (x-min)/(max-min);
	}
	 
	public static void main(String[] args) throws IOException {
		String input;
		String output;
		int Lx;	// width of the system in pixels
		int Ly;	// height
		// parameters for optional looping over time steps
		int start;
		int incr;
		int end;
		int log;	// apply a logarithmic transform?
		if(!(args.length == 4 || args.length == 5 || args.length == 7 || args.length == 8)) {
			System.out.println("Wrong syntax!");
			return;
		}
		
		input = args[0];
		output = args[1];
		Lx = Integer.parseInt(args[2]);
		Ly = Integer.parseInt(args[3]);
		if(args.length == 7 || args.length == 8) {
			start = Integer.parseInt(args[4]);
			incr = Integer.parseInt(args[5]);
			end = Integer.parseInt(args[6]);
		} else {
			start = 0;
			incr = 1;
			end = 0;
		}
		if(args.length == 5 || args.length == 8) log = Integer.parseInt(args[args.length-1]);
		else log = 0;	// log-transform not applied
		
		boolean real;	// real-valued data?
		
		// loop over different time steps
		for(int m = start; m <= end; m += incr) {
			input = args[0];
			output = args[1];
			if(args.length == 7 || args.length == 8) {	// loop over different time steps
				String[] pieces = input.split("#");	// split 'input' at "#"
				input = "";	// input emptied
				if(pieces.length >= 1) input += pieces[0];	// append piece before "#"
				input += m;	// append time step ("#" can be used as a wild card!)
				if(pieces.length >= 2) input += pieces[1];	// append piece after "#" 
				pieces = output.split("#");	// repeat for 'output'
				output = "";
				if(pieces.length >= 1) output += pieces[0];
				output += m;
				if(pieces.length >= 2) output += pieces[1];
			}
			if(!input.contains(".")) input += ".dat";	// no extension specified? -> use .dat
			if(!output.contains(".")) output += ".png";	// similarly for output
			File f = new File(input);
			if(!f.exists()) {	// check if input file is found
				System.out.println("File "+input+" not found!");
				continue;
			}
			
			// checks if the data is real or complex and determines the extrema
			BufferedReader reader = new BufferedReader(new FileReader(input));
			int parts = reader.readLine().split(" ").length;	// splits the line read at " "
			if(parts == 1) real = true;	// single column -> real
			else if(parts == 2) real = false;	// double -> complex
			else {
				System.out.println("Invalid data!");
				reader.close();
				return;
			}
			reader.close();
			reader = new BufferedReader(new FileReader(input));	// start again
			double[][][] data = new double[Ly][Lx][2];	// data array
			double min = 1.0e100;	// initially huge just to be safe...
			double max = -1.0e100;	// initially -huge ...
			if(real) {
				double a;
				for(int j = 0; j < Ly; j++) {
					for(int i = 0; i < Lx; i++) {
						a = Double.parseDouble(reader.readLine());
						data[j][i][0] = a;
						if(a < min) min = a;	// check if current minimum
						else if(a > max) max = a;	// maximum
					}
				}
			} else {
				min = 0.0;	// minimum of amplitude fixed at 0
				String[] temp;
				double a, b, A;
				for(int j = 0; j < Ly; j++) {
					for(int i = 0; i < Lx; i++) {
						temp = reader.readLine().split(" ");
						a = Double.parseDouble(temp[0]);
						data[j][i][0] = a;
						b = Double.parseDouble(temp[1]);
						data[j][i][1] = b;
						A = Math.sqrt(a*a+b*b);	// amplitude
						if(A > max) max = A;
					}
				}
			}
			reader.close();

			// pixel manipulations and writing output
			BufferedImage image = new BufferedImage(Lx, Ly, BufferedImage.TYPE_INT_RGB);
			double a, b, A;
			float arg;
			if(real) {	// real data
				if(log > 0) {	// logarithmic transform applied
					for(int j = 0; j < Ly; j++) {
						for(int i = 0; i < Lx; i++) {
							a = data[j][i][0];
							// Ly-j-1: image coordinates defined upside down
							// RGB & HSB: see Wikipedia for RGB and HSV/HSL color models
							// mapping: data is first shifted from between [min, max] to [1, max-min+1],
							// 			([0, max-min] -> log(0) = -inf) and then a logarithm is taken 
							image.setRGB(i, Ly-j-1, Color.HSBtoRGB(0, 0, (float)map(Math.log(a-min+1.0), 0.0, Math.log(max-min+1.0))));
						}
					}
				} else {	// log-transform not applied
					for(int j = 0; j < Ly; j++) {
						for(int i = 0; i < Lx; i++) {
							image.setRGB(i, Ly-j-1, Color.HSBtoRGB(0, 0, (float)map(data[j][i][0], min, max)));
						}
					}
				}
			} else {	// complex data
				double div2pi = 0.5/Math.PI;	// phase between [-pi, pi] -> [0, 1]
				if(log != 0) {
					for(int j = 0; j < Ly; j++) {
						for(int i = 0; i < Lx; i++) {
							a = data[j][i][0];
							b = data[j][i][1];
							A = Math.sqrt(a*a+b*b);	// amplitude -> brightness
							arg = (float)(div2pi*Math.atan2(b, a)+0.5);	// phase -> hue
							image.setRGB(i, Ly-j-1, Color.HSBtoRGB(arg, 1, (float)map(Math.log(A+1.0), 0.0, Math.log(max+1.0))));
						}
					}
				} else {	// log-transform not applied
					for(int j = 0; j < Ly; j++) {
						for(int i = 0; i < Lx; i++) {
							a = data[j][i][0];
							b = data[j][i][1];
							A = Math.sqrt(a*a+b*b);
							arg = (float)(div2pi*Math.atan2(b, a)+0.5);
							image.setRGB(i, Ly-j-1, Color.HSBtoRGB(arg, 1, (float)map(A, 0.0, max)));
						}
					}
				}
			}
			ImageIO.write(image, "png", new File(output));
		}
	}
}
