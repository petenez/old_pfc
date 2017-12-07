/*
Author: Petri Hirvonen, petri.hirvonen@aalto.fi
December 7, 2017
*/

#include <complex.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define pi 3.14159265358979323846264338327

// contains stuff related to output
struct Output {
	char name[128];		// run name
	int T_print;		// interval for printing output
	int T_write;		// interval for saving state
};

// contains stuff related to data arrays
struct Arrays {
	int W;				// system width
	int H;				// system height
	ptrdiff_t lH;		// local system height
	ptrdiff_t lW;		// local system width
	ptrdiff_t lh0;		// local vertical start index
	ptrdiff_t lw0;		// local horizontal start index
	
	double *A;		// operator for linear part, e^{-k^2 \hat{\mathcal{L}} \Delta t}
	double *B;		// operator for nonlinear part ...
	double *p;		// array for \psi_C
	double *q;		// another array
	
	fftw_plan p_P;	// FFTW plan for F(p)
	fftw_plan q_Q;	// F(q)
	fftw_plan Q_q;	// F^-1(Q)
};

// contains model parameters
struct Model {
	double alpha;
	double beta;
	double gamma;
	double delta;
	double c;
};

// contains stuff related to relaxation
struct Relaxation {
	time_t t0;			// start time
	int id;				// rank of process
	int ID;				// number of processes
	int t;				// time step count
	double f;			// average free energy density
	double p;			// average density
	double d;			// sampling step size for calculation box optimization
	int T;				// total number of iterations
	double dx;			// x-discretization
	double dy;			// y-
	double dt;			// time step
	int T_optimize;		// interval for calculation box optimization
};

// seeds processes' random number generators
void seed_rngs(struct Relaxation *relaxation, FILE *input) {
	int seed;								// seed for next process
	// seed for 0th process from input file
	if(fscanf(input, " %d", &seed) != 1) {
		printf("Invalid input for random number generator seed!\n");
		exit(1);
	}
	if(relaxation->id == 0) {
		if(seed == 0) srand(time(NULL));	// random seed
		else srand(seed);					// user-specified seed
		seed = rand();						// sample new seed for next process
	}
	else {
		MPI_Recv(&seed, 1, MPI_INT, relaxation->id-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);	// receive from previous
		srand(seed);																			// seed
		seed = rand();																			// new seed for next
	}
	if(relaxation->id != relaxation->ID-1) MPI_Send(&seed, 1, MPI_INT, relaxation->id+1, 0, MPI_COMM_WORLD);	// send to next
}

// configures arrays and FFTW plans
void configure_arrays(struct Arrays *arrays, FILE *input) {
	if(fscanf(input, " %d %d", &arrays->W, &arrays->H) != 2) {
		printf("Invalid input for array dimensions!\n");
		exit(1);
	}

	int lWHh = fftw_mpi_local_size_2d_transposed(arrays->H, arrays->W/2+1, MPI_COMM_WORLD, &arrays->lH, &arrays->lh0, &arrays->lW, &arrays->lw0);		// local data size halved
	int lWHp = 2*lWHh;	// local data size padded

	arrays->A = (double*)fftw_malloc(lWHh*sizeof(double));	// allocate arrays
	arrays->B = (double*)fftw_malloc(lWHh*sizeof(double));
	
	arrays->p = (double*)fftw_malloc(lWHp*sizeof(double));
	arrays->q = (double*)fftw_malloc(lWHp*sizeof(double));
	
	arrays->p_P = fftw_mpi_plan_dft_r2c_2d(arrays->H, arrays->W, arrays->p, (fftw_complex*)arrays->p, MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);	// set up FFTW plans
	arrays->q_Q = fftw_mpi_plan_dft_r2c_2d(arrays->H, arrays->W, arrays->q, (fftw_complex*)arrays->q, MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_OUT);
	arrays->Q_q = fftw_mpi_plan_dft_c2r_2d(arrays->H, arrays->W, (fftw_complex*)arrays->q, arrays->q, MPI_COMM_WORLD, FFTW_MEASURE|FFTW_MPI_TRANSPOSED_IN);
	
}

// configures output
void configure_output(struct Output *output, FILE *input) {
	// read intervals for printing and saving from input file
	if(fscanf(input, " %d %d", &output->T_print, &output->T_write) != 2) {
		printf("Invalid input for output configuration!\n");
		exit(1);
	}
}

// configures model
void configure_model(struct Model *model, FILE *input) {
	// read model parameters from input file
	if(fscanf(input, " %lf %lf %lf %lf %lf", &model->alpha, &model->beta, &model->gamma, &model->delta, &model->c) != 5) {
		printf("Invalid input for model parameters!\n");
		exit(1);
	}
}

// configures relaxation
void configure_relaxation(struct Relaxation *relaxation, FILE *input) {
	// read from input file
	if(fscanf(input, " %d %lf %lf %lf %d", &relaxation->T, &relaxation->dx, &relaxation->dy, &relaxation->dt, &relaxation->T_optimize) != 5) {
		printf("Invalid input for relaxation settings!\n");
		exit(1);
	}
}

// initializes the density field with random noise
// p and A are the average density and noise amplitude
void random_state(struct Arrays *arrays, double p, double A) {
	int Wp = 2*(arrays->W/2+1);
	int w, h, k;
	for(h = 0; h < arrays->lH; h++) {
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			// average + amplitude * random from U(-0.5, 0.5)
			arrays->q[k+w] = p+A*((double)rand()/RAND_MAX-0.5);
		}
	}
}

// returns the one-mode approximation for close-packed lattice
// x0, y0 give rotation (theta) axis; u, v give translation; l0 dimensionless lattice constant
double OMA(double x0, double y0, double u, double v, double l0, double theta) {
	double qx, qy;
	qx = 2.0*pi/l0;						// wave numbers
	qy = qx/sqrt(3.0);
	double cosine = cos(theta);
	double sine = sin(theta);
	double x = cosine*x0-sine*y0+u;		// rotation matrix applied
	double y = sine*x0+cosine*y0+v;
	return cos(qx*x)*cos(qy*y)+0.5*cos(2.0*qy*y);
}

// initializes the density field with random crystallites
// dx, dy and l0 control the lenght scale
// p and A are the average density and noise amplitude
// N is number of crystallites and R their radius
void polycrystalline_state(struct Arrays *arrays, double dx, double dy, double l0, double p, double A, int N, double R) {
	double crystallites[3*N];	// array for crystallite coordinates and orientations
	int w, h, gh, k, n, id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	if(id == 0) {						   // 0th process samples crystallites
		for(n = 0; n < N; n++) {
			crystallites[3*n] = (double)rand()/RAND_MAX*arrays->W;
			crystallites[3*n+1] = (double)rand()/RAND_MAX*arrays->H;
			crystallites[3*n+2] = 2.0*pi*rand()/RAND_MAX;
		}
	}
	// broadcasted to other processes
	MPI_Bcast(crystallites, 3*N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	int Wp = 2*(arrays->W/2+1);
	// x-, y- and squared distances from crystallite center, Gaussian envelope
	double u, v, r2, G;
	double div = 0.5/(R*R*l0*l0);
	int min;	// index for closest crystallite
	double r2min;	// squared distance for closest crystallite
	double umin = 0.0;	// u and v for closest crystallite
	double vmin = 0.0;
	for(h = 0; h < arrays->lH; h++) {
		gh = arrays->lh0+h;	// global vertical index
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			arrays->q[k+w] = p;	// (1) set to average density
			min = -1;
			r2min = 1.0e100;
			for(n = 0; n < N; n++) {
				u = w-crystallites[3*n];
				// closest periodic image searched
				if(u > 0.5*arrays->W) u -= arrays->W;
				else if(u < -0.5*arrays->W) u += arrays->W;
				u *= dx;
				v = gh-crystallites[3*n+1];
				if(v > 0.5*arrays->H) v -= arrays->H;
				else if(v < -0.5*arrays->H) v += arrays->H;
				v *= dy;
				r2 = u*u+v*v;
				if(r2 < r2min) {	// closest crystallite so far?
					min = n;
					r2min = r2;
					umin = u;
					vmin = v;
				}
			}
			G = exp(-r2min*div);
			// (2) superimpose crystallites
			arrays->q[k+w] += A*G*OMA(umin, vmin, 0.0, 0.0, l0, crystallites[3*min+2]);
		}
	}
}

// initializes the density field with data from a data file
void read_state(struct Arrays *arrays, FILE *file, double p, double A) {
	int Wp = 2*(arrays->W/2+1);
	int w, h, k;
	double dummy;
	double p0 = 0.0;
	double Amin = 1.0e100;
	double Amax = -1.0e100;
	for(h = 0; h < arrays->H; h++) {
		// data corresponds to local chunk?
		if(h >= arrays->lh0 && h < arrays->lh0+arrays->lH) {
			k = Wp*(h-arrays->lh0);
			for(w = 0; w < arrays->W; w++) {
				if(fscanf(file, "%lf", &arrays->q[k+w]) != 1) {
					printf("Invalid data for initialization!\n");
					exit(2);
				}
				p0 += arrays->q[k+w];
				if(arrays->q[k+w] < Amin) Amin = arrays->q[k+w];
				if(arrays->q[k+w] > Amax) Amax = arrays->q[k+w];
			}
		}
		else {
			for(w = 0; w < arrays->W; w++) {
				// fast-forward over other processes' chunks
				if(fscanf(file, "%lf", &dummy) != 1) {
					printf("Invalid data for initialization!\n");
					exit(2);
				}
			}
		}
	}
	p0 /= arrays->W*arrays->lH;
	for(h = 0; h < arrays->lH; h++) {
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			arrays->q[k+w] = (arrays->q[k+w] - p0)/(Amax - Amin)*A + p;
		}
	}
}

// initializes the density field
void initialize_system(struct Arrays *arrays, FILE *input) {
	int init;
	if(fscanf(input, " %d", &init) != 1) {
		printf("Invalid initialization type!\n");
		exit(1);
	}
	if(init == 0) {	// random initialization
		double p, A;
		if(fscanf(input, " %lf %lf", &p, &A) != 2) {
			printf("Invalid initialization parameters!\n");
			exit(1);
		}
		random_state(arrays, p, A);
	}
	else if(init == 1) {	// polycrystalline initialization
		int N;
		double dx, dy, l0, p, A, R;
		if(fscanf(input, " %lf %lf %lf %lf %lf %d %lf", &dx, &dy, &l0, &p, &A, &N, &R) != 7) {
			printf("Invalid initialization parameters!\n");
			exit(1);
		}
		polycrystalline_state(arrays, dx, dy, l0, p, A, N, R);
	}
	else if(init == 2) {			// initial state read from a data file
		char filename[128];
		double p, A;
		if(fscanf(input, " %s %lf %lf", filename, &p, &A) != 3) {
			printf("Invalid data file name!\n");
			exit(1);
		}
		FILE *file = fopen(filename, "r");
		read_state(arrays, file, p, A);
	}
}

// saves state into a file
void write_state(struct Arrays *arrays, struct Relaxation *relaxation, struct Output *output) {
	char filename[128];			// filename
	sprintf(filename, "%s-t:%d.dat", output->name, relaxation->t);
	FILE *file;					// file stream
	int Wp = 2*(arrays->W/2+1);
	int i, w, h, k;
	for(i = 0; i < relaxation->ID; i++) {	// go through ranks
		MPI_Barrier(MPI_COMM_WORLD);		// makes every process wait until previous has completed writing output
		if(relaxation->id == i) {			// ith process continues, others go back waiting
			if(relaxation->id == 0) file = fopen(filename, "w");	// 0th process overwrites possible previous data
			else file = fopen(filename, "a");						// others append to end
			for(h = 0; h < arrays->lH; h++) {
				k = Wp*h;
				for(w = 0; w < arrays->W; w++) {
					fprintf(file, "%e\n", arrays->q[k+w]);
				}
			}
			fclose(file);		// close file stream
		}
	}
}

// prints output, format:
// time step count, elapsed wall-clock time (s), dx (dimensionless*), dy (*), average free energy density (*), average density (*)
void print(struct Relaxation *relaxation, struct Output *output) {
	if(relaxation->id == 0) {			// only 0th process
		printf("%d %d %.9lf %.9lf %.9lf %.9lf\n", relaxation->t, (int)(time(NULL)-relaxation->t0), relaxation->dx, relaxation->dy, relaxation->f, relaxation->p);
		char filename[128];
		sprintf(filename, "%s.out", output->name);
		FILE *file;
		file = fopen(filename, "a");	// need only append, empty file was already generated
		fprintf(file, "%d %d %.9lf %.9lf %.9lf %.9lf\n", relaxation->t, (int)(time(NULL)-relaxation->t0), relaxation->dx, relaxation->dy, relaxation->f, relaxation->p);
		fclose(file);
	}
}

// updates operators for linear and nonlinear parts
void update_AB(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation) {
	int w, h, gw, k;
	double ky, kx2, k2;									// (squared) k-vector (components)
	double divWH = 1.0/arrays->W/arrays->H;				// 1/(W*H)
	double dkx = 2.0*pi/(relaxation->dx*arrays->W);		// discretization in k-space
	double dky = 2.0*pi/(relaxation->dy*arrays->H);
	double d2, l, expl;									// (1+\nabla^2)^2, variables for intermediate results
	double c = model->c;
	for(w = 0; w < arrays->lW; w++) {
		gw = arrays->lw0+w;								// global x-index
		kx2 = gw*dkx;									// kx
		kx2 *= kx2;										// kx^2
		k = arrays->H*w;
		for(h = 0; h < arrays->H; h++) {
			if(h < arrays->H/2) ky = h*dky;				// upper half negative ky-values
			else ky = (h-arrays->H)*dky;
			k2 = kx2+ky*ky;								// k^2
			d2 = (1.0-k2);
			d2 *= d2;
			l = model->alpha+model->beta*d2;			// alpha+beta*(1+\nabla^2)^2 in k-space
			expl = exp(-((1.0-c)+c*k2)*l*relaxation->dt);			// e^{-k^2 \hat{\mathcal{L}} \Delta t}
			arrays->A[k+h] = expl;
			if(l == 0.0) arrays->B[k+h] = -((1.0-c)+c*k2)*relaxation->dt;		// avoid divide by zero
			else arrays->B[k+h] = (expl-1.0)/l;
			arrays->B[k+h] *= divWH;					// descaling
		}
	}
}

// scales p arrays (in k-space) by 1/(W*H) ((I)FFTs cause scaling of data by sqrt(W*H))
void scale_P(struct Arrays *arrays) {
	fftw_complex *P;			// complex data pointers
	P = (fftw_complex*)&arrays->p[0];
	double divWH = 1.0/arrays->W/arrays->H;
	int i;
	int lA = arrays->lW*arrays->H;
	for(i = 0; i < lA; i++) P[i] *= divWH;
}

// scales q arrays (in k-space) by 1/(W*H)
void scale_Q(struct Arrays *arrays) {
	fftw_complex *Q;
	Q = (fftw_complex*)&arrays->q[0];
	double divWH = 1.0/arrays->W/arrays->H;
	int i;
	int lA = arrays->lW*arrays->H;
	for(i = 0; i < lA; i++) Q[i] *= divWH;
}

// computes average free energy density (f) and average density (p) (f & p -> fp)
void fp(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation) {
	int Wp = 2*(arrays->W/2+1);
	double dkx = 2.0*pi/(relaxation->dx*arrays->W);
	double dky = 2.0*pi/(relaxation->dy*arrays->H);
	memcpy(arrays->p, arrays->q, Wp*arrays->lH*sizeof(double));		// save current state to p array
	fftw_execute(arrays->q_Q);
	scale_Q(arrays);
	double kx2, ky, k2, d2;
	fftw_complex *Q = (fftw_complex*)&arrays->q[0];
	int w, h, gw, k;
	for(w = 0; w < arrays->lW; w++) {
		gw = arrays->lw0+w;
		kx2 = gw*dkx;
		kx2 *= kx2;
		k = arrays->H*w;
		for(h = 0; h < arrays->H; h++) {
			if(h < arrays->H/2) ky = h*dky;
			else ky = (h-arrays->H)*dky;
			k2 = kx2+ky*ky;
			d2 = (1.0-k2);
			d2 *= d2;
			Q[k+h] *= d2;		// compute (1-\nabla^2)^2 \psi_C in k-space
		}
	}
	fftw_execute(arrays->Q_q);
	relaxation->f = 0.0;		// reset variables for average (free energy) density
	relaxation->p = 0.0;
	double div3 = 1.0/3.0;
	double p, q, p2;				// variables to simplify expressions
	double alpha = model->alpha;	// variables for model parameters to simplify expressions
	double beta = model->beta;
	double gamma = model->gamma;
	double delta = model->delta;
	for(h = 0; h < arrays->lH; h++) {
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			p = arrays->p[k+w];
			q = arrays->q[k+w];
			p2 = p*p;
			// average free energy density
			relaxation->f += 0.5*( alpha*p2 + beta*p*q ) + div3*gamma*p*p2 + 0.25*delta*p2*p2;
			// restore q arrays
			arrays->q[k+w] = p;
			// average densities
			relaxation->p += p;
		}
	}
	// communication between processes, note that can't reduce as average because chunks may have different sizes!
	MPI_Allreduce(MPI_IN_PLACE, &relaxation->f, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	MPI_Allreduce(MPI_IN_PLACE, &relaxation->p, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	double divWH = 1.0/arrays->W/arrays->H;
	relaxation->f *= divWH;
	relaxation->p *= divWH;
	fftw_execute(arrays->p_P);
	scale_P(arrays);
}

// performs one iteration of the semi-implicit spectral method
void step(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation) {
	int w, h, k;
	int Wp = 2*(arrays->W/2+1);
	if(relaxation->t%100 == 0) {	// avoids numerical instability (I guess there's some accumulating numerical error, don't know really...)
		memcpy(arrays->p, arrays->q, Wp*arrays->lH*sizeof(double));
		fftw_execute(arrays->p_P);	// update k-space solution fully from real-space solution
		scale_P(arrays);
	}
	double q, q2;
	double gamma = model->gamma;
	double delta = model->delta;
	for(h = 0; h < arrays->lH; h++) {
		k = Wp*h;
		for(w = 0; w < arrays->W; w++) {
			q = arrays->q[k+w];
			q2 = q*q;
			// compute nonlinear part (C)
			arrays->q[k+w] = gamma*q2 + delta*q*q2;
		}
	}
	fftw_execute(arrays->q_Q);
	fftw_complex *P = (fftw_complex*)&arrays->p[0];
	fftw_complex *Q = (fftw_complex*)&arrays->q[0];
	for(w = 0; w < arrays->lW; w++) {
		k = arrays->H*w;
		for(h = 0; h < arrays->H; h++) {
			P[k+h] = arrays->A[k+h]*P[k+h] + arrays->B[k+h]*Q[k+h];
			// copy p's to q's
			Q[k+h] = P[k+h];
		}
	}
	fftw_execute(arrays->Q_q);
}

// optimizes calculation box size for system
// samples slightly different box sizes by varying dx and dy and interpolates optimum quadratically
void optimize(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation) {
	// sample sizes
	double dxs[] = {relaxation->dx, relaxation->dx-relaxation->d, relaxation->dx+relaxation->d, relaxation->dx, relaxation->dx};
	double dys[] = {relaxation->dy, relaxation->dy, relaxation->dy, relaxation->dy-relaxation->d, relaxation->dy+relaxation->d};
	double fs[5];	// average free energy densities of sizes sampled
	int i;
	for(i = 0; i < 5; i++) {		// sample sizes
		relaxation->dx = dxs[i];
		relaxation->dy = dys[i];
		fp(arrays, model, relaxation);				// compute energy (and density (not needed here though))
		fs[i] = relaxation->f;	// save average free energy density
	}
	// interpolate new dx and dy (minimum of surface f = a+b*dx+c*dy+d*dx^2+e*dy^2)
	relaxation->dx = (relaxation->d*(fs[1]-fs[2])+2.0*dxs[0]*(-2.0*fs[0]+fs[1]+fs[2]))/(2.0*(fs[1]-2.0*fs[0]+fs[2]));
	relaxation->dy = (relaxation->d*(fs[3]-fs[4])+2.0*dys[0]*(-2.0*fs[0]+fs[3]+fs[4]))/(2.0*(fs[3]-2.0*fs[0]+fs[4]));
	// check that change in discretization is acceptable
	double l0 = 4.0*pi/sqrt(3.0);						// approximate dimensionless lenght scale (lattice constant)
	double dw = arrays->W*(relaxation->dx-dxs[0]);		// change in horizontal system size
	double dh = arrays->H*(relaxation->dy-dys[0]);		// ... vertical ...
	double dr = sqrt(dw*dw+dh*dh);						// "change vector"
	double x = 0.25*l0;									// limit to 1/4 of lattice constant (to ensure stability)
	if(dr > x) {	// if the change in system dimensions exceeds 1/4 of the lattice constant ...
		x /= dr;
		relaxation->dx = x*relaxation->dx+(1.0-x)*dxs[0];	// ... truncate the change to 1/4 of the lattice constant
		relaxation->dy = x*relaxation->dy+(1.0-x)*dys[0];
	}
	// update A and B
	update_AB(arrays, model, relaxation);	// dx and dy changed -> need to update operators
	// update sampling step size (tries to keep it in the same ballpark with dr (for, hopefully, more accurate optimization))
	double ddx = relaxation->dx-dxs[0];
	double ddy = relaxation->dy-dys[0];
	double ddr = sqrt(ddx*ddx+ddy*ddy);					// discretization change
	if(ddr < relaxation->d) relaxation->d *= 0.5;		// if change vector < d, halve d
	else relaxation->d *= 2.0;							// otherwise double
	if(relaxation->d < 1.0e-6) relaxation->d *= 2.0;	// can cause numerical issues if d gets too small
}

// relaxes the system for T time steps
void relax(struct Arrays *arrays, struct Model *model, struct Relaxation *relaxation, struct Output *output) {
	for(relaxation->t = 0; relaxation->t <= relaxation->T; relaxation->t++) {
		if(relaxation->T_optimize > 0 && relaxation->t > 0 && relaxation->t%relaxation->T_optimize == 0) optimize(arrays, model, relaxation);	// optimize?
		if(relaxation->t%output->T_print == 0) {	// print output?
			fp(arrays, model, relaxation);
			print(relaxation, output);
		}
		if(relaxation->t%output->T_write == 0) write_state(arrays, relaxation, output);		// write out state?
		if(relaxation->t < relaxation->T) step(arrays, model, relaxation);	// perform iteration step
	}
}

// frees allocated arrays
void clear_arrays(struct Arrays *arrays) {
	fftw_free(arrays->p);
	fftw_free(arrays->q);
	fftw_free(arrays->A);
	fftw_free(arrays->B);
}

int main(int argc, char **argv) {
	// init MPI
	MPI_Init(&argc, &argv);
	fftw_mpi_init();
	// create structs
	struct Arrays arrays;
	struct Model model;
	struct Relaxation relaxation;
	struct Output output;
	relaxation.t0 = time(NULL);
	relaxation.t = 0;
	relaxation.d = 0.0001;
	MPI_Comm_rank(MPI_COMM_WORLD, &relaxation.id);
	MPI_Comm_size(MPI_COMM_WORLD, &relaxation.ID);
	strcpy(output.name, argv[1]);
	// input stream
	char filename[128];
	sprintf(filename, "%s.in", output.name);
	FILE *input = fopen(filename, "r");
	if(input == NULL) {
		printf("Error: Input file not found!\n");
		return 1;
	}
	// create empty output file
	sprintf(filename, "%s.out", output.name);
	FILE *out = fopen(filename, "w");
	fclose(out);
	// read input
	char label;
	char line[1024];
	while(fscanf(input, " %c", &label) != EOF) {
		if(label == '#' && fgets(line, 1024, input)) {}			// comment
		else if(label == '\n') {}								// empty line
		else if(label == 'S') {									// seed random number generators
			seed_rngs(&relaxation, input);
		}
		else if(label == 'O') {									// output
			configure_output(&output, input);
		}
		else if(label == 'A') {									// set up arrays and FFTW plans
			configure_arrays(&arrays, input);
		}
		else if(label == 'I') {									// initialization
			initialize_system(&arrays, input);
		}
		else if(label == 'M') {									// model
			configure_model(&model, input);
		}
		else if(label == 'R') {									// relaxation
			configure_relaxation(&relaxation, input);
			update_AB(&arrays, &model, &relaxation);
			relax(&arrays, &model, &relaxation, &output);
		}
		else {
			printf("Invalid label: %c!\n", label);				// bad input
			return 1;
		}
	}
	fclose(input);
	
	clear_arrays(&arrays);										// clean-up
	MPI_Finalize();

	return 0;

}
