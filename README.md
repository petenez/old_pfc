# pfc
High-performance phase field crystal (PFC) code for generating realistic model systems of polycrystalline graphene.

Background

Phase field crystal (PFC) models are a family of continuum methods for modeling the microstructure and its evolution in polycrystalline materials. The basics of phase field crystal (PFC) models are covered in the book by Provatas and Elder [1], for example. PFC comes with attractive multiscale characteristics - it gives simultaneous access to atomistic and mesoscopic length scales and to long time scales. This makes it ideal for modeling realistic microstructures and their evolution. Conventional atomistic techniques on the other hand have severe length and/or time scale limitations in this respect. We have exploited the standard one-mode PFC model for generating model systems of polycrystalline graphene for further atomistic calculations employing molecular dynamics (MD) and quantum-mechanical density functional theory [2]. Previously, Zhang et al. had used PFC similarly for generating model systems of defect-engineered graphene to initialize atomistic calculations [3].

Purpose of code

In short, this code is intended for generating highly-relaxed and realistic model systems of polycrystalline graphene. I give no promises of developing this code actively. The model systems are initialized with a crude initial guess and are then relaxed to equilibrium using PFC. The resulting PFC density fields can be mapped into atomic coordinates and studied further using atomistic methods. I have some codes for the mapping step, but they are messy and complicated, and not yet ready for publication. I describe in the last section how to approach this step.

Model

The code implements the standard one-mode model [4] with one controlled length scale. Its free energy is given by

F = int dr (psi/2*(alpha + beta*(1 + nabla^2)^2)*psi + gamma*psi^3/3 + delta*psi^4/4).

Here, psi is the density field that describes the system. Of the parameters, alpha is proportional to temperature and nabla^2 is the Laplacian. Besides the model parameters, average density also has a significant influence on the behavior of the model. Both conserved

dpsi/dt = nabla^2 dF/dpsi

and nonconserved dynamics

dpsi/dt = - dF/dpsi

can be used. Former conserves the average density of the initial state and is diffusive, whereas latter does not conserve the average density and follows the steepest descent path in energy. The semi-implicit spectral method [1] is used to solve these differential equations numerically. Periodic boundary conditions ensue from spectral treatment. For simplicity, the code is limited to planar and freestanding systems.

Implementation

The code is written in C and exploits MPI [5] and FFTW [6] for efficient and parallelizable computation. It requires little memory and scales well over multiple cores and nodes. Relatively large systems can be modeled even on a laptop. The code is simple, but not very user-friendly in the sense that few safeguards have been implemented against misuse. It does, however, give some error messages for invalid input. 

In order to run, the code requires an input file where the details of the calculation are specified. These include: a seed for the random number generator, frequency of output, dimensions of the 2D calculation grid, initialization type, model parameters and relaxation settings. The different options are explained in the comments in the sample input file.

The code produces two kinds of output: numerical data of elapsed time, spatial discretization, average free energy density and average density, as well as data files of the current state of the density field. The latter simply list the density values of the calculation grid row-by-row.

The calculations can be initialized (1) with random noise, (2) with randomly distributed and oriented crystallites or (3) with a state from a previous calculation loaded from a file. For the random initialization, the noise is sampled from a uniform distribution with user-defined mean and amplitude. In case of crystallites, their number and radii can be controlled. The average density as well as the amplitude and the length scale of the density oscillations in the initial state are also user-defined. The length scale is the lattice constant and is roughly ~4pi/sqrt(3) in dimensionless units. When reading the initial state from a data file, the user can specify an average density and an amplitude for density oscillations according to which the system will be normalized.

The user can control all of the model parameters and also choose between conserved and nonconserved dynamics. Recall that the average density influences strongly the behavior of the model. Conserved dynamics conserves the average density of the initial state, whereas nonconserved dynamics lets it seek an equilibrium value given by the model parameters.

In the relaxation settings, the user defines the number of iterations for the relaxation, spatial discretization and time step, and the interval between calculation box size optimizations. The latter tries to eliminate the mismatch between the system and the box. In practice, it varies the box size slightly, records the resulting changes in energy and uses quadratic interpolation to guess an optimal box size.

The code can be compiled with the command

mpicc pfc.c -lfftw3_mpi -lfftw3 -lm -Ofast -Wall -o pfc

and run by

mpirun -np 8 pfc case

Here -np 8 indicates that eight CPUs will be used for the computation. The text string "case" is the name for the study case - the input file must be "case.in", numerical output appears in "case.out" and data files begin with "case".

A plotter tool written in Java is also provided for visualization of the systems modeled. I typically compile it into a JAR file (also provided). For a 512-by-512 case "dummy" it can be run as

java -jar plotter.jar dummy-t:0[.dat] dummy-t:0[.png] 512 512

where the output image file is specified. Note that the file extensions can be omitted - the tool appends the text strings with ".dat" and ".png" by default. It can also be used in batch mode

java -jar plotter.jar dummy-t:# dummy-t:# 512 512 0 1000 10000

where the hashtags will be substituted with numbers going through 0, 1000, 2000, ..., 10000. The plotter can also be used to visualize complex-valued data (which the PFC code does not produce, however).

Practical considerations

To generate model systems of polycrystalline graphene for further atomistic calculations, the code gives three options: start with (1) random noise, (2) crystallites or (3) a polycrystalline state. Starting with a random state and waiting for coarsening of the emerging grain structure is feasible, but can take an unpractical amount of time. Starting with crystallites (crystallite radius << their separation) that grow in a constant density phase is much faster. The fastest is to begin with a polycrystalline state (crystalline radius >~ their separation), but here the grain structure is a bit more constrained compared to the crystallite growth alternative.

The choice of model parameters and average density strongly influences the growth and coarsening of the system: For a system that is effectively close to melting, coarsening is faster, but the defects are fuzzier and harder to resolve for atomic coordinates. On the other hand, under some conditions, a more rapidly crystallizing metastable stripe phase can grow out of the honeycomb crystallites and may take a significant amount of time to transform into honeycomb structure. In such cases, it is usually most straightforward to first relax the model systems closer to melting and then use this relaxed state as the starting point of a second final relaxation further from melting. The calculation box size optimization can be employed for the latter stage, for example, to relieve residual stress.

The dynamics also play a role: conserved dynamics is diffusive, whereas nonconserved dynamics follows the steepest descent path in energy. Latter is somewhat similar to quenching as it is more likely to eventually get stuck in a metastable state. However, temperature is built in in the model so, strictly speaking, this is not quenching. Furthermore, PFC displays Peierls barriers whereby conserved dynamics can also, in principle, get stuck, but in practice coarsening rarely stops.

The two sample input files demonstrate a quick two-stage process where crystallites are first grown close to melting and the resulting polycrystalline system is then relaxed further from melting. The first step exploits conserved and the second nonconserved dynamics. Calculation box size optimization is also employed for the second stage. The shell script provided compiles the PFC code, runs the two relaxation steps and plots their results.

Lastly, the relaxed density field needs to be mapped into atomic coordinates. One can do this naively by associating all local maxima with atom positions. However, this results in a large number of missing atoms along grain boundaries. A much more robust approach is to locate all triplets of local minima whose members are closest neighbors to each other and to place an atom at the center.

References

[1] Provatas and Elder, Phase-Field Methods in Materials Science and Engineering (Wiley-VCH, 2010).

[2] Hirvonen et al., Physical Review B 94, 035414 (2016); Fan et al., Nano Letters 17, 5919-5924 (2017); Azizi et al., Carbon 125, 384-390 (2017).

[3] Zhang et al., Extreme Mechanics Letters 1, 3-8 (2014).

[4] Elder et al., Physical Review Letters 88, 245701 (2002).

[5] MPI Forum, website. Available (visited December 7, 2017): http://mpi-forum.org/ .

[6] FFTW, website. Available (visited December 7, 2017): http://www.fftw.org/ .
