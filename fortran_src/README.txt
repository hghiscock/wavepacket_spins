
		#################################################################################
		##            Programme calcphis.py calculates the singlet yield               ##
		##     of a radical pair reaction as a function of applied field direction     ## 
		#################################################################################

/---------------------------------------------------------------------------------------------------------------------/
/---------------------------------------------------------------------------------------------------------------------/

Input directory:

  #input.dat:
	Contains parameters for the calculation. 
	b0 is the external field strength in milli Tesla 
	g_e is the gyromagnetic ratio of an electron
	k is the recombination rate
	theta range is domain of polar angles to use, with 'grid' number of points
	T and nt are the total time and number of time steps respectively
	radical A/B is file containing radical data in input directory
	n_spins is number of spins to include in that radical
	outputfile is the name of the data file into which the results go (in data directory)
	number of cores is for diagonalisation step

	Example:
	b0	g_e	k
	0.05	1.76d8	1.0d6
	theta range	grid
	0.0	1.0	100
	T	nt	
	1.0d-5	20000
	radical A	n_spins
	h1_hf.dat	1
	radical B	n_spins
	h2_hf.dat	1
	outputfile
	test.dat
	number of cores
	10

  #radical.dat:
	Contains information for the spin active nuclei in the radical, template:
	m1	m2	m3	m4...
	nucleus 1
	a[0,0]	a[1,0]	a[2,0]
	a[0,1]	a[1,1]	a[2,1]
	a[0,2]	a[1,2]	a[2,2]
	nucleus 2
	...
	where mi is the spin multiplicity 2S+1 of nucleus i, and a[i,j] are the components of the hyperfine 
	interaction tensor

  #noise.dat
	Contains information about the time dependent field
	
	Example
	number of frequencies
	1
	frequency	phase	strength (nT)
	13170245.0	0.5	50
	
/---------------------------------------------------------------------------------------------------------------------/
/---------------------------------------------------------------------------------------------------------------------/

Fortran files:

  #procedures.f90
	Contains miscellaneous subroutines

	#sprsin
          Converts a square matrix a(n,n) into a sparse matrix. Only elements with
          magnitude > thresh are retained. Output is two linear arrays with dimension
          nmax, sa contains array values and ija contains indices.

	#sprsax
          Multiply a sparse matrix in storage arrays sa and ija by a square matrix x

	#ran2
        Generates a random number, initial seed idum usually = -1

	#gasdev
          Generate normally distributed random number      

	#ampdist
          Use amplitude distribution as found in experiment
          modify param.f90 file to include distribution

        #kron
	  Calculate the kronecker product of mat1 and mat2, result mat3

	#y_rotate	
          Rotate a matrix A about the y axis

	#t_rotatevec
          Rotate a matrix vector about the y axis

  #param.f90
	Reads in and sets up parameters for calculation
	
        #r_readin(d, pd, t)
          Read in radical data

        #p_readin(d)
          Read in parameters for calculation

        #initialise(d)
          Set up calculation and read in/generate RF parameters
        
	#write_log_file(d, da, db)
          write the log file of the calculation

  #hamiltonians.f90
	Builds Hamiltonians for the calculation

        #pauli_matrices(c)
          build pauli matrices for S/I = 1/2 and 1

        #build(h, c, d)
          build static cartesian and hyperfine hamiltonians

        #noise(h, d)
          build zero time hamiltonian

        #get_time(c, d)
          populate time array        

        #time_step(h, d)
          calculate hamiltonian for next time step

        #integrate_sy(c, d)
          integrate singlet yield for that time step

        #calc_ps(c, d, rab_a, rab_b)
          calculate singlet probability from spin correlation tensors

        #rms_amp(c, d)
          calculate rms amplitude of RF components

  #sc_tensor.f90
	Propagates states over time and calculates spin correlation tensors

        #initialize(sc, h, d)
          initialise state vectors

        #create(sc, h, d)
          allocate memory for various things

        #evolve(sc, h, d)
          evolve state vectors for one time step

        #evolve_t_ind(sc, h, d)
          evolve state vectors for one time step with time independent
          Hamiltonian

        #calc_rab(sc, d)
          calculate spin correlation tensor from state vectors

        #ex_rab(sc, d, h)
          extract spin correlation tensor

        #calc_sc(sc, h, d)
          calculate spin correlation tensor for given time step

  #calcphis.py
	Carries out the calculation
	Initialises the arrays and builds Hamiltonians
	Calculates spin correlation tensors for both radicls
	Calculates singlet yield
	Writes to file 

  #Makefile
	Compiles above files together to give calcphis executable.

	Alter CC = --fcompiler=gfortran to change fortran compiler

/---------------------------------------------------------------------------------------------------------------------/
/---------------------------------------------------------------------------------------------------------------------/

Calculation details:
	
  Calculates singlet yield by numerically propagating each initial state and then calculating spin
  correlation tensors at each time step and combining to give the singlet probability and therefore
  the contribution to the singlet yield

/---------------------------------------------------------------------------------------------------------------------/
/---------------------------------------------------------------------------------------------------------------------/
