        module param
        use procedures
        implicit none

        private :: p_readin
        private :: r_readin
        private :: initialise
        private :: write_log_file

        real(8), parameter, private :: pi = 3.14159265359
        
!-------------------------------------------------------------------------------------------------!

        type radical_dat
        integer :: nspins                               
        integer :: m
        integer, allocatable :: m_array(:)             !2S+1 for nuclei
        real(8), allocatable :: a_tensor(:,:,:)        !hyperfine tensors
        real(8) :: alpha                               !relative orientation of radicals
        contains
        procedure       :: read_in => r_readin
        end type radical_dat

!-------------------------------------------------------------------------------------------------!

        type param_dat
        integer :: ncores                              !number of cores
        real(8) :: a(4), b(4)                          !integrator coefficient arrays
        real(8) :: omega_0                             !external field strength
        real(8) :: gamma_e                             

        real(8) :: w_min, w_max                        !bounding rf frequencies
        integer :: nfourier                            !number of frequency components
        real(8), allocatable :: w_rf(:)                !frequencies
        real(8), allocatable :: phase_rf(:)            !phases
        real(8), allocatable :: coeff_rf(:)            !amplitudes
        real(8), allocatable :: theta_rf(:)            !direction
        real(8), allocatable :: phi_rf(:)

        integer :: t_ind                               !if time-independent hamiltonian
        logical :: iden                                !if identical radical
        logical :: rms                                 !if want to calculate rms amplitude

        integer :: nt                                  !time steps
        real(8) :: time                                !overall time
        real(8) :: tau                                 !step size
        real(8) :: k                                   !rate constant
        real(8), allocatable :: t_array(:)             !time array

        real(8) :: theta_min, theta_max                !zeeman polar angle bounds
        real(8) :: phi_min, phi_max                    !zeeman azimuthal angle bounds
        real(8) :: d_theta, d_phi                      !angle step sizes
        integer :: ntheta, nphi                        !number of steps
        real(8) :: phi
        real(8), allocatable :: theta(:)

        character (len = 10) :: rada, radb             !radical input files
        integer :: na, nb                              !number of spins

        character (len = 200) :: outputfile, logfile   !data files to write to
        integer :: start, finish                       !timing the code
        contains 
        procedure       :: read_in => p_readin
        procedure       :: initialise
        procedure       :: write_log_file
        end type param_dat

!-------------------------------------------------------------------------------------------------!

        contains

!-------------------------------------------------------------------------------------------------!

        subroutine r_readin(d, pd, t)

        !Read in radical data

        class(radical_dat) :: d
        class(param_dat) :: pd
        integer, intent(in) :: t
        real(8) :: a1
        integer :: i
        character(len=80) :: testc

        select case(t)
        case(1)
                open(1, file = 'input/' // trim(pd%rada))
                d%nspins = pd%na
        case(2)
                open(1, file = 'input/'// trim(pd%radb))
                d%nspins = pd%nb
        end select

        allocate( d%m_array(d%nspins) )
        read(1,'(a)') testc
        do i = 1, d%nspins
          read(testc(2*i-1:2*i),*) d%m_array(i)
        end do

        d%m = 2*product(d%m_array)
        allocate( d%a_tensor(d%nspins+1,3,3) )
        d%a_tensor = 0.00d0
        do i = 1, d%nspins
                read(1,*)
                read(1,*) d%a_tensor(i,1,1), d%a_tensor(i,1,2), &
                d%a_tensor(i,1,3)
                read(1,*) d%a_tensor(i,2,1), d%a_tensor(i,2,2), &
                d%a_tensor(i,2,3)
                read(1,*) d%a_tensor(i,3,1), d%a_tensor(i,3,2), &
                d%a_tensor(i,3,3)
        end do

        close(1)

        end subroutine r_readin

!-------------------------------------------------------------------------------------------------!

        subroutine p_readin(d)
        
        !Read in parameters for calculation

        class(param_dat) :: d
        real(8) :: rtthree, b0
        real(8) :: tmin, tmax
        real(8) :: pmin, pmax

        call system_clock(d%start)
        
        open(1, file = 'input/input.dat')

        read(1,*)
        read(1,*) b0, d%gamma_e, d%k
        read(1,*)
        read(1,*) tmin, tmax, d%ntheta
        read(1,*)
        read(1,*) d%time, d%nt
        read(1,*)
        read(1,*) d%rada, d%na
        read(1,*)
        read(1,*) d%radb, d%nb
        read(1,*)
        read(1,*) d%outputfile
        read(1,*)
        read(1,*) d%ncores

        d%theta_min = tmin*pi
        d%theta_max = tmax*pi
        d%d_theta = (d%theta_max - d%theta_min)/dble(d%ntheta)

        d%omega_0 = b0*d%gamma_e
        d%tau = d%time/dble(d%nt)
        rtthree = 1.00d0/DSQRT(3.00d0)
        d%a = (/ 0.50d0*(1.00d0 - rtthree)*d%tau, rtthree*d%tau, &
        -rtthree*d%tau, 0.50d0*(1.00d0 + rtthree)*d%tau /)
        d%b = (/ 0.00d0, 0.50d0*(0.50d0 + rtthree)*d%tau, 0.50d0*d%tau,&
         0.50d0*(0.50d0 - rtthree)*d%tau /)

        close(1)

        end subroutine p_readin
 
!-------------------------------------------------------------------------------------------------!

        subroutine initialise(d)

        !Set up calculation and read in/generate RF parameters

        class(param_dat) :: d

        integer :: idum, i, j
        real(8) :: dw, t_tmp

        idum = -1
        open(57, file = 'input/noise.dat')
        read(57,*)
        read(57,*) d%nfourier
        read(57,*)

        if (d%nfourier == 0) then
                d%t_ind = 1
        else
                d%t_ind = 0
        end if

        allocate( d%w_rf(d%nfourier) )
        allocate( d%phase_rf(d%nfourier) )
        allocate( d%coeff_rf(d%nfourier) )
        allocate( d%theta_rf(d%nfourier) )
        allocate( d%phi_rf(d%nfourier) )

        do i = 1, d%nfourier
            read(57,*) d%w_rf(i), d%phase_rf(i), d%coeff_rf(i)
            d%coeff_rf(i) = d%coeff_rf(i)*d%gamma_e*1.0d-6
            d%phi_rf(i) = 0.00d0
            d%theta_rf(i) = 0.00d0
            d%phase_rf(i) = d%phase_rf(i)*pi
        end do

        allocate( d%t_array(3*d%nt) )        

        t_tmp = 0.00d0
        do i = 1, d%nt
                t_tmp = t_tmp + d%a(1)
                do j = 2, 4
                        d%t_array(3*i-4+j) = t_tmp
                        t_tmp = t_tmp + d%a(j)
                end do
        end do

        allocate( d%theta(d%ntheta) )
        
        do i = 1, d%ntheta
                d%theta(i) = d%theta_min + dble(i-1)*d%d_theta
        end do

        d%phi = 0.00d0*pi

        close(57)

        end subroutine initialise

!-------------------------------------------------------------------------------------------------!

        subroutine write_log_file(d, da, db)

        !write the log file of the calculation

        class(param_dat) :: d
        class(radical_dat) :: da, db

        open(10, file = 'log_files/nowrunning.dat')
        write(10, *) 'theta values computed', d%ntheta
        write(10, *) 'domain', d%theta_min, 'to', d%theta_max
        write(10, *) 'noise frequencies', d%nfourier
        write(10, *) 'frequency range (MHz)', d%w_min/(2.00d0*pi*1d6), &
        'to', d%w_max/(2.00d0*pi*1d6)
        write(10, *) 'timesteps', d%nt
        write(10, *) 'time ', d%time
        write(10, *) 'rate constant', d%k
        write(10, *) 'output ', 'file ', trim(adjustl(d%outputfile))
        write(10, *) 'Electron A coupled to ', da%nspins, 'nuclei ', &
        '2S+1 ', da%m_array(1:da%nspins)
        write(10, *) 'Electron B coupled to ', db%nspins, 'nuclei ', &
        '2S+1 ', db%m_array(1:db%nspins)

        call system_clock(d%finish)
        write(10, *) 'CPU time', (d%finish-d%start)

        end subroutine write_log_file

!-------------------------------------------------------------------------------------------------!

        end module param
        

