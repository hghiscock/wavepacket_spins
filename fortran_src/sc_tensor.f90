        module sc_tensor
        use hamiltonians
        use param
        implicit none

        private :: evolve
        private :: evolve_t_ind
        private :: calc_rab
        private :: ex_rab
        private :: initialize
        private :: create
        private :: calc_sc

        real(8), parameter, private :: pi = 3.14159265359

!-------------------------------------------------------------------------------------------------!

        type spin_corr
        integer :: m, mprime                            !spin multiplicity
        integer :: i                                    !time step
        real(8), allocatable :: r_ab(:,:,:)             !spin correlation tensor
        complex(8), allocatable :: p(:,:,:), q(:,:,:)   !state vector and time derivative
        complex(8), allocatable :: h_ze(:,:,:)          !zeeman H
        complex(8), allocatable :: rtmp(:,:,:)  
        complex :: xj

        contains
        procedure :: evolve
        procedure :: evolve_t_ind
        procedure :: calc_rab
        procedure :: ex_rab
        procedure :: initialize
        procedure :: calc_sc
        procedure :: create

        end type spin_corr

!-------------------------------------------------------------------------------------------------!

        contains

!-------------------------------------------------------------------------------------------------!

        subroutine initialize(sc, h, d)

        !initialise state vectors

        class(spin_corr) :: sc
        class(hamiltonian) :: h
        class(param_dat) :: d

        complex(8), allocatable :: htmp(:,:)
        real(8) :: rtthree
        integer :: i, j

        allocate( htmp(h%m,h%m) )

        sc%p = (0.00d0, 0.00d0)
        sc%q = (0.00d0, 0.00d0)
        sc%xj = (0.00d0, 1.00d0)
        sc%rtmp = 0.00d0
        sc%r_ab = 0.00d0
        sc%mprime = h%m/2

        do i = 1, d%ntheta
                sc%h_ze(i,:,:) = d%omega_0*(h%hx*sin(d%theta(i)) &
                + h%hz*cos(d%theta(i))) + h%h_hf
                htmp = sc%h_ze(i,:,:) + h%h_td0

                do j = 1, h%m
                        sc%q(i,j,j) = (1.00d0, 0.00d0)
                end do
                sc%p(i,:,:) = -1.00d0*sc%xj*matmul(htmp,sc%q(i,:,:))
        end do

        deallocate( htmp )

        end subroutine initialize

!-------------------------------------------------------------------------------------------------!

        subroutine create(sc, h, d)

        !allocate memory for various things

        class(spin_corr) :: sc
        class(hamiltonian) :: h
        class(param_dat) :: d
        
        allocate( sc%p(d%ntheta,h%m,h%m), sc%q(d%ntheta,h%m,h%m) )
        allocate( sc%h_ze(d%ntheta,h%m,h%m) )
        allocate( sc%r_ab(d%ntheta,3,3), sc%rtmp(d%ntheta,3,3) )

        h%i = 1

        call sc%initialize(h, d)

        end subroutine create

!-------------------------------------------------------------------------------------------------!

        subroutine evolve(sc, h, d)

        !evolve state vectors for one time step

        class(spin_corr) :: sc
        class(hamiltonian) :: h
        class(param_dat) :: d

        complex(8) :: htmp(h%m,h%m)
        complex(8) :: w(h%m,h%m)
        integer :: i, j

        sc%q = sc%q + d%a(1)*sc%p
        do j = 2, 4
                call h%time_step(d)
!$omp parallel do default(shared) private(i,htmp,w)
                do i = 1, d%ntheta
                        htmp = sc%h_ze(i,:,:) + h%h_td(:,:)
                        w = matmul(htmp,htmp) + sc%xj*h%dhdt(:,:)
                        sc%p(i,:,:) = sc%p(i,:,:) - d%b(j)*matmul(w,sc%q(i,:,:))
                        sc%q(i,:,:) = sc%q(i,:,:) + d%a(j)*sc%p(i,:,:)
                end do
!$omp end parallel do
        end do

        end subroutine evolve

!-------------------------------------------------------------------------------------------------!

        subroutine evolve_t_ind(sc, h, d)

        !evolve state vectors for one time step with time independent
        !Hamiltonian

        class(spin_corr) :: sc
        class(hamiltonian) :: h
        class(param_dat) :: d

        complex(8) :: w(h%m,h%m)
        integer :: i, j

        sc%q = sc%q + d%a(1)*sc%p
        do j = 2,4
!$omp parallel do default(shared) private(i,w)
                do i = 1, d%ntheta
                        w = matmul(sc%h_ze(i,:,:),sc%h_ze(i,:,:))
                        sc%p(i,:,:) = sc%p(i,:,:) - d%b(j)*matmul(w,sc%q(i,:,:))
                        sc%q(i,:,:) = sc%q(i,:,:) + d%a(j)*sc%p(i,:,:)
                end do
!$omp end parallel do
        end do

        end subroutine evolve_t_ind

!-------------------------------------------------------------------------------------------------!

        subroutine calc_rab(sc, d)

        !calculate spin correlation tensor from state vectors

        class(spin_corr) :: sc
        class(param_dat) :: d
        integer :: i, j, jprime, k, kprime

        complex(8), allocatable :: rtmp(:,:,:)

        allocate( rtmp(d%ntheta,3,3) )

        rtmp = (0.00d0, 0.00d0)
        sc%rtmp = (0.00d0, 0.00d0)
!$omp parallel do default(shared) private(i, j,jprime,k,kprime) reduction(+:rtmp)
        do j = 1, sc%mprime
                jprime = j+sc%mprime
                do k = 1, sc%mprime
                        kprime = k+sc%mprime
                        do i = 1, d%ntheta
                                !Rxx/Rxy
                                rtmp(i,1,1) = rtmp(i,1,1) + &
                                conjg(sc%q(i,j,kprime))*sc%q(i,jprime,k)
                                rtmp(i,1,1) = rtmp(i,1,1) + &
                                conjg(sc%q(i,j,k))*sc%q(i,jprime,kprime)
                                !Ryy/Ryx
                                rtmp(i,2,2) = rtmp(i,2,2) - &
                                conjg(sc%q(i,jprime,k))*sc%q(i,j,kprime)
                                rtmp(i,2,2) = rtmp(i,2,2) + &
                                conjg(sc%q(i,jprime,kprime))*sc%q(i,j,k)
                                !Rzz
                                rtmp(i,3,3) = rtmp(i,3,3) + &
                                conjg(sc%q(i,j,k))*sc%q(i,j,k)
                                rtmp(i,3,3) = rtmp(i,3,3) + &
                                conjg(sc%q(i,jprime,kprime))*sc%q(i,jprime,kprime)
                                rtmp(i,3,3) = rtmp(i,3,3) - &
                                conjg(sc%q(i,j,kprime))*sc%q(i,j,kprime)
                                rtmp(i,3,3) = rtmp(i,3,3) - &
                                conjg(sc%q(i,jprime,k))*sc%q(i,jprime,k)
                                !Rxz/Ryz
                                rtmp(i,1,3) = rtmp(i,1,3) + &
                                conjg(sc%q(i,j,kprime))*sc%q(i,j,k)
                                rtmp(i,1,3) = rtmp(i,1,3) - &
                                conjg(sc%q(i,jprime,kprime))*sc%q(i,jprime,k)
                                !Rzx/Rzy
                                rtmp(i,3,1) = rtmp(i,3,1) + &
                                conjg(sc%q(i,j,k))*sc%q(i,jprime,k)
                                rtmp(i,3,1) = rtmp(i,3,1) - &
                                conjg(sc%q(i,j,kprime))*sc%q(i,jprime,kprime)
                        end do
                end do
        end do
!$end parallel do

        sc%rtmp = rtmp

        end subroutine calc_rab


!-------------------------------------------------------------------------------------------------!

        subroutine ex_rab(sc, d, h)

        !extract spin correlation tensor

        class(spin_corr) :: sc 
        class(param_dat) :: d     
        class(hamiltonian) :: h  

        integer :: i

!$omp parallel do default(shared) private(i)
        do i = 1, d%ntheta
                sc%r_ab(i,1,1) = real(sc%rtmp(i,1,1))
                sc%r_ab(i,2,2) = real(sc%rtmp(i,2,2))
                sc%r_ab(i,3,3) = reaL(sc%rtmp(i,3,3))/2.00d0
                sc%r_ab(i,1,2) = aimag(sc%rtmp(i,1,1))
                sc%r_ab(i,2,1) = aimag(sc%rtmp(i,2,2))
                sc%r_ab(i,1,3) = real(sc%rtmp(i,1,3))
                sc%r_ab(i,3,1) = real(sc%rtmp(i,3,1))
                sc%r_ab(i,2,3) = aimag(sc%rtmp(i,1,3))
                sc%r_ab(i,3,2) = aimag(sc%rtmp(i,3,1))
        end do
!$omp end parallel do

        sc%r_ab = sc%r_ab/dble(h%m)


        end subroutine ex_rab
        
!-------------------------------------------------------------------------------------------------!

        subroutine calc_sc(sc, h, d)

        !calculate spin correlation tensor for given time step

        class(spin_corr) :: sc
        class(hamiltonian) :: h
        class(param_dat) :: d

        select case(d%t_ind)
        case(0)
            call sc%evolve(h, d)
            call sc%calc_rab(d)     
            call sc%ex_rab(d, h)
        case(1)
            call sc%evolve_t_ind(h, d)
            call sc%calc_rab(d)
            call sc%ex_rab(d, h)
        end select

        end subroutine calc_sc
        
!-------------------------------------------------------------------------------------------------!

        end module sc_tensor
