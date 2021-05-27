        module hamiltonians
        use param
        use procedures
        implicit none

        private :: pauli_matrices
        private :: build
        private :: noise
        private :: get_time
        private :: time_step
        private :: rms_amp
        private :: integrate_sy
        private :: calc_ps

        real(8), parameter, private :: pi = 3.14159265359

!-------------------------------------------------------------------------------------------------!
        
        type hamiltonian
        integer :: m                                 !spin multiplicity
        integer :: nspins                            !number of spins
        integer :: i                                 !time step
        complex(8), allocatable :: hx(:,:)           !cartesian hamiltonians
        complex(8), allocatable :: hy(:,:)
        complex(8), allocatable :: hz(:,:)    
        complex(8), allocatable :: h_hf(:,:)         !hyperfine hamiltonian
        complex(8), allocatable :: h_td0(:,:)        !zero time hamiltonian        
        complex(8), allocatable :: h_td(:,:)         !time-dependent H
        complex(8), allocatable :: dhdt(:,:)         !time derivative of H
        complex(8), allocatable :: hfourier(:,:,:)   !RF hamiltonians
        
        contains
        procedure       :: build
        procedure       :: noise
        procedure       :: time_step
        end type hamiltonian

!-------------------------------------------------------------------------------------------------!

        type calc_dat
        integer :: j
        complex(8), allocatable :: s2(:,:,:)         !Pauli matrices
        complex(8), allocatable :: s3(:,:,:)
        real(8), allocatable :: p_s(:)               !Singlet projection
        real(8), allocatable :: phis(:)              !singlet yield
        real(8), allocatable :: t_array(:)           !time array

        real(8) :: rms_a

        contains
        procedure       :: pauli_matrices
        procedure       :: get_time
        procedure       :: rms_amp
        procedure       :: integrate_sy
        procedure       :: calc_ps
        end type calc_dat

!-------------------------------------------------------------------------------------------------!

        contains

!-------------------------------------------------------------------------------------------------!

        subroutine pauli_matrices(c)

        !build pauli matrices for S/I = 1/2 and 1

        class(calc_dat) :: c
        
        allocate( c%s2(2,2,4), c%s3(3,3,4) )

        c%s2(:,1,1) = (/ (0.0, 0.0), (0.5, 0.0) /)                      !Sx   
        c%s2(:,2,1) = (/ (0.5, 0.0), (0.0, 0.0) /)                      
                                                                      
        c%s2(:,1,2) = (/ (0.0, 0.0), (0.0, -0.5) /)                     !Sy   
        c%s2(:,2,2) = (/ (0.0, 0.5), (0.0, 0.0) /)                      
                                                                        
        c%s2(:,1,3) = (/ (0.5, 0.0), (0.0, 0.0) /)                      !Sz   
        c%s2(:,2,3) = (/ (0.0, 0.0), (-0.5, 0.0) /)                     
                                                                       
        c%s2(:,1,4) = (/ (1.0, 0.0), (0.0, 0.0) /)                      !I2
        c%s2(:,2,4) = (/ (0.0, 0.0), (1.0, 0.0) /)                      

        c%s3(:,1,1) = (/ (0.0, 0.0), (1.0, 0.0), (0.0,0.0) /)           !Sx   
        c%s3(:,2,1) = (/ (1.0, 0.0), (0.0, 0.0), (1.0,0.0) /)          
        c%s3(:,3,1) = (/ (0.0, 0.0), (1.0, 0.0), (0.0,0.0) /)          
                                                                       
        c%s3(:,1,2) = (/ (0.0, 0.0), (0.0, -1.0), (0.0,0.0) /)          !Sy  
        c%s3(:,2,2) = (/ (0.0, 1.0), (0.0, 0.0), (0.0,-1.0) /)        
        c%s3(:,3,2) = (/ (0.0, 0.0), (0.0, 1.0), (0.0,0.0) /)         
                                                                      
        c%s3(:,1,3) = (/ (1.0, 0.0), (0.0, 0.0), (0.0,0.0) /)           !Sz   
        c%s3(:,2,3) = (/ (0.0, 0.0), (0.0, 0.0), (0.0,0.0) /)           
        c%s3(:,3,3) = (/ (0.0, 0.0), (0.0, 0.0), (-1.0,0.0) /)          
                                                                       
        c%s3(:,1,4) = (/ (1.0, 0.0), (0.0, 0.0), (0.0,0.0) /)           !I3  
        c%s3(:,2,4) = (/ (0.0, 0.0), (1.0, 0.0), (0.0,0.0) /)           
        c%s3(:,3,4) = (/ (0.0, 0.0), (0.0, 0.0), (1.0,0.0) /)           
                                                                        
        c%s3(:,:,1:2) = c%s3(:,:,1:2)/sqrt(2.00d0)

        end subroutine pauli_matrices

!-------------------------------------------------------------------------------------------------!

        subroutine build(h, c, d)

        !build static cartesian and hyperfine hamiltonians

        class(hamiltonian) :: h
        class(calc_dat) :: c
        class(radical_dat) :: d

        complex(8), allocatable :: hxtmp(:,:), hytmp(:,:), hztmp(:,:)
        complex(8), allocatable :: htmp(:,:,:), htmp1(:,:,:), htmp2(:,:)
        complex(8), allocatable :: hf1(:,:,:)
        real(8) :: a
        integer :: dim_tmp, i, j, k, l

        a = 1.76d8
        h%m = d%m
        h%nspins = d%nspins

        allocate( h%hx(h%m,h%m), h%hy(h%m,h%m), h%hz(h%m,h%m) )
        allocate( h%h_hf(h%m,h%m) )
        allocate( hxtmp(h%m,h%m), hytmp(h%m,h%m), hztmp(h%m,h%m) )
        allocate( htmp(h%m,h%m,3), htmp1(h%m,h%m,3), htmp2(h%m,h%m) )
        allocate( hf1(h%m,h%m,h%nspins) )
        

        !Build zeeman hamiltonians
        h%hx(1:2,1:2) = c%s2(:,:,1)
        h%hy(1:2,1:2) = c%s2(:,:,2)
        h%hz(1:2,1:2) = c%s2(:,:,3)
        dim_tmp = 2
        do i = 1, h%nspins
                hxtmp = h%hx
                hytmp = h%hy
                hztmp = h%hz
                if (d%m_array(i) == 2) then
                        call kron(hxtmp, c%s2(:,:,4), h%hx, dim_tmp, 2)
                        call kron(hytmp, c%s2(:,:,4), h%hy, dim_tmp, 2)
                        call kron(hztmp, c%s2(:,:,4), h%hz, dim_tmp, 2)
                        dim_tmp = dim_tmp*2
                else if (d%m_array(i) == 3) then
                        call kron(hxtmp, c%s3(:,:,4), h%hx, dim_tmp, 3)
                        call kron(hytmp, c%s3(:,:,4), h%hy, dim_tmp, 3)
                        call kron(hztmp, c%s3(:,:,4), h%hz, dim_tmp, 3)
                        dim_tmp = dim_tmp*3
                end if
        end do

        !Build hyperfine hamiltonian
        h%h_hf = 0.00d0
        hf1 = 0.00d0
        do k = 1, h%nspins
                dim_tmp = 2
                do i = 1, 3
                        htmp(1:2,1:2,i) = c%s2(:,:,i)
                end do
                do l = 1, k-1
                        if (d%m_array(l) == 2) then
                                htmp1 = htmp
                                do i = 1, 3
                                        call kron(htmp1(:,:,i), &
                                        c%s2(:,:,4), htmp(:,:,i), &
                                        dim_tmp, 2)
                                end do
                                dim_tmp = dim_tmp*2
                         else if (d%m_array(l) == 3) then
                                htmp1 = htmp
                                do i = 1,3
                                        call kron(htmp1(:,:,i), &
                                        c%s3(:,:,4), htmp(:,:,i), &
                                        dim_tmp, 3)
                                end do
                                dim_tmp = dim_tmp*3
                        end if
                end do
                if (d%m_array(k) == 2) then
                        do i = 1, 3
                                do j = 1, 3
                                        htmp2 = 0.00d0
                                        call kron(htmp(:,:,i), &
                                        c%s2(:,:,j), htmp2, dim_tmp, 2)
                                        hf1(:,:,k) = hf1(:,:,k) + &
                                        d%a_tensor(k,i,j)*htmp2*a
                                end do
                        end do
                        dim_tmp = dim_tmp*2
                else if (d%m_array(k) == 3) then
                        do i = 1, 3
                                do j = 1, 3
                                        htmp2 = 0.00d0
                                        call kron(htmp(:,:,i), &
                                        c%s3(:,:,j), htmp2, dim_tmp, 3)
                                        hf1(:,:,k) = hf1(:,:,k) + &
                                        d%a_tensor(k,i,j)*htmp2*a
                                end do
                        end do
                        dim_tmp = dim_tmp*3
                end if
                do l = k+1, h%nspins
                        if (d%m_array(l) == 2) then
                                htmp2 = hf1(:,:,k)
                                call kron(htmp2, c%s2(:,:,4), &
                                hf1(:,:,k), dim_tmp, 2)
                                dim_tmp = dim_tmp*2
                        else if (d%m_array(l) == 3) then
                                htmp2 = hf1(:,:,k)
                                call kron(htmp2, c%s3(:,:,4), &
                                hf1(:,:,k), dim_tmp, 3)
                                dim_tmp = dim_tmp*3
                        end if
                end do
        end do

        do k = 1, h%nspins
                h%h_hf = h%h_hf + hf1(:,:,k)
        end do

        end subroutine build

!-------------------------------------------------------------------------------------------------!

        subroutine noise(h, d)

        !build zero time hamiltonian

        class(hamiltonian) :: h
        class(param_dat) :: d

        integer :: i

        allocate( h%h_td(h%m, h%m), h%dhdt(h%m, h%m) )
        allocate( h%h_td0(h%m, h%m) )
        allocate( h%hfourier(d%nfourier, h%m, h%m) )

        h%h_td0 = 0.00d0

        do i = 1, d%nfourier
                h%hfourier(i,:,:) = d%coeff_rf(i)*((h%hx*cos(d%phi_rf(i)) &
                + h%hy*sin(d%phi_rf(i)))*sin(d%theta_rf(i)) + h%hz*cos(d%theta_rf(i)))
                h%h_td0 = h%h_td0 + h%hfourier(i,:,:)*sin(d%phase_rf(i))
        end do

        end subroutine noise

!-------------------------------------------------------------------------------------------------!

        subroutine get_time(c, d)

        !populate time array        

        class(calc_dat) :: c
        class(param_dat) :: d

        integer :: i
        real(8) :: third

        allocate( c%phis(d%ntheta) )
        allocate( c%p_s(d%ntheta) )
        allocate( c%t_array(d%nt) )

        c%t_array(1) = 0.00d0
        do i = 2, d%nt
                c%t_array(i) = c%t_array(i-1) + d%tau
        end do

        third = 1.00d0/3.00d0
        c%t_array(1) = third*exp(-1.00d0*d%k*c%t_array(1))*d%tau*d%k
        c%t_array(d%nt) = third*exp(-1.00d0*d%k*c%t_array(d%nt))*d%tau*d%k
        do i = 2, d%nt-1
                c%t_array(i) = (mod(i-1,2)+1.00d0)*2.00d0*third*&
                      exp(-1.00d0*d%k*c%t_array(i))*d%tau*d%k
        end do
        c%j = 0

        end subroutine get_time

!-------------------------------------------------------------------------------------------------!

        subroutine time_step(h, d)

        !calculate hamiltonian for next time step

        class(hamiltonian) :: h
        class(param_dat) :: d

        complex(8), dimension(h%m,h%m) :: hi, dhdti
        integer :: i, j, k
        real(8) :: arg

        hi = (0.00d0, 0.00d0)
        dhdti = (0.00d0, 0.00d0)
!$omp parallel do default(shared) private(i,arg) reduction(+:hi,dhdti)
        do i = 1, d%nfourier
                arg = d%w_rf(i)*d%t_array(h%i) + d%phase_rf(i)
                hi = hi + &
                sin(arg)*h%hfourier(i,:,:)
                dhdti = dhdti + &
                d%w_rf(i)*cos(arg)*h%hfourier(i,:,:)
        end do
!$omp end parallel do
        h%i = h%i + 1
        h%h_td = hi
        h%dhdt = dhdti

        end subroutine time_step

!-------------------------------------------------------------------------------------------------!

        subroutine integrate_sy(c, d)

        !integrate singlet yield for that time step

        class(calc_dat) :: c
        class(param_dat) :: d

        integer :: i, j
        real(8) :: third, arg

        do j = 1, d%ntheta
                c%phis(j) = c%phis(j) + c%p_s(j)*c%t_array(c%j)
        end do

        end subroutine integrate_sy
        
!-------------------------------------------------------------------------------------------------!

        subroutine calc_ps(c, d, rab_a, rab_b)

        !calculate singlet probability from spin correlation tensors

        class(calc_dat) :: c
        class(param_dat) :: d
        real(8), intent(in) :: rab_a(:,:,:), rab_b(:,:,:)

        integer :: i, j, k, l

        c%p_s = 0.25d0
        do j = 1, 3
                do k = 1, 3
                        do l = 1, d%ntheta
                                c%p_s(l) = c%p_s(l) + rab_a(l,j,k)* &
                                rab_b(l,j,k)
                        end do
                end do
        end do

        c%j = c%j+1
        
        end subroutine calc_ps

!-------------------------------------------------------------------------------------------------!

        subroutine rms_amp(c, d)

        !calculate rms amplitude of RF components

        class(calc_dat) :: c
        class(param_dat) :: d

        real(8) :: noisefield(3), gamma_e, rms_tmp
        real(8), allocatable :: n_vec(:,:)
        integer :: i, j, k

        if (.false.) then
                allocate( n_vec(d%nfourier,3) )

                gamma_e = 1.76d2

                d%coeff_rf = d%coeff_rf/gamma_e
                do j = 1, d%nfourier
                        n_vec(j,:) = (/1.00d0, 1.00d0, 1.00d0/)
                        n_vec(j,:) = n_vec(j,:)*d%coeff_rf(j)
                        call y_rotatevec(d%theta_rf(j),n_vec(j,:))
                end do

                c%rms_a = 0.00d0
                do i = 1, d%nt
                        noisefield = 0.00d0
                        do j = 1, d%nfourier
                                do k = 1, 3
                                        noisefield(k) = noisefield(k) + &
                                        n_vec(j,k)*sin(d%w_rf(j)* &
                                        c%t_array(i) + d%phi_rf(j))
                                end do
                        end do
                        rms_tmp = 0.00d0
                        do k = 1, 3
                                rms_tmp = rms_tmp + noisefield(k)*noisefield(k)
                        end do
                        c%rms_a = c%rms_a + sqrt(rms_tmp)
                end do

                c%rms_a = c%rms_a/(dble(d%nt)*sqrt(3.00d0))
                write(10, *) 'rms amplitude (nT)', c%rms_a
        
        else
                c%rms_a = 0.00d0
        end if

        open(1, file = 'data/'//d%outputfile)
        do i = 1, d%ntheta
                write(1,*) d%theta(i), c%phis(i)
        end do


        end subroutine rms_amp

!-------------------------------------------------------------------------------------------------!
        
        end module hamiltonians
