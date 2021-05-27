        program s_yield
        use param
        use hamiltonians
        use procedures
        use sc_tensor
        implicit none

!------------------------------------------------------------------------------!
        type(hamiltonian) :: ha, hb
        type(spin_corr) :: sca, scb
        type(radical_dat) :: da, db
        type(param_dat) :: p_d
        type(calc_dat) :: c_d

        integer :: i, j, k, t
        integer :: cr

!------------------------------------------------------------------------------!

        !time code
        call system_clock(count_rate = cr)

        !read in parameters
        call p_d%read_in()                                     
        call da%read_in(p_d,1)                                      
        call db%read_in(p_d,2)                                      
        call p_d%initialise() 
        call c_d%pauli_matrices() 
        call c_d%get_time(p_d)
        call omp_set_num_threads(p_d%ncores)

!------------------------------------------------------------------------------!

        !build hamiltonians
        call ha%build(c_d, da)
        call hb%build(c_d, db)
        call ha%noise(p_d)
        call hb%noise(p_d)

        !initialise arrays
        call sca%create(ha, p_d)
        call scb%create(hb, p_d)
        call system_clock(p_d%start)
       
!------------------------------------------------------------------------------!

        !propagate over nt time steps
        do i = 1, p_d%nt
                call sca%calc_sc(ha, p_d)
                call scb%calc_sc(hb, p_d)

                call c_d%calc_ps(p_d, sca%r_ab, scb%r_ab)
                call c_d%integrate_sy(p_d)
                
                if (mod(i,p_d%nt/100) == 0) then
                  call system_clock(p_d%finish)
                  print *, i/(p_d%nt/100), '% done'
                  print *, 'Time', (p_d%finish - p_d%start)/dble(cr)
                  call system_clock(p_d%start)
                end if
                call flush()

        end do

!------------------------------------------------------------------------------!

        !write output
        call p_d%write_log_file(da, db)
        call c_d%rms_amp(p_d)
        
        end program s_yield
