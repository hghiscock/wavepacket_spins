        module procedures
        implicit none
        
        contains

!-----------------------------------------------------------------------------!

        subroutine sprsin(a, n, sa, ija)

        !Converts a square matrix a(n,n) into a sparse matrix. Only elements with
        !magnitude > thresh are retained. Output is two linear arrays with dimension
        !nmax, sa contains array values and ija contains indices.

        integer, intent(in) :: n
        integer, allocatable, intent(out) :: ija(:)
        complex(8), intent(in) :: a(n,n)
        complex(8), allocatable, intent(out) :: sa(:)

        integer :: i, j, k, nmax
        real(8) :: thresh
        logical :: mask(n,n)

        mask = abs(a) .ne. 0
        nmax = count(mask)+1

        allocate( ija(nmax), sa(nmax) )

        thresh = 1.00d0
        ija = 0
        sa = (0.00d0, 0.00d0)
        do j = 1, n
            sa(j) = a(j,j)
        end do
        ija(1) = n+2
        k = n+1
        do i = 1, n
            do j = 1, n
                if (abs(a(i,j)) > thresh) then
                    if (i .ne. j) then
                        k = k+1
                        sa(k) = a(i,j)
                        ija(k) = j
                    end if
                end if
            end do
            ija(i+1) = k+1
        end do

        end subroutine sprsin

!-----------------------------------------------------------------------------!

        function sprsax(sa, ija, x, n)

        !Multiply a sparse matrix in storage arrays sa and ija by a square matrix x

        integer, intent(in) :: n, ija(:)
        complex(8), intent(in) :: sa(:)
        complex(8), intent(in) :: x(n,n)
        complex(8) :: sprsax(n,n)

        integer :: i, j, k

        sprsax = (0.00d0, 0.00d0)
        do j = 1, n
            do i = 1, n
                sprsax(i,j) = sa(i)*x(i,j)
                do k = ija(i), ija(i+1)-1
                    sprsax(i,j) = sprsax(i,j) + sa(k)*x(ija(k),j)
                end do
            end do
        end do  

        end function sprsax

!-----------------------------------------------------------------------------!

        real(8) function ran2(idum)

        !Generates a random number, initial seed idum usually = -1

        integer :: idum, IM1, IM2, IMM1, IA1, IA2, IQ1, IQ2, IR1, IR2, &
        NTAB, NDIV
        real(8) :: AM, EPS, RNMX
        parameter (IM1=2147483563, IM2=2147483399, AM=1./IM1, &
        IMM1=IM1-1, IA1=40014, &
        IA2=40692, IQ1=53668, IQ2=52774, IR1=12211, IR2=3791, NTAB=32, &
        NDIV=1+IMM1/NTAB, EPS=1.20e-7, RNMX=1.-EPS)
        integer :: idum2, j, k, iv(NTAB), iy
        save iv, iy, idum2
        data idum2 /123456789/, iv /NTAB*0/, iy /0/

        if (idum.le.0) then
                idum=max(-idum,1)
                idum2=idum
                do j = NTAB+8, 1, -1
                        k = idum/IQ1
                        idum = IA1*(idum-k*IQ1)-k*IR1
                        if (idum.lt.0) idum=idum+IM1
                        if (j.le.NTAB) iv(j)=idum
                end do
                iy = iv(1)
        end if

        k = idum/IQ1
        idum = IA1*(idum-k*IQ1)-k*IR1
        if (idum.lt.0) idum=idum+IM1
        k = idum2/IQ2
        idum2=IA2*(idum2-k*IQ2)-k*IR2
        if (idum2.lt.0) idum2=idum2+IM2
        j = 1+iy/NDIV
        iy = iv(j)-idum2
        iv(j) = idum
        if(iy.lt.1) iy = iy+IMM1
        ran2 = min(AM*iy, RNMX)
        return

        end function ran2

!-----------------------------------------------------------------------------!

      real(8) function gasdev(idum)

      !Generate normally distributed random number      

      integer :: idum

      integer :: iset
      real(8) :: fac, gset, rsq, v1, v2
      save iset, gset
      data iset/0/

      if (idum.lt.0) iset = 0

      if (iset.eq.0) then
      1   v1 = 2.*ran2(idum)-1
          v2 = 2.*ran2(idum)-1
          rsq = v1**2+v2**2
          if(rsq.ge.1..or.rsq.eq.0.) goto 1
          fac = sqrt(-2.*log(rsq)/rsq)
          gset = v1*fac
          gasdev = v2*fac
          iset = 1
      else
          gasdev = gset
          iset = 0
      end if

      return

      end function gasdev


!-----------------------------------------------------------------------------!

        real(8) function ampdist(freq, nfreq)

        !Use amplitude distribution as found in experiment
        !modify param.f90 file to include distribution

        real(8), intent(in) :: freq
        integer, intent(in) :: nfreq
        real(8) :: a, C, gamma_e
        integer :: idum

        ampdist = 50.00d0*gamma_e

        end function ampdist

!-----------------------------------------------------------------------------!

        subroutine kron(mat1, mat2, mat3, dim1, dim2)

        !Calculate the kronecker product of mat1 and mat2, result mat3

        implicit none

        complex(8), intent(in)  :: mat1(:,:), mat2(:,:)
        complex(8), intent(out) :: mat3(:,:)
        integer, intent(in) :: dim1, dim2
        integer :: I, J

        forall(I=1:dim1, J=1:dim1)
            mat3(dim2*(I-1)+1:dim2*I,dim2*(J-1)+1:dim2*J)=mat1(I,J)*mat2
        end forall

        end subroutine kron

!-----------------------------------------------------------------------------!

        subroutine y_rotate(theta, a)

        !Rotate a matrix A about the y axis

        real(8), intent(inout) :: a(:,:)
        real(8), intent(in) :: theta

        real(8), dimension(3,3) :: r_y

        r_y(:,1) = (/ cos(theta), 0.00d0, sin(theta) /)
        r_y(:,2) = (/ 0.00d0, 1.00d0, 0.00d0 /)
        r_y(:,3) = (/ -sin(theta), 0.00d0, cos(theta) /)

        a = matmul(r_y, matmul(a, transpose(r_y)))

        end subroutine y_rotate

!-----------------------------------------------------------------------------!

        subroutine y_rotatevec(theta, v)

        !Rotate a matrix vector about the y axis

        real(8), intent(inout) :: v(:)
        real(8), intent(in) :: theta
        real(8), dimension(3,3) :: r_y

        r_y(:,1) = (/ cos(theta), 0.00d0, sin(theta) /)
        r_y(:,2) = (/ 0.00d0, 1.00d0, 0.00d0 /)
        r_y(:,3) = (/ -sin(theta), 0.00d0, cos(theta) /)
        
        v = matmul(r_y, v)

        end subroutine y_rotatevec

!-----------------------------------------------------------------------------!

        end module procedures
