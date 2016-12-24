c trapp.f -- Parallel Trapezoidal Rule, first version
c
c Input: None.
c Output: Estimate of the integral from a to b of f(x)
c         using the trapezoidal rule and n trapezoids.
c Algorithm:
c   1. Each process calculates "its" interval of
c       integration.
c   2. Each process estimates the integral of f(x)
c       over its interval using the trapezoidal rule.
c   3a. Each process != 0 sends its integral to 0.
c   3b. Process 0 sums the calculations received from
c       the individual processes and prints the result.
c  Note:  f(x), a, b, and n are all hardwired.
c
      program trapezoidal
c
      IMPLICIT NONE
      include 'mpif.h'
c
      integer   my_rank  ! My process rank.
      integer   p        ! The number of processes.
      real      a        ! Left endpoint.
      real      b        ! Right endpoint.
      integer   n        ! Number of trapezoids.
      real      h        ! Trapezoid base length.
      real      local_a  ! Left endpoint for my process.
      real      local_b  ! Right endpoint my process.
      integer   local_n  ! Number of trapezoids for my
                         ! calculation.
      real      integral ! Integral over my interval.
      real      total    ! Total integral.
      integer   source   ! Process sending integal.
      integer   dest     ! All messages go to 0.
      integer   tag
      integer   status(MPI_STATUS_SIZE)
      integer   ierr
      
      real      Trap
      real   time_start
      real   time_end
      
      data a, b, n, dest, tag /0.0, 1.0, 1024, 0, 50/

      call CPU_TIME(time_start)

C Let the system do what it needs to start up MPI.
      call MPI_INIT(ierr)

C Get my process rank.
      call MPI_COMM_RANK(MPI_COMM_WORLD, my_rank, ierr)
C Find out how many processes are being used.
      call MPI_COMM_SIZE(MPI_COMM_WORLD, p, ierr)
      h = (b-a)/n        ! h is the same for all processes.
      local_n = n/p      ! So is the number of trapezoids.
C Length of each process' interval of integration = local_n*h.
C So my interval starts at :
      local_a = a + my_rank*local_n*h
      local_b = local_a + local_n*h
      integral = Trap(local_a, local_b, local_n, h)
C Add up the integals calculated by each process.
      if (my_rank .EQ. 0) then
         total = integral
         do source = 1, p-1
             call MPI_RECV(integral, 1, MPI_REAL, source, tag,
     +              MPI_COMM_WORLD, status, ierr)
             total = total + integral
         enddo
      else
         call MPI_SEND(integral, 1, MPI_REAL, dest,
     +          tag, MPI_COMM_WORLD, ierr)
      endif
c   Print the result.
      if (my_rank .EQ. 0) then
            write(6,200) n
 200        format(' ','With n = ',I4,' trapezoids, our estimate')
            write(6,300) a, b, total
 300        format(' ','of the integral from ',f6.2,' to ',f6.2,
     +             ' = ',f11.5)
      endif
C       Shut down MPI.
      call MPI_FINALIZE(ierr)

      call CPU_TIME(time_end)
      write(6,400) time_end - time_start
400   format(' ', 'The elapsed time is ',f11.5)
      end program trapezoidal
c      
c      
c
      real function Trap(local_a, local_b, local_n, h)

      IMPLICIT NONE

      real     local_a
      real     local_b
      integer  local_n
      integer  i
      real     h
      real     x
      real     f
      real     integ
      integ = (f(local_a) + f(local_b))/2.0
      x = local_a

      do i = 1, local_n-1
          x = x + h
          integ = integ + f(x)
      enddo

      integ = integ*h
      Trap = integ
      return
      end
c
c
      real function f(x)
      IMPLICIT NONE
      real x
      f = x*x
      end
