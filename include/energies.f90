module energies
   use types
   use parameters
   use potentials
   use omp_lib

   implicit none

   public energy, denergy

contains
   ! This configuration calculates the energy of a given configuration
   subroutine energy(x, y, ener)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), intent(out) :: ener

      ! Local variables
      integer :: i, j
      real(dp) :: rij, xij, yij, uij

      ener = 0.0_dp

      do i = 1, (np - 1)
         do j = (i + 1), np
            uij = 0.0_dp

            xij = x(j) - x(i)
            yij = y(j) - y(i)

            ! Minimum image convention
            xij = xij - boxl*nint(xij / boxl)
            yij = yij - boxl*nint(yij / boxl)

            rij = norm2([xij, yij])

            if (rij < rc) then
               call hardsphere(rij, uij)
               ener = ener + uij
            end if
         end do
      end do
   end subroutine energy

   ! This subroutine calculates the difference in energy when a particle is displaced
   subroutine denergy(x, y, no, dener)
      real(dp), intent(in) :: x(:), y(:)
      real(dp), intent(out) :: dener
      integer, intent(in) :: no
      ! Local variables
      integer :: i
      real(dp) :: rij, xij, yij, uij

      dener = 0.0_dp ! initializing

      do i = 1, np
         if ( i == no ) cycle

         xij = x(no)-x(i)
         yij = y(no)-y(i)

         ! Minimum image convention
         xij = xij-boxl*nint(xij/boxl)
         yij = yij-boxl*nint(yij/boxl)

         rij = norm2([xij, yij])

         if (rij < rc) then
            call hardsphere(rij, uij)
            dener = dener + uij
         end if
      end do
   end subroutine denergy
end module energies
