program main
   use types, only: dp
   use parameters
   use utils
   use energies
   use movement
   use, intrinsic :: iso_fortran_env, only: output_unit
   implicit none

   ! Local scalar variables
   real(dp) :: d, volratio, rng, ener
   real(dp) :: rhoaverage, rhoave, current_volume, rhoprom
   integer :: rngint, i, j, k, l, nattemp, nacc
   integer :: thermsteps, eqsteps, u, vacc, vattemp
   integer :: v, avevolfreq, accsize, w
   ! For writing to files
   character(len=8) :: i_char
   character(len=14) :: j_char, k_char, l_char, m_char
   ! Local arrays
   real(dp), allocatable :: x(:), y(:)
   real(dp), allocatable :: rhoacc(:), volacc(:), volsqacc(:)

   ! Initialize the RNG
   call random_init(.false., .true.)

   ! Read an input file that contains all the necessary information
   call parse_input('input.in', eqsteps, thermsteps, avevolfreq)

   ! Update the simulation parameters with this information
   boxl = sqrt(np / rho)
   rc = boxl / 2.0_dp
   d = sqrt(1.0_dp / rho)

   ! Initialization of variables
   nattemp = 0
   nacc = 1
   rhoaverage = 0.0_dp
   rhoave = 0.0_dp
   rhoprom = 0.0_dp
   vacc = 1
   vattemp = 0
   j = 0
   current_volume = 0.0_dp

   ! Write information to screen
   write(unit=output_unit, fmt='(a,f8.4)') 'Cut-off radius = ', rc
   write(unit=output_unit, fmt='(a,f8.4)') 'Mean interparticle distance: ', d
   write(unit=output_unit, fmt='(a,f8.4)') 'Reduced pressure = ', pressure
   write(unit=output_unit, fmt='(a,f8.4)') 'Reference reduced density = ', rho
   write(unit=output_unit, fmt='(a,f8.4)') 'Reduced temperature = ', ktemp
   write(unit=output_unit, fmt='(a,f8.4)') 'Shoulder length = ', lambdasw

   ! Allocate memory for arrays
   allocate(x(np), y(np))
   accsize = eqsteps / avevolfreq
   allocate(rhoacc(accsize), volacc(accsize), volsqacc(accsize))
   rhoacc = 0.0_dp
   volacc = 0.0_dp
   volsqacc = 0.0_dp

   ! Either read a configuration file or generate a new one
   if (from_file) then
      write(unit=output_unit, fmt='(a)') 'Reading from positions file...'
      open(newunit=u, file = 'configuration.dat', status = 'unknown')
      do i = 1, np
         read(u, *) x(i), y(i)
      end do
      close(u)
   else
      ! initial configuration as a simple lattice
      call iniconfig(x, y, d)
      ! write the initial configuration to file
      open(newunit=u, file = 'init_conf.dat', status = 'unknown')
      write(unit=u, fmt='(i3)') np
      do i = 1, np
         write(unit=u, fmt='(2f12.8)') x(i), y(i)
      end do
      close(u)
   end if

   ! Initial configuration energy, regardless of how it was created
   call energy(x, y, ener)
   write(unit=output_unit, fmt='(a,f8.4)') 'E/N for the initial configuration:', ener/np

   ! MC cycle to thermalize the system
   do i = 1, thermsteps
      ! Loop for the total number of particles
      do k = 1, np
         call random_number(rng)
         rngint = 1 + floor((np + 1) * rng)

         if (rngint <= np) then
            call mcmove(x, y, ener, nattemp, nacc)
            call adjust(nattemp, nacc, del, 0.4_dp, 2)
         else
            call mcvolume(x, y, rhoave, ener, vattemp, vacc)
            call adjust(vattemp, vacc, dispvol, 0.2_dp, 1)
         end if
      end do

      ! Print information to screen
      if (mod(i, 10000) == 0) then
         write(unit=output_unit, fmt='(a)') 'MC Step, Particle disp, Energy / N, Disp ratio'
         write(unit=output_unit, fmt='(i10, 3f8.4)') i, del, ener/real(np, dp), &
         & real(nacc, dp) / real(nattemp, dp)
         write(unit=output_unit, fmt='(a)') 'MC Step, Density average, box size, Vol ratio, Vol disp'
         volratio = real(vacc, dp) / real(vattemp, dp)
         write(unit=output_unit, fmt='(i10, 4f8.4)') i, rhoave / real(vacc, dp), &
         & boxl, volratio, dispvol
      end if
   end do

   ! Reset accumulation variables
   nattemp = 0
   nacc = 1
   vattemp = 0
   vacc = 1
   ! Start accumulating results
   write(unit=output_unit, fmt='(a)') 'Averaging starts...'
   ! Open the necessary files for saving thermodynamical quantities
   open(newunit=u, file='energy.dat', status='unknown')
   open(newunit=v, file='density.dat', status='unknown')
   ! Opne file for saving snapshots
   open(newunit=w, file='trajectory.dat', status='unknown')
   ! Production cyle
   do i = 1, eqsteps
      do k = 1, np
         call random_number(rng)
         rngint = 1 + floor((np + 1) * rng)

         if (rngint <= np) then
            call mcmove(x, y, ener, nattemp, nacc)
         else
            call mcvolume(x, y, rhoave, ener, vattemp, vacc)
         end if
      end do

      if (mod(i, avevolfreq) == 0) then
         ! Save the value for the energy
         write(unit=u, fmt='(i12.2, f15.10)') i, ener/real(np, dp)

         ! Update the accumulation index
         j = j + 1

         ! Accumulate the results for the density
         current_volume = real(np, dp) / rho
         rhoaverage = rhoaverage + rho
         rhoprom = rhoaverage / real(j, dp)
         rhoacc(j) = rho
         volacc(j) = current_volume
         volsqacc(j) = current_volume**2.0_dp

         ! Save all results to file
         write(unit=v, fmt='(2f17.10)') rhoprom, current_volume

         ! Save the positions and simulation information
         write(unit=w, fmt='(i3)') np
         write (i_char, '(i8)') 0
         write(unit=w, fmt='(3a2)') adjustl(i_char), adjustl(i_char), adjustl(i_char)
         write (j_char, '(f10.6)') boxl
         write (k_char, '(f8.6)') 0.0_dp
         write(unit=w, fmt='(a10, 2a9)') adjustl(j_char), adjustl(k_char), adjustl(k_char)
         write(unit=w, fmt='(a9, a10, a9)') adjustl(k_char), adjustl(j_char), adjustl(k_char)
         write (j_char, '(f8.6)') 1.0_dp
         write(unit=w, fmt='(3a9)') adjustl(k_char), adjustl(k_char), adjustl(j_char)
         do l = 1, np
            write (i_char, '(i8)') 1
            write (j_char, '(f12.6)') x(l)
            write (k_char, '(f14.6)') y(l)
            write (l_char, '(f8.5)') 0.0_dp
            write (m_char, '(f8.5)') 1.0_dp
            write(unit=w, fmt='(a11, a11, a8, a8, a)') adjustl(j_char), adjustl(k_char), &
            & adjustl(l_char), adjustl(m_char), adjustl(i_char)
         end do
         write(unit=w, fmt='(a)')
      end if
   end do

   ! Close off files that were opened for saving information
   close(u)
   close(v)

   ! Do some averaging for the density
   call calc_variable(rhoacc, 'Density', 'average_density.dat')
   ! Do some averaging for the volume
   call calc_variable(volacc, 'Volume', 'average_volume.dat')
   ! Do some averaging for the squared volume
   call calc_variable(volsqacc, 'Squared Volume', 'average_sqvolume.dat')

   ! write the final configuration to file
   open(newunit=u, file = 'configuration.dat', status = 'unknown')
   do i = 1, np
      write(unit=u, fmt='(2f14.8)') x(i), y(i)
   end do
   close(u)
end program main
