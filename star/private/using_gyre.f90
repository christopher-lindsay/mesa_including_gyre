! ***********************************************************************
!
!   Copyright (C) 2010-2019  The MESA Team
!
!   MESA is free software; you can use it and/or modify
!   it under the combined terms and restrictions of the MESA MANIFESTO
!   and the GNU General Library Public License as published
!   by the Free Software Foundation; either version 2 of the License,
!   or (at your option) any later version.
!
!   You should have received a copy of the MESA MANIFESTO along with
!   this software; if not, it is available at the mesa website:
!   http://mesa.sourceforge.net/
!
!   MESA is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!   See the GNU Library General Public License for more details.
!
!   You should have received a copy of the GNU Library General Public License
!   along with this software; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
!
! ***********************************************************************

      module using_gyre

      use star_private_def
      use const_def
      use utils_lib
      use gyre_lib
      use astero_lib
      use astero_def 

      implicit none

      private
      public :: get_Delta_nu, get_fundamental

      logical, parameter :: dbg = .false.



      contains

       

       subroutine get_Delta_nu(id, Delta_nu) 

         integer, intent(in) :: id
         real(dp), intent(out) :: Delta_nu 
         integer :: num_freqs
         integer :: ierr, n,i, num_sel
         real(dp), allocatable :: global_data(:), sel_ord(:), sel_nu(:) 
         real(dp), allocatable :: point_data(:,:)  
         real(dp) :: hold, nu_max, nu_range
         logical, allocatable :: use_n(:)
         type(star_info), pointer :: s
         ierr = 0 
         call star_ptr(id, s, ierr) 
         if (ierr /= 0) return 


        if (GYRE_IS_ENABLED) then         
         
         ! first try local directory
         filename = gyre_inlist_name

         if(level==1) then ! First pass either the user set the file or we load the defaults
            if (len_trim(filename) == 0) filename = 'gyre.in'

            exists=.false.
            inquire(file=filename,exist=exists)

            if(.not.exists) filename = trim(mesa_dir) // '/star/defaults/gyre_delta_nu.defaults'
         else
            ! User had include '' in their profile_columns.list file, so dont try to load the local one, jump to the defaults
            if (len_trim(filename) == 0) filename =trim(mesa_dir) // '/star/defaults/gyre_delta_nu.defaults'
         end if

         call init_gyre(filename)
         
         num_results = 0 

         call astero_gyre_get_modes(id, 0, .TRUE., ierr) 
         if (ierr /=0) then 
         print *, 'Failed when calling do_gyre_get_modes' 
         return 
         end if

         nu_max = s% nu_max 

         nu_range = 0.5 * nu_max**0.9 ! Mosser et al. 2010 

         num_sel = 0 

         allocate(use_n(num_results)) 
         use_n(:) = .false. 

         do n = 1, num_results 
         if (cyclic_freq(n) > nu_max-nu_range .and. cyclic_freq(n) < nu_max+nu_range) then 
              num_sel = num_sel + 1 
              use_n(n) = .true. 
         end if 
         end do 

         allocate(sel_ord(num_sel)) 
         allocate(sel_nu(num_sel)) 

         i = 1 
         do n =1, num_results 
         if (use_n(n)) then 
             sel_ord(i) = order(n) 
             sel_nu(i) = cyclic_freq(n) 
             i = i+1 
         endif 
         enddo 

         num_freqs = num_sel  
         call get_lsq_fit(num_sel, sel_ord, sel_nu, Delta_nu, hold) ! Or any built in function that will get a slope of a line fit to data points 
         
         deallocate(use_n) 
         deallocate(sel_ord) 
         deallocate(sel_nu) 
         
         
         ! Clear out frequency results     
         num_results = 0 
         el = 0 
         order = 0 
         em = 0 
         cyclic_freq = 0 
         growth_rate = 0 
         inertia = 0 

         return 
       endif

       end subroutine get_Delta_nu



       subroutine run_gyre_fund(id, ierr)

           integer, intent(in)  :: id
           integer, intent(out) :: ierr

           real(dp), allocatable :: global_data(:)
           real(dp), allocatable :: point_data(:,:)
           integer               :: ipar(0)
           real(dp)              :: rpar(0)

           type (star_info), pointer :: s
           ierr = 0
           call star_ptr(id, s, ierr)
           if (ierr /= 0) return

           ! Pass model data to GYRE
           if (allocated(xi_r_radial)) deallocate(xi_r_fund)
           allocate(xi_r_radial(s%nz))
           xi_r_fund(:)=0
           if (allocated(xi_r_dipole)) deallocate(xi_r_1o)
           allocate(xi_r_dipole(s%nz))
           xi_r_1o(:)=0

           ! We need a different name for GYRE.in for fundamental mode
           call star_get_pulse_data(id, 'GYRE', .FALSE., .TRUE., .FALSE., &
                global_data, point_data, ierr)
           if (ierr /= 0) then
              print *,'Failed when calling star_get_pulse_data'
              return
           end if

           call gyre_set_model(global_data, point_data, 101)

           ! Run GYRE to get modes
           call gyre_get_modes(0, process_fund_mode, ipar, rpar)

           ! show that fundamental gyre has run
           gyre_fund_has_run = .true.

        contains

           subroutine process_radial_modes(md, ipar, rpar, retcode)
              type(mode_t), intent(in) :: md
              integer, intent(inout)   :: ipar(:)
              real(dp), intent(inout)  :: rpar(:)
              integer, intent(out)     :: retcode
              integer :: k

              type (star_info), pointer :: s
              ierr = 0
              call star_ptr(id, s, ierr)
              if (ierr /= 0) return

            if (md%n_p >= 1 .and. md%n_p <= 100) then
                ! Print out degree, radial order, mode inertia, and frequency
                ! print *, 'Found mode: l, n_p, n_g, E, nu = ', &
                !     md%l, md%n_p, md%n_g, md%E_norm(), REAL(md%freq('HZ'))

                if (md%l == 0) then ! radial modes
                    frequencies(md%l+1, md%n_p) = (md%freq('UHZ') - s% nu_max) / s% delta_nu

                    if (md%n_p == 1) then ! store the fundamental eigenfunction
                       if (allocated(xi_r_fund)) deallocate(xi_r_fund)
                       allocate(xi_r_fund(md%n_k))
                       do k = 1, md%n_k
                          xi_r_radial(k) = md%xi_r(k)
                       end do
                       xi_r_radial = xi_r_radial(md%n_k:1:-1)

                    else if (md%n_p == 2) then ! store the 1o eigenfunction
                       if (allocated(xi_r_1o)) deallocate(xi_r_1o)
                       allocate(xi_r_fund(md%n_k))
                       do k = 1, md%n_k
                          xi_r_radial(k) = md%xi_r(k)
                       end do
                       xi_r_radial = xi_r_radial(md%n_k:1:-1) 
                    end if




                else if (inertias(md%n_p) > 0 .and. md%E_norm() > inertias(md%n_p)) then
                    write (*,*) 'Skipping mode: inertia higher than already seen'
                else ! non-radial modes

                    ! choose the mode with the lowest inertia
                    inertias(md%n_p) = md%E_norm()
                    frequencies(md%l+1, md%n_p) = (md%freq('UHZ') - s% nu_max) / s% delta_nu
                    ng_array(md%n_p) = md%n_g


                    if (md%n_p == s% x_integer_ctrl(1) - 1) then ! store the eigenfunction 
                       if (allocated(xi_r_dipole)) deallocate(xi_r_dipole)
                       allocate(xi_r_dipole(md%n_k))

                       write(*, *) "nk is", md%n_k 
                       write(*, *) "nz is", s%nz

                       do k = 1, md%n_k
                          xi_r_dipole(k) = md%xi_r(k)
                       end do
                       xi_r_dipole = xi_r_dipole(md%n_k:1:-1)
                    end if


                end if
            end if

            retcode = 0
         end subroutine process_radial_modes


        end subroutine run_gyre_fund


       
       subroutine get_lsq_fit(n, x, y, m, b)
         integer, intent(in) :: n 
         real(dp), intent(in) :: x(n), y(n) 
         real(dp), intent(out) :: m, b 
         integer :: i
         real(dp) :: sumX, sumX2, sumY, sumXY 

         sumX=0 
         sumX2=0
         sumY=0 
         sumXY=0

         do i = 1, n
         sumX = sumX + x(i)
         sumX2 = sumX2 + pow(x(i),2) 
         sumY = sumY + y(i) 
         sumXY = sumXY + x(i)*y(i) 
         end do 

         m = (n * sumXY - sumX * sumY)/(n*sumX2 - sumX*sumX) 
         b = (sumY - m*sumX)/n 

         return 
       end subroutine get_lsq_fit


  


      end module using_gyre

