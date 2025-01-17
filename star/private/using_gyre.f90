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

      implicit none

      private
      !public :: get_Delta_nu, get_nth_radial_mode

      logical, parameter :: dbg = .false.
      real(dp) :: radial_freqs(100), radial_orders(100), radial_inertias(100) 
      integer :: num_radial_modes         
      logical :: add_atmosphere, keep_surface_point, add_center_point 
      logical :: gyre_radial_modes_calculated = .false. 
      logical :: GYRE_IS_ENABLED = .true.
      contains


       subroutine get_Delta_nu(id, Delta_nu) 
       
         integer, intent(in) :: id
         real(dp), intent(out) :: Delta_nu 
         integer :: num_freqs, num_results, num_sel, level
         integer :: ierr, n,i, num_seli
         character(len=1000) :: filename
         real(dp), allocatable :: global_data(:), sel_ord(:), sel_nu(:) 
         real(dp), allocatable :: point_data(:,:)  
         real(dp) :: hold, nu_max, nu_range
         logical, allocatable :: use_n(:)
         type(star_info), pointer :: s
         ierr = 0 
         call star_ptr(id, s, ierr) 
         if (ierr /= 0) return 

       
        if (GYRE_IS_ENABLED) then         

         if (.not. gyre_radial_modes_calculated) call get_radial_modes(id,s,ierr)


           num_results = 0 
           nu_max = s% nu_max
           nu_range = 0.5 * nu_max**0.9 ! Mosser et al. 2010         

           num_sel = 0 

           allocate(use_n(num_results)) 
           use_n(:) = .false. 

           do n = 1, num_results 
           if (radial_freqs(n) > nu_max-nu_range .and. radial_freqs(n) < nu_max+nu_range) then 
                num_sel = num_sel + 1 
                use_n(n) = .true. 
           end if 
           end do 

           allocate(sel_ord(num_sel)) 
           allocate(sel_nu(num_sel))  
       
           i = 1 
           do n =1, num_results 
           if (use_n(n)) then 
               sel_ord(i) = radial_orders(n) 
               sel_nu(i) = radial_freqs(n) 
               i = i+1 
           endif 
           enddo 

           num_freqs = num_sel  
           call get_lsq_fit(num_sel, sel_ord, sel_nu, Delta_nu, hold) ! Or any built in function that will get a slope of a line fit to data points 
         
            deallocate(use_n) 
            deallocate(sel_ord) 
            deallocate(sel_nu) 
            
          
          return 
        
        else 
         write(*,*) "GYRE is not enabled"
         return 

       endif 

      end subroutine get_Delta_nu
      


        subroutine get_radial_modes(id,s,ierr) 
            integer, intent(in) :: id 
            integer, intent(out) :: ierr 
            type(star_info), pointer :: s 
            integer :: numl0, level
            logical :: exists
            character(len=1000) :: filename
            real(dp), allocatable :: global_data(:) 
            real(dp), allocatable :: point_data(:,:) 
            integer               :: ipar(0), num_modes, i ,itemp, istart 
            real(dp)              :: rpar(0), temp_freqs(100), temp_ens(100)  
            real(dp) :: temp 
            call star_ptr(id, s, ierr) 
            if (ierr /= 0) return
 
            add_atmosphere = s% add_atmosphere_to_pulse_data
            keep_surface_point = s% add_center_point_to_pulse_data
            add_center_point = s% keep_surface_point_for_pulse_data


            ! load the gyre inlist and initialize gyre
            ! first try local directory
            filename = s% job% gyre_inlist_name
       
            if(level==1) then ! First pass either the user set the file or we load the defaults
               if (len_trim(filename) == 0) filename = 'gyre.in'
       
               exists=.false.
               inquire(file=filename,exist=exists)
       
               if(.not.exists) filename = trim(mesa_dir) // '/star/defaults/gyre_delta_nu.defaults'
            else
               ! User had include '' in their gyre.in file, so dont try to load the local one, jump to the defaults
               if (len_trim(filename) == 0) filename =trim(mesa_dir) // '/star/defaults/gyre_delta_nu.defaults'
            end if
       
       
            ! initialize gyre and set constants using the subroutines in gyre_lib 
            call gyre_init(filename)
       
            call gyre_set_constant('G_GRAVITY', standard_cgrav)
            call gyre_set_constant('C_LIGHT', clight)
            call gyre_set_constant('A_RADIATION', crad)
       
            call gyre_set_constant('M_SUN', Msun)
            call gyre_set_constant('R_SUN', Rsun)
            call gyre_set_constant('L_SUN', Lsun)
       
            call gyre_set_constant('GYRE_DIR', TRIM(mesa_dir)//'/gyre/gyre')




            ! Pass model to GYRE 
            ! Logcials are add_center, keep_surface,add_atmosphere
            call star_get_pulse_data(id, 'GYRE', add_center_point, & 
                 keep_surface_point , add_atmosphere, & 
                    global_data, point_data, ierr) 
            if (ierr /= 0) then 
                    print *, 'Failed when calling star_get_pulse_data' 
                    return 
            end if 
            call gyre_set_model(global_data, point_data, 101) 
            !Zero out previous data 
            radial_freqs(:) = 0d0
            radial_orders(:) = 0  
            radial_inertias(:) = 0d0 
            ! Run GYRE to get modes
            ! Get l = 0 modes 
            num_modes = 0 
            temp_freqs(:) = 0 
            temp_ens(:) = 0 
            call gyre_get_modes(0, process_mode, ipar, rpar) 
            
            num_radial_modes = num_modes 
       
        contains 
          subroutine process_mode (md, ipar, rpar, retcode) 
            type(mode_t), intent(in) :: md
            integer, intent(inout) :: ipar(:) 
            real(dp), intent(inout) :: rpar(:) 
            integer, intent(out) :: retcode 
            integer :: k 
            type (star_info), pointer :: s
            ierr = 0 
            call star_ptr(id, s, ierr) 
            if (ierr/=0) return 
            if (md%n_p >= 1 .and. md%n_p <=50) then 
              num_modes = num_modes + 1 
              radial_freqs(num_modes) = REAL(md%freq('UHZ')) 
              radial_orders(num_modes) = md% n_p
              radial_inertias(num_modes) = md% E_norm()
            end if 
            retcode = 0 
          end subroutine process_mode 
        end subroutine get_radial_modes

      !subroutine run_gyre_radial(id, ierr)

      !    integer, intent(in)  :: id
      !    integer, intent(out) :: ierr

      !    real(dp), allocatable :: global_data(:)
      !    real(dp), allocatable :: point_data(:,:)
      !    integer               :: ipar(0)
      !    real(dp)              :: rpar(0)

      !    type (star_info), pointer :: s
      !    ierr = 0
      !    call star_ptr(id, s, ierr)
      !    if (ierr /= 0) return

      !    ! Pass model data to GYRE
      !    if (allocated(xi_r_fund)) deallocate(xi_r_fund)
      !    allocate(xi_r_fund(s%nz))
      !    xi_r_fund(:) = 0
      !    
      !    if (allocated(xi_r_1o)) deallocate(xi_r_1o)
      !    allocate(xi_r_1o(s%nz))
      !    xi_r_1o(:) = 0

      !    ! We need a different name for GYRE.in for fundamental mode
      !    call star_get_pulse_data(id, 'GYRE', .FALSE., .TRUE., .FALSE., &
      !         global_data, point_data, ierr)
      !    if (ierr /= 0) then
      !       print *,'Failed when calling star_get_pulse_data'
      !       return
      !    end if

      !    call gyre_set_model(global_data, point_data, 101)

      !    ! Run GYRE to get modes
      !    call gyre_get_modes(0, process_radial_modes, ipar, rpar)

      !    ! show that fundamental gyre has run
      !    gyre_fund_has_run = .true.

      ! contains

      !    subroutine process_radial_modes(md, ipar, rpar, retcode)
      !       type(mode_t), intent(in) :: md
      !       integer, intent(inout)   :: ipar(:)
      !       real(dp), intent(inout)  :: rpar(:)
      !       integer, intent(out)     :: retcode
      !       integer :: k

      !       type (star_info), pointer :: s
      !       ierr = 0
      !       call star_ptr(id, s, ierr)
      !       if (ierr /= 0) return

      !     if (md%n_p >= 1 .and. md%n_p <= 100) then

      !         if (md%l == 0) then ! radial modes
      !             frequencies(md%l+1, md%n_p) = md%freq('UHZ')

      !             if (md%n_p == 1) then ! store the fundamental eigenfunction
      !                if (allocated(xi_r_fund)) deallocate(xi_r_fund)
      !                allocate(xi_r_fund(md%n_k))
      !                do k = 1, md%n_k
      !                   xi_r_fund(k) = md%xi_r(k)
      !                end do
      !                xi_r_fund = xi_r_fund(md%n_k:1:-1)

      !             else if (md%n_p == 2) then ! store the 1o eigenfunction
      !                if (allocated(xi_r_1o)) deallocate(xi_r_1o)
      !                allocate(xi_r_1o(md%n_k))
      !                do k = 1, md%n_k
      !                   xi_r_1o(k) = md%xi_r(k)
      !                end do
      !                xi_r_1o = xi_r_1o(md%n_k:1:-1) 
      !             end if

      !         end if
      !     end if

      !     retcode = 0
      !  end subroutine process_radial_modes


      ! end subroutine run_gyre_radial


       
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

