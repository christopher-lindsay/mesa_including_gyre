    subroutine get_Delta_nu(id, Delta_nu, num_freqs) 
     use gyre_lib 
     use astero_lib 
     use astero_def 
         integer, intent(in) :: id
         real(dp), intent(out) :: Delta_nu 
         integer, intent(out) :: num_freqs
         integer :: ierr, n,i, num_sel
         real(dp), allocatable :: global_data(:), sel_ord(:), sel_nu(:) 
         real(dp), allocatable :: point_data(:,:)  
         real(dp) :: hold, nu_max, nu_range
         logical, allocatable :: use_n(:)
         type(star_info), pointer :: s
    
         ierr = 0 
         call star_ptr(id, s, ierr) 
         if (ierr /= 0) return 
         !   call star_get_pulse_data(s%id, 'GYRE', .TRUE., .FALSE., .TRUE., global_data, point_data,ierr) !logicals are add_center,
         if (ierr /=0) then                                                                             !keep_surface, add_atmosphere 
         print *, 'Failed when calling star_get_pulse_data' 
         return 
         end if 
         
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
       end subroutine get_Delta_nu


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