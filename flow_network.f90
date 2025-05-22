module flow_network
  ! Sparse Poiseuille flow solver on a general fracture network.
  ! Currently uses simple Gauss–Seidel iterations; suitable for a few
  ! hundred elements which is the typical size of the step-wise HF mesh.
  use fracture_types
  implicit none
  private
  public :: solve_hf_pressure
contains
  !------------------------------------------------------------------
  subroutine solve_hf_pressure(num_hf, aperture_hf, mu, dt, ds_hint, p_old_hf, &
                               inlet_type, inlet_val, leak_hf, p_new_hf)
    ! Arguments identical to the legacy routine so the drivers do not
    ! need to change, but inside we build a global solution that spans
    ! *all* active fracture elements (HF + NF) stored in frac_list.
    integer, intent(in)  :: num_hf
    real*8 , intent(in)  :: aperture_hf(num_hf)      ! ignored; we read from frac_list
    real*8 , intent(in)  :: mu, dt, ds_hint          ! mu fluid viscosity, ds for HF tips
    real*8 , intent(in)  :: p_old_hf(num_hf)
    integer, intent(in)  :: inlet_type               ! 0 = Dirichlet P, 1 = rate Q
    real*8 , intent(in)  :: inlet_val
    real*8 , intent(in)  :: leak_hf(num_hf)          ! leak only for HF faces here
    real*8 , intent(out) :: p_new_hf(num_hf)

    if (.not.allocated(frac_list)) then
       p_new_hf = p_old_hf
       return
    end if
    integer :: nTot, i, iter, j, neigh
    real*8, parameter :: tol = 1.0d-4
    integer, parameter :: maxIter = 200
    real*8, allocatable :: p_old(:), p_new(:), rhs(:), vol(:), Tsum(:)

    nTot = size(frac_list)
    allocate(p_old(nTot), p_new(nTot), rhs(nTot), vol(nTot), Tsum(nTot))

    !----------------------------------------------------------------
    ! 1.  Gather old pressures and element volumes (approx aperture*area)
    !----------------------------------------------------------------
    real*8 :: area_i, dx1,dy1,dz1, dx2,dy2,dz2, face_area

    do i = 1, nTot
       p_old(i) = frac_list(i)%pressure

       ! crude area estimate: half crossed-diagonal magnitude
       dx1 = frac_list(i)%x(3) - frac_list(i)%x(1)
       dy1 = frac_list(i)%y(3) - frac_list(i)%y(1)
       dz1 = frac_list(i)%z(3) - frac_list(i)%z(1)
       dx2 = frac_list(i)%x(4) - frac_list(i)%x(2)
       dy2 = frac_list(i)%y(4) - frac_list(i)%y(2)
       dz2 = frac_list(i)%z(4) - frac_list(i)%z(2)
       face_area = 0.5d0 * sqrt( (dy1*dz2 - dz1*dy2)**2 + &
                                 (dz1*dx2 - dx1*dz2)**2 + &
                                 (dx1*dy2 - dy1*dx2)**2 )
       vol(i) = max(frac_list(i)%aperture,1.d-9)*face_area   ! [m^3]
    end do

    ! RHS for backward Euler discretisation: vol/dt * p_old - leak
    rhs = (vol/dt) * p_old

    ! assign leak for HF faces
    do i = 1, num_hf
       rhs(i) = rhs(i) - leak_hf(i)
    end do

    ! assign Carter leak for active NF faces
    do i = num_hf+1, nTot
       if (frac_list(i)%is_active .and. .not.frac_list(i)%is_hf) then
          rhs(i) = rhs(i) - frac_list(i)%leak
       end if
    end do

    !------------------------------------------------------------
    ! 2A. If inlet is rate-controlled, we perform an outer Newton
    !     iteration on the inlet pressure Pin so that the flow rate
    !     between the inlet face (index 1) and its first active
    !     neighbour matches inlet_val.
    !------------------------------------------------------------
    integer, parameter :: maxNewton = 12
    real*8, parameter :: tolQ = 1.0d-6
    real*8 :: Pin, Qcalc, Tij_inlet, dx_in,dy_in,dz_in, ds_in
    integer :: neigh_inlet, nIt

    if (inlet_type == 1) then
       ! choose first active neighbour of face 1
       neigh_inlet = -1
       do j = 1, size(frac_list(1)%neighbors)
          if (frac_list(frac_list(1)%neighbors(j))%is_active) then
             neigh_inlet = frac_list(1)%neighbors(j)
             exit
          end if
       end do
       if (neigh_inlet < 0) then
          ! fallback: single face line, treat as Dirichlet
          inlet_type = 0
       end if
    end if

    Pin = p_old_hf(1)

    do nIt = 1, merge(maxNewton,1,inlet_type==1)

       ! Set initial guess pressure field
       p_new = p_old
       if (inlet_type == 1) p_new(1) = Pin

       ! === Gauss-Seidel loop (same as before) ===
       do iter = 1, maxIter
          real*8 :: maxErr
          maxErr = 0.0d0
          do i = 1, nTot
             if (.not.frac_list(i)%is_active) cycle
             Tsum(i) = 0.0d0
             do j = 1, size(frac_list(i)%neighbors)
                neigh = frac_list(i)%neighbors(j)
                if (.not.frac_list(neigh)%is_active) cycle
                dx = centroid_x(i) - centroid_x(neigh)
                dy = centroid_y(i) - centroid_y(neigh)
                dz = centroid_z(i) - centroid_z(neigh)
                ds_ij = sqrt(dx*dx + dy*dy + dz*dz)
                if (ds_ij<=0.d0) cycle
                w_avg = 0.5d0*(frac_list(i)%aperture + frac_list(neigh)%aperture)
                Tij   = (w_avg**3)/(12.d0*mu*ds_ij)
                Tsum(i) = Tsum(i) + Tij
             end do
          end do

          do i = 1, nTot
             if (.not.frac_list(i)%is_active) cycle

             if (i==1) then
                if (inlet_type==0) then
                   p_new(i) = inlet_val
                   cycle
                else
                   p_new(i) = Pin
                   cycle
                end if
             end if

             sumTpj = 0.0d0
             do j = 1, size(frac_list(i)%neighbors)
                neigh = frac_list(i)%neighbors(j)
                if (.not.frac_list(neigh)%is_active) cycle
                dx = centroid_x(i) - centroid_x(neigh)
                dy = centroid_y(i) - centroid_y(neigh)
                dz = centroid_z(i) - centroid_z(neigh)
                ds_ij = sqrt(dx*dx + dy*dy + dz*dz)
                if (ds_ij<=0) cycle
                w_avg = 0.5d0*(frac_list(i)%aperture + frac_list(neigh)%aperture)
                Tij   = (w_avg**3)/(12.d0*mu*ds_ij)
                sumTpj = sumTpj + Tij * p_new(neigh)
             end do
             denom = (vol(i)/dt) + Tsum(i)
             if (denom>0) p_new(i) = (rhs(i) + sumTpj)/denom
          end do

          maxErr = maxval(abs(p_new - p_old))
          if (maxErr < tol) exit
          p_old = p_new
       end do   ! Gauss-Seidel loop

       if (inlet_type /= 1) exit   ! Dirichlet case done

       ! Compute flow rate between inlet face and neighbour
       dx_in = centroid_x(1) - centroid_x(neigh_inlet)
       dy_in = centroid_y(1) - centroid_y(neigh_inlet)
       dz_in = centroid_z(1) - centroid_z(neigh_inlet)
       ds_in = sqrt(dx_in*dx_in + dy_in*dy_in + dz_in*dz_in)
       w_avg = 0.5d0*(frac_list(1)%aperture + frac_list(neigh_inlet)%aperture)
       Tij_inlet = (w_avg**3)/(12.d0*mu*ds_in)
       Qcalc = Tij_inlet * (Pin - p_new(neigh_inlet))

       if (abs(Qcalc - inlet_val) < tolQ*max(1.d0, inlet_val)) exit

       Pin = Pin + (inlet_val - Qcalc)/Tij_inlet   ! Newton update
    end do   ! Newton outer

    !----------------------------------------------------------------
    ! 3. Copy back to frac_list and driver array for HF faces
    !----------------------------------------------------------------
    do i = 1, nTot
       if (frac_list(i)%is_active) frac_list(i)%pressure = p_new(i)
    end do

    do i = 1, num_hf
       p_new_hf(i) = frac_list(i)%pressure
    end do

    deallocate(p_old,p_new,rhs,vol,Tsum)
  contains
    pure function centroid_x(idx) result(val)
      integer,intent(in)::idx
      real*8 :: val
      val = sum(frac_list(idx)%x)/4.d0
    end function centroid_x
    pure function centroid_y(idx) result(val)
      integer,intent(in)::idx
      real*8 :: val
      val = sum(frac_list(idx)%y)/4.d0
    end function centroid_y
    pure function centroid_z(idx) result(val)
      integer,intent(in)::idx
      real*8 :: val
      val = sum(frac_list(idx)%z)/4.d0
    end function centroid_z
  end subroutine solve_hf_pressure
end module flow_network 