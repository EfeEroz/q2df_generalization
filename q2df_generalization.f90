program q2df_generalization
    use precision
    implicit none
    
    real(WP) :: Z1_CFD_glob, Z2_CFD_glob, chi11_CFD_glob, chi12_CFD_glob, chi22_CFD_glob, leftover_ox1_CFD_glob, &
        leftover_ox2_CFD_glob, leftover_ox3_CFD_glob
    real(WP), dimension(3) :: OMIX_glob, FMIX_glob
    real(WP) :: chi_xi_glob
    real(WP) :: chi_eta_glob, x_prime_glob, eta_opt_glob, Z_opt_glob
    
    print *, "Starting"
    Z1_CFD_glob = 0.334_WP
    Z2_CFD_glob = 0.129_WP
    chi11_CFD_glob = 3.743_WP
    chi12_CFD_glob = -1.348_WP
    chi22_CFD_glob = 2.1893_WP
    
    leftover_ox1_CFD_glob = -1.0_WP / 92.14_WP * 9.0_WP * 31.999_WP
    leftover_ox2_CFD_glob = 1.0_WP * 0.23292_WP
    leftover_ox3_CFD_glob = -1.0_WP / 100.21_WP * 11.0_WP * 31.999_WP

    call get_optimal_bcs(Z1_CFD_glob, Z2_CFD_glob, chi11_CFD_glob, chi12_CFD_glob, chi22_CFD_glob, &
        leftover_ox1_CFD_glob, leftover_ox2_CFD_glob, leftover_ox3_CFD_glob, OMIX_glob, FMIX_glob, &
        chi_xi_glob, chi_eta_glob, x_prime_glob, Z_opt_glob, eta_opt_glob)

    print *, "Optimal 1D domain (OMIX, FMIX, chi):"
    print *, OMIX_glob
    print *, FMIX_glob
    print *, chi_xi_glob
    print *, Z_opt_glob

    contains

subroutine get_optimal_bcs(Z1_CFD, Z2_CFD, chi11_CFD, chi12_CFD, chi22_CFD, leftover_ox1_CFD, &
    leftover_ox2_CFD, leftover_ox3_CFD, OMIX, FMIX, chi_xi, chi_eta, x_prime, Z_opt, eta_opt)
    ! Here, the real numbers leftover_ox# represent the kg's of oxygen mass leftover if complete combustion is
    ! done on 1 kg of stream #. As expected, negative indicates stream # is fuel-rich.
    use precision
    implicit none

    real(WP), intent(in) :: Z1_CFD, Z2_CFD, chi11_CFD, chi12_CFD, chi22_CFD, leftover_ox1_CFD, &
        leftover_ox2_CFD, leftover_ox3_CFD
    real(WP), dimension(3), intent(out) :: OMIX, FMIX
    real(WP), intent(out) :: chi_xi
    real(WP), intent(out) :: chi_eta, x_prime, Z_opt, eta_opt ! Outputs not necessary for CFD, but instructive to dump

    ! Making arrays of the mixture fractions, dissipation rates, and leftover oxygen to make the data storage
    ! "symmetric" with regard to Z1, Z2, and Z3 rather than just storing Z1, Z2, and their dissipation rates
    real(WP), dimension(3) :: Z
    real(WP) :: chi13_CFD, chi23_CFD, chi33_CFD
    real(WP), dimension(3,3) :: chi
    real(WP), dimension(3) :: leftover_ox_CFD

    real(WP) :: Z1, Z2, chi11, chi12, chi22
    real(WP), dimension(3) :: leftover_ox
    integer :: i, j, i_opt, j_opt
    real(WP) :: x_opt, chi_ratio_opt
    real(WP), allocatable, dimension(:) :: left_endpoints, right_endpoints
    real(WP), allocatable, dimension(:) :: x_extrema
    real(WP), allocatable, dimension(:) :: x_opt_candidates
    integer :: m, n
    real(WP) :: x_opt_candidate, chi_ratio, a, b, leftover_ox_left, leftover_ox_right, TOL
    real(WP), dimension(3) :: OMIX_shuffled, FMIX_shuffled
    integer :: k
    real(WP) :: x_prime_left, x_prime_right
    logical :: complement_xi

    Z = [Z1_CFD, Z2_CFD, 1.0_WP - Z1_CFD - Z2_CFD]

    chi13_CFD = -chi11_CFD - chi12_CFD
    chi23_CFD = -chi22_CFD - chi12_CFD
    chi33_CFD = chi11_CFD + 2.0_WP*chi12_CFD + chi22_CFD
    chi(1,1) = chi11_CFD; chi(1,2) = chi12_CFD; chi(1, 3) = chi13_CFD
    chi(2,1) = chi12_CFD; chi(2,2) = chi22_CFD; chi(2, 3) = chi23_CFD
    chi(3,1) = chi13_CFD; chi(3,2) = chi23_CFD; chi(3, 3) = chi33_CFD
    
    chi_ratio_opt = 100000.0_WP
    leftover_ox_CFD = [leftover_ox1_CFD, leftover_ox2_CFD, leftover_ox3_CFD]
    ! Permutes the stream identities by changing i and j, thereby giving every possible one-dimensional 
    ! domain consideration
    do i = 1, 3
        do j = 1, 3
            if (i /= j) then
                Z1 = Z(i)
                Z2 = Z(j)
                chi11 = chi(i, i)
                chi12 = chi(i, j)
                chi22 = chi(j, j)
                leftover_ox = [leftover_ox_CFD(i), leftover_ox_CFD(j), leftover_ox_CFD(1+2+3-i-j)]
                
                call get_valid_x_ranges(Z1, Z2, leftover_ox, left_endpoints, right_endpoints)
                call get_x_extrema(Z1, Z2, chi11, chi12, chi22, x_extrema)

                allocate(x_opt_candidates(0))
                
                do m = 1, size(x_extrema)
                    do n = 1, size(left_endpoints)
                        if (left_endpoints(n) < x_extrema(m) .and. x_extrema(m) < right_endpoints(n)) then
                            call append(x_opt_candidates, x_extrema(m))
                        end if
                    end do
                end do

                deallocate(x_extrema)

                do n = 1, size(left_endpoints)
                    call append(x_opt_candidates, left_endpoints(n))
                    call append(x_opt_candidates, right_endpoints(n))
                end do

                deallocate(left_endpoints)
                deallocate(right_endpoints)
                
                do k = 1, size(x_opt_candidates)
                    x_opt_candidate = x_opt_candidates(k)
                    chi_ratio = chi_ratio_val(Z1, Z2, chi11, chi12, chi22, x_opt_candidate)
                    if (chi_ratio < chi_ratio_opt) then
                        i_opt = i
                        j_opt = j
                        x_opt = x_opt_candidate
                        chi_ratio_opt = chi_ratio
                    end if
                end do

                deallocate(x_opt_candidates)
            end if
        end do
    end do

    ! At this point, the locally-optimal boundary conditions are determined by i_opt, j_opt, and x_opt
    Z1 = Z(i_opt)
    Z2 = Z(j_opt)
    chi11 = chi(i_opt, i_opt)
    chi12 = chi(i_opt, j_opt)
    chi22 = chi(j_opt, j_opt)
    leftover_ox = [leftover_ox_CFD(i_opt), leftover_ox_CFD(j_opt), leftover_ox_CFD(1+2+3-i_opt-j_opt)]

    b = b_val(Z1, Z2, x_opt)
    a = x_opt * b
    leftover_ox_left = a*leftover_ox(1) + (1-a)*leftover_ox(3)
    leftover_ox_right = b*leftover_ox(1) + (1-b)*leftover_ox(2)

    TOL = 0.00000001_WP
    ! Assigns OMIX and FMIX by determining which endpoint is fuel-lean or fuel-rich
    if (leftover_ox_left >= -TOL .and. leftover_ox_right <= TOL) then
        OMIX_shuffled = [a, 0.0_WP, 1.0_WP-a]
        FMIX_shuffled = [b, 1.0_WP-b, 0.0_WP]
        complement_xi = .false.
    else if (leftover_ox_left <= TOL .and. leftover_ox_right >= -TOL) then
        OMIX_shuffled = [b, 1.0_WP-b, 0.0_WP]
        FMIX_shuffled = [a, 0.0_WP, 1.0_WP-a]
        complement_xi = .true.
    else
        print *, "Both boundary conditions are fuel-lean or both are fuel-rich, so can't assign OMIX and FMIX identities."
    end if

    OMIX(i_opt) = OMIX_shuffled(1)
    FMIX(i_opt) = FMIX_shuffled(1)
    OMIX(j_opt) = OMIX_shuffled(2)
    FMIX(j_opt) = FMIX_shuffled(2)
    OMIX(1+2+3-i_opt-j_opt) = OMIX_shuffled(3)
    FMIX(1+2+3-i_opt-j_opt) = FMIX_shuffled(3)
    chi_xi = chi_xi_val(Z1, Z2, chi11, chi12, chi22, x_opt)

    chi_eta = chi_eta_val(Z1, Z2, chi11, chi12, chi22, x_opt)
    eta_opt = eta_val(Z1, Z2, x_opt)
    Z_opt = xi_val(Z1, Z2, x_opt)
    if (complement_xi) then
        Z_opt = 1 - Z_opt
    end if

    if (i_opt == 1 .and. j_opt == 2) then
        x_prime_left = 3
        x_prime_right = 4
    else if (i_opt == 1 .and. j_opt == 3) then
        x_prime_left = 5
        x_prime_right = 4
    else if (i_opt == 2 .and. j_opt == 1) then
        x_prime_left = 3
        x_prime_right = 2
    else if (i_opt == 2 .and. j_opt == 3) then
        x_prime_left = 1
        x_prime_right = 2
    else if (i_opt == 3 .and. j_opt == 1) then
        x_prime_left = 5
        x_prime_right = 6
    else if (i_opt == 3 .and. j_opt == 2) then
        x_prime_left = 1
        x_prime_right = 0
    else
        print *, "ERROR, stream permutation invalid."
    end if
    x_prime = x_prime_left + x_opt * (x_prime_right - x_prime_left)

end subroutine get_optimal_bcs

subroutine get_valid_x_ranges(Z1, Z2, leftover_ox, left_endpoints, right_endpoints)
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2
    real(WP), dimension(3), intent(in) :: leftover_ox
    real(WP), dimension(:), allocatable, intent(out) :: left_endpoints, right_endpoints
    real(WP) :: stream1_stream3_a, stream1_stream3_b, stream1_stream2_a, stream1_stream2_b
    real(WP), allocatable, dimension(:) :: x_boundaries
    real(WP) :: temp
    integer :: i, j
    logical :: x1_valid
    real(WP) :: x1_left_leftover_ox, x1_right_leftover_ox
    real(WP), allocatable, dimension(:) :: x_boundaries_valid
    allocate(left_endpoints(0))
    allocate(right_endpoints(0))

    ! stream1_stream3 variables denote constant-eta line intersecting stoichiometric stream1-stream3 mixture if it exists
    ! stream1_stream2 variables denote constant-eta line intersecting stoichiometric stream1-stream2 mixture if it exists
    stream1_stream3_a = leftover_ox(3) / (leftover_ox(3) - leftover_ox(1))
    stream1_stream2_b = leftover_ox(2) / (leftover_ox(2) - leftover_ox(1))
    stream1_stream3_b = 1.0_WP + Z2 / ((1.0_WP-Z1-Z2)/(1.0_WP-stream1_stream3_a) - 1.0_WP)
    stream1_stream2_a = 1.0_WP - (1.0_WP-Z1-Z2)/(Z2/(stream1_stream2_b-1.0_WP) + 1.0_WP)

    allocate(x_boundaries(2))
    x_boundaries = [0.0_WP, 1.0_WP]
    if (0.0_WP <= stream1_stream3_a .and. stream1_stream3_a <= stream1_stream3_b .and. stream1_stream3_b <= 1.0_WP) then
        call append(x_boundaries, stream1_stream3_a / stream1_stream3_b)
    end if
    if (0.0_WP <= stream1_stream2_a .and. stream1_stream2_a <= stream1_stream2_b .and. stream1_stream2_b <= 1.0_WP) then
        call append(x_boundaries, stream1_stream2_a / stream1_stream2_b)
    end if

    do i = 1, size(x_boundaries) - 1
        do j = 1, size(x_boundaries) - 1
            if (x_boundaries(j) < x_boundaries(j+1)) then
                temp = x_boundaries(j)
                x_boundaries(j) = x_boundaries(j + 1)
                x_boundaries(j + 1) = temp
            end if
        end do
    end do

    x1_left_leftover_ox = leftover_ox(1) * Z1 + leftover_ox(3) * (1.0_WP - Z1)
    x1_right_leftover_ox = leftover_ox(1) * Z1 + leftover_ox(2) * (1.0_WP - Z1)
    x1_valid = (x1_left_leftover_ox >= 0.0_WP .and. x1_right_leftover_ox <= 0.0_WP) .or. &
        (x1_left_leftover_ox <= 0.0_WP .and. x1_right_leftover_ox >= 0.0_WP)
    allocate(x_boundaries_valid(1))
    if (x1_valid) then
        x_boundaries_valid = [1.0_WP]
    else
        x_boundaries_valid = [0.0_WP]
    end if

    do i = 1, size(x_boundaries) - 2
        ! The -x + 1.0 is equivalent to switching between 0.0 (False) and 1.0 (True) and was done to circumvent
        ! the need to implement resizable arrays for logical types
        call append(x_boundaries_valid, -x_boundaries_valid(size(x_boundaries_valid)) + 1.0_WP)
    end do

    do i = 1, size(x_boundaries_valid)
        if (x_boundaries_valid(i) > 0.5_WP) then   ! Again, using 1.0 for True and 0.0 for False
            call append(left_endpoints, x_boundaries(i + 1))
            call append(right_endpoints, x_boundaries(i))
        end if
    end do
    deallocate(x_boundaries)
    deallocate(x_boundaries_valid)

end subroutine get_valid_x_ranges

subroutine get_x_extrema(Z1, Z2, chi11, chi12, chi22, x_extrema)
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, chi11, chi12, chi22
    real(WP), allocatable, dimension(:), intent(out) :: x_extrema
    integer, parameter :: GRID = 1000
    real(WP), dimension(GRID) :: dchiratio_dx_arr
    real(WP) :: x, step_size
    integer :: i

    x = 0.000001_WP
    step_size = (1.0_WP - 2.0_WP*x) / (GRID - 1.0_WP)
    do i = 1, GRID
        dchiratio_dx_arr(i) = dchiratio_dx_val(Z1, Z2, chi11, chi12, chi22, x + step_size * (i - 1))
    end do

    allocate(x_extrema(0))
    do i = 1, GRID
        if (dchiratio_dx_arr(i) * dchiratio_dx_arr(i + 1) <= 0) then
            call append(x_extrema, x + step_size * (i - 0.5_WP))
        end if
    end do
end subroutine get_x_extrema

real(WP) function chi_ratio_val(Z1, Z2, chi11, chi12, chi22, x)
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, chi11, chi12, chi22, x
    real(WP) :: chi_eta, chi_xi

    chi_eta = chi_eta_val(Z1, Z2, chi11, chi12, chi22, x)
    chi_xi = chi_xi_val(Z1, Z2, chi11, chi12, chi22, x)
    chi_ratio_val = chi_eta / chi_xi
end function chi_ratio_val

real(WP) function chi_eta_val(Z1, Z2, chi11, chi12, chi22, x)
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, chi11, chi12, chi22, x
    real(WP) :: eta, xi
    eta = eta_val(Z1, Z2, x)
    xi = xi_val(Z1, Z2, x)

    chi_eta_val = ((1.0_WP-eta)**2.0_WP*chi11 - 2.0_WP*eta*(1.0_WP-eta)*(1.0_WP-x)*chi12 + (eta*(1.0_WP-x))**2.0_WP*chi22) / &
        (((1.0_WP-eta)*(x+(1.0_WP-x)*xi) + eta*xi*(1.0_WP-x))**2.0_WP)
end function chi_eta_val

real(WP) function chi_xi_val(Z1, Z2, chi11, chi12, chi22, x)
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, chi11, chi12, chi22, x
    real(WP) :: eta, xi
    eta = eta_val(Z1, Z2, x)
    xi = xi_val(Z1, Z2, x)

    chi_xi_val = (xi**2.0_WP*chi11 + 2.0_WP*xi*(x+(1.0_WP-x)*xi)*chi12 + (x+(1.0_WP-x)*xi)**2.0_WP*chi22) / &
        (((1.0_WP-eta)*(x+(1.0_WP-x)*xi) + eta*xi*(1.0_WP-x))**2.0_WP)
end function chi_xi_val

real(WP) function b_val(Z1, Z2, x)
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, x
    real(WP) :: Z3
    Z3 = 1.0_WP - Z1 - Z2
    b_val = (Z1+Z2+x*(Z1+Z3) - sqrt((Z1+Z2+x*(Z1+Z3))**2.0_WP - 4.0_WP*Z1*x)) / (2.0_WP*x)
end function b_val

real(WP) function eta_val(Z1, Z2, x)
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, x
    real(WP) :: b
    b = b_val(Z1, Z2, x)
    eta_val = b
end function eta_val

real(WP) function xi_val(Z1, Z2, x)
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, x
    real(WP) :: b
    b = b_val(Z1, Z2, x)
    xi_val = Z2 / (1.0_WP - b)
end function xi_val

real(WP) function db_dx_val(Z1, Z2, x)
    ! At the point (Z1, Z2), find the derivative of b with respect to x
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, x
    real(WP) :: Z3, b
    Z3 = 1.0_WP - Z1 - Z2
    b = b_val(Z1, Z2, x)
    
    db_dx_val = 1.0_WP/(4.0_WP*x**2.0_WP) * ((Z1+Z3-0.5_WP/sqrt((Z1+Z2+x*(Z1+Z3))**2.0_WP-4.0_WP*Z1*x)*(2.0_WP*(Z1+Z2+x*(Z1+ &
        Z3))*(Z1+Z3)-4.0_WP*Z1)) * 2.0_WP*x - 2.0_WP * (Z1+Z2+x*(Z1+Z3) - sqrt((Z1+Z2+x*(Z1+Z3))**2.0_WP - 4.0_WP*Z1*x)))
end function db_dx_val

real(WP) function deta_dx_val(Z1, Z2, x)
    ! At the point (Z1, Z2), find the derivative of eta with respect to x
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, x
    real(WP) :: db_dx
    db_dx = db_dx_val(Z1, Z2, x)

    deta_dx_val = db_dx
end function deta_dx_val

real(WP) function dxi_dx_val(Z1, Z2, x)
    ! At the point (Z1, Z2), find the derivative of xi with respect to x
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, x
    real(WP) :: b, db_dx
    b = b_val(Z1, Z2, x)
    db_dx = db_dx_val(Z1, Z2, x)

    dxi_dx_val = Z2/((1.0_WP-b)**2.0_WP) * db_dx
end function dxi_dx_val

real(WP) function dchiratio_dx_val(Z1, Z2, chi11, chi12, chi22, x)
    ! At the point (Z1, Z2), find the derivative of chi_eta / chi_xi with respect to x
    use precision
    implicit none

    real(WP), intent(in) :: Z1, Z2, chi11, chi12, chi22, x
    real(WP) :: eta, xi, deta, dxi
    eta = eta_val(Z1, Z2, x)
    xi = xi_val(Z1, Z2, x)
    deta = deta_dx_val(Z1, Z2, x)
    dxi = dxi_dx_val(Z1, Z2, x)

    dchiratio_dx_val = (-2.0_WP*(1.0_WP-eta)*deta*chi11 + 2.0_WP*deta*(2.0_WP*eta-1.0_WP)*(1.0_WP-x)*chi12 + 2.0_WP*eta* &
        (1.0_WP-eta)*chi12 + 2.0_WP*eta*(1.0_WP-x)*(deta*(1.0_WP-x)-eta)*chi22) * (xi**2.0_WP*chi11 + 2.0_WP*xi*(x+ &
        (1.0_WP-x)*xi)*chi12 + (x+(1.0_WP-x)*xi)**2.0_WP*chi22) - (2.0_WP*xi*dxi*chi11 + 2.0_WP*dxi*(x+(1.0_WP-x)*xi)* &
        chi12 + 2.0_WP*xi*(1.0_WP+(1.0_WP-x)*dxi-xi)*chi12 + 2.0_WP*(x+(1.0_WP-x)*xi)*(1.0_WP+(1.0_WP-x)*dxi-xi)*chi22) * &
        ((1.0_WP-eta)**2.0_WP*chi11 - 2.0_WP*eta*(1.0_WP-eta)*(1.0_WP-x)*chi12 + (eta*(1.0_WP-x))**2.0_WP*chi22)
end function dchiratio_dx_val

subroutine append(arr, val)
    use precision
    implicit none

    real(WP), allocatable, dimension(:), intent(inout) :: arr
    real(WP), intent(in) :: val
    real(WP), allocatable, dimension(:) :: temp
    
    allocate(temp(size(arr) + 1))
    if (size(arr) > 0) then
        temp(1:size(arr)) = arr
    end if
    temp(size(arr) + 1) = val
    call move_alloc(temp, arr)
end subroutine append

end program q2df_generalization