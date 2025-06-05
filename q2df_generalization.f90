program q2df_generalization
    use precision
    implicit none
    
    call test_get_optimal_bcs()

    contains

! IMPORTANT: Do not forget to use the mkl library erfc^-1 function to find chi_ZZ at reference mixture fraction
! given that it is chi_xi at Z_opt. Z_stoic is also output by the program for convenience.
subroutine get_optimal_bcs(Z1_CFD, Z2_CFD, chi11_CFD, chi12_CFD, chi22_CFD, leftover_ox1_CFD, &
    leftover_ox2_CFD, leftover_ox3_CFD, OMIX, FMIX, chi_xi, chi_eta, x_prime, Z_opt, eta_opt, Z_stoic)
    ! Here, the real numbers leftover_ox# represent the kg's of oxygen mass leftover if complete combustion is
    ! done on 1 kg of stream #. As expected, negative indicates stream # is fuel-rich.
    use precision
    use, intrinsic :: ieee_arithmetic
    implicit none

    real(WP), intent(in) :: Z1_CFD, Z2_CFD, chi11_CFD, chi12_CFD, chi22_CFD, leftover_ox1_CFD, &
        leftover_ox2_CFD, leftover_ox3_CFD
    real(WP), dimension(3), intent(out) :: OMIX, FMIX
    real(WP), intent(out) :: chi_xi
    real(WP), intent(out) :: chi_eta, x_prime, Z_opt, eta_opt ! Outputs not necessary for CFD, but instructive to dump
    real(WP), intent(out) :: Z_stoic

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

    real(WP) :: q2df2_chi_ratio, q2df2_xi, q2df2_eta, q2df2_chi_ratio_check, generalized_chi_ratio

    Z = [Z1_CFD, Z2_CFD, 1.0_WP - Z1_CFD - Z2_CFD]

    chi13_CFD = -chi11_CFD - chi12_CFD
    chi23_CFD = -chi22_CFD - chi12_CFD
    chi33_CFD = chi11_CFD + 2.0_WP*chi12_CFD + chi22_CFD
    chi(1,1) = chi11_CFD; chi(1,2) = chi12_CFD; chi(1, 3) = chi13_CFD
    chi(2,1) = chi12_CFD; chi(2,2) = chi22_CFD; chi(2, 3) = chi23_CFD
    chi(3,1) = chi13_CFD; chi(3,2) = chi23_CFD; chi(3, 3) = chi33_CFD
    
    chi_ratio_opt = huge(1.0_WP)
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

                    if (ieee_is_nan(chi_ratio)) then
                        x_opt_candidate = abs(x_opt_candidate - 0.0000000001_WP)
                    end if
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
    Z_stoic = leftover_ox_left / (leftover_ox_left - leftover_ox_right)
    if (complement_xi) then
        Z_opt = 1 - Z_opt
        Z_stoic = 1 - Z_stoic
    end if

    ! Add testing with Q2DF2 (make sure chi ratio is computed correctly for Q2DF2 and that the generalized choice
    ! is at least as good as Q2DF2, as quantified by psi (chi ratio)
    q2df2_chi_ratio = chi_ratio_val(Z(1), Z(3), chi(1, 1), chi(1, 3), chi(3, 3), 0.0000000001_WP)
    q2df2_xi = Z(1) + Z(3)
    q2df2_eta = Z(3) / (Z(1) + Z(3))
    q2df2_chi_ratio_check = (chi(1, 1) + 2.0_WP*(1.0_WP-q2df2_eta)*chi(1, 2) + (1.0_WP-q2df2_eta)**2.0_WP*chi(2, 2)) / &
        (q2df2_xi**2.0_WP * chi(2, 2))
    if (abs(q2df2_chi_ratio - q2df2_chi_ratio_check) > 0.001_WP) then
        print *, "ALERT", q2df2_chi_ratio, q2df2_chi_ratio_check
    end if

    ! Continued testing:
    generalized_chi_ratio = chi_eta / chi_xi
    if (ieee_is_nan(q2df2_chi_ratio) .or. ieee_is_nan(generalized_chi_ratio)) then
        print *, "NaN ISSUE:", q2df2_chi_ratio, generalized_chi_ratio
    end if
    if (q2df2_chi_ratio < generalized_chi_ratio) then
        print *, ">>>>   ISSUE", q2df2_chi_ratio, "<", generalized_chi_ratio

        print *, "Optimal i, j, and x", i_opt, j_opt, x_opt
        Z1 = Z(i_opt)
        Z2 = Z(j_opt)
        leftover_ox = [leftover_ox_CFD(i_opt), leftover_ox_CFD(j_opt), leftover_ox_CFD(1+2+3-i_opt-j_opt)]
        call get_valid_x_ranges(Z1, Z2, leftover_ox, left_endpoints, right_endpoints)
        print *, "Optimal left endpoints", left_endpoints
        print *, "Optimal right endpoints", right_endpoints
        deallocate(left_endpoints)
        deallocate(right_endpoints)

        call get_valid_x_ranges(Z(3), Z(1), [leftover_ox_CFD(3), leftover_ox_CFD(1), leftover_ox_CFD(2)], &
            left_endpoints, right_endpoints)
        print *, "Q2DF2-1 left endpoints", left_endpoints
        print *, "Q2DF2-1 right endpoints", right_endpoints
        deallocate(left_endpoints)
        deallocate(right_endpoints)

        call get_valid_x_ranges(Z(1), Z(3), [leftover_ox_CFD(1), leftover_ox_CFD(3), leftover_ox_CFD(2)], &
            left_endpoints, right_endpoints)
        print *, "Q2DF2-2 left endpoints", left_endpoints
        print *, "Q2DF2-2 right endpoints", right_endpoints
        deallocate(left_endpoints)
        deallocate(right_endpoints)
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
    do i = 1, GRID - 1
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

subroutine test_get_optimal_bcs()
    use precision
    implicit none
    
    integer, parameter :: N = 150
    integer :: i
    real(WP) :: leftover_ox1_CFD, leftover_ox2_CFD, leftover_ox3_CFD
    real(WP), dimension(N) :: Z1_CFD_list, Z2_CFD_list, chi11_CFD_list, chi12_CFD_list, chi22_CFD_list
    real(WP), dimension(N, 3) :: OMIX_truth, FMIX_truth
    real(WP), dimension(N) :: chi_xi_truth
    real(WP), dimension(N) :: chi_eta_truth, x_prime_truth, eta_opt_truth, Z_opt_truth, Z_stoic_truth

    real(WP), dimension(3) :: OMIX, FMIX
    real(WP) :: chi_xi
    real(WP) :: chi_eta, x_prime, eta_opt, Z_opt, Z_stoic

    integer :: num_failures
    num_failures = 0

    leftover_ox1_CFD = -1.0_WP / 92.14_WP * 9.0_WP * 31.999_WP
    leftover_ox2_CFD = 1.0_WP * 0.23292_WP
    leftover_ox3_CFD = -1.0_WP / 100.21_WP * 11.0_WP * 31.999_WP

    Z1_CFD_list = [0.2199014833190605_WP, 0.33724889867020524_WP, &
        0.4986614756676147_WP, 0.5473975725738289_WP, 0.6518446856741693_WP, &
        0.10011579482263833_WP, 0.2540913890650932_WP, 0.22912771734103338_WP, &
        0.5532287089156661_WP, 0.42617771347983435_WP, 0.05017197579677181_WP, &
        0.39596483660459614_WP, 0.46850268483227786_WP, 0.21152852794938165_WP, &
        0.31978147817863894_WP, 0.08521825767419386_WP, 0.05860632441878657_WP, &
        0.26087143946165775_WP, 0.6852275764235453_WP, 0.3353462937557347_WP, &
        0.5603224618444705_WP, 0.28192800758686465_WP, 0.15240056342074515_WP, &
        0.38160982109391983_WP, 0.15648022727515504_WP, 0.4127767238789505_WP, &
        0.19185443818191086_WP, 0.7227385333773537_WP, 0.08481392176514818_WP, &
        0.3557030916482899_WP, 0.44986806407635177_WP, 0.5291482372882572_WP, &
        0.18575446369421877_WP, 0.410379911989409_WP, 0.4262601527358208_WP, &
        0.10301395132197028_WP, 0.26732406560547023_WP, 0.015577285366665146_WP, &
        0.440566052479085_WP, 0.2182344295966475_WP, 0.1402416982133526_WP, &
        0.5785432965280617_WP, 0.3403689992201045_WP, 0.6238387205390453_WP, &
        0.3248215986019519_WP, 0.23609184978666217_WP, 0.5447778257083787_WP, &
        0.37332424851416274_WP, 0.3790366126450805_WP, 0.2984766060199765_WP, &
        0.09141904812902701_WP, 0.5794287269844878_WP, 0.4702399633396789_WP, &
        0.38082351796535896_WP, 0.366452490710751_WP, 0.41042525537703367_WP, &
        0.32471802425049406_WP, 0.39981158709912334_WP, 0.6527941297948753_WP, &
        0.4043399535575942_WP, 0.3737644840093439_WP, 0.33539588601218445_WP, &
        0.4693843619898315_WP, 0.3417361669300392_WP, 0.42444859120006884_WP, &
        0.5936092401607916_WP, 0.4046750123424271_WP, 0.3995551672280777_WP, &
        0.41840635583207475_WP, 0.3207757646093318_WP, 0.4645072359392396_WP, &
        0.23103225424586576_WP, 0.14031228757077374_WP, 0.5679987279397191_WP, &
        0.042928147514109914_WP, 0.08382232557046158_WP, 0.6448766846548963_WP, &
        0.28958802854562904_WP, 0.2050861958148722_WP, 0.10083073983905352_WP, &
        0.28072608914963615_WP, 0.3704117950187408_WP, 0.27863002074775045_WP, &
        0.586435512639438_WP, 0.8626431115515417_WP, 0.207390174164105_WP, &
        0.3232471329590064_WP, 0.2953569342936011_WP, 0.6167907921731448_WP, &
        0.517376967819739_WP, 0.5638173904655239_WP, 0.13958116275434468_WP, &
        0.18671133826199568_WP, 0.782938069412775_WP, 0.3925356832399642_WP, &
        0.47031330206222344_WP, 0.8011651267471145_WP, 0.3060994102202625_WP, &
        0.00399882969264699_WP, 0.5222969889557657_WP, 0.10830361514460254_WP, &
        0.5023802478647574_WP, 0.3042369767885765_WP, 0.3466747177091465_WP, &
        0.35820610982018725_WP, 0.4936346904206096_WP, 0.3809677511562454_WP, &
        0.33348140630186907_WP, 0.332920218958037_WP, 0.7100720807137153_WP, &
        0.14749782530644379_WP, 0.6786469761448588_WP, 0.41175940237229086_WP, &
        0.37900417894494526_WP, 0.21329443289138078_WP, 0.39248584947045156_WP, &
        0.39643645494672664_WP, 0.5613989593303148_WP, 0.2125595542874605_WP, &
        0.48534364232014165_WP, 0.07030791070963247_WP, 0.24845062668257015_WP, &
        0.5388951988688745_WP, 0.36928628524833174_WP, 0.42107987829336113_WP, &
        0.13772788309397188_WP, 0.5818216418492889_WP, 0.5031793325623714_WP, &
        0.16871974416159205_WP, 0.5877516691573587_WP, 0.7595470831297896_WP, &
        0.3320904697213568_WP, 0.1352489077352564_WP, 0.3322297593342785_WP, &
        0.4076134570189294_WP, 0.23770575038778705_WP, 0.31066432093050494_WP, &
        0.5380878837333941_WP, 0.3459591689840704_WP, 0.2796180286556356_WP, &
        0.4083006099648495_WP, 0.37744115692403224_WP, 0.2529412733776454_WP, &
        0.1921897322956072_WP, 0.5034849512181531_WP, 0.10785131678733441_WP, &
        0.3830587212862429_WP, 0.46171878326383614_WP, 0.2_WP, 0.2_WP]
    Z2_CFD_list = [0.3705680937512448_WP, 0.34392899696950824_WP, &
        0.32705817515950003_WP, 0.31998335176127896_WP, 0.046115196195361974_WP, &
        0.4497745859137792_WP, 0.4504820151801578_WP, 0.31941769907896295_WP, &
        0.0814056258180373_WP, 0.18621384531728674_WP, 0.4739657505431522_WP, &
        0.3064783260060796_WP, 0.38796853211730536_WP, 0.24532477722231166_WP, &
        0.24313528433168716_WP, 0.6421625137595265_WP, 0.3626880748917589_WP, &
        0.2717566926176325_WP, 0.048905252402952896_WP, 0.4642142604341024_WP, &
        0.05617292409388961_WP, 0.13177021107482748_WP, 0.10494491750435105_WP, &
        0.3870044922330624_WP, 0.46322163807163014_WP, 0.04916554287202971_WP, &
        0.4114918693296534_WP, 0.1398867526165806_WP, 0.6609007393246171_WP, &
        0.36441405001385996_WP, 0.018435103788960176_WP, 0.28617649604162654_WP, &
        0.4142045497766184_WP, 0.3868160343638349_WP, 0.06911246898010968_WP, &
        0.6207405109209677_WP, 0.2826328581986067_WP, 0.5664602650861379_WP, &
        0.2819561026635625_WP, 0.5350636983204851_WP, 0.7730029852720429_WP, &
        0.26556542417113727_WP, 0.36132427507712656_WP, 0.2685228765249786_WP, &
        0.41756888548365234_WP, 0.38940371020191705_WP, 0.4517058857701718_WP, &
        0.2925467696602508_WP, 0.5158959100072742_WP, 0.2064114849399907_WP, &
        0.784052168928812_WP, 0.1856240320522414_WP, 0.37477895221933083_WP, &
        0.34126810812432_WP, 0.27506489629970954_WP, 0.4338615394967612_WP, &
        0.3334863868284931_WP, 0.3989873628897666_WP, 0.17275641071008677_WP, &
        0.18672046802426834_WP, 0.17866973403018435_WP, 0.24402278507822608_WP, &
        0.23494653535259788_WP, 0.3088041275444134_WP, 0.2463194502096872_WP, &
        0.04882008621365231_WP, 0.43651196458612196_WP, 0.0895003733545358_WP, &
        0.37416858022650495_WP, 0.47002297238415236_WP, 0.10479788953145505_WP, &
        0.4031677857690975_WP, 0.28612020597500043_WP, 0.008492308474438865_WP, &
        0.5518617637366549_WP, 0.4065115450893835_WP, 0.1611668592886759_WP, &
        0.3201205215964496_WP, 0.6251392296144408_WP, 0.2016286212834643_WP, &
        0.35270984329647836_WP, 0.2900585557345272_WP, 0.462715639482829_WP, &
        0.13559452113671347_WP, 0.0895590273985805_WP, 0.5911628581781774_WP, &
        0.29806084528983595_WP, 0.3659500303594625_WP, 0.18463588635622183_WP, &
        0.07541279585211762_WP, 0.37395397940961217_WP, 0.42719626898798646_WP, &
        0.4875616991428157_WP, 0.046855476157732195_WP, 0.17739139049486743_WP, &
        0.4332643553915542_WP, 0.08214465617291684_WP, 0.48798773495581793_WP, &
        0.9276180424799187_WP, 0.46793765336655935_WP, 0.46948012081118123_WP, &
        0.13873536490036467_WP, 0.4236691547151846_WP, 0.06622778204991378_WP, &
        0.3426569173980422_WP, 0.283108888342331_WP, 0.3939151715604325_WP, &
        0.09586823519507531_WP, 0.39693912354456784_WP, 0.12999319943765003_WP, &
        0.10066456031830967_WP, 0.2814714237195239_WP, 0.5180677339749454_WP, &
        0.25726316947696887_WP, 0.3824543742217693_WP, 0.42488768478570216_WP, &
        0.20740996824395966_WP, 0.2339370305712483_WP, 0.38607983379980637_WP, &
        0.41240799076257556_WP, 0.4345804727585336_WP, 0.4543281382020148_WP, &
        0.4445272766526614_WP, 0.2074102748829358_WP, 0.284625854669468_WP, &
        0.3007142481803916_WP, 0.37814076642593225_WP, 0.10710190737014635_WP, &
        0.5117412450144577_WP, 0.35927797873800077_WP, 0.10434478839910205_WP, &
        0.3981091361962083_WP, 0.28573847311088313_WP, 0.4703936787525612_WP, &
        0.5422076537450176_WP, 0.34866512803915783_WP, 0.2847647685621477_WP, &
        0.21141064098174234_WP, 0.6385804230981139_WP, 0.2794595354383969_WP, &
        0.18056745029495191_WP, 0.3035192193450436_WP, 0.2240367308446115_WP, &
        0.4511447049753618_WP, 0.3882567126350668_WP, 0.23555933663973871_WP, &
        0.49771856107555434_WP, 0.46056076397270657_WP, 0.4_WP, 0.4_WP]
    chi11_CFD_list = [0.37448988013761686_WP, 1.2080334384442448_WP, &
        0.5504055865908986_WP, 1.0925378386908877_WP, 1.0989881299301363_WP, &
        1.9749569160344527_WP, 0.5024357055548725_WP, 0.6823016980747854_WP, &
        0.0273243485975192_WP, 1.8482880046591774_WP, 0.8464867171963026_WP, &
        0.729763741342554_WP, 0.32913223002243175_WP, 1.5601136768866235_WP, &
        0.3967115478546368_WP, 0.6252017135087902_WP, 0.46485487718235663_WP, &
        1.5170183562340402_WP, 0.8723349175607615_WP, 1.8883864428777344_WP, &
        0.1672949127537655_WP, 0.834374688562945_WP, 1.410442075366921_WP, &
        0.25396961381471295_WP, 1.3932056861325763_WP, 0.40150820183539215_WP, &
        1.068154220090624_WP, 1.7562223826988685_WP, 1.2726279788886494_WP, &
        0.35558205073943383_WP, 1.5897733657280848_WP, 0.8628043291847629_WP, &
        0.09657572886132093_WP, 0.06437340550152304_WP, 0.9155817521092293_WP, &
        0.7674240652837501_WP, 1.8635666321390718_WP, 0.5247805688594531_WP, &
        1.3702097619538702_WP, 1.6981280138795412_WP, 0.6404995369312549_WP, &
        0.701071945978126_WP, 1.5119059331790106_WP, 0.5275943002401089_WP, &
        0.7180401350831824_WP, 1.1904368800158878_WP, 1.7647003400074834_WP, &
        0.2925069867300407_WP, 1.4652238128030857_WP, 0.4076734744295445_WP, &
        1.2258662323573903_WP, 0.7350400795391094_WP, 0.07141439700788532_WP, &
        0.31286561789027023_WP, 1.9125013592185598_WP, 1.4058174561693906_WP, &
        0.44441661252136155_WP, 0.5222624226660373_WP, 1.3951101074366286_WP, &
        1.6447689543571302_WP, 1.4554268516222_WP, 0.3825984765893611_WP, &
        1.7627177830302416_WP, 0.4217835541174666_WP, 1.053895106439828_WP, &
        1.7726365914389606_WP, 1.3237436348992115_WP, 1.3051290571271017_WP, &
        0.8602153714801837_WP, 0.5645751340876433_WP, 0.20816495254153233_WP, &
        1.9213858921410247_WP, 1.0981524961503244_WP, 0.36669991435611315_WP, &
        1.6823613653769602_WP, 0.07059586738277801_WP, 0.956955922331739_WP, &
        0.30525397377881736_WP, 1.5117377454715273_WP, 1.2623391616877935_WP, &
        1.3278824921452705_WP, 1.1992303352251623_WP, 0.8635962753437723_WP, &
        1.793024279328282_WP, 0.712881958640577_WP, 0.24396886160691245_WP, &
        0.6068370527734206_WP, 0.6047549920909814_WP, 1.5294785365014805_WP, &
        1.9409091944775263_WP, 0.89826244736222_WP, 1.6888735150082952_WP, &
        1.7108055505649407_WP, 0.589531422394793_WP, 1.8716786147932314_WP, &
        0.8664413120458836_WP, 1.1462600779351972_WP, 1.2102165492924912_WP, &
        1.991649132992543_WP, 1.9282751744859152_WP, 1.5336567625132957_WP, &
        1.7614800730224924_WP, 0.16665827603304684_WP, 1.9874975043479999_WP, &
        0.4639688996857916_WP, 0.21090616711623134_WP, 0.4029628949014801_WP, &
        1.939942672444102_WP, 1.6330542930873646_WP, 0.32546523127861837_WP, &
        1.2635767067696342_WP, 1.5662429024762834_WP, 0.05282729445019818_WP, &
        1.2975436691021063_WP, 0.943665859368576_WP, 1.7367828948941442_WP, &
        0.878955812488238_WP, 0.7918423109744743_WP, 1.6322182974200652_WP, &
        0.19944507611308904_WP, 1.988104889740077_WP, 1.9224147089722698_WP, &
        1.4315311807753552_WP, 0.42573785209207826_WP, 0.7158613273175176_WP, &
        0.21790630516364184_WP, 0.821630496987328_WP, 0.9291607418439691_WP, &
        1.4396935838593845_WP, 1.2930480756487517_WP, 1.1252820095911835_WP, &
        0.045619568993587656_WP, 0.8370699005664852_WP, 1.5899711517601285_WP, &
        1.2687672851711982_WP, 1.0836998083729703_WP, 0.2689098074149967_WP, &
        0.14144037324212322_WP, 1.7026429054443015_WP, 1.1162043416209115_WP, &
        0.08768801746026589_WP, 1.916932459697617_WP, 1.3293598094721137_WP, &
        0.32691044324662166_WP, 1.2176379983119752_WP, 1.5441898518498265_WP, &
        0.37831050354988904_WP, 0.2685390249994337_WP, 1.3_WP, 2.6_WP]
    chi12_CFD_list = [-0.5320037866604121_WP, 0.23563586287195093_WP, &
        0.2945280199871081_WP, 0.478713175030655_WP, -0.6923906372141617_WP, &
        0.2676786392517665_WP, -0.22443454591710185_WP, 0.24440230493603893_WP, &
        -0.09347096424153109_WP, 0.8661307950603732_WP, 0.6624171190445489_WP, &
        -0.04887573540833268_WP, -0.057337822486441814_WP, 0.30925022269070035_WP, &
        -0.1604091744958509_WP, 0.4120984945549112_WP, -0.28127129242758414_WP, &
        1.0458068601154276_WP, 0.35438080247724835_WP, -0.8725486001797165_WP, &
        0.20552544948667334_WP, -0.994208845838252_WP, -1.3753994835604468_WP, &
        -0.13042237092158096_WP, 0.007107620815502758_WP, -0.1224633564874176_WP, &
        0.440479577838399_WP, 1.120208378032193_WP, -0.80663702806795_WP, &
        0.4387220050661861_WP, 0.42374144689981924_WP, -0.30303341667467576_WP, &
        -0.2509403468688064_WP, -0.13207927714366266_WP, -0.1302705180215593_WP, &
        1.0086455296692498_WP, 0.06130944592639698_WP, 0.26853266846917656_WP, &
        -0.8167952437943864_WP, -1.1224999668605493_WP, 0.14625478817439308_WP, &
        -0.1522014148272574_WP, -1.194010427956042_WP, -0.5450983601742062_WP, &
        1.1381741716023244_WP, 0.8185083306916844_WP, 0.9029168287906786_WP, &
        0.3205044306725915_WP, -1.1663052832946548_WP, -0.321513349589935_WP, &
        -0.2311043871134415_WP, -0.7130532011711956_WP, -0.05214470234861428_WP, &
        -0.3090089215338304_WP, -0.12038731386336754_WP, 1.422446084034911_WP, &
        0.34628940193636926_WP, 0.35551191346501765_WP, -0.8087339077372864_WP, &
        0.9170578212528193_WP, -0.6663547469505713_WP, -0.21719495544653233_WP, &
        -0.03230668135477632_WP, 0.14523454601698482_WP, 0.03437743783244551_WP, &
        -0.08119054340141996_WP, 0.3715922386193339_WP, 0.3946647601593958_WP, &
        0.6561867102068776_WP, 0.6750836374977776_WP, 0.21321738972189763_WP, &
        0.05924700990519216_WP, 0.39046700208212015_WP, 0.45195126837380206_WP, &
        0.9197628809966887_WP, 0.02583081965901865_WP, 0.6323430068550735_WP, &
        0.01995865536266095_WP, 0.6077609082730555_WP, 0.13844779773203397_WP, &
        0.5649853464925103_WP, -0.847124480048162_WP, -0.012177059803286427_WP, &
        -1.599243963243535_WP, -0.18156026271803466_WP, -0.34582851903870293_WP, &
        0.4947595669752256_WP, 0.4155260274056104_WP, 0.460549990730005_WP, &
        0.9059306098166788_WP, 0.7622571582395441_WP, -0.18437593873550218_WP, &
        -0.8453218273531848_WP, 0.02187584573339664_WP, 0.09964869772989343_WP, &
        0.5524246961495529_WP, -0.007132188757407287_WP, 0.7891452591607631_WP, &
        -0.5876832475550082_WP, 0.38801546874723725_WP, -1.0152104995956592_WP, &
        0.737265270268444_WP, 0.03390790782393477_WP, 0.06893043805817389_WP, &
        -0.19498521573118963_WP, -0.010588859924052202_WP, 0.3688654218725612_WP, &
        1.3679766934908717_WP, -0.02742112031635921_WP, 0.5226459226822942_WP, &
        -0.0004606786115146644_WP, -0.009814961749722073_WP, -0.10747335304396946_WP, &
        -0.019922285233516468_WP, -0.3592477321150176_WP, 0.950191885475488_WP, &
        0.08117813545016149_WP, -0.4307292947884317_WP, 0.16681136863838764_WP, &
        0.520543163143804_WP, 1.0971597652873715_WP, 0.05252284659858297_WP, &
        -1.2146986661836985_WP, -0.39760285257397565_WP, 0.06652641864863493_WP, &
        -0.248299582484595_WP, 0.1794433990425347_WP, -0.8725307734560586_WP, &
        -0.22926706174969946_WP, -0.4556978592929579_WP, 0.6466101743347573_WP, &
        -0.17109814321494915_WP, -0.48400981477323174_WP, -0.26857577879559064_WP, &
        -0.023294509931201413_WP, -0.25264020026066714_WP, 7.978640186638575e-05_WP, &
        0.3703198907514363_WP, -0.3222928792748565_WP, 0.36577224625176097_WP, &
        0.21076664846239274_WP, 0.2963424642545056_WP, 0.6477186306682372_WP, &
        -0.28049935564988376_WP, 0.8604013472607417_WP, 1.1128460826811306_WP, &
        -0.2717489494596904_WP, -0.22010359530426848_WP, -0.8_WP, -1.6_WP]
    chi22_CFD_list = [0.9767321738539616_WP, 0.2640920929278583_WP, &
        0.2299821167604854_WP, 1.0120684731457443_WP, 1.708256001518653_WP, &
        0.7938850940011903_WP, 0.42306362958535493_WP, 0.23568383253036562_WP, &
        1.5960817452382452_WP, 0.7069485633681638_WP, 0.5942785990647463_WP, &
        0.09409702609616755_WP, 0.3634071367874687_WP, 1.9145055695548616_WP, &
        1.9394571621156174_WP, 0.40815156154820187_WP, 0.2830523779120533_WP, &
        1.3078704865042496_WP, 0.315940641705168_WP, 0.6385183027054233_WP, &
        0.5427153842709038_WP, 1.357977209296319_WP, 1.540279139499389_WP, &
        0.1345394993070026_WP, 0.5279132176174415_WP, 0.18249120545801456_WP, &
        0.6106636155704146_WP, 0.9904161597982089_WP, 1.9810419461654414_WP, &
        0.9163101187783578_WP, 0.13140705966415278_WP, 0.2724687795662817_WP, &
        0.830356065406284_WP, 0.8950795907034639_WP, 1.8873531752438442_WP, &
        1.3497359251170455_WP, 1.3313246913747685_WP, 0.9906787344183767_WP, &
        1.2530021881821847_WP, 1.2753393253474221_WP, 0.7930912536617669_WP, &
        0.31046533720711_WP, 1.8283370624450053_WP, 0.625673633237817_WP, &
        1.949154596462041_WP, 0.5932859052820307_WP, 0.7861454939316261_WP, &
        1.02306371341613_WP, 1.7744547268038224_WP, 0.9727149037656684_WP, &
        0.7919226533241905_WP, 1.5848924101211843_WP, 0.4079812716625326_WP, &
        0.969344589596074_WP, 0.912155691675336_WP, 1.9565635873338854_WP, &
        0.6810648357324167_WP, 0.9247913432983428_WP, 0.9403531793858524_WP, &
        1.1461975820128807_WP, 0.4940694348629191_WP, 0.6539174775158023_WP, &
        0.9955305789590359_WP, 1.6012321147979975_WP, 1.034700391688032_WP, &
        1.6320217939066162_WP, 0.16019390392001243_WP, 1.1383713729218887_WP, &
        0.6519252895384495_WP, 1.8830039833941365_WP, 0.44714725551420287_WP, &
        0.48593463333132636_WP, 1.3552940652182421_WP, 1.1858084893097764_WP, &
        1.3994261643728207_WP, 0.5763485789042675_WP, 1.196715601680969_WP, &
        1.2022363144875858_WP, 1.5734165755739793_WP, 1.5896383588150451_WP, &
        0.28582500348348594_WP, 0.6001123225065945_WP, 1.363753615741378_WP, &
        1.4498873280402937_WP, 0.19231815963567866_WP, 1.7997687608046742_WP, &
        0.4264648136643914_WP, 1.0412127869525951_WP, 1.7444439879385827_WP, &
        1.6251203183255543_WP, 1.8252587949813641_WP, 0.17712428563628158_WP, &
        0.5539123303798772_WP, 0.023696424513304004_WP, 0.06183607949490222_WP, &
        1.4989353929430698_WP, 0.8247493447365684_WP, 1.3653621000088398_WP, &
        1.748453430610656_WP, 1.6942489525958737_WP, 1.7626205955701821_WP, &
        0.36625544985602065_WP, 0.21888729954507213_WP, 0.7312330492739738_WP, &
        0.13415107971212503_WP, 0.17347336508952105_WP, 1.2443335777766495_WP, &
        1.6994405241074764_WP, 0.014359867080624422_WP, 1.4711411695479462_WP, &
        1.3561534152412895_WP, 1.0688090294016515_WP, 0.7389032278557499_WP, &
        0.10591065375928133_WP, 0.802176261652614_WP, 0.5903916287151363_WP, &
        0.009192193776444535_WP, 1.647965254558636_WP, 1.7556462747029364_WP, &
        1.7654989894979967_WP, 1.4307919175052428_WP, 0.6228890876233206_WP, &
        1.6822588168068204_WP, 1.0712154765271242_WP, 0.020938521588204972_WP, &
        0.5255475434993422_WP, 0.9875342742838447_WP, 1.6529814530496385_WP, &
        0.47535889442163093_WP, 1.6221954269794843_WP, 0.8032661205650753_WP, &
        0.8929282339645246_WP, 0.9072547064615017_WP, 0.0904612789554422_WP, &
        1.8678838297615477_WP, 0.8344595484872188_WP, 0.012376579391412124_WP, &
        1.5456250659930137_WP, 1.2620283818320337_WP, 0.3445464048215645_WP, &
        1.6243113582108568_WP, 1.5577998625247516_WP, 0.5931734472070427_WP, &
        0.950257383851105_WP, 0.9336767732120987_WP, 0.8270085962827145_WP, &
        1.571592874474108_WP, 1.6546377880046177_WP, 1.4_WP, 2.8_WP]
    OMIX_truth = reshape([0.0_WP, 0.06935237628171348_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171348_WP, &
        0.06935237628171337_WP, 0.06935237628171342_WP, 0.06935237628171348_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.05017197579677174_WP, 0.0_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171337_WP, 0.05860632441878655_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.0_WP, 0.06935237628171348_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.06935237628171337_WP, 0.0_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171337_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, 0.0_WP, &
        0.04617693959125235_WP, 0.06935237628171342_WP, 0.06935237628171348_WP, &
        0.06935237628171337_WP, 0.06935237628171342_WP, 0.06935237628171287_WP, 0.0_WP, &
        0.0_WP, 0.06935237628171342_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.06935237628171342_WP, &
        0.06935237628171348_WP, 0.0_WP, 0.06935237628171342_WP, 0.0_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, 0.0_WP, &
        0.06935237628171342_WP, 0.0_WP, 0.06935237628171348_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171348_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171337_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171337_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.04292814751410995_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, 0.0_WP, &
        0.06935237628171342_WP, 0.0_WP, 0.0_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171348_WP, 0.0_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.0693523762817134_WP, 0.04920946739103449_WP, 0.0_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171337_WP, 0.0_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171348_WP, 0.0_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.0_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171337_WP, &
        0.06935237628171342_WP, 0.0_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.0_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.0_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.0_WP, 0.0_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.0_WP, 0.06935237628171348_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171342_WP, 0.06935237628171342_WP, &
        0.06935237628171342_WP, 0.06935237628171348_WP, 0.06935237628171345_WP, &
        0.06935237628171345_WP, 0.0_WP, 0.0_WP, 0.937812275603084_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9498280242032282_WP, 0.937812275603084_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9413936755812135_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9378122756030841_WP, 0.9306476237182866_WP, 0.937812275603084_WP, &
        0.937812275603084_WP, 0.937812275603084_WP, 0.9306476237182866_WP, &
        0.937812275603084_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.937812275603084_WP, 0.9538230604087476_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182871_WP, 0.937812275603084_WP, 0.937812275603084_WP, &
        0.9306476237182866_WP, 0.937812275603084_WP, 0.937812275603084_WP, &
        0.937812275603084_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9378122756030851_WP, 0.9306476237182866_WP, 0.9378122756030839_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.937812275603084_WP, 0.9306476237182866_WP, 0.937812275603084_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9570718524858901_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.937812275603084_WP, &
        0.9306476237182866_WP, 0.937812275603084_WP, 0.9522021389501223_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.937812275603084_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9507905326089655_WP, 0.9378122756030834_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.937812275603084_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.937812275603084_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9601183998643827_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9834224755215359_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.937812275603084_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9470296478953595_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.937812275603084_WP, 0.9498211107639469_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9378122756030839_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.9306476237182866_WP, 0.9306476237182866_WP, 0.9306476237182866_WP, &
        0.937812275603084_WP, 0.937812275603084_WP, 0.062187724396916_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.062187724396916_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.06218772439691594_WP, 0.0_WP, 0.062187724396916_WP, 0.062187724396916_WP, &
        0.062187724396916_WP, 0.0_WP, 0.062187724396916_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.062187724396916_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.062187724396916_WP, 0.06218772439691597_WP, 0.0_WP, &
        0.06218772439691601_WP, 0.062187724396916_WP, 0.062187724396916_WP, 0.0_WP, &
        0.0_WP, 0.06218772439691489_WP, 0.0_WP, 0.06218772439691606_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.062187724396916_WP, &
        0.0_WP, 0.06218772439691599_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.062187724396916_WP, 0.0_WP, &
        0.06218772439691595_WP, 0.047797861049877766_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.06218772439691597_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.062187724396916555_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.062187724396915944_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.062187724396916_WP, 0.0_WP, &
        0.0_WP, 0.03988160013561731_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.016577524478464112_WP, 0.0_WP, 0.0_WP, &
        0.062187724396916_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.05297035210464057_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.06218772439691599_WP, 0.05017888923605307_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.06218772439691611_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.062187724396915944_WP, 0.062187724396915944_WP], &
        shape(OMIX_truth))
    FMIX_truth = reshape([0.36355826481408804_WP, 0.4942873461237324_WP, &
        0.7312849102882517_WP, 0.7978895552023769_WP, 0.682212986817348_WP, &
        0.12888971797439097_WP, 0.42740993017098217_WP, 0.3126234196184367_WP, &
        0.5996115488764736_WP, 0.5154345689237414_WP, 0.05017197579677174_WP, &
        0.588184248149486_WP, 0.7538605725770947_WP, 0.2624232747282652_WP, &
        0.40834447175190736_WP, 0.12053541591318306_WP, 0.05860632441878655_WP, &
        0.3398626422940706_WP, 0.7193866835408825_WP, 0.6640481184236258_WP, &
        0.591860510184272_WP, 0.3280170576532143_WP, 0.17160369871027903_WP, &
        0.6497336920748618_WP, 0.24282439105331538_WP, 0.4356141360551886_WP, &
        0.2889517004619201_WP, 0.8383234985178208_WP, 0.12269590458794256_WP, &
        0.5399913937855471_WP, 0.45755798944414344_WP, 0.7615323284989437_WP, &
        0.2928924591315496_WP, 0.652945664979061_WP, 0.45489134380595303_WP, &
        0.1704377129139125_WP, 0.35366979007535615_WP, 0.0_WP, 0.6299677723197672_WP, &
        0.5081654873019352_WP, 0.48784411601301564_WP, 0.8070919250810602_WP, &
        0.553701422080125_WP, 0.8741265152621437_WP, 0.532735102885292_WP, &
        0.3560543435874866_WP, 0.5774273763235112_WP, 0.5126846469898145_WP, &
        0.8425015411279363_WP, 0.36377833561301487_WP, 0.2094405921651169_WP, &
        0.7065152111322048_WP, 0.7405272407826968_WP, 0.5611744878333519_WP, &
        0.4911077153653768_WP, 0.7082967318632838_WP, 0.46732769480575503_WP, &
        0.6478065169086281_WP, 0.8002008434363034_WP, 0.4884194722746056_WP, &
        0.4617326814041486_WP, 0.42994632249944975_WP, 0.6044799233510052_WP, &
        0.47700048881478146_WP, 0.5522631512550497_WP, 0.6226333567869313_WP, &
        0.7008939553228397_WP, 0.43468964940064514_WP, 0.6531052784208797_WP, &
        0.5773290514696111_WP, 0.5146512193007922_WP, 0.3546087660148411_WP, &
        0.17181298695929903_WP, 0.5725908596207185_WP, 0.04292814751410995_WP, &
        0.09504498640716907_WP, 0.7654195877559891_WP, 0.40506521478655966_WP, &
        0.4828282710350445_WP, 0.1095368764735284_WP, 0.4097254237981287_WP, &
        0.5362790173821596_WP, 0.48557471004224356_WP, 0.6855575304097905_WP, &
        0.9522021389501223_WP, 0.44776272520611926_WP, 0.44287670532834184_WP, &
        0.4418182896074945_WP, 0.7522803187869928_WP, 0.5568828274474208_WP, &
        0.89596984383598_WP, 0.19917277031439165_WP, 0.38889497131344786_WP, &
        0.8207698124569848_WP, 0.4686451579146322_WP, 0.8195853586593806_WP, &
        0.8720128539434517_WP, 0.5670890588919855_WP, 0.0_WP, 0.6195916597294462_WP, &
        0.14795695630983619_WP, 0.5782425408235468_WP, 0.5005241294718743_WP, &
        0.367921850988942_WP, 0.5644410024275783_WP, 0.6791341636365786_WP, &
        0.6096664694116283_WP, 0.36381467077862556_WP, 0.5772456387362621_WP, &
        0.814098489356488_WP, 0.15697570216181822_WP, 0.9601183998643827_WP, &
        0.8417125739003475_WP, 0.49730510449445575_WP, 0.31371753958852855_WP, &
        0.6639494881083599_WP, 0.49023756977463584_WP, 0.7266152165156362_WP, &
        0.3140885181868398_WP, 0.8163837853196826_WP, 0.07114500827317544_WP, &
        0.41928003243670403_WP, 0.9834224755215359_WP, 0.4553014427974507_WP, &
        0.5760448433567242_WP, 0.20273630414595223_WP, 0.9325603128607829_WP, &
        0.5595984160010883_WP, 0.2901081572719364_WP, 0.9470296478953595_WP, &
        0.8467042579111438_WP, 0.5285052601561206_WP, 0.16444553839727216_WP, &
        0.6665741344561484_WP, 0.9498211107639469_WP, 0.33856608755311923_WP, &
        0.4170568958697838_WP, 0.6758667750582312_WP, 0.46042475952990347_WP, &
        0.3698542332002919_WP, 0.4898959120494432_WP, 0.526550757984087_WP, &
        0.31114962643241995_WP, 0.3077623817470957_WP, 0.8142476925332582_WP, &
        0.12089827075767522_WP, 0.7437124267289192_WP, 0.8461340639945406_WP, &
        0.34875078838669316_WP, 0.34875078838669316_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.4609643874258464_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.4225726236764888_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.9255684619763515_WP, 0.3804083402705538_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.5395752404700965_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, 0.0_WP, &
        0.636441735185912_WP, 0.5057126538762676_WP, 0.2687150897117483_WP, &
        0.20211044479762308_WP, 0.317787013182652_WP, 0.871110282025609_WP, &
        0.5725900698290178_WP, 0.6873765803815632_WP, 0.4003884511235264_WP, &
        0.4845654310762586_WP, 0.9498280242032282_WP, 0.411815751850514_WP, &
        0.2461394274229053_WP, 0.7375767252717348_WP, 0.5916555282480926_WP, &
        0.8794645840868169_WP, 0.9413936755812135_WP, 0.6601373577059294_WP, &
        0.28061331645911747_WP, 0.33595188157637423_WP, 0.408139489815728_WP, &
        0.6719829423467857_WP, 0.828396301289721_WP, 0.3502663079251383_WP, &
        0.7571756089466846_WP, 0.5643858639448114_WP, 0.7110482995380799_WP, &
        0.16167650148217916_WP, 0.8773040954120574_WP, 0.46000860621445294_WP, &
        0.5424420105558565_WP, 0.23846767150105624_WP, 0.7071075408684504_WP, &
        0.34705433502093896_WP, 0.545108656194047_WP, 0.8295622870860875_WP, &
        0.6463302099246439_WP, 0.5390356125741536_WP, 0.37003222768023275_WP, &
        0.4918345126980647_WP, 0.5121558839869844_WP, 0.19290807491893985_WP, &
        0.44629857791987504_WP, 0.12587348473785634_WP, 0.467264897114708_WP, &
        0.6439456564125134_WP, 0.0_WP, 0.48731535301018547_WP, 0.15749845887206368_WP, &
        0.6362216643869851_WP, 0.7905594078348831_WP, 0.29348478886779517_WP, &
        0.2594727592173032_WP, 0.4388255121666481_WP, 0.5088922846346232_WP, &
        0.2917032681367162_WP, 0.532672305194245_WP, 0.3521934830913719_WP, &
        0.19979915656369665_WP, 0.5115805277253944_WP, 0.5382673185958514_WP, &
        0.5700536775005502_WP, 0.3955200766489948_WP, 0.5229995111852186_WP, &
        0.44773684874495034_WP, 0.3773666432130687_WP, 0.2991060446771603_WP, &
        0.5653103505993549_WP, 0.34689472157912027_WP, 0.4226709485303889_WP, &
        0.48534878069920784_WP, 0.6453912339851589_WP, 0.828187013040701_WP, &
        0.42740914037928146_WP, 0.9570718524858901_WP, 0.9049550135928309_WP, &
        0.23458041224401094_WP, 0.5949347852134403_WP, 0.5171717289649556_WP, &
        0.8904631235264716_WP, 0.5902745762018713_WP, 0.46372098261784045_WP, &
        0.5144252899577564_WP, 0.31444246959020955_WP, 0.047797861049877766_WP, &
        0.5522372747938807_WP, 0.5571232946716582_WP, 0.5581817103925055_WP, &
        0.24771968121300725_WP, 0.4431171725525792_WP, 0.10403015616402_WP, &
        0.8008272296856084_WP, 0.6111050286865521_WP, 0.17923018754301523_WP, &
        0.5313548420853678_WP, 0.18041464134061935_WP, 0.12798714605654826_WP, &
        0.43291094110801454_WP, 0.07443153802364855_WP, 0.0_WP, 0.8520430436901638_WP, &
        0.42175745917645324_WP, 0.4994758705281257_WP, 0.6320781490110581_WP, &
        0.4355589975724217_WP, 0.32086583636342136_WP, 0.3903335305883717_WP, &
        0.6361853292213744_WP, 0.4227543612637379_WP, 0.185901510643512_WP, &
        0.8430242978381818_WP, 0.03988160013561731_WP, 0.15828742609965252_WP, &
        0.5026948955055442_WP, 0.6862824604114715_WP, 0.33605051189164015_WP, &
        0.5097624302253642_WP, 0.2733847834843638_WP, 0.6859114818131602_WP, &
        0.1836162146803174_WP, 0.9288549917268245_WP, 0.5807199675632959_WP, &
        0.016577524478464112_WP, 0.5446985572025493_WP, 0.42395515664327577_WP, &
        0.7972636958540478_WP, 0.06743968713921711_WP, 0.4404015839989117_WP, &
        0.7098918427280636_WP, 0.05297035210464057_WP, 0.1532957420888562_WP, &
        0.47149473984387935_WP, 0.8355544616027278_WP, 0.3334258655438515_WP, &
        0.05017888923605307_WP, 0.6614339124468808_WP, 0.5829431041302162_WP, &
        0.32413322494176877_WP, 0.0_WP, 0.6301457667997081_WP, 0.5101040879505567_WP, &
        0.473449242015913_WP, 0.68885037356758_WP, 0.6922376182529043_WP, &
        0.18575230746674176_WP, 0.8791017292423248_WP, 0.25628757327108076_WP, &
        0.1538659360054594_WP, 0.6512492116133068_WP, 0.6512492116133068_WP], &
        shape(FMIX_truth))
    chi_xi_truth = [1.132767773077786_WP, 0.380258937331631_WP, &
        0.32231588082852164_WP, 1.3066749627115932_WP, 1.9776570610390156_WP, &
        1.440161370765379_WP, 0.47742345670617875_WP, 0.35845548947936606_WP, &
        1.8636791942422235_WP, 0.9118666466943984_WP, 1.6251318242617525_WP, &
        0.10931395782490425_WP, 0.4530590013957429_WP, 2.369873025118244_WP, &
        2.3143107845353157_WP, 1.313111204636009_WP, 0.15269532742593664_WP, &
        1.7878640292672934_WP, 0.37192234868749646_WP, 0.9083243673581205_WP, &
        0.6354924137077123_WP, 1.5618954507082945_WP, 1.7749079095951896_WP, &
        0.1647563280549461_WP, 0.7519996339650575_WP, 0.20803103551783703_WP, &
        0.971010398224703_WP, 1.209928734025957_WP, 1.8202926661762868_WP, &
        1.2309182203977063_WP, 0.1551220353491466_WP, 0.3361589011119739_WP, &
        0.9041758742123089_WP, 1.107842812814551_WP, 2.197591808967397_WP, &
        3.0756719635466254_WP, 1.639633592510386_WP, 14.567247508513784_WP, &
        1.4253126854138518_WP, 1.624219476249524_WP, 1.7740105668728232_WP, &
        0.3330145389147089_WP, 2.1062844527975386_WP, 0.6986640791402295_WP, &
        2.767604901056922_WP, 1.0133749625138821_WP, 263.87394827333515_WP, &
        1.2967370767243962_WP, 1.793364998969004_WP, 1.1223800433561055_WP, &
        1.409600215125818_WP, 1.851099242515194_WP, 0.5066671160570511_WP, &
        1.1531132218397124_WP, 1.0995692175472214_WP, 2.8251984098634746_WP, &
        0.9090899875305495_WP, 1.2445399351466786_WP, 1.0736576829419413_WP, &
        1.4428886234643905_WP, 0.5883358639300641_WP, 0.7607424125502945_WP, &
        1.2004244186282103_WP, 1.9891823952573993_WP, 1.256110940211172_WP, &
        1.8970502025089544_WP, 0.2951291741443223_WP, 1.347593416745002_WP, &
        0.9346012843479854_WP, 2.659689359291771_WP, 0.5331631508740942_WP, &
        0.6686046677315067_WP, 1.7905202616886762_WP, 1.3719438334601741_WP, &
        3.2964274161676066_WP, 0.7117842897910027_WP, 1.4493973255384396_WP, &
        1.4809946482255596_WP, 2.7531776582913454_WP, 1.9402824162719996_WP, &
        0.49921050854881327_WP, 0.7505100552313926_WP, 1.7777134910678325_WP, &
        1.688851087365818_WP, 0.21516738723968015_WP, 2.3014747854924_WP, &
        0.6072399746565132_WP, 1.3981488501394388_WP, 2.107150930638721_WP, &
        1.921143725328399_WP, 2.403354399731514_WP, 0.22278556895518137_WP, &
        0.7912887415808697_WP, 0.02779841849336256_WP, 0.08400577345805631_WP, &
        2.0427984074141037_WP, 0.9648774401629001_WP, 2.066971162055553_WP, &
        723.9149993485896_WP, 57.666649480234_WP, 1.5457982629993905_WP, &
        0.4692759062184804_WP, 0.2901537274112154_WP, 0.8547990207190364_WP, &
        0.17481984151708096_WP, 0.21136323279175298_WP, 1.649364537858087_WP, &
        2.0586800999325883_WP, 0.03930251523249754_WP, 1.7531774544188121_WP, &
        1.5846668328358011_WP, 0.7297184131804725_WP, 0.9764498728431649_WP, &
        0.12975627116140157_WP, 0.8966057390655601_WP, 0.9690419804519858_WP, &
        0.020169420432177565_WP, 1.9576514016417363_WP, 2.2726322724489827_WP, &
        2.3481089076835024_WP, 3.3118870956066986_WP, 0.8631138136499218_WP, &
        1.4469790739829962_WP, 1.2419938325355253_WP, 0.035955981397891516_WP, &
        0.6089505582901727_WP, 1.274426390807465_WP, 1.908596058830712_WP, &
        0.5740061245606682_WP, 1.1434507897167672_WP, 0.9592031998790499_WP, &
        1.0900663349426187_WP, 0.9250252017969107_WP, 0.21888524439590384_WP, &
        0.8522583283581401_WP, 0.9753484678345014_WP, 0.01658704228973489_WP, &
        1.8808398348178168_WP, 41.820309884841855_WP, 0.4912091066244235_WP, &
        1.9491941128624581_WP, 1.9613668481762625_WP, 0.8181042591045876_WP, &
        1.0909958141858958_WP, 1.3150243210327155_WP, 1.4463092959111938_WP, &
        2.0331818496814593_WP, 2.1305406253890107_WP, 1.6406789539363185_WP, &
        3.281357907872637_WP]
    chi_eta_truth = [0.2831924404418469_WP, 3.407916517622378_WP, &
        2.5061864608009254_WP, 5.6072410009639855_WP, 1.0234385602837304_WP, &
        3.483074861699091_WP, 1.2322792388487118_WP, 1.6250828086322269_WP, &
        0.5221621639634288_WP, 4.28848662574346_WP, 0.8464867171963026_WP, &
        1.4228124930789279_WP, 1.2540724501214253_WP, 2.82590581675392_WP, &
        0.9132325631956881_WP, 1.4309665303150847_WP, 0.46485487718235663_WP, &
        3.9714534739655387_WP, 1.6895335607300108_WP, 2.995545668982516_WP, &
        0.6403872760759186_WP, 0.40856608789820786_WP, 1.2112815898066858_WP, &
        0.338513593398175_WP, 3.703795036168991_WP, 0.3622568823872017_WP, &
        3.209078856136006_WP, 5.903162472436505_WP, 2.657399592838161_WP, &
        2.577354702433898_WP, 2.0391371353075285_WP, 0.954115496882195_WP, &
        0.05871094530685169_WP, 0.6944830433791418_WP, 1.3016661395835953_WP, &
        3.0327780598467524_WP, 3.721933250412852_WP, 4.722663460405124_WP, &
        1.5589375268250998_WP, 3.8996349671864756_WP, 16.97401707030845_WP, &
        1.0575447103223115_WP, 1.7693484561770363_WP, 0.07695067003153801_WP, &
        6.92493257026316_WP, 4.302276279374929_WP, 363.1557337012671_WP, &
        1.656246414345011_WP, 1.9270866838040923_WP, 0.4629775234811379_WP, &
        7.545083859573488_WP, 0.7706960448298249_WP, 0.5641114978879843_WP, &
        0.5995962220864525_WP, 3.7638366361870417_WP, 14.250448276612632_WP, &
        1.915205399907995_WP, 3.813732616617387_WP, 0.9356599526924778_WP, &
        4.069262857133191_WP, 1.3657306353567547_WP, 0.5357332599926765_WP, &
        3.5586976834473742_WP, 1.7711570633070162_WP, 2.423792678731428_WP, &
        2.4967861470021933_WP, 6.382454095436403_WP, 2.1471891759147326_WP, &
        5.148009832777838_WP, 6.822227123198178_WP, 0.641961330553597_WP, &
        5.119690663028974_WP, 1.861908215201699_WP, 1.2231976163499538_WP, &
        1.6823613653769602_WP, 0.09505068392788996_WP, 3.72696301605096_WP, &
        0.989431158141284_WP, 16.19805402671838_WP, 1.526935502024379_WP, &
        4.0640215531783985_WP, 0.835862176038681_WP, 3.8304775938130295_WP, &
        0.29908051185041784_WP, 0.5420795928401864_WP, 1.4535929017286104_WP, &
        2.1294212046665977_WP, 2.6639269091303093_WP, 4.83513645032439_WP, &
        3.9093403699251157_WP, 10.20020165844678_WP, 3.5642397264888745_WP, &
        4.383944288652058_WP, 0.7087748081731517_WP, 2.8888335892589865_WP, &
        9.29435292952162_WP, 2.0977387003031054_WP, 9.56875417987804_WP, &
        2.777289988839199_WP, 107.9619421261674_WP, 2.711644166229236_WP, &
        3.6280991748249143_WP, 0.7319805869248714_WP, 2.3949972510024042_WP, &
        0.6263861676729577_WP, 0.5439368983029066_WP, 3.510064128144979_WP, &
        3.5902210728287494_WP, 4.241942658760962_WP, 2.8271078641501677_WP, &
        1.4586376741662748_WP, 2.6154220083784905_WP, 1.8837297828857968_WP, &
        2.349002948931458_WP, 1.8827561304951321_WP, 10.15411831324852_WP, &
        1.5121167101194146_WP, 1.7606429044665168_WP, 4.322783677070838_WP, &
        6.823351830769706_WP, 2.0436864802156505_WP, 6.485792925337634_WP, &
        0.6843926652147787_WP, 0.44109685085343947_WP, 1.57194041935644_WP, &
        0.2881564765076766_WP, 5.644405936988331_WP, 0.5893767736019121_WP, &
        4.481751339899176_WP, 2.0038477840423203_WP, 3.496654185998911_WP, &
        0.2609257329573239_WP, 1.1416953304562616_WP, 3.929161800872756_WP, &
        3.0900620950703432_WP, 2.1792212137398694_WP, 0.5137434912049087_WP, &
        2.099721600519302_WP, 20.989824122383684_WP, 2.550176724234632_WP, &
        0.9047009600091841_WP, 5.338764816636616_WP, 2.672373150515512_WP, &
        0.6901559290067871_WP, 9.17607639752359_WP, 2.136839414694431_WP, &
        3.529996727291321_WP, 4.022878039891345_WP, 2.5935654845324105_WP, &
        5.187130969064821_WP]
    x_prime_truth = [5.097711575088253_WP, 4.85969218749866_WP, &
        4.905163671086312_WP, 4.91308023042848_WP, 4.898342046220411_WP, &
        4.461924679705692_WP, 4.837738032305498_WP, 4.77816001000066_WP, &
        4.88433782436035_WP, 4.865448729163576_WP, 4.0_WP, 5.151008610325061_WP, &
        4.908003709433125_WP, 4.735723226708734_WP, 4.830162078638721_WP, &
        4.424630713253064_WP, 4.0_WP, 4.795939983831158_WP, 4.903595134760689_WP, &
        5.185109022474037_WP, 4.882823105971167_WP, 5.092543605615547_WP, &
        5.075070016971463_WP, 5.177544122828414_WP, 4.714392874698959_WP, &
        5.110186537916189_WP, 4.759986267009863_WP, 4.917272536909284_WP, &
        4.434762093204139_WP, 4.871567626669887_WP, 4.848429318508972_WP, &
        5.260780524275973_WP, 4.8423416576576574_WP, 4.893785379088262_WP, &
        4.847540786990008_WP, 4.593092543334333_WP, 4.803906417149917_WP, &
        2.504684291156173_WP, 5.16806029244203_WP, 5.126440342821352_WP, &
        4.857839063739239_WP, 5.322369731920488_WP, 5.139341076744548_WP, &
        5.494049438024441_WP, 4.869818271958988_WP, 4.805219687582121_WP, &
        1.4505940417603749_WP, 4.8647270272497725_WP, 5.39484655813319_WP, &
        4.809355397250786_WP, 4.66886850555198_WP, 4.901838806597561_WP, &
        4.906347298975239_WP, 4.876415664315965_WP, 4.858783777750027_WP, &
        4.902085703403895_WP, 4.851597974927337_WP, 4.892942762273114_WP, &
        5.3112511857731_WP, 4.858006528775902_WP, 5.115533160287609_WP, &
        4.838695268101049_WP, 4.885269347082281_WP, 4.854607326600366_WP, &
        4.874421503364644_WP, 4.888614422074006_WP, 4.9010515417417615_WP, &
        4.840455422903818_WP, 4.893811337049062_WP, 4.879873745994291_WP, &
        4.865243928935142_WP, 4.804425657433379_WP, 4.596349626945591_WP, &
        4.878879700720944_WP, 4.0_WP, 4.270320519752504_WP, 4.909392995174011_WP, &
        4.828787134145159_WP, 4.856362229715667_WP, 4.36685818954794_WP, &
        4.830734505955668_WP, 5.134105910079479_WP, 4.8571746533593565_WP, &
        5.197771390353093_WP, 0.0_WP, 4.845113555958039_WP, 4.843404777340238_WP, &
        4.843029639304147_WP, 4.907810460343373_WP, 4.875463252117859_WP, &
        4.92259518915861_WP, 4.6517979030354315_WP, 5.101762743681845_WP, &
        4.915503256541434_WP, 4.8520151651831505_WP, 4.915381142978012_WP, &
        4.9204686307454235_WP, 4.877704612363323_WP, 2.0265274734734735_WP, &
        1.4056337821190523_WP, 4.531266538516223_WP, 4.880063517666929_WP, &
        4.861440493678316_WP, 4.811502426139409_WP, 5.142776810359831_WP, &
        4.897881184609606_WP, 4.8862453820879415_WP, 4.809374437448365_WP, &
        5.147101319572478_WP, 4.914810827942288_WP, 4.5581967442947215_WP, 0.0_WP, &
        4.917605631147523_WP, 4.860543606621101_WP, 4.778933698215676_WP, &
        4.895545704117788_WP, 4.858533126472548_WP, 4.904554192225314_WP, &
        4.779194805712515_WP, 4.915049297243751_WP, 4.025196876562005_WP, &
        4.834591750342456_WP, 0.0_WP, 4.847678110010809_WP, 4.87960594199996_WP, &
        5.078001450110304_WP, 4.925632288523019_WP, 4.876067597229262_WP, &
        4.760943032647286_WP, 0.0_WP, 4.918091381218741_WP, 4.868776374598_WP, &
        4.57826538221931_WP, 5.186511398254846_WP, 6.0_WP, 4.795158526410793_WP, &
        4.833710035804402_WP, 4.897387504696117_WP, 1.5753552757913187_WP, &
        4.812487271859462_WP, 4.858434466228585_WP, 4.868289286018254_WP, &
        4.777109241374016_WP, 4.774656097057685_WP, 4.91482643805112_WP, &
        4.426357582725718_WP, 4.906748396572118_WP, 4.9180361845328555_WP, &
        5.095489903539171_WP, 5.095489903539171_WP]
    eta_opt_truth = [0.636441735185912_WP, 0.4942873461237324_WP, &
        0.7312849102882517_WP, 0.7978895552023769_WP, 0.682212986817348_WP, &
        0.12888971797439097_WP, 0.42740993017098217_WP, 0.3126234196184367_WP, &
        0.5996115488764736_WP, 0.5154345689237414_WP, 0.05017197579677174_WP, &
        0.411815751850514_WP, 0.7538605725770947_WP, 0.2624232747282652_WP, &
        0.40834447175190736_WP, 0.12053541591318306_WP, 0.05860632441878655_WP, &
        0.3398626422940706_WP, 0.7193866835408825_WP, 0.33595188157637423_WP, &
        0.591860510184272_WP, 0.6719829423467857_WP, 0.828396301289721_WP, &
        0.3502663079251383_WP, 0.24282439105331538_WP, 0.5643858639448114_WP, &
        0.2889517004619201_WP, 0.8383234985178208_WP, 0.12269590458794256_WP, &
        0.5399913937855471_WP, 0.45755798944414344_WP, 0.23846767150105624_WP, &
        0.2928924591315496_WP, 0.652945664979061_WP, 0.45489134380595303_WP, &
        0.1704377129139125_WP, 0.35366979007535615_WP, 0.9306476237182871_WP, &
        0.37003222768023275_WP, 0.4918345126980647_WP, 0.48784411601301564_WP, &
        0.19290807491893985_WP, 0.44629857791987504_WP, 0.12587348473785634_WP, &
        0.532735102885292_WP, 0.3560543435874866_WP, 0.9378122756030851_WP, &
        0.5126846469898145_WP, 0.15749845887206368_WP, 0.36377833561301487_WP, &
        0.2094405921651169_WP, 0.7065152111322048_WP, 0.7405272407826968_WP, &
        0.5611744878333519_WP, 0.4911077153653768_WP, 0.7082967318632838_WP, &
        0.46732769480575503_WP, 0.6478065169086281_WP, 0.19979915656369665_WP, &
        0.4884194722746056_WP, 0.5382673185958514_WP, 0.42994632249944975_WP, &
        0.6044799233510052_WP, 0.47700048881478146_WP, 0.5522631512550497_WP, &
        0.6226333567869313_WP, 0.7008939553228397_WP, 0.43468964940064514_WP, &
        0.6531052784208797_WP, 0.5773290514696111_WP, 0.5146512193007922_WP, &
        0.3546087660148411_WP, 0.17181298695929903_WP, 0.5725908596207185_WP, &
        0.04292814751410995_WP, 0.09504498640716907_WP, 0.7654195877559891_WP, &
        0.40506521478655966_WP, 0.4828282710350445_WP, 0.1095368764735284_WP, &
        0.4097254237981287_WP, 0.46372098261784045_WP, 0.48557471004224356_WP, &
        0.31444246959020955_WP, 0.047797861049877766_WP, 0.44776272520611926_WP, &
        0.44287670532834184_WP, 0.4418182896074945_WP, 0.7522803187869928_WP, &
        0.5568828274474208_WP, 0.89596984383598_WP, 0.19917277031439165_WP, &
        0.6111050286865521_WP, 0.8207698124569848_WP, 0.4686451579146322_WP, &
        0.8195853586593806_WP, 0.8720128539434517_WP, 0.5670890588919855_WP, &
        0.9507905326089655_WP, 0.9378122756030834_WP, 0.14795695630983619_WP, &
        0.5782425408235468_WP, 0.5005241294718743_WP, 0.367921850988942_WP, &
        0.4355589975724217_WP, 0.6791341636365786_WP, 0.6096664694116283_WP, &
        0.36381467077862556_WP, 0.4227543612637379_WP, 0.814098489356488_WP, &
        0.15697570216181822_WP, 0.03988160013561731_WP, 0.8417125739003475_WP, &
        0.49730510449445575_WP, 0.31371753958852855_WP, 0.6639494881083599_WP, &
        0.49023756977463584_WP, 0.7266152165156362_WP, 0.3140885181868398_WP, &
        0.8163837853196826_WP, 0.07114500827317544_WP, 0.41928003243670403_WP, &
        0.016577524478464112_WP, 0.4553014427974507_WP, 0.5760448433567242_WP, &
        0.7972636958540478_WP, 0.9325603128607829_WP, 0.5595984160010883_WP, &
        0.2901081572719364_WP, 0.05297035210464057_WP, 0.8467042579111438_WP, &
        0.5285052601561206_WP, 0.16444553839727216_WP, 0.3334258655438515_WP, &
        0.05017888923605307_WP, 0.33856608755311923_WP, 0.4170568958697838_WP, &
        0.6758667750582312_WP, 0.9378122756030839_WP, 0.3698542332002919_WP, &
        0.4898959120494432_WP, 0.526550757984087_WP, 0.31114962643241995_WP, &
        0.3077623817470957_WP, 0.8142476925332582_WP, 0.12089827075767522_WP, &
        0.7437124267289192_WP, 0.8461340639945406_WP, 0.6512492116133068_WP, &
        0.6512492116133068_WP]
    Z_opt_truth = [0.6048589857571001_WP, 0.6304412237196907_WP, &
        0.6485692685134896_WP, 0.6561713116691522_WP, 0.9504482738470732_WP, &
        0.5167079628734655_WP, 0.5159478155864052_WP, 0.6567791171025946_WP, &
        0.9125279818662285_WP, 0.7999093958105298_WP, 0.5009983507901415_WP, &
        0.6731986411576979_WP, 0.5831198380250244_WP, 0.7363934845262408_WP, &
        0.7387461396395444_WP, 0.30998318010650877_WP, 0.614732832501943_WP, &
        0.7079918481585301_WP, 0.9474503011058494_WP, 0.5050030027218636_WP, &
        0.9396410384959057_WP, 0.8594919105851018_WP, 0.8880960292006596_WP, &
        0.5873326652882749_WP, 0.5022588289422737_WP, 0.9475742169823781_WP, &
        0.5578435286971573_WP, 0.8496888091137218_WP, 0.289848571595691_WP, &
        0.6084296131785207_WP, 0.9801911020679289_WP, 0.6948467156099088_WP, &
        0.565742780847512_WP, 0.5843582205492991_WP, 0.9257372315592673_WP, &
        0.33300156245940177_WP, 0.6963051846955969_WP, 0.775389306007491_WP, &
        0.6993469695390334_WP, 0.42945543341667314_WP, 0.16939240420170432_WP, &
        0.7168245382580872_WP, 0.6147157757721099_WP, 0.7136709728476329_WP, &
        0.5513136499341094_WP, 0.5815777096747929_WP, 0.9434568710215755_WP, &
        0.685652483061831_WP, 0.44989426623199813_WP, 0.7782066168983495_WP, &
        0.15751982926015617_WP, 0.8005431622867054_WP, 0.5972923127209542_WP, &
        0.633300403475135_WP, 0.7044371153060898_WP, 0.5338068583216047_WP, &
        0.6416620229511898_WP, 0.5712798778815312_WP, 0.8157878551984283_WP, &
        0.7993650192988728_WP, 0.8094824106292635_WP, 0.7377925018459048_WP, &
        0.7475451187272166_WP, 0.6681836178653472_WP, 0.7353246879570287_WP, &
        0.9475418139267402_WP, 0.5309589220868658_WP, 0.9038300092607037_WP, &
        0.5979481699726886_WP, 0.4949506554304258_WP, 0.8873925137070161_WP, &
        0.566787927574251_WP, 0.6925579578317271_WP, 0.9908748399952833_WP, &
        0.42338523246373416_WP, 0.5631949894577528_WP, 0.8268228971081947_WP, &
        0.6560239198619046_WP, 0.3282750488130244_WP, 0.7833459022031537_WP, &
        0.6210060238618886_WP, 0.6907072307749464_WP, 0.5028025348261178_WP, &
        0.8554140048449289_WP, 0.9059453620873753_WP, 0.3647833582637219_WP, &
        0.6797274954628142_WP, 0.6067791707269882_WP, 0.8016049451471953_WP, &
        0.9189674008398419_WP, 0.5981787629612958_WP, 0.5409688284797022_WP, &
        0.4801073606876419_WP, 0.9496528278119626_WP, 0.8093893048519027_WP, &
        0.5344485449170326_WP, 0.9117338785600525_WP, 0.4756471488036213_WP, &
        0.9187386105833867_WP, 0.8429696893980713_WP, 0.4955339606032282_WP, &
        0.8509259988801509_WP, 0.5447587852613137_WP, 0.9288368869569459_WP, &
        0.6346209936549525_WP, 0.6957936805219522_WP, 0.5767300517175411_WP, &
        0.8969876108294939_WP, 0.5767392538242199_WP, 0.8603196353542725_WP, &
        0.8918338609020277_WP, 0.7068367570507118_WP, 0.4433255716002688_WP, &
        0.7235654366696753_WP, 0.5890449140204963_WP, 0.5434494496551591_WP, &
        0.7771337260656412_WP, 0.7486298523639032_WP, 0.5851492831870428_WP, &
        0.5568591373877356_WP, 0.5330343497550432_WP, 0.5118150773471024_WP, &
        0.547979339787901_WP, 0.7771333965757589_WP, 0.6941636690240706_WP, &
        0.6793449435421294_WP, 0.5936799742580133_WP, 0.8849167991830958_WP, &
        0.45012351402149464_WP, 0.6206264718993268_WP, 0.8878793801866641_WP, &
        0.5722235505146267_WP, 0.6929681376402697_WP, 0.4984138179999163_WP, &
        0.42914760727005047_WP, 0.6253521535400157_WP, 0.6940144031911826_WP, &
        0.7728349209799973_WP, 0.7513913225198774_WP, 0.69971498522518_WP, &
        0.805976563316718_WP, 0.6738623603503435_WP, 0.7592679279086313_WP, &
        0.5152357417806868_WP, 0.5828101821354946_WP, 0.7468866511488089_WP, &
        0.4651911761328296_WP, 0.5051179928525552_WP, 0.5734754061064401_WP, &
        0.5734754061064401_WP]
    Z_stoic_truth = [-1.6463158024381813e-17_WP, -5.849863305819067e-17_WP, &
        -8.594240447965067e-18_WP, -8.66337333164654e-18_WP, -5.980805317468996e-17_WP, &
        4.809433157684395e-17_WP, -8.292339052923502e-18_WP, -5.728621573144398e-17_WP, &
        -8.460765213780366e-18_WP, -8.37758738573573e-18_WP, 0.018107396470563804_WP, &
        -1.6898753580871714e-17_WP, -8.617549112552957e-18_WP, &
        -8.137142028115456e-18_WP, -8.274102999690523e-18_WP, 4.8049474900873435e-17_WP, &
        0.010235780751963344_WP, -8.209256620157286e-18_WP, -8.582006568166591e-18_WP, &
        5.967893108599089e-17_WP, -5.917126039177966e-17_WP, -1.6396285698571213e-17_WP, &
        -1.6108330136910382e-17_WP, -1.7022163249898935e-17_WP, &
        4.8714547412522704e-17_WP, -1.6600422224750248e-17_WP, &
        -8.161703240469812e-18_WP, -8.705887314655204e-18_WP, 4.8061067169188126e-17_WP, &
        -8.401683283438663e-18_WP, -8.321340209280807e-18_WP, &
        -1.7250997538977416e-17_WP, 0.02238547634417851_WP, -8.514326459457057e-18_WP, &
        -5.823136783643968e-17_WP, 4.8318664627076393e-17_WP, -8.222248860551893e-18_WP, &
        9.992007221626409e-16_WP, -1.6982335586553113e-17_WP, 1.6740962465130275e-17_WP, &
        -8.35067915028676e-18_WP, -2.6019037683127538e-17_WP, -1.683039273865218e-17_WP, &
        -1.748776395155398e-17_WP, -8.394548776359884e-18_WP, -5.757147789149611e-17_WP, &
        2.4424906541753444e-15_WP, -8.37489768187216e-18_WP, -9.581334493840124e-17_WP, &
        -8.231786898225685e-18_WP, -8.088528072290238e-18_WP, -8.56881114661237e-18_WP, &
        -8.603767604280193e-18_WP, -8.422580406070728e-18_WP, -8.353853034592022e-18_WP, &
        -8.570635084653719e-18_WP, -8.330781820868213e-18_WP, -8.509135919781499e-18_WP, &
        -1.733158465407318e-17_WP, -8.351238515113582e-18_WP, -8.325372083649709e-18_WP, &
        -5.806339815234538e-17_WP, -8.46562635616401e-18_WP, -8.340150907052806e-18_WP, &
        -8.413776718518887e-18_WP, -5.938661523864116e-17_WP, -8.563061136735905e-18_WP, &
        -8.299323306339173e-18_WP, 5.1086926622582913e-17_WP, -8.43858688283744e-18_WP, &
        -8.376821014022313e-18_WP, 4.9338803466524754e-17_WP, -8.054354338431491e-18_WP, &
        -8.433885823747393e-18_WP, 0.024757175521933298_WP, -7.985520841422192e-18_WP, &
        -8.629532441810166e-18_WP, -8.270974481541365e-18_WP, -8.345805887760902e-18_WP, &
        -7.998424663380739e-18_WP, -8.275421181313268e-18_WP, &
        -1.6796063396373584e-17_WP, -8.348473550599255e-18_WP, 3.418965043775032e-17_WP, &
        0.016853270759143335_WP, -8.31189535794413e-18_WP, -8.307192139207555e-18_WP, &
        -8.306174023381261e-18_WP, -8.61591343550356e-18_WP, -8.418338290342682e-18_WP, &
        -8.7672257417074e-18_WP, -5.65542180471359e-17_WP, 1.6511164147218813e-17_WP, &
        -8.687379472537803e-18_WP, -8.332056679605199e-18_WP, -8.686133469673512e-18_WP, &
        -8.74162971995407e-18_WP, 2.528530112372974e-17_WP, 0.5959910750594837_WP, &
        -1.1102230246251565e-15_WP, -8.032837249782922e-18_WP, &
        -8.439493815913692e-18_WP, -8.36302413303543e-18_WP, 4.941421775131033e-17_WP, &
        4.212906038615433e-17_WP, -8.540875505249165e-18_WP, -8.470811346728179e-18_WP, &
        -5.762274855830446e-17_WP, -1.687700815696857e-17_WP, -8.680366139921226e-18_WP, &
        -8.040958231788112e-18_WP, 0.025909319446478563_WP, -8.709469694295349e-18_WP, &
        -8.359886705428539e-18_WP, -8.184766770567444e-18_WP, -8.525461644264765e-18_WP, &
        -8.353006573105146e-18_WP, -8.589434864577792e-18_WP, -8.185113240444484e-18_WP, &
        -8.682767317529918e-18_WP, 4.778598367269689e-17_WP, -8.284553009166778e-18_WP, &
        0.05172238317307287_WP, -8.319162485339444e-18_WP, -8.437312221300571e-18_WP, &
        -1.6164835922832837e-17_WP, -8.806610171604891e-18_WP, &
        -8.421022032400447e-18_WP, -8.16277731345076e-18_WP, 0.010854260363507562_WP, &
        -8.71475146360714e-18_WP, -8.390395511423465e-18_WP, -8.047696931829223e-18_WP, &
        -8.52812192656292e-18_WP, 0.014099890717023961_WP, -8.208038697499568e-18_WP, &
        -8.28242644527101e-18_WP, -8.537554096908828e-18_WP, -4.440892098500626e-16_WP, &
        -5.766271376901939e-17_WP, -8.352674261706364e-18_WP, -8.388477781810242e-18_WP, &
        -8.182369310713652e-18_WP, -8.179209052409164e-18_WP, -8.680522868294186e-18_WP, &
        -5.605999167502276e-17_WP, -3.44282233831814e-17_WP, -3.485659123814957e-17_WP, &
        4.1088077343559854e-17_WP, 4.1088077343559854e-17_WP]

    do i = 1, N
        print *, "Checking ", i
        call get_optimal_bcs(Z1_CFD_list(i), Z2_CFD_list(i), chi11_CFD_list(i), chi12_CFD_list(i), &
            chi22_CFD_list(i), leftover_ox1_CFD, leftover_ox2_CFD, leftover_ox3_CFD, OMIX, FMIX, chi_xi, &
            chi_eta, x_prime, Z_opt, eta_opt, Z_stoic)
        num_failures = num_failures + is_failure(OMIX_truth(i, 1), OMIX(1), "OMIX1")
        num_failures = num_failures + is_failure(OMIX_truth(i, 2), OMIX(2), "OMXI2")
        num_failures = num_failures + is_failure(OMIX_truth(i, 3), OMIX(3), "OMIX3")
        num_failures = num_failures + is_failure(FMIX_truth(i, 1), FMIX(1), "FMIX1")
        num_failures = num_failures + is_failure(FMIX_truth(i, 2), FMIX(2), "FMIX2")
        num_failures = num_failures + is_failure(FMIX_truth(i, 3), FMIX(3), "FMIX3")
        num_failures = num_failures + is_failure(x_prime_truth(i), x_prime, "x_prime")
        num_failures = num_failures + is_failure(chi_xi_truth(i), chi_xi, "chi_xi")
        num_failures = num_failures + is_failure(chi_eta_truth(i), chi_eta, "chi_eta")
        num_failures = num_failures + is_failure(eta_opt_truth(i), eta_opt, "eta")
        num_failures = num_failures + is_failure(Z_opt_truth(i), Z_opt, "Z_opt")
        num_failures = num_failures + is_failure(Z_stoic_truth(i), Z_stoic, "Z_stoic")
    end do

    print *, "Number of failures? ", num_failures
end subroutine test_get_optimal_bcs

integer function is_failure(truth, actual_value, label)
    real(WP), intent(in) :: truth, actual_value
    character(len=*), intent(in) :: label

    if (abs(truth - actual_value) > 0.00000001_WP) then
        print *, ">>>>   ERROR ", label, ":", truth, "versus", actual_value
        is_failure = 1
    else
        ! print *, "    Checked ", label, ":", truth, "versus", actual_value
        is_failure = 0
    end if
end function is_failure

end program q2df_generalization