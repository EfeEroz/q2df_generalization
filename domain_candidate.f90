module domain_candidate
    use precision
    implicit none
    real(WP) :: cand_Z1, cand_Z2, cand_chi11, cand_chi12, cand_chi22, x

    real(WP) :: b
    real(WP) ::  cand_xi, cand_eta
    real(WP) :: cand_chi_xi, cand_chi_eta, cand_chi_ratio

    real(WP), parameter :: TOL = 0.000001_WP

    contains

subroutine set_dom_cand_inputs(Z1_in, Z2_in, chi11_in, chi12_in, chi22_in, x_in)
    use precision
    implicit none

    real(WP), intent(in) :: Z1_in, Z2_in, chi11_in, chi12_in, chi22_in, x_in
    cand_Z1 = Z1_in
    cand_Z2 = Z2_in
    cand_chi11 = chi11_in
    cand_chi12 = chi12_in
    cand_chi22 = chi22_in
    x = x_in
end subroutine set_dom_cand_inputs

subroutine compute_intermediate_vars()
    use precision
    implicit none
    real(WP) :: cand_Z3, cand_chi_xi_numer, cand_chi_eta_numer, cand_chi_denom
    integer :: attempt
    logical :: success
    
    success = .false.
    do attempt = 1, 5
        cand_Z3 = 1.0_WP - cand_Z1 - cand_Z2
        if ((cand_Z1+cand_Z2+x*(cand_Z1+cand_Z3))**2 - 4.0_WP*cand_Z1*x < -TOL .or. eq(2.0_WP*x, 0.0_WP)) then
            call perturb_x()
            cycle
        end if
        b = (cand_Z1+cand_Z2+x*(cand_Z1+cand_Z3) - sqrt((cand_Z1+cand_Z2+x*(cand_Z1+cand_Z3))**2 - &
            4.0_WP*cand_Z1*x)) / (2.0_WP*x)
        cand_eta = b
        if (eq(b, 1.0_WP)) then
            call perturb_x()
            cycle
        end if
        cand_xi = cand_Z2 / (1.0_WP - b)
        cand_chi_eta_numer = (1.0_WP-cand_eta)**2*cand_chi11 - 2.0_WP*cand_eta*(1.0_WP-cand_eta)*(1.0_WP-x)*cand_chi12 + &
            (cand_eta*(1.0_WP-x))**2*cand_chi22
        cand_chi_xi_numer = cand_xi**2*cand_chi11 + 2.0_WP*cand_xi*(x+(1.0_WP-x)*cand_xi)*cand_chi12 + &
            (x+(1.0_WP-x)*cand_xi)**2*cand_chi22
        cand_chi_denom = ((1.0_WP-cand_eta)*(x+(1.0_WP-x)*cand_xi) + cand_eta*cand_xi*(1.0_WP-x))**2
        if (eq(cand_chi_denom, 0.0_WP) .or. eq(cand_chi_xi_numer, 0.0_WP)) then
            call perturb_x()
            cycle
        end if
        cand_chi_eta = cand_chi_eta_numer / cand_chi_denom
        cand_chi_xi = cand_chi_xi_numer / cand_chi_denom
        cand_chi_ratio = cand_chi_eta_numer / cand_chi_xi_numer
        success = .true.
    end do

    if (.not. success) cand_chi_ratio = -1.0_WP   ! So candidates with NaN issues are not selected

end subroutine compute_intermediate_vars

real(WP) function get_chi_ratio()
    use precision
    implicit none
    call compute_intermediate_vars()
    get_chi_ratio = cand_chi_ratio
end function get_chi_ratio

subroutine get_final_vars(i_opt, j_opt, leftover_ox, OMIX, FMIX, chi_xi, chi_eta, x_prime, Z_val, eta, xi_stoic)
    use precision
    implicit none
    integer, intent(in) :: i_opt, j_opt
    real(WP), dimension(3), intent(in) :: leftover_ox
    real(WP), dimension(3), intent(out) :: OMIX, FMIX
    real(WP), intent(out) :: chi_xi, chi_eta, x_prime, Z_val, eta, xi_stoic

    real(WP) :: a
    real(WP) :: leftover_ox_left, leftover_ox_right
    real(WP), dimension(3) :: OMIX_shuffled, FMIX_shuffled
    real(WP) :: x_prime_left, x_prime_right
    logical :: complement_xi

    call compute_intermediate_vars()
    chi_xi = cand_chi_xi
    chi_eta = cand_chi_eta
    Z_val = cand_xi
    eta = cand_eta
    a = x*b
    
    ! Computation of OMIX, FMIX, x_prime, xi_stoic below:
    leftover_ox_left = a*leftover_ox(1) + (1.0_WP-a)*leftover_ox(3)
    leftover_ox_right = b*leftover_ox(1) + (1.0_WP-b)*leftover_ox(2)

    ! Assigns OMIX and FMIX by determining which endpoint is fuel-lean or fuel-rich
    if (leftover_ox_left >= -TOL .and. leftover_ox_right <= TOL) then   ! Robust (but necessary, don't use 0.0_WP)
        OMIX_shuffled = [a, 0.0_WP, 1.0_WP-a]
        FMIX_shuffled = [b, 1.0_WP-b, 0.0_WP]
        complement_xi = .false.
    else if (leftover_ox_left <= TOL .and. leftover_ox_right >= -TOL) then   ! Robust (but necessary, don't use 0.0_WP)
        OMIX_shuffled = [b, 1.0_WP-b, 0.0_WP]
        FMIX_shuffled = [a, 0.0_WP, 1.0_WP-a]
        complement_xi = .true.
    else
        print *, "ERROR, Both boundary conditions are fuel-lean or both are fuel-rich, so can't assign OMIX and FMIX identities."
    end if

    OMIX(i_opt) = OMIX_shuffled(1)
    FMIX(i_opt) = FMIX_shuffled(1)
    OMIX(j_opt) = OMIX_shuffled(2)
    FMIX(j_opt) = FMIX_shuffled(2)
    OMIX(1+2+3-i_opt-j_opt) = OMIX_shuffled(3)
    FMIX(1+2+3-i_opt-j_opt) = FMIX_shuffled(3)
    
    ! print *, "VALUES FOR COMPARING", cand_Z1, ",", cand_Z2, ",", cand_chi11, ",", cand_chi12, ",", cand_chi22, ",", x
    xi_stoic = leftover_ox_left / (leftover_ox_left - leftover_ox_right)
    if (complement_xi) then
        Z_val = 1.0_WP - Z_val
        xi_stoic = 1.0_WP - xi_stoic
    end if

    if (i_opt == 1 .and. j_opt == 2) then
        x_prime_left = 3.0_WP
        x_prime_right = 4.0_WP
    else if (i_opt == 1 .and. j_opt == 3) then
        x_prime_left = 5.0_WP
        x_prime_right = 4.0_WP
    else if (i_opt == 2 .and. j_opt == 1) then
        x_prime_left = 3.0_WP
        x_prime_right = 2.0_WP
    else if (i_opt == 2 .and. j_opt == 3) then
        x_prime_left = 1.0_WP
        x_prime_right = 2.0_WP
    else if (i_opt == 3 .and. j_opt == 1) then
        x_prime_left = 5.0_WP
        x_prime_right = 6.0_WP
    else if (i_opt == 3 .and. j_opt == 2) then
        x_prime_left = 1.0_WP
        x_prime_right = 0.0_WP
    else
        print *, "ERROR, stream permutation invalid."
    end if
    x_prime = x_prime_left + x * (x_prime_right - x_prime_left)

    if (.not. (eq(sum(OMIX), 1.0_WP) .and. eq(sum(FMIX), 1.0_WP))) print *, "OMIX/FMIX sum issue", OMIX, FMIX
    if (min(minval(OMIX), minval(FMIX)) < -TOL) print *, "OMIX/FMIX negative values", OMIX, FMIX
    if (min(chi_xi, chi_eta) < -TOL) print *, "Negative chi", chi_xi, chi_eta
    if (.not. (in_range(x_prime/6.0_WP) .and. in_range(x) .and. in_range(Z_val) .and. in_range(eta) &
        .and. in_range(xi_stoic))) then
        print *, "Not in range", x_prime, x, Z_val, eta, xi_stoic
    end if
    ! Robust:
    ! chi_xi = max(chi_xi, 0.0_WP)
    ! chi_eta = max(chi_eta, 0.0_WP)
    ! Z_val = max(min(Z_val, 1.0_WP), 0.0_WP)
    ! eta = max(min(eta, 1.0_WP), 0.0_WP)
    ! xi_stoic = max(min(xi_stoic, 1.0_WP), 0.0_WP)    
end subroutine get_final_vars

logical function eq(num1, num2)
    use precision
    implicit none
    real(WP), intent(in) :: num1, num2

    eq = abs(num1 - num2) < TOL
end function eq

logical function in_range(num)
    use precision
    implicit none
    real(WP), intent(in) :: num

    in_range = -TOL < num .and. num < 1.0_WP + TOL
end function in_range

subroutine perturb_x()
    use precision
    implicit none
    ! x = mod(x + 0.0000001_WP, 1.0_WP)   ! Robust
end subroutine perturb_x

end module domain_candidate