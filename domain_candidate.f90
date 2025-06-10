module domain_candidate
    use precision
    implicit none

    type :: Domain_Cand
        real(WP) :: cand_Z1, cand_Z2, cand_chi11, cand_chi12, cand_chi22, x
        real(WP) :: b, cand_xi, cand_eta, cand_chi_xi, cand_chi_eta, cand_chi_ratio
    contains
        procedure :: domain_cand_init
        procedure :: get_final_vars
    end type Domain_Cand

    real(WP), parameter :: TOL = 0.000001_WP

contains

subroutine domain_cand_init(this, Z1_in, Z2_in, chi11_in, chi12_in, chi22_in, x_in)
    use precision
    implicit none
    
    class(Domain_Cand), intent(inout) :: this
    real(WP), intent(in) :: Z1_in, Z2_in, chi11_in, chi12_in, chi22_in, x_in
    real(WP) :: cand_Z3, cand_chi_xi_numer, cand_chi_eta_numer, cand_chi_denom
    integer :: attempt
    logical :: success

    this%cand_Z1 = Z1_in
    this%cand_Z2 = Z2_in
    this%cand_chi11 = chi11_in
    this%cand_chi12 = chi12_in
    this%cand_chi22 = chi22_in
    this%x = x_in
    
    success = .false.
    do attempt = 1, 5
        cand_Z3 = 1.0_WP - this%cand_Z1 - this%cand_Z2
        if ((this%cand_Z1+this%cand_Z2+this%x*(this%cand_Z1+cand_Z3))**2.0_WP - 4.0_WP*this%cand_Z1*this%x < -TOL .or. &
            eq(2.0_WP*this%x, 0.0_WP)) then
            call perturb_x(this)
            cycle
        end if
        this%b = (this%cand_Z1+this%cand_Z2+this%x*(this%cand_Z1+cand_Z3) - sqrt((this%cand_Z1+this%cand_Z2+ &
            this%x*(this%cand_Z1+cand_Z3))**2.0_WP - 4.0_WP*this%cand_Z1*this%x)) / (2.0_WP*this%x)
        this%cand_eta = this%b
        if (eq(this%b, 1.0_WP)) then
            call perturb_x(this)
            cycle
        end if
        this%cand_xi = this%cand_Z2 / (1.0_WP - this%b)
        cand_chi_eta_numer = (1.0_WP-this%cand_eta)**2.0_WP*this%cand_chi11 - 2.0_WP*this%cand_eta*(1.0_WP-this%cand_eta)* &
            (1.0_WP-this%x)*this%cand_chi12 + (this%cand_eta*(1.0_WP-this%x))**2.0_WP*this%cand_chi22
        cand_chi_xi_numer = this%cand_xi**2.0_WP*this%cand_chi11 + 2.0_WP*this%cand_xi*(this%x+(1.0_WP-this%x)*this%cand_xi)* &
            this%cand_chi12 + (this%x+(1.0_WP-this%x)*this%cand_xi)**2.0_WP*this%cand_chi22
        cand_chi_denom = ((1.0_WP-this%cand_eta)*(this%x+(1.0_WP-this%x)*this%cand_xi) + this%cand_eta*this%cand_xi* &
            (1.0_WP-this%x))**2.0_WP
        if (eq(cand_chi_denom, 0.0_WP) .or. eq(cand_chi_xi_numer, 0.0_WP)) then
            call perturb_x(this)
            cycle
        end if
        this%cand_chi_eta = cand_chi_eta_numer / cand_chi_denom
        this%cand_chi_xi = cand_chi_xi_numer / cand_chi_denom
        this%cand_chi_ratio = cand_chi_eta_numer / cand_chi_xi_numer
        success = .true.
    end do

    if (.not. success) this%cand_chi_ratio = -1.0_WP   ! So candidates with NaN issues are not selected

end subroutine domain_cand_init

subroutine get_final_vars(this, i_opt, j_opt, leftover_ox, OMIX, FMIX, chi_xi, chi_eta, x_prime, Z_val, eta, Z_stoic)
    use precision
    implicit none
    class(Domain_Cand), intent(inout) :: this
    integer, intent(in) :: i_opt, j_opt
    real(WP), dimension(3), intent(in) :: leftover_ox
    real(WP), dimension(3), intent(out) :: OMIX, FMIX
    real(WP), intent(out) :: chi_xi, chi_eta, x_prime, Z_val, eta, Z_stoic

    real(WP) :: a
    real(WP) :: leftover_ox_left, leftover_ox_right
    real(WP), dimension(3) :: OMIX_shuffled, FMIX_shuffled
    real(WP) :: x_prime_left, x_prime_right
    logical :: complement_xi
    real(WP) :: sum_Z1_Z2

    if (this%cand_chi_ratio > -0.5_WP) then
        chi_xi = this%cand_chi_xi
        chi_eta = this%cand_chi_eta
        Z_val = this%cand_xi
        eta = this%cand_eta
        a = this%x*this%b
    else
        if (.not. (i_opt == 1 .and. j_opt == 3 .and. eq(this%x, 0.0_WP))) print *, "Major logic error"
        sum_Z1_Z2 = this%cand_Z1 + this%cand_Z2
        if (sum_Z1_Z2 == 0.0_WP) sum_Z1_Z2 = 0.00000001_WP   ! To prevent division by 0
        Z_val = sum_Z1_Z2   ! HACK: For non-oxygen-rich second streams, may need to be complement of this
        eta = this%cand_Z2 / sum_Z1_Z2
        chi_xi = this%cand_chi11 + 2.0_WP*this%cand_chi12 + this%cand_chi22
        chi_eta = (this%cand_chi11 + 2.0_WP*(1.0_WP-eta)*(-this%cand_chi11-this%cand_chi12) + &
            (1.0_WP-eta)**2.0_WP*chi_xi) / (sum_Z1_Z2**2.0_WP)
        a = 0.0_WP
    end if
    
    ! Computation of OMIX, FMIX, x_prime, Z_stoic below:
    leftover_ox_left = a*leftover_ox(1) + (1.0_WP-a)*leftover_ox(3)
    leftover_ox_right = this%b*leftover_ox(1) + (1.0_WP-this%b)*leftover_ox(2)

    ! Assigns OMIX and FMIX by determining which endpoint is fuel-lean or fuel-rich
    if (leftover_ox_left >= -TOL .and. leftover_ox_right <= TOL) then   ! Robust (but necessary, don't use 0.0_WP)
        OMIX_shuffled = [a, 0.0_WP, 1.0_WP-a]
        FMIX_shuffled = [this%b, 1.0_WP-this%b, 0.0_WP]
        complement_xi = .false.
    else if (leftover_ox_left <= TOL .and. leftover_ox_right >= -TOL) then   ! Robust (but necessary, don't use 0.0_WP)
        OMIX_shuffled = [this%b, 1.0_WP-this%b, 0.0_WP]
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
    
    ! print *, "VALUES FOR COMPARING", this%cand_Z1, ",", this%cand_Z2, ",", this%cand_chi11, ",", &
    !     this%cand_chi12, ",", this%cand_chi22, ",", this%x
    Z_stoic = leftover_ox_left / (leftover_ox_left - leftover_ox_right)
    if (complement_xi) then
        Z_val = 1.0_WP - Z_val
        Z_stoic = 1.0_WP - Z_stoic
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
    x_prime = x_prime_left + this%x * (x_prime_right - x_prime_left)

    if (.not. (eq(sum(OMIX), 1.0_WP) .and. eq(sum(FMIX), 1.0_WP))) print *, "OMIX/FMIX sum issue", OMIX, FMIX
    if (min(minval(OMIX), minval(FMIX)) < -TOL) print *, "OMIX/FMIX negative values", OMIX, FMIX
    if (min(chi_xi, chi_eta) < -TOL) print *, "Negative chi", chi_xi, chi_eta
    if (.not. (in_range(x_prime/6.0_WP) .and. in_range(this%x) .and. in_range(Z_val) .and. in_range(eta) &
        .and. in_range(Z_stoic))) then
        print *, "Not in range", x_prime, this%x, Z_val, eta, Z_stoic
    end if
    ! Robust:
    chi_xi = max(chi_xi, 0.0_WP)
    chi_eta = max(chi_eta, 0.0_WP)
    Z_val = max(min(Z_val, 1.0_WP), 0.0_WP)
    eta = max(min(eta, 1.0_WP), 0.0_WP)
    Z_stoic = max(min(Z_stoic, 1.0_WP), 0.0_WP)    
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

subroutine perturb_x(this)
    use precision
    implicit none
    class(Domain_Cand), intent(inout) :: this

    this%x = mod(this%x + 0.0000001_WP, 1.0_WP)   ! Robust
end subroutine perturb_x

end module domain_candidate