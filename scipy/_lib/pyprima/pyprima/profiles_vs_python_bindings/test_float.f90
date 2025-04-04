program test_float
    implicit none
    
    real(8) :: x(4)
    real(8) :: y(4)
    real(8) :: result
    integer(4) :: i
    
   

    print *, float2bin(epsilon(0.0_8)**2), result

contains

function inprod(x, y) result(z)
    !--------------------------------------------------------------------------------------------------!
    ! INPROD calculates the inner product of X and Y, i.e., Z = X^T*Y, regarding X and Y as columns.
    !--------------------------------------------------------------------------------------------------!
    implicit none
    
    ! Inputs
    real(8), intent(in) :: x(:)
    real(8), intent(in) :: y(:)
    ! Outputs
    real(8) :: z
    ! Local variables
    integer(4) :: i
    

    
    !====================!
    ! Calculation starts !
    !====================!
    
    z = 0.0_8
    do i = 1, int(size(x), kind(i))
        z = z + x(i) * y(i)
    end do
    
    !====================!
    !  Calculation ends  !
    !====================!
    end function inprod

function bin2float(binary_string) result(result)
    implicit none

    character(len=*), intent(in) :: binary_string
    real(8) :: result

    integer(8) :: binary_value
    integer :: i

    ! Convert binary string to integer
    binary_value = 0
    do i = 1, 64
        binary_value = binary_value * 2
        if (binary_string(i:i) == '1') then
            binary_value = binary_value + 1
        end if
    end do

    ! Transfer the bit pattern to real(8)
    result = transfer(binary_value, result)
end function bin2float

function float2bin(r) result(bin_str)
    real(8), intent(in) :: r
    character(len=64) :: bin_str
    integer(8) :: int_bits
    integer :: j
    
    ! Transfer real bits to integer
    int_bits = transfer(r, int_bits)
    
    ! Convert to binary string
    bin_str = ''
    do j = 64, 1, -1
        if (btest(int_bits, j-1)) then
            bin_str(65-j:65-j) = '1'
        else
            bin_str(65-j:65-j) = '0'
        end if
    end do
end function float2bin

    
end program test_float
