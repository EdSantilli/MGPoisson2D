! ------------------------------------------------------------------------------
! Defines the precision of all numbers that we will use.
! ------------------------------------------------------------------------------
module precision
    integer, parameter  :: dp = kind(0.d0)

    real(dp), parameter :: zero      =  0.0_dp
    real(dp), parameter :: one       =  1.0_dp
    real(dp), parameter :: two       =  2.0_dp
    real(dp), parameter :: three     =  3.0_dp
    real(dp), parameter :: four      =  4.0_dp
    real(dp), parameter :: five      =  5.0_dp
    real(dp), parameter :: six       =  6.0_dp
    real(dp), parameter :: seven     =  7.0_dp
    real(dp), parameter :: eight     =  8.0_dp
    real(dp), parameter :: nine      =  9.0_dp
    real(dp), parameter :: ten       = 10.0_dp
    real(dp), parameter :: eleven    = 11.0_dp
    real(dp), parameter :: twelve    = 12.0_dp
    real(dp), parameter :: thirteen  = 13.0_dp
    real(dp), parameter :: fourteen  = 14.0_dp
    real(dp), parameter :: fifteen   = 15.0_dp
    real(dp), parameter :: sixteen   = 16.0_dp
    real(dp), parameter :: seventeen = 17.0_dp
    real(dp), parameter :: eighteen  = 18.0_dp
    real(dp), parameter :: nineteen  = 19.0_dp
    real(dp), parameter :: twenty    = 20.0_dp

    real(dp), parameter :: half        = one / two
    real(dp), parameter :: third       = one / three
    real(dp), parameter :: fourth      = one / four
    real(dp), parameter :: fifth       = one / five
    real(dp), parameter :: sixth       = one / six
    real(dp), parameter :: seventh     = one / seven
    real(dp), parameter :: eighth      = one / eight
    real(dp), parameter :: ninth       = one / nine
    real(dp), parameter :: tenth       = one / ten
    real(dp), parameter :: eleventh    = one / eleven
    real(dp), parameter :: twelfth     = one / twelve
    real(dp), parameter :: thirteenth  = one / thirteen
    real(dp), parameter :: fourteenth  = one / fourteen
    real(dp), parameter :: fifteenth   = one / fifteen
    real(dp), parameter :: sixteenth   = one / sixteen
    real(dp), parameter :: seventeenth = one / seventeen
    real(dp), parameter :: eighteenth  = one / eighteen
    real(dp), parameter :: nineteenth  = one / nineteen
    real(dp), parameter :: twentieth   = one / twenty

    real(dp), parameter :: threehalves  = three / two
    real(dp), parameter :: threefourths = three / four

    real(dp), parameter :: pi        = four * atan(one)
    real(dp), parameter :: bogus_val = 1.2345E300_dp

end module precision
