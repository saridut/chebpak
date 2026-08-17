!> Predefined constants, including math constants.
! TODO: Replace the names with equivalent macro names in glibc math.h
! 
module constants_m

use, intrinsic ::  iso_fortran_env, only: int32, int64, real64

implicit none

!> Default integer precision
integer, parameter :: ip = int32

!> Default long integer precision
integer, parameter :: ip_long = int64

!> Default real precision
integer, parameter :: rp = real64

!> Size of a char in bytes
integer, parameter :: sizeof_char = 1

!> Size of a default int in bytes
integer, parameter :: sizeof_int  = 4

!> Size of a default long int in bytes
integer, parameter :: sizeof_long_int  = 8

!> Size of a default real in bytes
integer, parameter :: sizeof_real = 8

!> `One third` : \( \frac{1}{3} \)
real (rp), parameter :: math_third = 0.333333333333333_rp

!> `pi` : \( \pi \)
real (rp), parameter :: math_pi = 3.1415926535897931_rp

!> `pi/2` : \( \pi / 2 \)
real (rp), parameter :: math_pi_2 = 1.5707963267948966_rp

!> `pi/4` : \( \pi / 4 \)
real (rp), parameter :: math_pi_4 = 0.78539816339744828_rp

!> Reciprocal of `pi` : \( 1/\pi \)
real (rp), parameter :: math_1_pi = 0.31830988618379069_rp

!> Two times reciprocal of `pi` : \( 2/\pi \)
real (rp), parameter :: math_2_pi = 0.63661977236758138_rp

!>Two times the reciprocal of the square root of `pi` : \( 2/\sqrt{\pi} \)
real (rp), parameter :: math_2_sqrtpi = 1.1283791670955126_rp

!>Square root of `two` : \(\sqrt{2}\)
real (rp), parameter :: math_sqrt2 = 1.4142135623730951_rp

!>Cube root of `two` : \(\sqrt[3]{2}\)
real (rp), parameter :: math_cbrt2 = 1.2599210498948732_rp

!>Reciprocal of the square root of two : \(1/\sqrt{2}\)
real (rp), parameter :: math_sqrt1_2 = 0.70710678118654746_rp

!>Square root of three : \(\sqrt{3}\)
real (rp), parameter :: math_sqrt3 = 1.7320508075688772_rp

!>Square root of `e` : \(\sqrt{e}\)
real (rp), parameter :: math_sqrt_e = 1.6487212707001282_rp

!>Square root of `pi` : \(\sqrt{\pi}\)
real (rp), parameter :: math_sqrt_pi = 1.7724538509055159_rp

!>The base of natural logarithm : \( e \)
real (rp), parameter :: math_e = 2.7182818284590451_rp

!>The logarithm of `e` to base two : \( \log_2 e \)
real (rp), parameter :: math_log2e = 1.4426950408889634_rp

!>The logarithm of `e` to base 10 : \( \log_{10} e \)
real (rp), parameter :: math_log10e = 0.43429448190325182_rp

!>The natural logarithm of two : \( \ln 2\)
real (rp), parameter :: math_ln2 = 0.69314718055994529_rp

!>The natural logarithm of 10 : \( \ln 10 \)
real (rp), parameter :: math_ln10 = 2.3025850929940459_rp

end module constants_m
