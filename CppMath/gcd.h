/** \file gcd.h
 * \brief GCD of 2 integers
 * cf. Chapter 3 Greatest Common Divisor
 * Scheinerman, C++ for Mathematicians Introduction
 * */
#ifndef __GCD_H__
#define __GCD_H__

/**
 * Calculate the greatest common divisor of two integers.
 * Note: gcd(0,0) will return 0 and print an error message.
 * @param a the first integer
 * @param b the second integer
 * @return the greatest common divisor of a and b
 * */

long gcd(long a, long b);

#endif // __GCD_H__  this line and previous are the mechanism to prevent double inclusion


