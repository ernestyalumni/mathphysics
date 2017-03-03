/* gcd-iter.cpp
 * cf. Chapter 3 Greatest Common Divisor
 * Scheinerman, C++ for Mathematicians Introduction
 * pp. 40 Program 3.8: An iterative procedure for gcd
 * */
#include "gcd.h"
#include <iostream>
using namespace std;

long gcd(long a, long b) {
	// Make sure a and b are both nonnegative
	if (a<0) a = -a;
	if (b<0) b = -b;
	
	// if a and b are both zero, print an error and return 0
	if ( (a==0) && (b==0) ) {
		cerr << "WARNING: gcd called with both arguments equal to zero."
			<< endl;
		return 0;
	}
	
	long new_a, new_b; // place to hold new versions of a and b
	
	/*
	 * We use the fact that gcd(a,b) = gcd(b,c) where c = a%b.  
	 * Note that if b is zero, gcd(a,b) = gcd(a,0) = a.  If a is zero,
	 * and b is not, we get a%b equal to zero, so new_b will be zero,
	 * hence b will be zero and the loop will exit with a == 0, which
	 * is what we want.
	 * */
	 
	 while (b != 0) {
		 new_a = b;
		 new_b = a%b;
		 
		 a = new_a;
		 b = new_b;
	}
	
	return a;
}
