/* gcd-recursive.cpp
 * cf. Chapter 3 Greatest Common Divisor
 * Scheinerman, C++ for Mathematicians Introduction
 * pp. 38 Program 3.7: A recursive procedure for gcd
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
	
	// If b is zero, the answer is a
	if (b==0) return a;
	// If a is zero, the answer is b
	if (a==0) return b;
	
	long c = a%b;
	
	// "The previous (sic original) call to gcd goes "on hold" pending this result of gcd(b,c)
	return gcd(b,c);
}
