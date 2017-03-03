/* gcd-prob2.cpp
 * cf. Chapter 3 Greatest Common Divisor
 * Scheinerman, C++ for Mathematicians Introduction
 * 
 * 
 * pp. 43  3.5 An exhaustive approach to the GCD problem
 * Problem 3.10: A slightly better program to calculate p_n
 * 
 * 
 * Compiling tip - I did this : 
 * g++ gcd-iter.cpp gcd-prob2.cpp -o gcd-prob2
 */
#include <iostream>
#include "gcd.h"
using namespace std;

/**
 * Find the probability that two integers in {1,...,n} are relatively 
 * prime.
 * */
 
int main() {
	long n;
	cout << "Enter n --> ";
	cin >> n;
	
	long count = 0;
	
	for (long a=1; a<=n; a++) {
		for (long b=a+1; b<=n; b++) {
			if (gcd(a,b) == 1) {
				count++;
			}
		}
	}
	// b=a+1, so inner loop runs from a+1 up to n; Thus program calculates gcd(5,10), but not gcd(10,5)
	// we need to double count at the end to correct.  Also, we do not calculate gcd(t,t) for any t,
	// so we add 1 to (the doubled) count at the end to correct for that.
	count = 2*count + 1;
	
	
	cout << double(count) / (double(n) * double(n)) << endl;
	return 0;
}
 
