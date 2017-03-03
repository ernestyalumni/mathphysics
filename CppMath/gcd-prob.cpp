/* gcd-prob.cpp
 * cf. Chapter 3 Greatest Common Divisor
 * Scheinerman, C++ for Mathematicians Introduction
 * 
 * 
 * pp. 43  3.5 An exhaustive approach to the GCD problem
 * Problem 3.9: A program to calculate p_n
 * 
 * Compiling tip - I did this : 
 * g++ gcd-iter.cpp gcd-prob.cpp -o gcd-prob
 */
#include <iostream>
#include "gcd.h"
using namespace std;

/**
 * Find the probability that 2 integers in {1,...,n} are relatively
 * prime.
 * */

int main() {
	long n ;
	cout << "Enter n --> ";
	cin >> n;
	
	long count = 0;
	
	for (long a=1; a<=n; a++) {
		for (long b=1; b<=n; b++) {
			if (gcd(a,b) == 1) {
				count++;
			}
		}
	}

	cout << double(count) / (double(n) * double(n)) << endl;
	return 0;
}

