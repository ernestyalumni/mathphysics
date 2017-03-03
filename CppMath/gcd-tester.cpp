/* gcd-tester.cpp
 * cf. Chapter 3 Greatest Common Divisor
 * Scheinerman, C++ for Mathematicians Introduction
 * 
 * Compiler tip - I did this:
 * g++ gcd.cpp gcd-tester.cpp -o gcd-tester
 * or this
 * g++ gcd-recursive.cpp gcd-tester.cpp -o gcd-recursive-tester
 * or this
 * g++ gcd-iter.cpp gcd-tester.cpp -o gcd-iter-tester
 * 
 * */
#include "gcd.h"
#include <iostream>
using namespace std;

/**
 * A program to test the gcd procedure
 * */
int main() {
	long a,b;
	
	cout << "Enter the first number  --> ";
	cin >> a;
	cout << "Enter the second number --> ";
	cin >> b;
	
	cout << "The gcd of " << a << " and " << b << " is "
		<< gcd(a,b) << endl;
	return 0;
}
	
