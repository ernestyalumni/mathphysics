/* gcd-tester2.cpp
 * cf. Chapter 3 Greatest Common Divisor
 * Scheinerman, C++ for Mathematicians Introduction
 * 
 * Compiler tip - I did this:
 * g++ gcd.cpp gcd-tester2.cpp -o gcd-tester
 * or this
 * g++ gcd-recursive.cpp gcd-tester2.cpp -o gcd-recursive-tester
 * or this
 * g++ gcd-iter.cpp gcd-tester2.cpp -o gcd-iter-tester
 * or this :
 * g++ -std=c++14 gcd-tester2.cpp gcd.cpp -o gcd_naive
 * */
#include "gcd.h"
#include <iostream>


/**
 * A program to test the gcd procedure
 * */
int main() {
	long a,b;
	
	std::cout << "Enter the first number  --> ";
	std::cin >> a;
	std::cout << "Enter the second number --> ";
	std::cin >> b;
	
	std::cout << "The gcd of " << a << " and " << b << " is "
		  << gcd(a,b) << std::endl;
	return 0;
}
	
