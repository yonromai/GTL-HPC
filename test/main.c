#include <stdlib.h>
#include <stdio.h>

typedef double (*bin_operator) (const double a, const double b);

double add(const double d1, const double d2) {
	return d1+d2;
}

double diff(const double d1, const double d2) {
	return d1-d2;
}


void aff(bin_operator op, double d1, double d2) {
	printf("%lf", op(d1,d2));
}

int main (void) {
	
	aff(add,3.5,2.5);
	aff(diff,3.5,2.5);
	return 0;
}
