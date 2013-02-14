#!/usr/bin/env python
from generic_scan import *

def op_plus(a,b):
	return a + b
	
def op_cross(x,f):
	return 0 if f==1 else x
	
def op_companion(x,f):
	return (x or f)
	
x = [5,2,3]
mat = [[0, 3, 1],[2, 2, 0],[0, 0, 1]]
val = [3, 1, 2, 2, 1]
ind = [1, 2, 0, 1, 2]
ptr = [0, 2, 4, 5]

f = [1 if i in ptr else 0 for i in range(len(val))]

a = [val[i]*x[ind[i]] for i in range(len(val))]

print a
print f

s = Scan(f, a, op_plus, op_cross, op_companion)
c = s.scan()
print c
print [c[ptr[i]-1] for i in range(1,len(ptr))]

