#!/usr/bin/env python
from generic_scan import *

def op_plus(a,b):
	return a + b
	
def op_cross(x,f):
	return 0 if f==1 else x
	
def op_companion(x,f):
	return (x or f)
	
x = [3,2,-1,5]
mat = [[2, 0, 0, 1],[4, 3, 0, 0],[0, 0, 7, 0],[1, 0, 1, 0]]
val = [2,1,4,3,7,1,1]
ind = [0,3,0,1,2,0,2]
ptr = [0,2,4,5,7]

f = [1 if i in ptr else 0 for i in range(len(val))]

a = [val[i]*x[ind[i]] for i in range(len(val))]

print a
print f

s = Scan(f, a, op_plus, op_cross, op_companion)
c = s.scan()
print c
print [c[ptr[i]-1] for i in range(1,len(ptr))]

