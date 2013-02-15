#!/usr/bin/env python
from generic_scan import *

def op_plus(a,b):
	return a+b

a = [1]*61
print a
s = Scan.simple_scan(a, op_plus)
print s.scan()
