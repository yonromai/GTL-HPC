#!/usr/bin/env python
import math
import generic_scan as scan
import pdb

def bin_to_dec(bin_array):
	num = 0
	for bit in reversed(bin_array):
		num = (num << 1) | bit;
	return num
	
def dec_to_bin(dec_int):
	bin_array = []
	while dec_int > 0:
		bin_array.append(dec_int & 1)
		dec_int = dec_int >> 1
	return bin_array

def get_x(a,b,i):
	return a[i]^b[i]
	
def get_y(a,b,i):
	return a[i] & b[i]

def get_z(a,b,i):
	return [get_x(a,b,i), get_y(a,b,i)]

def op_plus(a,b):
	return a | b
	
def op_cross(a,b):
	return a & b

def is_pow_2(l):
	r = math.log(l,2)
	return math.floor(r) == math.ceil(r)

def get_carries(a,b):
	if len(a) != len(b):
		print "get_carries: Not same size"
		throw
	while not is_pow_2(len(a)):
		a.append(0)
		b.append(0)
	z=[]
	for i in range(len(a)):
		z.append(get_z(a,b,i))

	scan.up_sweep(z,op_plus,op_cross)
	scan.down_sweep(z,op_plus,op_cross)
	return z

for i in range(1,128):
	for j in range(1,128):
		expected = i + j
		a = dec_to_bin(max(i,j))
		b = dec_to_bin(min(i,j))
		for k in range(len(a) - len(b)):
			b.append(0)
			
		print "a:", a
		print "b:", b
		z = [[0,0]]
		z.extend(get_carries(a,b))
		t = []
		for k in range(len(a)):
			t.append(a[k]^b[k])
		t.append(0)
		s = []
		
		#print 'a',[x for x in a]
		for k in range(len(t)):
			s.append(t[k]^z[k][1])
		
		actual = bin_to_dec(s)
		if actual != expected:
			print "You suck bitch (i:",i,"j:",j,"). Was given", actual, "expected", expected
			print 'a',[x for x in a]
			print 'b',[x for x in b]
			print 'c',[x[1] for x in z]
			print 's',[x for x in s]
			throw
		else :
			print "ok"
		 
print bin_to_dec(dec_to_bin(255))

