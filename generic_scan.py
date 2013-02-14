#!/usr/bin/env python
import math

class Scan:
	def __init__(self, a, b, plus, cross = None, companion = None):
		self.a = a
		self.b = b
		self.c = [[a[i],b[i]] for i in range(len(a))]
		self.plus = plus
		self.cross = cross
		self.companion = companion if companion else cross
		
	def __op_point(self, z1, z2):
		if self.cross:
			return [self.companion(z1[0],z2[0]), self.plus(self.cross(z1[1],z2[0]), z2[1])]
		else:
			return self.plus(z1,z2)
	
	def __up_sweep(self):
		for d in range(int(math.log(len(self.c),2))):
			p = int(math.pow(2,d))
			i = p-1
			while i+p < len(self.c):
				self.c[i+p] = self.__op_point(self.c[i],self.c[i+p])
				i += p*2

	def __down_sweep(self):
		for d in reversed(range(1,int(math.log(len(self.c),2)))):
			p = int(math.pow(2,d))
			i = p-1
			p /= 2
			while i+p < len(self.c):
				self.c[i+p] = self.__op_point(self.c[i],self.c[i+p])
				i += p*2
		
	def scan(self):
		if len(self.a) != len(self.b):
			return None
		self.__up_sweep()
		self.__down_sweep()
		return [v[1] for v in self.c]
	

	
