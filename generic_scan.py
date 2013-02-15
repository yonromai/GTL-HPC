#!/usr/bin/env python
import math

class Scan:
  def __init__(self, a, b, plus, cross, companion):
    self.a = a
    self.simple_rec = (not b)
    self.plus = plus
    self.true_size = len(a)
    self.size = 1 << int(math.ceil(math.log(self.true_size, 2)))
    
    if self.simple_rec:
      self.c = [a[i] if i < self.true_size else 0 for i in range(self.size)]
    else:
      self.b = b
      self.c = [[a[i],b[i] if b else 0] if i < self.true_size else [0,0] for i in range(self.size)]
      self.cross = cross
      self.companion = companion
  
  @classmethod
  def simple_scan(cls, a, plus):
    return cls(a, None, plus, None, None)
    
  @classmethod  
  def general_scan(cls, a, b, plus, cross, companion = None):
    comp = cross if not companion else companion
    return cls(a, b, plus, cross, comp)
  
  @staticmethod
  def get_pow_2(v):
    return math.pow(2,math.floor(math.log(v,2)))
  
  def __op_point(self, z1, z2):
    if self.simple_rec:
      return self.plus(z1,z2)
    else:
      return [self.companion(z1[0], z2[0]), self.plus(self.cross(z1[1], z2[0]), z2[1])] 
  
  def __up_sweep(self):
    for d in range(int(math.log(self.size,2))):
      p = int(math.pow(2,d))
      i = p-1
      while i+p < self.size:
        self.c[i+p] = self.__op_point(self.c[i],self.c[i+p])
        i += p*2

  def __down_sweep(self):
    for d in reversed(range(1,int(math.log(self.size,2)))):
      p = int(math.pow(2,d))
      i = p-1
      p /= 2
      while i+p < self.size:
        self.c[i+p] = self.__op_point(self.c[i],self.c[i+p])
        i += p*2
    
  def scan(self):
    if (not self.simple_rec) and (len(self.a) != len(self.b)):
      return None
    self.__up_sweep()
    self.__down_sweep()
    return self.c[:self.true_size] if self.simple_rec else [v[1] for v in self.c[:self.true_size]]
  

  
