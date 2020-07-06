#!/usr/bin/env python
import sys,os
import numpy as np
import dss14fpi
import csv

class DSS:
  
  def __init__(self):

    root='./'

    self.io=1
    dss14fpi.root.root=root.ljust(255)
    
    self.pi={0:{},1:{},-1:{}}
  def _get_f(self,z,Q2,storage,func,ic):
    if (z,Q2) not in storage[ic]:
      storage[ic][(z,Q2)]=np.array(func(0,1,ic,self.io,z,Q2))/z
    return storage[ic][(z,Q2)]

  def get_f(self,z,Q2,hadron):
    """
    out: U,UB,D,DB,S,SB,C,B,GL
    """
    if   hadron=='pi+':return self._get_f(z,Q2,self.pi,dss14fpi.fdssh, 1)
    elif hadron=='pi-':return self._get_f(z,Q2,self.pi,dss14fpi.fdssh,-1)
    elif hadron=='pi0':return self._get_f(z,Q2,self.pi,dss14fpi.fdssh, 0)

if __name__=='__main__':

  dss=DSS()
  Q2 = 10.
  with open('output.csv','w') as file:
      file.write('z U UB D DB S SB C B GL')
      file.write('\n')
      for i in range(1,20):
          z = i/20.
          file.write(str(z)+' ')
          for item in dss.get_f(z,Q2,'pi+'):
              file.write(str(item)+' ')
          file.write('\n')

