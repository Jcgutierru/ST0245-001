# -*- coding: utf-8 -*-
"""
Created on Tue Jun  8 00:22:00 2021

@author: User
"""
# -*- coding: utf-8 -*-

import numpy as np

class Point:
    def __init__(self,x,y,z):
        self.x = x
        self.y = y
        self.z = z

class Rectangulo:
    def __init__(self,x,y,z,w,h,d):
        self.x = x
        self.y = y
        self.z = z
        self.w = w
        self.h = h
        self.d = d

    def contains(self, point):
        return (point.x > self.x - self.w/2 and point.x < self.x + self.w/2 and point.y > self.y - self.h/2 and point.y < self.y + self.h/2 and point.z > self.z - self.d/2 and point.z < self.z + self.d/2)


class OcTree:
    def __init__(self,boundary,n):
        self.boundary = boundary
        self.capacity = n
        self.points = []
        self.divide = False


    def insert(self,point):
        if not self.boundary.contains(point):
            return

        if len(self.points) < self.capacity and self.divide == False:
            self.points.append(point)

        else:
            if not self.divide:
                self.subdivide()
                self.divide = True
            
            for i in self.points:
                self.NE1.insert(i)
                self.NO1.insert(i)
                self.SE1.insert(i)
                self.SO1.insert(i)
                self.NE2.insert(i)
                self.NO2.insert(i)
                self.SE2.insert(i)
                self.SO2.insert(i)
            
            self.remove_points()
            
            self.NE1.insert(point)
            self.NO1.insert(point)
            self.SE1.insert(point)
            self.SO1.insert(point)
            self.NE2.insert(point)
            self.NO2.insert(point)
            self.SE2.insert(point)
            self.SO2.insert(point)
            
            
    def remove_points(self):
        self.points = []
            
    def subdivide(self):    
        ne1_boundary = Rectangulo(self.boundary.x + self.boundary.w/4, self.boundary.y + self.boundary.h/4, self.boundary.z + self.boundary.d/4, self.boundary.w/2,self.boundary.h/2,self.boundary.d/2)
        self.NE1 = OcTree(ne1_boundary,self.capacity)

        no1_boundary = Rectangulo(self.boundary.x - self.boundary.w/4, self.boundary.y + self.boundary.h/4 , self.boundary.z + self.boundary.d/4,self.boundary.w/2,self.boundary.h/2,self.boundary.d/2)
        self.NO1 = OcTree(no1_boundary,self.capacity)

        se1_boundary = Rectangulo(self.boundary.x + self.boundary.w/4, self.boundary.y - self.boundary.h/4 , self.boundary.z + self.boundary.d/4,self.boundary.w/2,self.boundary.h/2,self.boundary.d/2)
        self.SE1 = OcTree(se1_boundary,self.capacity)

        so1_boundary = Rectangulo(self.boundary.x - self.boundary.w/4, self.boundary.y - self.boundary.h/4 , self.boundary.z + self.boundary.d/4,self.boundary.w/2,self.boundary.h/2,self.boundary.d/2)
        self.SO1 = OcTree(so1_boundary,self.capacity)
        
        ne2_boundary = Rectangulo(self.boundary.x + self.boundary.w/4, self.boundary.y + self.boundary.h/4, self.boundary.z - self.boundary.d/4, self.boundary.w/2,self.boundary.h/2,self.boundary.d/2)
        self.NE2 = OcTree(ne2_boundary,self.capacity)

        no2_boundary = Rectangulo(self.boundary.x - self.boundary.w/4, self.boundary.y + self.boundary.h/4 , self.boundary.z - self.boundary.d/4,self.boundary.w/2,self.boundary.h/2,self.boundary.d/2)
        self.NO2 = OcTree(no2_boundary,self.capacity)

        se2_boundary = Rectangulo(self.boundary.x + self.boundary.w/4, self.boundary.y - self.boundary.h/4 , self.boundary.z - self.boundary.d/4,self.boundary.w/2,self.boundary.h/2,self.boundary.d/2)
        self.SE2 = OcTree(se2_boundary,self.capacity)

        so2_boundary = Rectangulo(self.boundary.x - self.boundary.w/4, self.boundary.y - self.boundary.h/4 , self.boundary.z - self.boundary.d/4,self.boundary.w/2,self.boundary.h/2,self.boundary.d/2)
        self.SO2 = OcTree(so2_boundary,self.capacity)

    def numpts(self):
        if self.divide == True:
            self.numpts_hijos(self.NE1)
            self.numpts_hijos(self.NE2)
            self.numpts_hijos(self.NO1)
            self.numpts_hijos(self.NO2)
            self.numpts_hijos(self.SE1)
            self.numpts_hijos(self.SE2)
            self.numpts_hijos(self.SO1)
            self.numpts_hijos(self.SO2)

        else:
            return len(self.points)

    def numpts_hijos(self,hijo,cont=0):
        if self.divide == True:
            self.numpts_hijos(self.NE1,cont)
            self.numpts_hijos(self.NE2,cont)
            self.numpts_hijos(self.NO1,cont)
            self.numpts_hijos(self.NO2,cont)
            self.numpts_hijos(self.SE1,cont)
            self.numpts_hijos(self.SE2,cont)            
            self.numpts_hijos(self.SO1,cont)
            self.numpts_hijos(self.SO2,cont)

        else:
            cont = cont + len(self.points)
            return cont
        