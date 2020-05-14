# -*- coding:utf-8 -*-
import sys
import math
import numpy as np
import random

LAT_CONST = 2.855324 
POISSON = 0.337253707481338
STRAIN_Y = 0.000
STRAIN_X = -POISSON * STRAIN_Y
STRAIN_Z = -POISSON * STRAIN_Y

DIST_X = LAT_CONST 
DIST_Y = LAT_CONST 
DIST_Z = LAT_CONST 

num =  20 
N_X = 304 
N_X_D = 190
N_Y = num
N_Z = num 
n_atoms = N_X *N_Y *N_Z *2
n_arrange = N_X_D *N_Y *N_Z *2
RATE = 0.000
n_vacancy = math.floor(n_arrange * RATE)
cut_off = 5.6

print(n_vacancy)

vert = [0.,0.,0.]
surf1 = [0.5*LAT_CONST,0.5*LAT_CONST,0.5*LAT_CONST]
#surf2 = [0.5*LAT_CONST,0.5*LAT_CONST,0.]
#surf3 = [0.5*LAT_CONST,0.,0.5*LAT_CONST]

size_x = (DIST_X * N_X) * (1.0 + STRAIN_X)
size_y = (DIST_Y * N_Y) * (1.0 + STRAIN_Y)
size_z = (DIST_Z * N_Z) * (1.0 + STRAIN_Z)

perfect =np.zeros((n_atoms, 8)) # atom_number x y z grouping_num TF inv_y inv_z

def make_perfect():

    count = 0
    # make perfect lattice
    for ix in range(N_X):
        for iy in range(N_Y):
            for iz in range(N_Z): 
                if ix==0:                   
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,2 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=1      
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,2 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=1
                elif ix==N_X_D-1:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,3 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1 
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,3 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                elif ix < N_X_D:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,1 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,1 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                else:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,4 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,4 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1

def center(x):
    if x%2==0:
        center = x/2
    else:
        center = (x-1)/2
    return center

def arrange_lat(R,Vac_Cu):
    del_num = 0   
    for j in range(n_atoms): # any big value is ok
        if del_num >= n_vacancy:
            break
        p = random.randint(1,n_arrange)        
        if R == 0:
            if perfect[p,4] ==1 and perfect[p,1] >= (cut_off) and perfect[p,1] <= (N_X_D*LAT_CONST -(cut_off)):
                for k in range(n_arrange):
                    if perfect[k,5] ==0 and Vac_Cu ==1 :          # if Vac [k,5] ==0
                        if neighbor(p,k,R,cut_off) == 1:
                            break
                    if perfect[k,4] ==5 and Vac_Cu ==2 :          # if Cu [k,4] ==5
                        if neighbor(p,k,R,cut_off) == 1:
                            break
                    if k == n_arrange-1:
                        if Vac_Cu ==1:    
                            perfect[p,5] =0           # if Vac[p,5] ==0
                            del_num +=1 
                            print(del_num)
                        if Vac_Cu ==2:
                            perfect[p,4] =5           # if Cu [p,4] ==5
                            del_num +=1
                            print(del_num)      
        else:
            delta = 0.05
            if perfect[p,4] ==1 and perfect[p,1] >= N_X_D*LAT_CONST*(0.5-delta) and perfect[p,1] <= N_X_D*LAT_CONST*(0.5+delta):
                for k in range(n_arrange):
                    if perfect[k,5] ==0 and Vac_Cu ==1 :          # if Vac [k,5] ==0
                        if neighbor(p,k,R,cut_off) == 1:
                            break
                    if perfect[k,4] ==5 and Vac_Cu ==2 :          # if Cu [k,4] ==5
                        if neighbor(p,k,R,cut_off) == 1:
                            break
                    if k == n_arrange-1:
                        for v in range(n_arrange):
                            if neighbor(p,v,R,0) ==1:
                                if Vac_Cu ==1:
                                    perfect[v,5] = 0  # if Vac [v,5] ==0
                                    del_num +=1 
                                    print(del_num) 
                    if k == n_arrange-1:
                        for v in range(n_arrange):
                            if neighbor(p,v,R,0) ==1:
                                if Vac_Cu ==2:        
                                    perfect[v,4] = 5  # if Cu [v,4] ==5
                                    del_num +=1 
                                    print(del_num) 
    
    name(Vac_Cu+3)
                   
def name(type_num):
    f = open("Fe.lat", "w")
    f.write("# Fe\n")
    f.write("\n")
    #f.write("%d atoms\n" % perfect_num)
    f.write("\n")
    f.write("%d atom types\n" % type_num)
    f.write("0 %12.6f xlo xhi\n" % size_x)
    f.write("0 %12.6f ylo yhi\n" % size_y)
    f.write("0 %12.6f zlo zhi\n" % size_z)
    f.write("\n")
    f.write("Atoms\n")
    f.write("\n")
    lat =0
    for i in range(n_atoms):
        if perfect[i,5] ==1 :
            lat +=1
            f.write("%12d %6d %12.6f %12.6f %12.6f\n" % (lat, perfect[i,4],perfect[i,1],perfect[i,2],perfect[i,3]))    
    f.close()

    with open("Fe.lat") as f:
        l = f.readlines()

    l.insert(2, "%d atoms" % lat)
    with open("Fe.lat", mode='w') as f:
        f.writelines(l)


def inverse_lat():
    f = open("inverse.lat", "w")
    lat =0
    for i in range(perfect_num): 
        if perfect[i,5] ==0 : # if Cu [i,4] ==5
            lat +=1
            f.write("%12d %6d %12.6f %12.6f %12.6f\n" % (lat, perfect[i,4],perfect[i,1],perfect[i,2],perfect[i,3]))    
    f.close()


def neighbor(p,k,R,cut):
     if distance(p,k,1)**2 + distance(p,k,2)**2 + distance(p,k,3)**2 <= (2*cut+R)**2:
         return 1
     else: 
         return 0

def distance(p,k,a):
    if a == 1:
        return abs(perfect[p,1]-perfect[k,1])
    elif a == 2:
        b1 = abs(perfect[p,2]-perfect[k,2])
        b2 = abs(perfect[p,2]-perfect[k,6])
        b3 = abs(perfect[p,6]-perfect[k,2])
        #theorically b4 is not needed
        return min(b1,b2,b3)
        
    else:
        b1 = abs(perfect[p,3]-perfect[k,3])
        b2 = abs(perfect[p,3]-perfect[k,7])
        b3 = abs(perfect[p,7]-perfect[k,3])
        #theorically b4 is not needed
        return min(b1,b2,b3)
            
            
            
if __name__ == "__main__":
    argvs = sys.argv
    make_perfect()
    arrange_lat(0,1)  #(R ,type) type= vacancy:1 , Cu: 2
                      #and RATE or num (=N_X,N_Y)are all parameters
                      

################moving test
#check num of arrage
#if you do not need, delete 'dump'
#
#Use OVITO to check molphorogy
#FFT configuration check
#arraving time should be defined consistnatly
#change x_detec, A2, Time to get beta
