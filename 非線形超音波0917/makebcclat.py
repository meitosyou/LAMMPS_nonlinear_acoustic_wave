# -*- coding:utf-8 -*-
import sys
import math
import numpy as np
import random

#格子定数・ひずみなどの設定
LAT_CONST = 2.855324 
POISSON = 0.293
STRAIN_X = 0.000
STRAIN_Y = -POISSON * STRAIN_X
STRAIN_Z = -POISSON * STRAIN_X

DIST_X = LAT_CONST 
DIST_Y = LAT_CONST 
DIST_Z = LAT_CONST 

#系の大きさの設定
num =  20
N_X =  304
N_X_D = 190
N_Y = num
N_Z = num 
n_atoms = N_X *N_Y *N_Z * 2
n_arrange = N_X_D *N_Y *N_Z * 2

#欠陥の比率を設定
RATE = 0.000
n_vacancy = math.floor(n_arrange * RATE)
cut_off = 5.6

#in.lammpsではx方向にも周期的境界条件を与えているため、非周期系にするためには真空領域を設定する必要がある。indent分左右に空間ができる。多めに設定しておくと良い。
indent = 3  #indent should not be zero

print(n_vacancy)

#系の設定
vert = [0.,0.,0.]
surf1 = [0.5*LAT_CONST,0.5*LAT_CONST,0.5*LAT_CONST]
#surf2 = [0.5*LAT_CONST,0.5*LAT_CONST,0.]
#surf3 = [0.5*LAT_CONST,0.,0.5*LAT_CONST]

size_x = (DIST_X * (N_X+indent*2)) * (1.0 + STRAIN_X)
size_y = (DIST_Y * N_Y) * (1.0 + STRAIN_Y)
size_z = (DIST_Z * N_Z) * (1.0 + STRAIN_Z)


#この関数で完全結晶の格子の原子リストを作成する。
#各原子にはIDとなる原子の番号、xyz座標 原子グループ番号 粒子の存在の有無(TF) 周期境涯越しのyz座標(距離を求めるときに使う)を割り当てている。
#原子グループ番号1:系の左端の原子(振動を導入するsourceにあたるところ)
#                2:detectorの原子
#                3:欠陥や析出物を導入する領域の原子
#                4:それ以外の原子
#                5:系の右端の原子(固定条件と自由条件を比較するために使用した)
#                6:銅原子

perfect =np.zeros((n_atoms, 8)) # atom_number x y z grouping_num TF inv_y inv_z

def make_perfect():

    count = 0
    # make perfect lattice
    for ix in range(indent,N_X+indent):
        for iy in range(N_Y):
            for iz in range(N_Z): 
                if ix==indent:                   
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,2 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=1      
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,2 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=1
                elif ix==N_X_D+indent-1:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,3 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count += 1 
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,3 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                elif ix < N_X_D+indent-1:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,1 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,1 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                elif ix == N_X+indent-1:
                    perfect[count] = [count+1 ,vert[0]+ix*LAT_CONST ,vert[1]+iy*LAT_CONST ,vert[2]+iz*LAT_CONST ,5 , 1,vert[1]+(iy-N_Y)*LAT_CONST,vert[2]+(iz-N_Z)*LAT_CONST]
                    count +=  1
                    perfect[count] = [count+1 ,surf1[0]+ix*LAT_CONST ,surf1[1]+iy*LAT_CONST ,surf1[2]+iz*LAT_CONST ,5 , 1,surf1[1]+(iy-N_Y)*LAT_CONST,surf1[2]+(iz-N_Z)*LAT_CONST]
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


#欠陥や析出物を導入する関数。欠陥や析出物同士がカットオフ長以上の距離となるような条件で指定の濃度に達するまで導入し続ける。
#具体的には欠陥の場合は粒子の存在を無にする。析出物の場合Cuを表すグループ6に書き換える。
def arrange_lat(R,Vac_Cu):
    del_num = 0   
    klist = [1] # absolutely not be void in this position
    for j in range(n_atoms): # any big value is ok
        if del_num >= n_vacancy:
            break
        p = random.randint(1,n_arrange)  
        if R == 0:
            if perfect[p,4] ==1 and perfect[p,1] >= (cut_off+indent*LAT_CONST) and perfect[p,1] <= ((N_X_D+indent)*LAT_CONST -(cut_off)):
                for k in range(len(klist)):
                    if perfect[klist[k],5] ==0 and Vac_Cu ==1 :          # if Vac [k,5] ==0
                        if neighbor(p,klist[k],R,cut_off) == 1:
                            break
                    if perfect[klist[k],4] ==6 and Vac_Cu ==2 :          # if Cu [k,4] ==5
                        if neighbor(p,klist[k],R,cut_off) == 1:
                            break
                    if k == len(klist)-1:
                        if Vac_Cu ==1:    
                            perfect[p,5] =0           # if Vac[p,5] ==0
                            del_num +=1
                            klist.append(p)
                            print(len(klist)-1) 
                            #print(del_num)
                        if Vac_Cu ==2:
                            perfect[p,4] =6           # if Cu [p,4] ==5
                            del_num +=1
                            #print(del_num)
                            klist.append(p)
                            print(len(klist)-1)      
        else:
            #delta = 0.05
            delta = 0.5
            if perfect[p,4] ==1 and perfect[p,1] >= (N_X_D*(0.5-delta)+indent)*LAT_CONST+R+cut_off and perfect[p,1] <= (N_X_D*(0.5+delta)+indent)*LAT_CONST-(R+cut_off):
                for k in range(len(klist)):
                    if perfect[klist[k],5] ==0 and Vac_Cu ==1 :          # if Vac [k,5] ==0
                        if neighbor(p,klist[k],R,cut_off) == 1:
                            break
                    if perfect[klist[k],4] ==6 and Vac_Cu ==2 :          # if Cu [k,4] ==5
                        if neighbor(p,klist[k],R,cut_off) == 1:
                            break
                    if k == len(klist)-1:
                        for v in range(n_arrange):
                            if neighbor(p,v,R,0) ==1:
                                if Vac_Cu ==1:
                                    perfect[v,5] = 0  # if Vac [v,5] ==0
                                    del_num +=1 
                                    klist.append(v)
                                    #print(klist[del_num-1])
                                    #print(del_num)
                                    print(len(klist)-1) 
                                if Vac_Cu ==2:        
                                    perfect[v,4] = 6  # if Cu [v,4] ==5
                                    del_num +=1 
                                    #print(del_num)
                                    klist.append(v)
                                    print(len(klist)-1)
    
    name(Vac_Cu+4)
                   
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

#現在は未使用の関数。はじめに非周期条件で出発させて、波動が協会から離れた段階でindentで空白となっていた場所に原子を埋めることにより周期的境界条件に変更させるための関数
def input_empty():
    f = open("in.lammps", "a")
    for jx in range(indent):
        for jy in range(num):
            for jz in range(num):
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (4, vert[0]+jx*LAT_CONST,vert[1]+jy*LAT_CONST,vert[2]+jz*LAT_CONST))
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (4, surf1[0]+jx*LAT_CONST,surf1[1]+jy*LAT_CONST,surf1[2]+jz*LAT_CONST))

    for jx in range(N_X+indent, N_X+2*indent):
        for jy in range(num):
            for jz in range(num):
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (4, vert[0]+jx*LAT_CONST,vert[1]+jy*LAT_CONST,vert[2]+jz*LAT_CONST))
                f.write("create_atoms %d single %12.6f %12.6f %12.6f\n" % (4, surf1[0]+jx*LAT_CONST,surf1[1]+jy*LAT_CONST,surf1[2]+jz*LAT_CONST))
    f.write("write_restart restart.equil\n")
    f.close()


#現在は未使用の関数。欠陥や析出物導入時に変化する部分のみを出力させる関数
def inverse_lat():
    f = open("inverse.lat", "w")
    lat =0
    for i in range(perfect_num): 
        if perfect[i,5] ==0 : # if Cu [i,4] ==5
            lat +=1
            f.write("%12d %6d %12.6f %12.6f %12.6f\n" % (lat, perfect[i,4],perfect[i,1],perfect[i,2],perfect[i,3]))    
    f.close()
    

#２つの空間座標間の距離がカットオフ長以内であるかどうか判定する関数。
def neighbor(p,k,R,cut):
     if distance(p,k,1)**2 + distance(p,k,2)**2 + distance(p,k,3)**2 <= (cut+R)**2:
         return 1
     else: 
         return 0

#２つの空間座標間のxyz各成分の距離を測る関数
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
            
            
#このプログラムを動作させる部分。無事に動けばFe.latが生成される。
if __name__ == "__main__":
    argvs = sys.argv
    make_perfect()
    arrange_lat(0,1)  #(R ,type) type= vacancy:1 , Cu: 2
    #input_empty()
                       #and RATE or num (=N_X,N_Y)are all parameters
                       #R unit is angstrome
                      

################moving test
#check num of arrage
#if you do not need, delete 'dump'
#
#Use OVITO to check molphorogy
#FFT configuration check
#arraving time should be defined consistnatly
#change x_detec, A2, Time to get beta
