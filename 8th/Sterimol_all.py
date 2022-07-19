#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np

def Sterimol(mol_name,txt_name,mol):
    x,y,z = np.loadtxt(mol_name,skiprows=2,usecols=(1,2,3),unpack=True)
    L=max(z) + 1.2

    #原子名の読み込み
    Atoms = np.loadtxt(mol_name,skiprows=2,usecols=(0,),unpack=True,dtype='S3')
    #全原子数
    number_atoms = len(x)

    #水素の数をnumber_Hに格納
    number_H = 0
    for i in range(number_atoms):
        if Atoms[i] == "H".encode():
            number_H = number_H + 1

    B5 = np.zeros(number_H-1)
    x_B1 = np.zeros(number_H-1)
    y_B1 = np.zeros(number_H-1)

    #xy平面上に射影した全原子の原点からの距離
    distxy = np.sqrt(x**2 + y**2)

    j = 0
    for i in range(number_atoms):
    #カルボキシル基上の水素以外の水素の、xy平面上に射影した水素の原点からの距離、x座標、y座標を代入（今回のxyzファイルでは、カルボキシル基の水素は上から4つ目に統一）
        if i != 3 and Atoms[i] == "H".encode():
            B5[j] = distxy[i]
            x_B1[j] = x[i]
            y_B1[j] = y[i]
            j = j + 1

    B5 = max(B5) + 1.20

    #x軸、y軸正負それぞれ4方向のもっとも原点から離れている4点を代入するための変数
    B1_minus_x = B1_minus_y = B1_plus_x = B1_plus_y = 0.0
    for i in range(j):
    #負の方向は、x_B1(y_B1)に格納されているi番目の値が、現在代入されている値より小さい場合、その値を代入
        if x_B1[i] <= 0 and B1_minus_x > x_B1[i]:
            B1_minus_x = x_B1[i]
        if y_B1[i] <= 0 and B1_minus_y > y_B1[i]:
            B1_minus_y = y_B1[i]
    #正の方向は、x_B1(y_B1)に格納されているi番目の値が、現在代入されている値より大きい場合、その値を代入
        if x_B1[i] > 0 and B1_plus_x < x_B1[i]:
            B1_plus_x = x_B1[i]
        if y_B1[i] > 0 and B1_plus_y < y_B1[i]:
            B1_plus_y = y_B1[i]

    #４点の原点からの距離のうち最大のものを抽出
    B1 = min(-B1_minus_x, B1_plus_x, -B1_minus_y, B1_plus_y )+1.20

    print (mol_name, "L:", round(L,2), "B1:", round(B1,2), "B5:", round(B5,2))
    if mol == 1:
        f = open(txt_name,'w')
        f.write('\t'+'\t'+"L"+'\t'+"B1"+'\t'+"B5"+'\n')
        f.write(str(mol_name)+'\t'+str(round(L,2))+'\t'+str(round(B1,2))+'\t'+str(round(B5,2))+'\n')
        f.close
    elif mol != 1:
        f = open(txt_name,'a')
        f.write(str(mol_name)+'\t'+str(round(L,2))+'\t'+str(round(B1,2))+'\t'+str(round(B5,2))+'\n')
        f.close

for mol in range(1,27):
    mol_name = "ChartonT" + str(mol) + ".xyz"
    txt_name = "Sterimol_all.txt"
    Sterimol(mol_name,txt_name,mol)
