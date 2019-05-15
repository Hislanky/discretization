#! /usr/bin/env python3
# coding: utf-8

import os
import numpy as np
import struct
import argparse as ag
import math
import time

###### DEFINITION FONCTIONS ######
# ARGUMENTS
def arguments():
    parser = ag.ArgumentParser()
    # fichier dat à import
    parser.add_argument("dat_file",help="""DATA file containing the 3D matrix from the discretized texturation""")
    return parser.parse_args()

# EXTRACTION DATA FROM DAT
def extract(data_file) :
    directory = os.path.dirname(__file__)
    path = os.path.join(directory, "data_texture",data_file)
    data_from_conv = open(path,"r")
    data = data_from_conv.readlines()
    data.remove(data[4])
    data.remove(data[2])
    data.remove(data[0])
    for l in range(len(data)) :
        data[l] = data[l].replace('\n','')
    data_from_conv.close()
    return data

def reshape_data(lines_data) :
    space_step = int(lines_data[0]) #in µm
    i_line = lines_data[1]
    i = i_line.split()
    iz = int(i[0])
    ix = int(i[1])
    iy = int(i[2])
    flatM3D = lines_data[2]
    list_data = [space_step, iz, ix, iy, flatM3D]
    return list_data

def reshape_M3D(flatM3D,iz,ix,iy) :
    matrix_1D = [int(i) for i in flatM3D]
    matrix_1D = np.array(matrix_1D)
    matrix_3D = matrix_1D.reshape(iz,ix,iy)
    return matrix_3D

#main
args = arguments()
data_name = args.dat_file
lines = extract(data_name)
shape = reshape_data(lines)

space_step = shape[0] # pas d'espace
indice_z = shape[1] # nb indice en z
indice_x = shape[2] # nb indice en x
indice_y = shape[3] # nb indice en y
flat_3Dmatrix = shape[4]
M3D = reshape_M3D(flat_3Dmatrix, indice_z, indice_x, indice_y)

print(space_step)
print(M3D)
