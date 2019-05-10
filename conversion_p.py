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
    # fichier stl à traiter
    parser.add_argument("STL_file",help="""STL file containing the textured surface to discretize""")
    # pas d'espace pour la discrétisation
    parser.add_argument("Space_step",help="""Step for space discretization in µm""")
    return parser.parse_args()

# EXTRACTION TRIANGLE NOMBRE FROM STL *SIMPLIFIED*
def stl(data_file):
    directory = os.path.dirname(__file__)
    path = os.path.join(directory, "texture_stl",data_file)
    texture = open(path,"r")
    stl = texture.readlines()
    texture.close()
    return stl

def min_axe(axe,triangles) :
    min_value = min(triangles[axe::3])
    return min_value

def max_axe(axe,triangles) :
    max_value = max(triangles[axe::3])
    return max_value

# DETERMINE INTERSECTION BETWEEN TRIANGLES AND Z PLAN
def inter_from_tri(Z,i,triangles) :
    z_values = triangles[2::3]
    z_tri = z_values[i*3:(i*3)+3:1]
    x_values = triangles[0::3]
    x_tri = x_values[i*3:(i*3)+3:1]
    y_values = triangles[1::3]
    y_tri = y_values[i*3:(i*3)+3:1]
    sommet_a = z_tri[0]
    if Z-sommet_a > 0 :
        sign_a = -1
    elif Z-sommet_a == 0 :
        sign_a = 3
    else :
        sign_a = 1
    sommet_b = z_tri[1]
    if Z-sommet_b > 0 :
        sign_b = -1
    elif Z-sommet_b == 0 :
        sign_b = 3
    else :
        sign_b = 1
    sommet_c = z_tri[2]
    if Z-sommet_c > 0 :
        sign_c = -1
    elif Z-sommet_c == 0 :
        sign_c = 3
    else :
        sign_c = 1
    dot_from_inter = list() # de la forme [X1,Y1,X2,Y2]
    if sign_a + sign_b == 0 :
        #récupérer le point d'intersection de Z avec AB
        A = [x_tri[0],y_tri[0],z_tri[0]]
        B = [x_tri[1],y_tri[1],z_tri[1]]
        U = [B[0]-A[0],B[1]-A[1],B[2]-A[2]]
        k = (Z-A[2])/U[2]
        Ix = A[0]+k*U[0]
        Iy = A[1]+k*U[1]
        dot_from_inter.extend([Ix,Iy])
    elif sign_a + sign_b == 6 :
        dot_from_inter.extend([x_tri[0],y_tri[0],x_tri[1],y_tri[1]])  
    if sign_b + sign_c == 0 :
        #récupérer le point d'intersection de Z avec BC
        B = [x_tri[1],y_tri[1],z_tri[1]]
        C = [x_tri[2],y_tri[2],z_tri[2]]
        U = [C[0]-B[0],C[1]-B[1],C[2]-B[2]]
        k = (Z-B[2])/U[2]
        Ix = B[0]+k*U[0]
        Iy = B[1]+k*U[1]
        dot_from_inter.extend([Ix,Iy])
    elif sign_b + sign_c == 6 :
        dot_from_inter.extend([x_tri[1],y_tri[1],x_tri[2],y_tri[2]])
    if sign_c + sign_a == 0 :
        #récupérer le point d'intersection de Z avec CA
        C = [x_tri[2],y_tri[2],z_tri[2]]
        A = [x_tri[0],y_tri[0],z_tri[0]]
        U = [A[0]-C[0],A[1]-C[1],A[2]-C[2]]
        k = (Z-C[2])/U[2]
        Ix = C[0]+k*U[0]
        Iy = C[1]+k*U[1]
        dot_from_inter.extend([Ix,Iy])
    elif sign_c + sign_a == 6 :
        dot_from_inter.extend([x_tri[0],y_tri[0],x_tri[2],y_tri[2]])
    if sign_a + sign_b + sign_c == 9 :
        dot_from_inter.extend([x_tri[0],y_tri[0],x_tri[1],y_tri[1]])
        dot_from_inter.extend([x_tri[0],y_tri[0],x_tri[2],y_tri[2]])
        dot_from_inter.extend([x_tri[1],y_tri[1],x_tri[2],y_tri[2]])
    else :
        if sign_a == 3 :
            dot_from_inter.extend([x_tri[0],y_tri[0],x_tri[0],y_tri[0]])
        if sign_b == 3 :
            dot_from_inter.extend([x_tri[1],y_tri[1],x_tri[1],y_tri[1]])
        if sign_c == 3 :
            dot_from_inter.extend([x_tri[2],y_tri[2],x_tri[2],y_tri[2]])
    return dot_from_inter

# DETERMINE INTERSECTION BETWEEN SEGMENTS AND LINE X 
def inter_from_seg(X,i,inter_plan) :
    x_values = inter_plan[0::2]
    x_tri = x_values[i*2:(i*2)+2:1]
    y_values = inter_plan[1::2]
    y_tri = y_values[i*2:(i*2)+2:1]
    x_p1 = x_tri[0]
    if X-x_p1 > 0 :
        sign_p1 = -1
    elif X-x_p1 == 0 :
        sign_p1 = 3
    else :
        sign_p1 = 1
    x_p2 = x_tri[1]
    if X-x_p2 > 0 :
        sign_p2 = -1
    elif X-x_p2 == 0 :
        sign_p2 = 3
    else :
        sign_p2 = 1
    dot_inter = list()
    if sign_p1 + sign_p2 == 0 :
        #récupérer le point d'intesection de X avec P1P2
        P1 = [x_tri[0],y_tri[0]]
        P2 = [x_tri[1],y_tri[1]]
        U = [P2[0]-P1[0],P2[1]-P1[1]]
        k = (X-P1[0])/U[0]
        YI = P1[1]+k*U[1]
        dot_inter.extend([YI])
    if sign_p1 + sign_p2 == 6 :
        dot_inter.extend([y_tri[0],y_tri[1]])
    if sign_p1 + sign_p2 == 4 or sign_p1 + sign_p2 == 2 :
        if sign_p1 == 3 :
            dot_inter.extend([y_tri[0],y_tri[0]])
        if sign_p2 == 3 :
            dot_inter.extend([y_tri[1],y_tri[1]])
    return dot_inter

# TEST NODES IN OR OUT SOLID DOMAIN
# WARNING : " pour l'instant simplifier en comblant les "troues"
def test_solid(Y, list_yi, delta_x) :
    d_min = min(list_yi)-((delta_x)/2)
    d_max = max(list_yi)+((delta_x)/2)
    test_value = 0
    if Y > d_min and Y < d_max :
        test_value = 1
    return test_value

# EXTRACT NAME FROME FILE
def extract_name(namefile) :
    name = list(namefile)
    i = 0
    while i < 4 :
        name.pop()
        i += 1
    name = ''.join(name)
    return name

# DETERMINE EXTERNAL NODE
def extern_node(matrix, i_x, i_y, i_z, min_x, min_y, min_z, pas) :
    Node = list()
    for z in range(i_z) :
        for x in range(i_x) :
            for y in range(i_y) :
                if matrix[z,x,y] == 1 :
                    if z != 0 and z != i_z and x != 0 and x != i_x and y != 0 and y != i_y :
                        # SELECT IF 26 NEIGHBOUR NODE        # COULD BE CHAGE TO 18 OR 6 NEIGHBOUR
                        nb_voisin1 = matrix[z-1,x-1,y+1]+matrix[z-1,x,y+1]+matrix[z-1,x+1,y+1]+matrix[z-1,x-1,y]+matrix[z-1,x,y]+matrix[z-1,x+1,y]+matrix[z-1,x-1,y-1]+matrix[z-1,x,y-1]+matrix[z-1,x+1,y-1]
                        nb_voisin2 = matrix[z,x-1,y+1]+matrix[z,x,y+1]+matrix[z,x+1,y+1]+matrix[z,x-1,y]+matrix[z,x+1,y]+matrix[z,x-1,y-1]+matrix[z,x,y-1]+matrix[z,x+1,y-1]
                        nb_voisin3 = matrix[z+1,x-1,y+1]+matrix[z+1,x,y+1]+matrix[z+1,x+1,y+1]+matrix[z+1,x-1,y]+matrix[z+1,x,y]+matrix[z+1,x+1,y]+matrix[z+1,x-1,y-1]+matrix[z+1,x,y-1]+matrix[z+1,x+1,y-1]
                        nb_voisin = nb_voisin1 + nb_voisin2 + nb_voisin3
                        if nb_voisin != 26 :
                            X = min_x + x*pas
                            Y = min_y + y*pas
                            Z = min_z + z*pas
                            Node.append([X,Y,Z])
                    elif z == 0 or z == i_z or x == 0 or x == i_x or y == 0 or y == i_y :
                        X = min_x + x*pas
                        Y = min_y + y*pas
                        Z = min_z + z*pas
                        Node.append([X,Y,Z])
    return Node

###### MAIN PROGRAM ######
if __name__ == "__main__":
    args = arguments()
    texture_file = args.STL_file
    delta_X = args.Space_step # LA VALEUR EN [µm]

    time_try1 = time.time()
    ###### DATA EXTRACTION ######
    try:
        delta_x = float(delta_X)
        delta_x = delta_x/1000 # LA VALEUR EN [mm]
        stl_raw_data = stl(texture_file)
        # Suppression des '\n'
        for l in range(len(stl_raw_data)) :
            stl_raw_data[l] = stl_raw_data[l].replace('\n','')
        # Nombre de triangle : première ligne
        nb_triangles = int(stl_raw_data[0])
        # Coordonné des triangles : A(xyz) B(xyz) C(xyz) ligne 3 5 7 9
        triangles = list()
        for i in range(nb_triangles) :
            triangles.extend(list(stl_raw_data[2+2*i]))
        for k in range(len(triangles)) :
            triangles[k] = int(triangles[k])
        # Normales des triangles : Nx Ny Nz ligne 2 4 6 8
        normal_triangles = list()
        for i in range(nb_triangles) :
            normal_triangles.extend(list(stl_raw_data[1+2*i]))
        count_ = 0
        for k in range(len(normal_triangles)) :
            if normal_triangles[k] == '-' :
                normal_triangles[k+1] = '-1'
                count_ += 1
        for l in range(count_) :
            normal_triangles.remove("-")
        for f in range(len(normal_triangles)) :
            normal_triangles[f] = int(normal_triangles[f])
        #PRINT TEST
        #print(stl_raw_data)
        #print(nb_triangles)
        #print(triangles)
        #print(normal_triangles)

    except ValueError :
        print('WARNING : Le programme attend une chaine de caractères pour le nom du fichier STL et un nombre pour le pas.')
        print('Il se peut que vous ayez inversé les deux arguments.')
        print('Rapel : python conversion.py STL_file Space_step')
    except FileNotFoundError :
        print('WARNING : Le fichier STL que vous avez indiqué est introuvable !')
        print('Vérifier que le fichier se trouve bien dans le dossier "texture_stl".')
    else:
        print('Le fichier STL à bien été trouvé.')
        print('Le fichier STL est composé de {} triangles.'.format(nb_triangles))
        print('Le pas de discrétisation spaciale choisit est de {} µm.'.format(delta_X))
        print('Temps extraction : {} secondes'.format(time.time()-time_try1))
        print('** Discrétisation en cours **')
    finally:
        pass

    time_try2 = time.time()
    ###### DISCRETIZATION ######
    try :
        # DETERMINATION OF EXTREMUMS VALUES ALONG THE 3 AXIS
        min_x=min_axe(0,triangles)
        max_x=max_axe(0,triangles)
        min_y=min_axe(1,triangles)
        max_y=max_axe(1,triangles)
        min_z=min_axe(2,triangles)
        max_z=max_axe(2,triangles)
        
        # DETERMINATION OF THE RANGE OF VALUES ALONG THE 3 AXIS
        nb_values_x = len([i for i in np.arange(min_x, max_x + delta_x, delta_x)])
        nb_values_y = len([i for i in np.arange(min_y, max_y + delta_x, delta_x)])
        nb_values_z = len([i for i in np.arange(min_z, max_z + delta_x, delta_x)])

        M3D = np.zeros((nb_values_z,nb_values_x,nb_values_y)) # martrice 3D remplis de 0
   
        Zi = min_z
        indice_z = 0
        # DECOUPAGE EN Z
        while Zi < max_z + delta_x :
            # DETERMINATION DES TRIANGLES COUPEES ET INTERSECTION
            inter_plan_z = list() 
            for i in range(nb_triangles) :
                seg = inter_from_tri(Zi,i,triangles)
                inter_plan_z.extend(seg)
            nb_seg = math.floor(len(inter_plan_z)/4)
            Xj = min_x
            indice_x = 0
            # DECOUPAGE EN X
            while Xj < max_x + delta_x :
                inter_line_x = list()
                # DETERMINATION DES SEGEMENTS COUPEES ET INTERSECTION
                for j in range(nb_seg) :
                    dot = inter_from_seg(Xj, j, inter_plan_z)
                    inter_line_x.extend(dot)
                Yk = min_y
                indice_y = 0
                while Yk < max_y + delta_x :
                    # TEST POUR PASSER A 1 LES NOEUDS SOLIDES
                    if inter_line_x != [] :
                        test = test_solid(Yk, inter_line_x, delta_x)
                        M3D[indice_z,indice_x,indice_y] = test
                    Yk += delta_x
                    indice_y += 1 
                Xj += delta_x
                indice_x += 1
            Zi += delta_x
            indice_z += 1
        M3D = M3D.astype(int)
        print(M3D)
    except :
        print('WARNING : erreur(s) lors de la discrétisation.')
    else :
        print('La discrétisation a été effectué avec succès.')
        print('Temps discrétisation : {} secondes'.format(time.time()-time_try2))
        print('** Ecriture en cours **')     
    finally:
        pass

    time_try3 = time.time()
    ###### DATA FILES WRITING ######
    try:
        name = extract_name(texture_file)
        # CREATION FICHIER DAT
        new_data = open('data_texture\{}_X{}.dat'.format(name, delta_X),'w') # il a un soucis avec le "\"
        new_data.write('Space step [µm]:'+'\n')
        new_data.write('{}'.format(delta_X)+'\n')
        new_data.write('Maximum indice on each axes [Z;X;Y] :'+'\n')
        new_data.write('{} {} {}'.format(indice_z, indice_x, indice_y)+'\n') # remplacer par Zi Xj Yk pour avoir la taille total du lattice en [mm]
        new_data.write('Values of the 3D matrix from list.flatten() : \n')
        flat_M3D = M3D.flatten()
        for i in range(len(flat_M3D)) :
            new_data.write('{}'.format(flat_M3D[i]))
        new_data.close()

        # CREATION FICHIER VTK (POUR VISUALISATION PARAVIEW)
        # Détermination des noeud de la peau extérieure
        Node_ext = extern_node(M3D, indice_x, indice_y, indice_z, min_x, min_y, min_z, delta_x)
        #print(Node_ext[1][1])
        vtk_data = open('data_texture\{}_X{}.vtk'.format(name, delta_X),'w') # il a un soucis avec le "\"
        # ecrire le fichier vtk pour qu'il soit lisible par paraview
        vtk_data.close()
    except :
        print('WARNING : erreur(s) lors de la création des fichiers .dat et .vtk')
    else:
        print('Les fichiers {}_X{}.dat'.format(name, delta_X)+' et {}_X{}.vtk'.format(name, delta_X)+' ont été crées dans le dossier "data_texture".')
        print('Temps écriture : {} secondes'.format(time.time()-time_try3))
    finally:
        print('- - - - - - -  FIN PROGRAMME  - - - - - - -')