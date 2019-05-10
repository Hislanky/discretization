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

# EXTRACTION TRIANGLE NOMBRE FROM STL
def stl_nb(data_file):
    directory = os.path.dirname(__file__)
    path = os.path.join(directory, "texture_stl",data_file)
    texture = open(path,"rb")
    texture.seek(80,0)
    nb_triangles = int.from_bytes(texture.read(4), 'little')
    texture.close()
    return nb_triangles

# EXTRACTION TRIANGLE COORDINATE FROM STL
def stl_analysis(data_file):
    directory = os.path.dirname(__file__)
    path = os.path.join(directory, "texture_stl",data_file)
    texture = open(path,"rb")
    texture.seek(84,0)
    # Liste des coordonnées des sommets des triangles [x,y,z,x,y,z...x,y,z]
    triangles = list()
    paring_range = 0
    while paring_range < nb_triangles :
        sommets = 0
        # passer/(stocker) les 12 octets de la direction normale
        texture.seek(12,1)
        # coordonnees des sommets {X;Y;Z}
        while sommets < 3 :
            x = struct.unpack('f', texture.read(4))
            x = x[0]
            y = struct.unpack('f', texture.read(4))                    
            y =y[0]
            z = struct.unpack('f', texture.read(4))
            z = z[0]
            triangles.extend([x,y,z])
            sommets += 1
        # passer les 2 octets de control
        texture.seek(2,1)
        paring_range += 1
    texture.close()
    return triangles

def min_axe(axe,triangles) :
    min_value = min(triangles[axe::3])
    return min_value

def max_axe(axe,triangles) :
    max_value = max(triangles[axe::3])
    return max_value

# DETERMINE INTERSECTION BETWEEN TRIANGLES AND Z PLAN
def inter_from_tri(Z,i,triangles) :
    z_values = triangles[2::3]
    z_tri = z_values[i:i+3:1]
    x_values = triangles[0::3]
    x_tri = x_values[i:i+3:1]
    y_values = triangles[1::3]
    y_tri = y_values[i:i+3:1]
    sommet_a = z_tri[0]
    if Z-sommet_a > 0 :
        sign_a = 1
    else :
        sign_a = 0
    sommet_b = z_tri[1]
    if Z-sommet_b > 0 :
        sign_b = 1
    else :
        sign_b = 0
    sommet_c = z_tri[2]
    if Z-sommet_c > 0 :
        sign_c = 1
    else :
        sign_c = 0
    dot_from_inter = list() # de la forme [X1,Y1,X2,Y2]
    if sign_a + sign_b == 1 :
        #récupérer le point d'intersection de Z avec AB
        A = [x_tri[0],y_tri[0],z_tri[0]]
        B = [x_tri[1],y_tri[1],z_tri[1]]
        U = [B[0]-A[0],B[1]-A[1],B[2]-A[2]]
        k = (Z-A[2])/U[2]
        Ix = A[0]+k*U[0]
        Iy = A[1]+k*U[1]
        dot_from_inter.extend([Ix,Iy])     
    if sign_b + sign_c == 1 :
        #récupérer le point d'intersection de Z avec BC
        B = [x_tri[1],y_tri[1],z_tri[1]]
        C = [x_tri[2],y_tri[2],z_tri[2]]
        U = [C[0]-B[0],C[1]-B[1],C[2]-B[2]]
        k = (Z-B[2])/U[2]
        Ix = B[0]+k*U[0]
        Iy = B[1]+k*U[1]
        dot_from_inter.extend([Ix,Iy])
    if sign_c + sign_a == 1 :
        #récupérer le point d'intersection de Z avec CA
        C = [x_tri[2],y_tri[2],z_tri[2]]
        A = [x_tri[0],y_tri[0],z_tri[0]]
        U = [A[0]-C[0],A[1]-C[1],A[2]-C[2]]
        k = (Z-C[2])/U[2]
        Ix = C[0]+k*U[0]
        Iy = C[1]+k*U[1]
        dot_from_inter.extend([Ix,Iy])
    return dot_from_inter

# DETERMINE INTERSECTION BETWEEN SEGMENTS AND LINE X 
def inter_from_seg(X,i,inter_plan) :
    x_values = inter_plan[0::2]
    x_tri = x_values[i:i+2:1]
    y_values = inter_plan[1::2]
    y_tri = y_values[i:i+2:1]
    x_p1 = x_tri[0]
    if X-x_p1 > 0 :
        sign_p1 = 1
    else :
        sign_p1 = 0
    x_p2 = x_tri[1]
    if X-x_p2 > 0 :
        sign_p2 = 1
    else :
        sign_p2 = 0
    dot_inter = list()
    if sign_p1 + sign_p2 == 1 :
        #récupérer le point d'intesection de X avec P1P2
        P1 = [x_tri[0],y_tri[0]]
        P2 = [x_tri[1],y_tri[1]]
        U = [P2[0]-P1[0],P2[1]-P1[1]]
        k = (X-P1[0])/U[0]
        YI = P1[1]+k*U[1]
        dot_inter.extend([YI])
    return dot_inter

# TEST NODES IN OR OUT SOLID DOMAIN
# WARNING : " pour l'instant simplifier en comblant les "troues"
def test_solid(Y, list_yi, delta_x) :
    d_min = min(list_yi)-(delta_x)/2
    d_max = max(list_yi)+(delta_x)/2
    test_value = 0
    if Y >= d_min and Y <= d_max :
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
        nb_triangles = stl_nb(texture_file)
        triangles = stl_analysis(texture_file)
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
            inter_plan_z = list() 
            # DETERMINATION DES TRIANGLES COUPEES ET INTERSECTION
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
                for i in range(nb_seg) :
                    dot = inter_from_seg(Xj,i,inter_plan_z)
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
        new_data = open('data_texture\{}_X{}.dat'.format(name, delta_X),'w')
        new_data.write('Space step [µm]:'+'\n')
        new_data.write('{}'.format(delta_X)+'\n')
        new_data.write('Z-X-Y|test values: 1 if in solid \n')
        for i in range(nb_values_z) :
            for j in range(nb_values_x) :
                for k in range(nb_values_y) :
                    new_data.write('{}-{}-{}|{}'.format(i,j,k,M3D[i,j,k])+'\n')
        new_data.close()

        # (TEST OPTIMIZATION) CREATION FICHIER DAT
        new_data1 = open('data_texture\{}_X{}_v2.dat'.format(name, delta_X),'w')
        new_data1.write('Space step [µm]:'+'\n')
        new_data1.write('{}'.format(delta_X)+'\n')
        new_data1.write('For each Node Z X Y Value \n')
        for i in range(nb_values_z) :
            new_data1.write('{} '.format(i)+'\n')
            for j in range(nb_values_x) :
                new_data1.write('{} '.format(j))
                for k in range(nb_values_y) :
                    new_data1.write('{} {}'.format(k,M3D[i,j,k])+'|')
                new_data1.write('\n')
        new_data1.close()

        # (TEST OPTI) AVEC FLATTEN (PUIS ###########DANS IMPORT)

        # CREATION FICHIER VTK (POUR VISUALISATION PARAVIEW)
        vtk_data = open('data_texture\{}_X{}.vtk'.format(name, delta_X),'w')
        # ecrir le fichier vtk pour qu'il soit lisible par paraview
        vtk_data.close()
    except :
        print('WARNING : erreur(s) lors de la création des fichier .dat et .vtk')
    else:
        print('Les fichiers {}_X{}.dat'.format(name, delta_X)+' et {}_X{}.vtk'.format(name, delta_X)+' ont été crées dans le dossier "data_texture".')
        print('Temps écriture : {} secondes'.format(time.time()-time_try3))
    finally:
        print('- - - - - - -  FIN PROGRAMME  - - - - - - -')