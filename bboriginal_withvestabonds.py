#!/usr/bin/env python
# coding: utf-8

# In[2]:



import winsound
import argparse
import os
import shutil
from pymatgen.core.composition import Composition
import csv
#import math
import pandas as pd
#import time
import itertools
import argparse
#import numpy as np

#issue with fraction import gcd
#from matminer.featurizers.structure import SiteStatsFingerprint
from pymatgen.core.structure import Structure
from pymatgen.analysis.molecule_structure_comparator import CovalentRadius
from pymatgen.io.cif import CifWriter
from pymatgen.io.cif import CifParser
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.core.periodic_table import Element
import time
import re
#from pymatgen.transformations.standard_transformations import AutoOxiStateDecorationTransformation


class Error(Exception):
    pass

class OxiStateError(Error):
    pass

class BuildingBlock:
    def __init__(self, center_site):
        self.center_atoms = [center_site]
        self.satellite_atoms = []
        self.coordination_number = 0
    def __str__(self):
        text = "Center atom(s): " + self.center_atoms[0].specie.symbol + f"    Coordinates: {self.center_atoms[0].x:.2f},{self.center_atoms[0].y:.2f},{self.center_atoms[0].z:.2f}"+ " \nCoordination Number: " + str(self.coordination_number)
        text += "\nSatellite atoms: \n"
        for atoms in self.satellite_atoms:
            text += atoms.specie.symbol + "\n"
        return text
    def export(self):
        bb = []
        bb.append(self.coordination_number)
        bb.append(self.center_atoms[0].specie.symbol)
        bb.append([self.center_atoms[0].x,self.center_atoms[0].y,self.center_atoms[0].z])
        for i,atom in enumerate(self.satellite_atoms):
            #atom_name = 'Satellite Atom'+str(i)
            bb.append(atom.specie.symbol)
            #coord_name = 'Satellite' + str(i) + ' Coords'
            bb.append([atom.x,atom.y,atom.z])
        #dict.update({'coord num': self.coordination_number})
        return bb

class StructureBB():
    def __init__(self):
        self.bb_list = []
    def set_bb_list(self, bb_list):
        self.bb_list = bb_list
        self.num_bb = len(self.bb_list)
    def __str__(self):
        text = "Number of building blocks: " + str(self.num_bb)
        return text
    def export(self):
        data = {}
        for i,bb in enumerate(self.bb_list):
            name = 'BuildingBlock'+str(i)
            data.update({name: bb.export()})
            #print(data)
        df = pd.DataFrame.from_dict(data, orient='index')
        #(data, columns = ['Center Atom(s)','Satellite Atoms'])
        return df

def readthreshold():
    df = pd.read_table(r"C:\Users\janic\Downloads\element_pair_bond_threshold.txt", delimiter='    ',  header=None)
    df = df.drop(df.columns[0], axis=1)
    df = df.drop(df.columns[2], axis=1)
    df.columns = ['Element1', 'Element2', 'Distance']
    for index in df.index:
        df.loc[index, 'Element1'] = df.loc[index, 'Element1'].lstrip()
        df.loc[index, 'Element2'] = df.loc[index, 'Element2'].lstrip()
    return df

# Search dataframe for specific element bond pair
# Input: dataframe of all element pairs, both elements in pair to search for
# Output: distance
def searchdf(df, elem, bond):
    
    row = df.loc[(df['Element1'] == elem) & (df['Element2'] == bond)]
    try:
        out = (row['Distance'])
        out = out.astype(str).values.flatten().tolist()[0]
        out = out.split()
        return out[0]
    except:
        return None

def findBBs():
    
    parser = argparse.ArgumentParser(description='Find BBs in crystal structures')
    parser.add_argument("--input", help="Path to directory with crystal structures")
    parser.add_argument("--output", help="Path to directory for output")
    args = parser.parse_args()
    
    duration = 1000
    frequency = 440
    start_time = time.time()
    stop_time = 0
    directory = args.input
    output_dir = args.output
    for path, dirs, files in os.walk(directory):
        for filename in files: 
            print(filename)
            fullname = os.path.join(path, filename)
            text_line = []
            text_line.append(filename)
            parser = ""
            if re.search(r"^.*\.cif$",filename):
                fullname = os.path.join(directory, filename)
                parser = CifParser(fullname)
            else:
                continue        
            structure = parser.get_structures()[0]
            formula = structure.formula
            neighbor_finder = CrystalNN()
            bb_list = []
            cn_count = {}
            distances_df = readthreshold()

            for i,x in enumerate(structure.sites):
                #cn needed for finding number of neighbors
                e1 = x.as_dict()['species'][0]['element']
                coordination_number = neighbor_finder.get_cn(structure,i)
                #print('coordination number: ' + str(coordination_number) + ' element: ' + x.specie.symbol)

                bb = BuildingBlock(x)
                neighbors = neighbor_finder.get_nn_info(structure, i)

                for site in neighbors:
                    e2 = site['site'].specie.symbol
                    #get the distance to compare the radius to
                    distance = bb.center_atoms[0].distance(site['site'])
                    #if the distance between the satellite atom and the center is equal to the sum of the two covalent radii, add it to the list
                    #if the radius is within 32%
                    radius = searchdf(distances_df, e1, e2)
                    distance = distance.item()
                    #a = f"{distance:.4f}, {covalent_radius:4f}"
                    #print(a)
                    #print('satellite atom: ' + site['site'].specie.symbol)
                    if(radius):
                        if (distance < float(radius)):
                        #metal bonded to metal: only bonds to itself?
                            #metal bonded to nometal
                            if (bb.center_atoms[0].specie.is_metal):
                                if (not site['site'].specie.is_metalloid and not site['site'].specie.is_metal):
                                    #print('metal bonded to nonmetal')
                                    #print(" used " + str(covalent_radius))
                                    bb.satellite_atoms.append(site['site'])
                                    #cn_count[bb.center_atoms[0]] = cn_count[bb.center_atoms[0]] + 1

                            #metal bonded to metalloid
                            if (bb.center_atoms[0].specie.is_metal):
                                if (site['site'].specie.is_metalloid):
                                    #print('metal bonded to metalloid')
                                    #print(" used " + str(covalent_radius))
                                    bb.satellite_atoms.append(site['site'])

                            #metalloid bonded to metalloid
                            #only bonds to itself??
                            if (bb.center_atoms[0].specie.is_metalloid):
                                if (site['site'].specie.is_metalloid):
                                    #print('metalloid bonded to metalloid')
                                    #print(' used ' + str(covalent_radius))
                                    bb.satellite_atoms.append(site['site'])

                            #metalloid bonded to nonmetal
                            if (bb.center_atoms[0].specie.is_metalloid):
                                if (not site['site'].specie.is_metalloid and not site['site'].specie.is_metal):
                                    #print('metalloid bonded to nonmetal')
                                    #print(" used " + str(covalent_radius))
                                    bb.satellite_atoms.append(site['site'])

                            #nonmetal bonded to nonmetal
                            #doesnt work lol
                            if (not bb.center_atoms[0].specie.is_metal and not bb.center_atoms[0].specie.is_metalloid):
                                if (not site['site'].specie.is_metal and not bb.center_atoms[0].specie.is_metalloid):
                                    if coordination_number >= 3:
                                        #print('nonmetal bonded to nonmetal')
                                        #print(" used " + str(covalent_radius))
                                        bb.satellite_atoms.append(site['site'])

                            else:
                                continue
                        #might readd this if the coord. num keeps being a problem
                        #bb.coordination_number += 1

                    bb.coordination_number = coordination_number
                if len(bb.satellite_atoms) >= 3:
                    bb_list.append(bb)

                #print(len(bb.satellite_atoms))
                #if (bb.coordination_number != len(bb.satellite_atoms)):
                #    print("uh oh")

                #print('next')

            # assumption 1:
            # if there are 4+ different coordination numbers assigned to the same center atom, none of them are valid
            # could help with structures with no bonds
            center_list = {}
            center = []

            for i, e in enumerate(bb_list):    
                if e.center_atoms[0].specie.symbol not in center:
                    center.append(e.center_atoms[0].specie.symbol)
            for i, e in enumerate(center):
                center_list[e] = []

            for i, e in enumerate(bb_list):
                if len(e.satellite_atoms) not in center_list[e.center_atoms[0].specie.symbol]:
                    center_list[e.center_atoms[0].specie.symbol].append(len(e.satellite_atoms))

            for i,e in enumerate(center_list):
                if len(center_list[e]) >= 4:
                    bb_list = [x for x in bb_list if e != x.center_atoms[0].specie.symbol]

            # assumption 2:
            # if a structure contains polyhedra where one element is bonded to itself, then that polyhedra is the only one
            contains_self_bonded_structure = False
            self_structure_center = ''
            for i, e in enumerate(bb_list):
                if all(f.specie.symbol == e.center_atoms[0].specie.symbol for f in e.satellite_atoms):
                    contains_self_bonded_structure = True
                    self_structure_center = e.center_atoms[0].specie.symbol
            if contains_self_bonded_structure == True:
                bb_list = [x for x in bb_list if self_structure_center == x.center_atoms[0].specie.symbol]
                    #print(len(bb.satellite_atoms))
                    #if (bb.coordination_number != len(bb.satellite_atoms)):
                    #    print("uh oh")

                    #print('next')
                        #print("\n")
                #if cn_count > 4:
                    #print('BAD FILE')
                    #continue        
                #print each building block
                #for bb in bb_list:
                    #df = pd.DataFrame(bb.export())
                    #print(bb)
            bbs = StructureBB()
            bbs.set_bb_list(bb_list)
            df = bbs.export()
        #print(fullname)

        #print(df)
        text_output = []
        text_line.append(bbs.num_bb)
        text_output.append(text_line)
        #time.sleep(4)

        csvfile = filename + " " + formula + '.csv'
        #result = df.to_csv(csvfile)

        #args.output
        result_fullname = os.path.join(output_dir, csvfile)
        with open(result_fullname, 'w') as f:
            f.write(df.to_csv())
            f.close()

    #make this so anyone can use it
    #with open(r"C:\Users\janic\Desktop\Output.csv", 'w') as f:
    #   f.truncate(0)
    #  wr = csv.writer(f)
    # wr.writerows(text_output)
    stop_time = time.time()
    print("Runtime: ", stop_time - start_time)


findBBs()


# In[ ]:




