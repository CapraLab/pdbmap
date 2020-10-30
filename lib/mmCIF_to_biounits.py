#!/usr/bin/env pythin
#*************************
# * Python3 code to re-implement
# * Assemblies.py located at http://mmcif.wwpdb.org/docs/sw-examples/python/html/assemblies.html
# *
# * Given a Biopython mmCIF dictionary, generate a new structure every assembly
# * listed in the pdbx_struct_assembly category table by performing
# * rotation and translation matrix operations and creating a new atom_site
# * category table for each assembly.
# *
# * Lines with superscriptions contain footnoted references or explanations.
# *************************/

# I do not think this is working quite yet

import copy
from os.path import splitext
import sys
import numpy as np
# from pdbx.reader.PdbxReader import PdbxReader
# from pdbx.writer.PdbxWriter import PdbxWriter
# from pdbx.reader.PdbxContainers import *
import logging
LOGGER = logging.getLogger(__name__)

def parseOperationExpression(expression: str) :
    LOGGER.debug('Parsing biounit operation %s',str)
    operations = []
    stops = [ "," , "-" , ")" ]

    currentOp = ""
    i = 1
	
    # Iterate over the operation expression
    while i in range(1, len(expression) - 1):
        pos = i

        # Read an operation
        while expression[pos] not in stops and pos < len(expression) - 1 : 
            pos += 1    
        currentOp = expression[i : pos]

        # Handle single operations
        if expression[pos] != "-" :
            operations.append(currentOp)
            i = pos

        # Handle ranges
        if expression[pos] == "-" :
            pos += 1
            i = pos
			
            # Read in the range's end value
            while expression[pos] not in stops :
                pos += 1
            end = int(expression[i : pos])
			
            # Add all the operations in [currentOp, end]
            for val in range((int(currentOp)), end + 1) :
                operations.append(str(val))
            i = pos
        i += 1
    return operations

def getValue(possible_list,list_element):
    assert type(possible_list) == list or list_element == 0,\
        "Failed request for element %d of scaler %s"%(list_element,possible_list)
    if type(possible_list) == list:
        return possible_list[list_element]
    return possible_list

def prepareOperation(mmcif_dict, op1index, op2index) :
    # Prepare matrices for operations 1 & 2
    op1 = [[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 1]]
    op2 = [[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 1]]
	
    # Fill the operation matrices for operations 1 & 2
    for i in range(3) :
        op1[i][3] = float(getValue(mmcif_dict["_pdbx_struct_oper_list.vector[" + str(i + 1) + "]"], op1index))
		
        if (op2index != -1) :
            op2[i][3] = float(getValue(mmcif_dict["_pdbx_struct_oper_list.vector[" + str(i + 1) + "]"], op2index))
        for j in range(3) :
            op1[i][j] = float(getValue(mmcif_dict["_pdbx_struct_oper_list.matrix[" + str(i + 1) + "][" + str(j + 1) + "]"], op1index))
            if (op2index != -1) :
                op2[i][j] = float(getValue(mmcif_dict["_pdbx_struct_oper_list.matrix[" + str(i + 1) + "][" + str(j + 1) + "]"], op2index))
    
    # Handles non-Cartesian product expressions
    if (op2index == -1) :
        return op1

    operation = [[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 0],[0, 0, 0, 1]]

    # Handles Cartesian product expressions (4x4 matrix multiplication)
    for row in range(4) :
        for col in range(4) :
            operation[row][col] = sum([ (op1[row][r] * op2[r][col]) for r in range(4)])
    return operation

def mmCIF_to_biounits(mmcif_dict):
    """Return a list of mmcif dictionaries for each biounit, and a dictionary of chains copies to new chains"""
    biounits_as_mmcif_dicts = []

    assembly_ids = mmcif_dict.get('_pdbx_struct_assembly_gen.assembly_id',None)
    if not assembly_ids: # We cannot create biounits if there are not assemblies to create
        return []
    assembly_count = len(assembly_ids) if type(assembly_ids)==list else 1

    atom_site_ids = mmcif_dict.get('_atom_site.id',None)
    if not atom_site_ids or (not type(atom_site_ids)== list) or len(atom_site_ids) < 2: # We cannot create biounits if there are no atoms to speak of
        return []

    atom_site_count = len(atom_site_ids)
    del atom_site_ids

    atom_site_Cartn_x = mmcif_dict.get('_atom_site.Cartn_x')
    atom_site_Cartn_y = mmcif_dict.get('_atom_site.Cartn_y')
    atom_site_Cartn_z = mmcif_dict.get('_atom_site.Cartn_z')

    assert len(atom_site_Cartn_x) == atom_site_count,"Cartn_x entries are missing in the structure"
    assert len(atom_site_Cartn_y) == atom_site_count,"Cartn_x entries are missing in the structure"
    assert len(atom_site_Cartn_z) == atom_site_count,"Cartn_x entries are missing in the structure"

    # Create a CIF dictionary for every assembly specified in pdbx_struct_assembly_gen
    
    for index in range(assembly_count):
        LOGGER.info("Processing biounit assembly %d of %d",index,assembly_count)
        # Ultimately we want to create transformed lists of atomic coordinates
        # setup numpy matrix to receive  maximum number of calculation results
        # Keep in mind that often many atoms are excluded from biounits
        new_Cartn = np.empty((atom_site_count,3))
    	
        # Lists to hold the individual operations
        oper = []
        oper2 = []

        # Keep track of the current atom and model number
        # Because often atoms are left out of biounits
        # I suppose, in theory, a model could be skipped too
        atom_index = 0
        model_index = 0


        # Retrieve the assembly_id attribute value for this assembly
        assemblyId = getValue(mmcif_dict['_pdbx_struct_assembly_gen.assembly_id'],index)

        # Retrieve the operation expression for this assembly from the oper_expression attribute	
        oper_expression = getValue(mmcif_dict["_pdbx_struct_assembly_gen.oper_expression"], index)

        # Count the number of left parentheses in the operation expression
        parenCount = oper_expression.count("(")

        # Handles one operation assemblies (e.g., "1")
        if parenCount == 0 : oper.append(oper_expression)
    	
        # Handles multiple operation assemblies, no Cartesian products (e.g., "(1-5)")
        if parenCount == 1 : oper.extend(parseOperationExpression(oper_expression))
    	
        # Handles Cartesian product expressions (e.g., "(X0)(1-60)")
        if parenCount == 2 :
            # Break the expression into two parenthesized expressions and parse them
            temp = oper_expression.find(")")
            oper.extend(parseOperationExpression(oper_expression[0:temp+1]))
            oper2.extend(parseOperationExpression(oper_expression[temp+1:]))

        # Retrieve the asym_id_list, which indicates which atoms to apply the operations to
        asym_id_list = getValue(mmcif_dict["_pdbx_struct_assembly_gen.asym_id_list"], index)

        temp = (1 > len(oper2)) and 1 or len(oper2)

        # Retrieve the pdbx_struct_oper_list category table, which details translation and rotation 
        # operations required to generate/transform assembly coordinates
        oper_list_ids = mmcif_dict["_pdbx_struct_oper_list.id"]
        oper_list_count = len(oper_list_ids) if type(oper_list_ids)==list else 1

        # Ultimately, we hvae to replace all the _atom_sites.* dictionary elements
        # in our created biounit...  In theory those keys could vary from cif to cif
        # The coordinates and model number are handled by other mechanisms besides simple copy 
        new_atom_site_lists = {}
        for key in mmcif_dict:
            if key.startswith("_atom_site.") and (".Cartn_" not in key) and (".pdbx_PDB_model_num" not in key):
                new_atom_site_lists[key] = []
        new_atom_site_pdbx_PDB_model_nums = []

        # For every operation in the first parenthesized list
        for op1 in oper :
            LOGGER.debug("op1=%s",str(op1))
            # Find the index of the current operation in the oper_list category table
            op1index = 0
            for row in range(oper_list_count):
                if getValue(mmcif_dict['_pdbx_struct_oper_list.id'],row) == op1:
                    op1index = row
                    break

            # For every operation in the second parenthesized list (if there is one)
            for i in range(temp) :		
                # Find the index of the second operation in the oper_list category table
                op2index = -1
                if (oper2) :
                    for row in range(oper_list_count):
                        if getValue(mmcif_dict['_pdbx_struct_oper_list.id'],row) == oper2[i]:
                            op2index = row
                            break

                # Prepare the operation matrix
                operation = prepareOperation(mmcif_dict, op1index, op2index)

                # Iterate over every atom in the atom_site reference table
                for r in range(atom_site_count):
                    # If the asym_id of the current atom is not in the asym_id list, skip to the next atom
                    if (asym_id_list.find(mmcif_dict["_atom_site.label_asym_id"][r]) == -1) :
                        continue
                    
                    # Add this row to the atom_site table for this assembly
                    # Simply copy the new mmcif_dict[_atom_site.*] components to the growing lists
                    for key in new_atom_site_lists:
                        new_atom_site_lists[key].append(mmcif_dict[key][r])

                    new_atom_site_pdbx_PDB_model_nums.append(str(model_index+1))

                    # Update the atom number and model number for this row
                    # atom_site.setValue(str(atomNum), "id", atomNum - 1)
                    # atom_site.setValue(str(modelNum), "pdbx_PDB_model_num", atomNum - 1) 

                    # Determine and set the new coordinates
                    coords = np.array([float(atom_site_Cartn_x[r]), float(atom_site_Cartn_y[r]), float(atom_site_Cartn_z[r]), 1.0])

                    for xyz in range(3):
                        new_Cartn[atom_index,xyz] = sum([(operation[xyz][b] * coords[b]) for b in range(4)])
                    atom_index += 1
                model_index += 1

        # Create the new dictionary, but with the new coordinates
        # First, shallow copy the dictionary
        new_mmcif_dict = mmcif_dict.copy()
        # Then, replace the _atom_site.* elements with the atom_sites lists of atoms retained in the biounit
        for key in new_atom_site_lists:
                new_mmcif_dict[key] = new_atom_site_lists[key] 
        new_mmcif_dict['_atom_site.id'] = [str(n) for n in range(1,atom_index+1)]
        new_mmcif_dict['_atom_site.Cartn_x'] = ["%.3f"%x for x in new_Cartn[0:atom_index,0]]
        new_mmcif_dict['_atom_site.Cartn_y'] = ["%.3f"%y for y in new_Cartn[0:atom_index,1]]
        new_mmcif_dict['_atom_site.Cartn_z'] = ["%.3f"%z for z in new_Cartn[0:atom_index,2]]
        new_mmcif_dict['_atom_site.pdbx_PDB_model_num'] = new_atom_site_pdbx_PDB_model_nums
        # for key in new_mmcif_dict:
        #   if '_atom_site.' in key:
        #      print(key,len(new_mmcif_dict[key]))
        biounits_as_mmcif_dicts.append(new_mmcif_dict)
    return biounits_as_mmcif_dicts        
