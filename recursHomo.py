import sys
import argparse
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from pulp import LpProblem, LpVariable, lpSum, LpMinimize, LpInteger, PULP_CBC_CMD, LpStatus
from rdkit.Chem import Descriptors, rdMolDescriptors
from itertools import combinations
from AaronTools.geometry import Geometry
from AaronTools.atoms import Atom
from AaronTools.theory import Theory, OptimizationJob, FrequencyJob
import importlib.resources

def get_Hf(inchi):
    try:
        with importlib.resources.open_text("hfrpkg.data", "ATcT_lib.txt", encoding="utf-8") as f:
            for line in f:
                parts = line.strip().split("\t")
                if len(parts) >= 6 and parts[3] == inchi:
                    return float(parts[5])
    except Exception:
        pass
    return None
def kekulize_and_make_aliphatic(mol):
    try:
        mol = Chem.Mol(mol)
        Chem.Kekulize(mol, clearAromaticFlags=True)
        for atom in mol.GetAtoms():
            atom.SetIsAromatic(False)
        for bond in mol.GetBonds():
            bond.SetIsAromatic(False)
        mol.UpdatePropertyCache(strict=False)
        mol = Chem.AddHs(mol)
        return mol 
    except Chem.KekulizeException:
        return None  
def get_fragments(input_mol,level):
    wildcard_reactants = []
    wildcard_products= []
    mol= kekulize_and_make_aliphatic(input_mol)
    if level==1:
        wildcard_reactants = [
            Chem.MolFromSmiles("[H][H]"),
        ]
        wildcard_products = [
            Chem.MolFromSmiles("[*]"),
        ]
    
    if level==2:
        wildcard_reactants = [
            Chem.MolFromSmiles("[*]"),
        ]
        wildcard_products = [
            Chem.MolFromSmiles("[*][*]"),
            Chem.MolFromSmiles("[*]=[*]"),
            Chem.MolFromSmiles("[*]#[*]"),
        ]

    if level ==3:
        wildcard_reactants = [
            Chem.MolFromSmiles("[*][*]"),
            Chem.MolFromSmiles("[*]=[*]"),
            Chem.MolFromSmiles("[*]#[*]"),
        ]
        wildcard_products = [
            Chem.MolFromSmiles("[*][*][*]"),
            Chem.MolFromSmiles("[*]=[*][*]"),
            Chem.MolFromSmiles("[*]#[*][*]"),
            Chem.MolFromSmiles("[*][*](*)[*]"),
            Chem.MolFromSmiles("[*][*](=[*])[*]"),
            Chem.MolFromSmiles("[*][*@]([*])([*])[*]"),
            Chem.MolFromSmiles("[*]=[*]=[*]"),    
        ]
    if level == 4:
        wildcard_reactants = [
            Chem.MolFromSmiles("[*]=[*]"),
            Chem.MolFromSmiles("[*][*]=[*][*]"),
            Chem.MolFromSmiles("[*][*]"),
            Chem.MolFromSmiles("[*]#[*]"),
            Chem.MolFromSmiles("[*]=[*][*]"),
            Chem.MolFromSmiles("[*]#[*][*]"),
        ]

        wildcard_products = [
            Chem.MolFromSmiles("[*][*][*]"),
            Chem.MolFromSmiles("[*]=[*][*]"),
            Chem.MolFromSmiles("[*][*]#[*]"),
            Chem.MolFromSmiles("[*]=[*]=[*]"),
            Chem.MolFromSmiles("[*][*]([*])[*]"),
            Chem.MolFromSmiles("[*][*](=[*])[*]"),
            Chem.MolFromSmiles("[*][*@]([*])([*])[*]"),
            Chem.MolFromSmiles("[*]=[*]=[*]=[*]"),
            Chem.MolFromSmiles("[*]=[*][*]#[*]"),
            Chem.MolFromSmiles("[*]#[*][*]#[*]"),
            Chem.MolFromSmiles("[*][*](=[*])[*]=[*]"),
            Chem.MolFromSmiles("[*]=[*][*]=[*]"),
            Chem.MolFromSmiles("[*][*](=[*])[*]#[*]"),
            Chem.MolFromSmarts("[*]=[*]-[!#6]-[*]=[*]"),
            Chem.MolFromSmiles("[*]([*])[*]=[*]"),
            Chem.MolFromSmarts("[*]#[*]-[!#6]-[*]=[*]"),
            Chem.MolFromSmarts("[*]#[*]-[!#6]-[*]#[*]"),
        ]
    wildcard_reactants = sorted(wildcard_reactants, key=lambda x: x.GetNumAtoms(), reverse=True)
    wildcard_products = sorted(wildcard_products, key=lambda x: x.GetNumAtoms(), reverse=True)
    return mol, wildcard_reactants, wildcard_products
def map_atoms(matched_structure, mol):
    p = Chem.AdjustQueryParameters.NoAdjustments()
    p.makeDummiesQueries = True
    query = Chem.AdjustQueryProperties(matched_structure, p) 
    matches = mol.GetSubstructMatches(query)
    if matches:
        return matches
    else:
        print(f"No atom mapping found for substructure: {Chem.MolToSmiles(matched_structure)}")
    return matched_structure
def match_substructure(mol, sub):
    p = Chem.AdjustQueryParameters()
    p.makeDummiesQueries = True
    p.adjustDegree = False
    p.adjustRingCount = False
    p.adjustAromatic = False
    p.adjustValence = False
    substructure = Chem.AdjustQueryProperties(sub, p)


    match = mol.HasSubstructMatch(substructure)
    return match
def append_unique(mol_list, mol):
    mol_smiles = Chem.MolToSmiles(mol)
    mol_list[:] = [m for m in mol_list if Chem.MolToSmiles(m) != mol_smiles]
    mol_list.append(mol)
def product_options(mol, wildcard_products):
    # Define wildcard libraries
    # Sort the wildcard libraries by size (larger first)


    # Search for matching substructures (reactants and products)
    #mol= Chem.RemoveHs(mol)
    matched_products = []
    for product in wildcard_products:
        if match_substructure(mol, product):
            matched_products.append(product)
    # Map atom IDs from reactants/products to the molecule
    mapped_products = []
    for product in matched_products:
        matches = map_atoms(product, mol)  # Get matched atom indices    
        if not matches:
            continue  # Skip if no matches found
        for match in matches:
            # Convert to an editable molecule
            rw_product = Chem.RWMol(product)
            for wild_idx, mol_idx in enumerate(match):  
                product_atom = rw_product.GetAtomWithIdx(wild_idx)  # Get atom in product
                mol_atom = mol.GetAtomWithIdx(mol_idx)  # Get corresponding atom in mol
                if product_atom.GetAtomicNum() == 0:  # If it's a wildcard atom
                    new_atomic_num = mol_atom.GetAtomicNum()  # Get atomic number from mol
                    product_atom.SetAtomicNum(new_atomic_num)  # Assign atomic number
                    product_atom.SetNoImplicit(False)
            try:
                if Chem.SanitizeMol(rw_product):
                    Chem.SanitizeMol(rw_product)
            except Exception as e:
                continue
            mapped_products.append(rw_product.GetMol()) # Convert back to immutable Mol
    #return mapped_products
    unique_products = []
    unique_smiles = set()    
    for prod in mapped_products:
            prod = Chem.AddHs(prod)
            smiles = Chem.MolToSmiles(prod)  # Convert to SMILES for easy comparison
            inchi = Chem.MolToInchi(prod)
            val = get_Hf(inchi)
            if smiles not in unique_smiles:
                if val is not None:
                    unique_smiles.add(smiles)
                    unique_products.append(prod)
    return unique_products
def reactant_options(unique_products, wildcard_reactants):
    matched_reactants = []
    mapped_reactants = []

    for rhs_product in unique_products:  # Iterate over RHS products
        for reactant in wildcard_reactants:  # Iterate over wildcard reactants
            if match_substructure(rhs_product, reactant):  # Check if reactant matches product
                matched_reactants.append((reactant, rhs_product))  # Store pair

    for reactant, rhs_product in matched_reactants:
        matches = map_atoms(reactant, rhs_product)  # Get atom mapping
        if not matches:
            continue  # Skip unmatched reactants

        rw_reactant = Chem.RWMol(reactant)  # Editable reactant molecule
        for match in matches:
            for wild_idx, prod_idx in enumerate(match):  
                reactant_atom = rw_reactant.GetAtomWithIdx(wild_idx)
                product_atom = rhs_product.GetAtomWithIdx(prod_idx)

                if reactant_atom.GetAtomicNum() == 0:  # If it's a wildcard atom
                    new_atomic_num = product_atom.GetAtomicNum()  # Assign atomic number from product
                    reactant_atom.SetAtomicNum(new_atomic_num)
                    reactant_atom.SetNoImplicit(False)
                    new_atomic_num = product_atom.GetAtomicNum()  # Get atomic number from mol
                    reactant_atom.SetAtomicNum(new_atomic_num)  # Assign atomic number
                    reactant_atom.SetNoImplicit(False)
            Chem.SanitizeMol(rw_reactant)
            mapped_reactants.append(rw_reactant.GetMol())
              # Convert back to immutable Mol


    # Ensure uniqueness of mapped reactants
    unique_reactants = []
    unique_smiles = set()

    for reactant in mapped_reactants:
        reactant = Chem.AddHs(reactant)
        smiles = Chem.MolToSmiles(reactant)
        inchi = Chem.MolToInchi(prod)
        val = get_Hf(inchi)
        if smiles not in unique_smiles:
            if val is not None:
                unique_smiles.add(smiles)
                unique_reactants.append(reactant)


    return unique_reactants
def Homodesmotic(input_mol, lhs_required=None, rhs_required=None, Substruct=None, Replacement=None):
    original_mol= input_mol
    x = kekulize_and_make_aliphatic(input_mol)
    input_mol, wildcard_reactants, wildcard_products = get_fragments(x,4)    
    def count_bonds(mol):
        bond_counts = {}
        for bond in mol.GetBonds():
            if bond.GetBeginAtom().GetSymbol()== 'H' or bond.GetEndAtom().GetSymbol()== 'H':
                continue 
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
    
            # Create a tuple where atom order doesn't matter by sorting the atoms based on their symbol and hybridization
            atom_pair = tuple(sorted(
                [(atom1.GetSymbol(), atom1.GetHybridization()), 
                 (atom2.GetSymbol(), atom2.GetHybridization())]
            ))
    
            bond_type = (atom_pair, bond.GetBondType())
            bond_counts[bond_type] = bond_counts.get(bond_type, 0) + 1
    
        return bond_counts
    
    def count_hydrogens(mol):
        #mol = Chem.AddHs(mol)
        hydrogen_counts = {}
        mol = Chem.MolFromSmiles(Chem.MolToSmiles(mol))
        for atom in mol.GetAtoms():
            key = (atom.GetSymbol(), atom.GetHybridization(), atom.GetNumExplicitHs())
            hydrogen_counts[key] = hydrogen_counts.get(key, 0) + 1
        return hydrogen_counts
    def Balance(input_mol):
        prob = LpProblem("Chemical_Balancing", LpMinimize)
        input_mol = kekulize_and_make_aliphatic(input_mol)
        dummy = input_mol
        if Substruct and Replacement:
            if isinstance(Substruct, list) and isinstance(Replacement, list):
                for sub, rep in zip(Substruct, Replacement):
                    dummy = Chem.ReplaceSubstructs(dummy, sub, rep, True)[0]
            elif isinstance(Substruct, Chem.Mol) and isinstance(Replacement, Chem.Mol):
                dummy = Chem.ReplaceSubstructs(dummy, Substruct, Replacement, True)[0]
        rhs_options = product_options(dummy, wildcard_products)
        lhs_addition_options = reactant_options(rhs_options, wildcard_reactants)
        lhs_addition_options.append(input_mol)
        input_index = lhs_addition_options.index(input_mol)
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_addition_options.append(mol)
            elif isinstance (lhs_required, Chem.Mol):
                    lhs_addition_options.append(lhs_required)
                
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    append_unique(rhs_options, mol)
            elif isinstance (rhs_required, Chem.Mol):
                append_unique(rhs_options, rhs_required)
        
        
        x_vars = [LpVariable(f"x_{i}", 0, None, LpInteger) for i in range(len(rhs_options))]
        y_vars = [LpVariable(f"y_{j}", 0, None, LpInteger) for j in range(len(lhs_addition_options))]
        prob += y_vars[input_index] == 1       
        if lhs_required:
            if isinstance (lhs_required, list):
                for mol in lhs_required:
                    lhs_r_idx = lhs_addition_options.index(mol) 
                    prob += y_vars[lhs_r_idx] == 1
            elif isinstance (lhs_required, Chem.Mol):
                lhs_r_idx = lhs_addition_options.index(lhs_required) 
                prob += y_vars[lhs_r_idx] == 1
        if rhs_required:
            if isinstance (rhs_required, list):
                for mol in rhs_required:
                    rhs_r_idx = rhs_options.index(mol) 
                    prob += x_vars[rhs_r_idx] == 1
            elif isinstance (rhs_required, Chem.Mol):
                rhs_r_idx = rhs_options.index(rhs_required) 
                prob += x_vars[rhs_r_idx] == 1
            
        if not rhs_options:
            print("Error: No RHS molecules found.")
            return [], []
    
        bond_types_rhs = [count_bonds(mol) for mol in rhs_options]
        bond_types_lhs = [count_bonds(mol) for mol in lhs_addition_options]
        hydrogen_rhs = [count_hydrogens(mol) for mol in rhs_options]
        hydrogen_lhs = [count_hydrogens(mol) for mol in lhs_addition_options]
        for k in set().union(*bond_types_rhs, *bond_types_lhs):
            prob += lpSum(x_vars[i] * bond_types_rhs[i].get(k, 0) for i in range(len(rhs_options))) == \
                    lpSum(y_vars[j] * bond_types_lhs[j].get(k, 0) for j in range(len(lhs_addition_options)))
    
        for t in set().union(*hydrogen_rhs, *hydrogen_lhs):
            prob += lpSum(x_vars[i] * hydrogen_rhs[i].get(t, 0) for i in range(len(rhs_options))) == \
                    lpSum(y_vars[j] * hydrogen_lhs[j].get(t, 0) for j in range(len(lhs_addition_options)))
        
        solver = PULP_CBC_CMD(msg=False)
        prob.solve(solver)
        #prob.solve()
        solution_x = [int(var.varValue) for var in x_vars]
        solution_y = [int(var.varValue) for var in y_vars]
    
        balanced_rhs = [(rhs_options[i], solution_x[i]) for i in range(len(rhs_options)) if solution_x[i] > 0]
        balanced_lhs = [(lhs_addition_options[j], solution_y[j]) for j in range(len(lhs_addition_options)) if solution_y[j] > 0]
        #print("Solver status:", LpStatus[prob.status])
        return balanced_rhs, balanced_lhs, LpStatus[prob.status]
    return Balance(x)
