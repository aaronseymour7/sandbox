import recurs
import homodes

def Combo(input_mol)
    missing_inchis = []
    
    reaction = recurs(input_mol)
    if reaction is None:
      reaction = homodes(input_mol)
    if reaction.len()!=3:
      print('Infeasible')
      return None
    rhs, lhs, status = reaction
    if status != 'Optimal':
      print('Infeasible')
      return None
    for mol, coeff in lhs[:-1]:
      inchi = Chem.MolToInchis(mol)
      Hf = get_Hf(inchi)
      if Hf is None:
        missing_inchis.append(inchi)
    for mol, coeff in rhs:
      inchi = Chem.MolToInchis(mol)
      Hf = get_Hf(inchi)
      if Hf is None:
        missing_inchis.append(inchi)
    
    for i in missing_inchis:
      feas = recurs(Chem.MolFromInchi(i))
      if feas is None:
        print("Uncalculable")
        return None
      print(f"Homodesmotic can solve for {i}")
  return rhs,lhs,status
