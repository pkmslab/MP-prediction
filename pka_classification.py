from rdkit import Chem
from pkasolver.query import calculate_microstate_pka_values

ACIDIC_PATTERNS = {
    "carboxylate": "[CX3](=O)[OX1H0-,OX2H1]",
    "phenol": "[OX2H][cX3]:[c]",
    "sulfonamide": "[SX4](=[OX1])(=[OX1])[NX3]",
    "heterocyclic_N_acidic": "[nX3;H1]",
    "hydroxamate": "[CX3](=O)[NX3][OX2H]",
    "carbon_acid": "[CX4;H1,H2,H3][CX3](=O)",
    "phosphate": "[PX4](=O)([OX2H])",
    "tetrazole": "c1nnn[nH]1",
    "thiol": "[SX2H]",
    "alcohol": "[CX4][OX2H]",
    "acidic_amide": "[CX3](=O)[NX3H2,NX3H1]",
    "acidic_aniline": "[NX3H2,NX3H1][cX3]:[c]",
    "carbamate": "[NX3][CX3](=O)[OX2]",
    "hydrazide": "[NX3][NX3][CX3](=O)",
    "imide": "[CX3](=O)[NX3][CX3](=O)",
    "sulfate": "[OX2][SX4](=O)(=O)[OX2]"
}

BASIC_PATTERNS = {
    "heterocyclic_N_basic": "[nX2]",
    "aliphatic_amine": "[NX3;H2,H1,H0;!$(NC=O);!$(N=*);!$(N-a)]",
    "guanidine": "[NX3][CX3](=[NX2])[NX3]",
    "amidine": "[NX3][CX3]=[NX2]",
    "aniline": "[NX3H2,NX3H1][cX3]:[c]",
    "basic_amide": "[NX3][CX3]=[NX2]"
}

def identify_functional_group(mol, atom_idx, patterns):
    for fg_name, smarts in patterns.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern is None:
            continue
        matches = mol.GetSubstructMatches(pattern)
        for match in matches:
            if atom_idx in match:
                return fg_name
    return None

def classify_site(state):
    idx = state.reaction_center_idx
    pka = state.pka
    prot_charge = state.protonated_mol.GetAtomWithIdx(idx).GetFormalCharge()
    deprot_charge = state.deprotonated_mol.GetAtomWithIdx(idx).GetFormalCharge()

    fg_acid = identify_functional_group(state.protonated_mol, idx, ACIDIC_PATTERNS)
    fg_base = identify_functional_group(state.deprotonated_mol, idx, BASIC_PATTERNS)

    if fg_acid is not None:
        return ("ACID", pka, fg_acid, "SMARTS")
    elif fg_base is not None:
        return ("BASE", pka, fg_base, "SMARTS")
    else:
        if prot_charge == 0 and deprot_charge < 0:
            return ("ACID", pka, "unknown", "charge_fallback")
        elif prot_charge > deprot_charge and deprot_charge == 0:
            return ("BASE", pka, "unknown", "charge_fallback")
        else:
            return ("OTHER", pka, "unknown", "charge_fallback")

def process_drug(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None

    try:
        states = calculate_microstate_pka_values(mol, only_dimorphite=False)
    except Exception:
        return None

    all_sites = [classify_site(s) for s in states]

    approach1_acid = [pka for (t, pka, fg, m) in all_sites if t == "ACID" and 3 <= pka <= 11]
    approach1_base = [pka for (t, pka, fg, m) in all_sites if t == "BASE" and 3 <= pka <= 11]
    a1_acid_final = min(approach1_acid) if approach1_acid else None
    a1_base_final = max(approach1_base) if approach1_base else None

    approach2_acid = [pka for (t, pka, fg, m) in all_sites
                       if t == "ACID" and 3 <= pka <= 11 and fg != "unknown"]
    approach2_base = [pka for (t, pka, fg, m) in all_sites
                       if t == "BASE" and 3 <= pka <= 11 and fg != "unknown"]
    a2_acid_final = min(approach2_acid) if approach2_acid else None
    a2_base_final = max(approach2_base) if approach2_base else None

    return {
        "Approach1_pKa_acid": a1_acid_final,
        "Approach1_pKa_base": a1_base_final,
        "Approach2_pKa_acid": a2_acid_final,
        "Approach2_pKa_base": a2_base_final,
        "All_sites_detail": str(all_sites)
    }
