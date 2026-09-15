"""
pkasolver site classification for the MP-prediction dataset.

Each ionisable site returned by pkasolver is classified as ACID or BASE by
matching the reaction centre atom against a sourced SMARTS pattern, falling
back to formal charge change when no pattern matches.

Two pKa selections are produced per molecule:

  v1  all sites classified ACID/BASE with 3 <= pKa <= 11
  v2  the same, restricted to sites whose functional group was identified by
      SMARTS (i.e. the charge fallback did not decide it)

`process_drug_full` is the entry point used by 02_pkasolver_logd.Rmd. It runs
`calculate_microstate_pka_values` exactly once per molecule and returns the raw
site list together with both selections.
"""

from rdkit import Chem
from pkasolver.query import calculate_microstate_pka_values


# ACIDIC FUNCTIONAL GROUP PATTERNS
# Sources:
#   [RDKit]     - rdkit.Chem.FunctionalGroups.BuildFuncGroupHierarchy()
#   [SMARTS-RX] - Journal of Cheminformatics (2025), MolecularAI/smartsrx GitHub
#   [Brenk]     - Brenk et al., ChemMedChem 2008, Supplementary Table S1
#   [Daylight]  - Daylight SMARTS Tutorial (official SMARTS specification)


ACIDIC_PATTERNS = {
    "carboxylate": "C(=O)[O;H,-]",  
    "phenol": "[O;H1;$(O-!@c)]",  
    "sulfonamide": "[N;$(NS(=O)(=O))]",  
    "heterocyclic_N_acidic": "[n;D2;H1;$(n1cccc1)]",  
    "hydroxamate": "C(=O)N[OH]",  
    "phosphate": "[O;$([O;D1]P(=O));!$(OP(=O)(=O))]", 
    "tetrazole": "[nH]1nnnc1", 
    "thiol": "[SH]", 
    "alcohol": "[O;H1;$(O-!@[C;!$(C=!@[O,N,S])])]",  
    "acidic_amide": "[CX3](=O)[NX3H2,NX3H1]",  
    "acidic_aniline": "[NX3H2,NX3H1][cX3]:[c]", 
    "carbamate": "[N;!$(NC(=O)OC(C)(C)C);!$(NC(=O)OC1c2ccccc2-c3c1cccc3);$(NC(=O)[O;D2])]",  
    "hydrazide": "C(=O)N[NH2]", 
    "imide": "[N;$(N(C(=O)[#6])C(=O)[#6])]", 
    "sulfate": "OS(=O)(=O)[O-]", 
    "carbonic_acid": "[CX3](=[OX1])([OX2H])[OX2H]",  
}


# BASIC FUNCTIONAL GROUP PATTERNS

BASIC_PATTERNS = {
    "heterocyclic_N_basic_6ring": "[n;D2;$(n1ccccc1)]", 
    "heterocyclic_N_basic_polyhet": "[n;D2;$(n1aaaaa1);!$(n1ccccc1)]", 
    "heterocyclic_N_basic_5ring_polyhet": "[n;H0;$(n1aaaa1);!$(n1cccc1)]",  
    "aliphatic_amine": "[N;$(N-[#6]);!$(N-[!#6;!#1]);!$(N-C=[O,N,S])]",  
    "guanidine": "[N;$(NC(~N)=N)]", 
    "amidine": "[N;!R;$(N([#6])=[C;$(C(~N)~N);!$(C(~N)(~N)~N)]),$([N;D1]=[C;$(C(~N)~N);!$(C(~N)(~N)~N)])]",
    "aniline": "c1cc([NH2])ccc1",  
}

PKA_MIN = 3.0
PKA_MAX = 11.0


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
    """Return (type, pka, functional_group, method) for one pkasolver state."""
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


def _empty_result():
    return {
        "pKa_pkasolver": None,
        "N_pka_pkasolver": 0,
        "Type_pkasolver": None,
        "pKa_acid_v1": None,
        "pKa_base_v1": None,
        "pKa_acid_v2": None,
        "pKa_base_v2": None,
        "FG_acid_v2": None,
        "FG_base_v2": None,
        "All_sites_detail": None,
    }


def _select(sites, site_type, require_named_fg):
    """Acid sites take the lowest pKa, base sites the highest."""
    hits = [
        (pka, fg)
        for (t, pka, fg, m) in sites
        if t == site_type
        and PKA_MIN <= pka <= PKA_MAX
        and (fg != "unknown" if require_named_fg else True)
    ]
    if not hits:
        return None, None
    chosen = min(hits, key=lambda x: x[0]) if site_type == "ACID" else max(hits, key=lambda x: x[0])
    return chosen[0], chosen[1]


def process_drug_full(smiles):
    """
    Classify every ionisable site of one molecule in a single pkasolver pass.

    Returns a dict with the raw site list (pKa/type strings and count) plus the
    v1 and v2 acid/base selections and the v2 functional groups. Returns the
    empty result - not None  when the molecule or pkasolver fails, so a failed
    row is never silently dropped from the dataset.
    """
    mol = Chem.MolFromSmiles(smiles) if smiles else None
    if mol is None:
        return _empty_result()

    try:
        states = calculate_microstate_pka_values(mol, only_dimorphite=False)
    except Exception:
        return _empty_result()

    sites = [classify_site(s) for s in states]
    if not sites:
        return _empty_result()

    acid_v1, _ = _select(sites, "ACID", require_named_fg=False)
    base_v1, _ = _select(sites, "BASE", require_named_fg=False)
    acid_v2, fg_acid_v2 = _select(sites, "ACID", require_named_fg=True)
    base_v2, fg_base_v2 = _select(sites, "BASE", require_named_fg=True)

    return {
        "pKa_pkasolver": ", ".join(str(round(pka, 3)) for (t, pka, fg, m) in sites),
        "N_pka_pkasolver": len(sites),
        "Type_pkasolver": ", ".join(t for (t, pka, fg, m) in sites),
        "pKa_acid_v1": acid_v1,
        "pKa_base_v1": base_v1,
        "pKa_acid_v2": acid_v2,
        "pKa_base_v2": base_v2,
        "FG_acid_v2": fg_acid_v2,
        "FG_base_v2": fg_base_v2,
        "All_sites_detail": str(sites),
    }


def process_drug(smiles):
    """Single molecule inspection helper - same numbers, printed per site."""
    result = process_drug_full(smiles)
    if result["N_pka_pkasolver"] == 0:
        return None
    return result
