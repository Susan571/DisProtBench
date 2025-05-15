from pymol import cmd
import sys


def get_interface(file:str, ligchain:str, dist:float=5.0) -> set:
    '''
    Returns which receptor residues are within 5 angstroms of the ligand chain
        All files should be built completely, and have a numeric residue names, and only protein chains
    file:str - PDB file of the complex
    ligchain: str Ligand chain e.g "B"
    dist:float = 5.0/8/10  Distance cutoff to use for determining if the receptor is the "interface"

    returns: A set of strings defining which residues of the receptor are considered the interface
        example {"44_A", "45_A"} #Residues 44 and 45 of chain A
    '''

    cmd.delete("all") #Clear PyMOL workspace
    cmd.remove("all") 
    cmd.load(file, "myfile") #Load the file
    cmd.remove("ele H")
    
    #Select all receptor residues within cutoff of lig (IDP) chain
    #Also select only the CA atom so there's one atom selected for each residue
    selection = "(byres( (not chain {0}) within {1} of chain {0} )) and name CA"
    selection = selection.format(ligchain, dist) 
    model= cmd.get_model(selection)
    
    #Return set of receptor residues (in the form of 'resi_chain')
    return {"{}_{}".format(at.resi, at.chain) for at in model.atom}


# Receptor recall (we used to call this finterface)
def get_finterface(truth_interface:set, target_interface:set) -> float:
    '''
    Calculates the receptor recall from the receptor interfaces of the ground truth (target)
        REQUIRES TRUTH AND PREDICTION TO HAVE IDENTICAL CHAIN AND RESIDUE NUMBERING
    and the prediction (target)
    '''
    made_interface = truth_interface.intersection(target_interface)
    return float(len(made_interface) / len(truth_interface))


# Receptor precision (we used to call this fpredistrue_iface)
def get_fpredistrue_iface (truth_interface:set, target_interface:set)-> float:
    '''
    Calculates the receptor precision from the receptor interfaces of the ground truth (target)
        REQUIRES TRUTH AND PREDICTION TO HAVE IDENTICAL CHAIN AND RESIDUE NUMBERING
    and the prediction (target)
    '''
    made_interface = target_interface.intersection(truth_interface)
    return float(len(made_interface) / len(target_interface)) if len(target_interface) > 0 else 0.0


# Ligand recall (used to call this flig_atinterface)
def get_flig_atinterface(file:str, chain:str , truth_interface:set, pred_contacting_resis:set, dist=5.0) -> float:
    '''
    Calculates the fraction of the predicted ligand that exists on the correct receptor interface
        REQUIRES TRUTH AND PREDICTION TO HAVE IDENTICAL CHAIN AND RESIDUE NUMBERING

    file:str Predicted complex pdb file path
    chain: str Ligand chain e.g "B"
    truth_interface:set The receptor interface residues of the true complex 
    pred_contacting_resis:list List of all residues within {dist} of the receptor
        You'll have to calculate this separately. We don't use the total residue count to account for some ligand residues not being at the interface
    dist:float = 5.0 Distance for determining a contact
    '''

    #Clear pymol workspace
    cmd.delete("all")
    cmd.remove("all")
    cmd.load(file, "mod")
    cmd.remove("ele H")

    #Create a pymol selection string that selects all the true receptor interface residues
    true_interface_sele = "None "
    for resi in truth_interface:
        true_interface_sele = true_interface_sele + " or (resi {} and chain {}) ".format(*resi.split("_"))

    # Selects all the predicted ligand residues that are within 5 angstroms of the true receptor interface residues
    selection = "(byres ( chain {} within {} of ({}))) and name CA".format(chain, dist, true_interface_sele)
    interface_resi_count = len(list(cmd.get_model(selection).atom))
    
    #Divide the count of ligand residues at the true interface by the length of all contacting ligand residues.
    idp_resi_count =  len(pred_contacting_resis) 
    return float(interface_resi_count) / float(idp_resi_count) if idp_resi_count > 0 else 0.0

#Ligand precision (used to call this flig_atpredinterface)
def get_flig_atpredinterface(file:str, chain:str , pred_interface:set, truth_contacting_resis:set, dist:float=5.0) -> float:
    '''
    Calculates the fraction of the predicted ligand that exists on the correct receptor interface
        REQUIRES TRUTH AND PREDICTION TO HAVE IDENTICAL CHAIN AND RESIDUE NUMBERING

    file:str True complex pdb file path
    chain: str Ligand chain e.g "B"
    pred_interface:set The receptor interface residues of the predicted complex 
    truth_contacting_resis:list List of all residues within {dist} of the receptor
        You'll have to calculate this separately. We don't use the total residue count to account for some ligand residues not being at the interface
    dist:float = 5.0 Distance for determining a contact
    '''

    #Clear pymol workspace
    cmd.delete("all")
    cmd.remove("all")
    cmd.load(file, "mod")
    cmd.remove("ele H")

    #Create a pymol selection string that selects all the predicted receptor interface residues
    pred_interface_sele = "None "
    for resi in pred_interface:
        pred_interface_sele = pred_interface_sele + " or (resi {} and chain {}) ".format(*resi.split("_"))

    # Selects all the true ligand residues that are within 5 angstroms of the predicted receptor interface residues
    selection = "(byres ( chain {} within {} of ({}))) and name CA".format(chain, dist, pred_interface_sele)
    interface_resi_count = len(list(cmd.get_model(selection).atom))
    
    idp_resi_count =  len(truth_contacting_resis)
    f_atinterface = float(interface_resi_count) / float(idp_resi_count) #True Contacting ligand residues at predicted receptor interface / all contacting residues
    return f_atinterface

def get_contacting_resis(file:str, ligchain:str, dist=5.0) -> set:
    '''
    Gets which ligand residues are within {dist} of the interface
        file:str complex pdb file path
        chain: str Ligand chain e.g "B"
        dist:float = 5.0 Distance for determining a contact
    '''
    cmd.delete("all")
    cmd.remove("all")
    cmd.load(file, "myfile")
    cmd.remove("ele H")
    selection = "(byres( (chain {0}) within {1} of (not chain {0}) )) and name CA"
    selection = selection.format(ligchain, dist)
    model= cmd.get_model(selection)

    return {"{}_{}".format(at.resi, at.chain) for at in model.atom}


def get_metrics(true_path:str, predicted_path:str, ligand_chain:str) -> dict:
    true_receptor_interface:set = get_interface(true_path, ligand_chain) # Residues of the receptor within 5 angstroms of the true ligand
    pred_receptor_interface:set = get_interface(predicted_path, ligand_chain) # Residues of the receptor within 5 angstroms of the predicted ligand

    #What fraction of the true interface is captured by the prediction interface?
    receptor_recall:float = get_finterface(true_receptor_interface, pred_receptor_interface)
    #What fraction of the predicted interface is actually the true interface?
    receptor_precision:float = get_fpredistrue_iface(true_receptor_interface, pred_receptor_interface)

    #These are used to account for the fact that some ligand residues aren't making contact with the receptor at all
    true_contacting_ligand_resis:set = get_contacting_resis(true_path, ligand_chain) #Residues of the true ligand making contact with the receptor
    pred_contacting_ligand_resis:set = get_contacting_resis(predicted_path, ligand_chain) #Residues of the predicted ligand making contact with the receptor

    #What fraction of the predicted (contacting) ligand residues are at the true interface
    ligand_recall:float = get_flig_atinterface(predicted_path, ligand_chain,
                                               true_receptor_interface, 
                                               pred_contacting_ligand_resis)
    
    #What fraction of the true (contacting) ligand residues are at the predicted interface
    ligand_precision:float = get_flig_atpredinterface(true_path, ligand_chain,
                                                       pred_receptor_interface,
                                                       true_contacting_ligand_resis)

    metrics = {"receptor_recall":receptor_recall, 
               "receptor_precision":receptor_precision, 
               "ligand_recall":ligand_recall,
               "ligand_precision":ligand_precision}
    
    return metrics


def main():
    true_path =sys.argv[1]
    predicted_path =sys.argv[2]
    ligand_chain = sys.argv[3]
    metrics:dict = get_metrics(true_path, predicted_path, ligand_chain)

    for metric, val in metrics.items():
        print(f"{metric.ljust(20)}:\t{val}")


if __name__ == "__main__":
    main()