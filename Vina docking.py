import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import vina
import meeko


def read_pdb(pdb_file):
    # pdb_file = open(pdb_file, "r").read()
    # prot = Chem.MolFromMolBlock(pdb_file, sanitize=False, removeHs=False)
    prot = Chem.MolFromPDBFile(
        pdb_file,
        removeHs=False,
        sanitize=False,
        proximityBonding=True)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETAROMATICITY)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETCONJUGATION)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SETHYBRIDIZATION)
    Chem.SanitizeMol(prot, Chem.SANITIZE_SYMMRINGS)
    # Chem.SanitizeMol(prot, Chem.SANITIZE_PROPERTIES)
    Chem.SanitizeMol(prot, Chem.SANITIZE_CLEANUP)
    # prot = Chem.AddHs(prot, addCoords=True)
    return prot


def parse_poses(poses):
    pmol = meeko.PDBQTMolecule(poses)
    return meeko.RDKitMolCreate.from_pdbqt_mol(pmol)[0]


ligs = {
    # "CCR2RAR": Chem.MolFromSmiles("Clc1ccc(c(c1)F)N1[C@H](C2CCCCC2)C(=C(C1=O)[O-])C(=O)C"),
    # "compound39": Chem.MolFromSmiles("[H]c1c([H])c(C([H])([H])c2c(C3([H])C([H])([H])C3([H])[H])nc3nc(N([H])[H])nn3c2[O-])c([H])c(Cl)c1Cl"),
    # "SD24": Chem.MolFromSmiles("c1cccc(c1C([O-])=O)Oc(ccc(Cl)c2)c2[N-]S(=O)(=O)c(c3)ccc(Cl)c3Cl"),
    # "SD24para": Chem.MolFromSmiles("c1cc(C([O-])=O)ccc1Oc(ccc(Cl)c2)c2[N-]S(=O)(=O)c(c3)ccc(Cl)c3Cl"),
    # "JNJ27141491": Chem.MolFromSmiles("FC1=CC([C@H](CC)N2C(C(OC)=O)=C(C3=CC=NO3)[N-]C2=S)=CC=C1F"),
    # "CCX140": Chem.MolFromSmiles("Cc1cnc(C(=O)c2ncnc3[nH]ccc23)c([N-]S(=O)(=O)c2ccc(Cl)c(C(F)(F)F)c2)c1"),
}

for lig_name, lig in ligs.items():
    lig = Chem.AddHs(lig)
    Chem.AllChem.EmbedMolecule(lig, randomSeed=42)
    meeko_prep = meeko.MoleculePreparation()
    meeko_prep.prepare(lig)
    lig_pdbqt = meeko_prep.write_pdbqt_string()
    v = vina.Vina(sf_name='vina', seed=42)
    v.set_receptor('5T1A_clean_mutations_reversed_withHs.pdbqt')
    v.set_ligand_from_string(lig_pdbqt)
    previous_ligand = next(Chem.SDMolSupplier('5t1a_C_VT5.sdf'))  # coords from 5t1a
    centroid = Chem.rdMolTransforms.ComputeCentroid(previous_ligand.GetConformer())
    # v.compute_vina_maps(center=[centroid.x, centroid.y,
    #                             centroid.z], box_size=[30, 30, 30])
    v.compute_vina_maps(
        # this box is used in the drugex project
        center=[5.1, 28.0, 187.6],
        box_size=[16.2, 17.8, 17.4]
    )
    v.dock(exhaustiveness=8)
    output_pdbqt = v.poses()
    with open(f"{lig_name}_poses.pdbqt", "w") as pout:
        pout.write(output_pdbqt)
    poses = parse_poses(output_pdbqt)
    f = Chem.SDWriter(f'{lig_name}_poses.sdf')
    for pose_id in range(poses.GetNumConformers()):
        # conf = poses.GetConformer(pose_id)
        # conf = Chem.Mol(conf, 0)
        f.write(poses, confId=pose_id)
    f.close()
    # combine into a complex
    for idx, pose in enumerate(Chem.SDMolSupplier(f'{lig_name}_poses.sdf')):
        protein = read_pdb('5T1A_clean_mutations_reversed_withHs.pdb')
        complex = Chem.CombineMols(protein, pose)
        Chem.MolToPDBFile(complex, lig_name + f"_complex_{idx}.pdb", flavor=4)
