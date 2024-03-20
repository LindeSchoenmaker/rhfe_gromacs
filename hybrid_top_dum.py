#! /usr/bin/python
import collections
from collections import Counter
from itertools import chain, compress


def process_atoms(line, in_f, out_file):
    dummies = {'A': [], 'B': []}
    out_file.write(line  + "\n")

    while True:
        line = in_f.readline()
        if line.strip() == "":
            out_file.write("\n")
            break
        elif line.split()[0] == ";":
            out_file.write(line + "\n")
        else:
            nr, typ, resnr, res, atom, cgnr, charge, mass, typB, chargeB, massB = line.strip().split()
            if typ.startswith('DUM_'):
                dummies['A'].append(nr)
            if typB.startswith('DUM_'):
                dummies['B'].append(nr)
            out_file.write(line)

    return dummies


def process_bonds(line, in_f, out_file, dummies, line_num):
    bridge = {'A': [], 'B': []}
    physical = []
    P2s = collections.defaultdict(list)
    P3s = collections.defaultdict(list)
    out_file.write(line + "\n")
    while True:
        line = in_f.readline()
        if line.strip() == "":
            break
        elif line.split()[0] == ";":
            out_file.write(line + "\n")
        else:
            ai, aj, funct, c1, k1, c2, k2 = line.strip().split()[:7]
            if ai in dummies['A'] and aj not in dummies['A']:
                bridge['A'].append([aj, ai])  #physical, dummy
                physical.append(aj)
            elif ai not in dummies['A'] and aj in dummies['A']:
                bridge['A'].append([ai, aj])  #physical, dummy
                physical.append(ai)
            elif ai in dummies['B'] and aj not in dummies['B']:
                bridge['B'].append([aj, ai])  #physical, dummy
            elif ai not in dummies['B'] and aj in dummies['B']:
                bridge['B'].append([ai, aj])  #physical, dummy
            out_file.write(line)

    in_f.seek(line_num)
    while True:
        line = in_f.readline()
        if line.strip() == "":
            break
        elif line.split()[0] == ";":
            continue
        else:
            ai, aj, funct, c1, k1, c2, k2 = line.strip().split()[:7]
            if ai in physical:
                if aj not in dummies['A'] and aj not in dummies['B']:
                    P2s[ai].append(aj)
            elif aj in physical:
                if ai not in dummies['A'] and ai not in dummies['B']:
                    P2s[aj].append(ai)

    in_f.seek(line_num)
    P2s_all = [item for row in P2s.values() for item in row]
    while True:
        line = in_f.readline()
        if line.strip() == "":
            break
        elif line.split()[0] == ";":
            continue
        else:
            ai, aj, funct, c1, k1, c2, k2 = line.strip().split()[:7]
            if aj in P2s_all and ai not in P2s.keys(): # if attached to P2s but not the P1 atom
                P3s[aj].append(ai)
            elif ai in P2s_all and aj not in P2s.keys(): # if attached to P2s but not the P1 atom
                P3s[ai].append(aj)

    return bridge, P2s, P3s

def process_other(line, in_f, out_file):
    out_file.write(line + "\n")
    while True:
        line = in_f.readline().strip()
        if line == "":
            out_file.write("\n")
            break
        elif line.split()[0] == ";":
            out_file.write(line + "\n")
        else:
            out_file.write(line + "\n")
    return


def process_angles(line,
                   in_f,
                   out_file,
                   fc='20.92',
                   bridge_atoms={
                       'A': [["C2", "DC10"]],
                       'B': [["C2", "C4"]]
                   },
                   P2s=[['C5', 'H8', 'H9']],
                   multiplicity_atoms=None):
    """lower force constant of angle dummy bridge atom, physical bridge atom and second physical shell"""
    out_file.write(line + "\n")
    dual_angle = '90.0'  #degree
    dual_fc = '627.6'  #kJoule mol-1 rad-2
    while True:
        line = in_f.readline()
        if line.strip() == "":
            break
        elif line.split()[0] == ";":
            out_file.write(line)
        else:
            ai, aj, ak, funct, c1, k1, c2, k2 = line.strip().split()[:8]
            comment = " ".join(line.split()[8:])
            atoms = [ai, aj, ak]
            for endstate in ['A', 'B']:
                for i in range(len(bridge_atoms[endstate])):
                    if bridge_atoms[endstate][i][0] in atoms and bridge_atoms[
                            endstate][i][1] in atoms:
                        if any(ext in atoms
                               for ext in P2s[bridge_atoms[endstate][i][0]]):
                            # for triple junction set force constant of angle low
                            if len(P2s[bridge_atoms[endstate][i][0]]) == 3:
                                if endstate == 'A':
                                    line = ai.rjust(6) + aj.rjust(
                                        7) + ak.rjust(6) + funct.rjust(
                                            8) + c1.rjust(15) + fc.rjust(
                                                10) + c2.rjust(20) + k2.rjust(
                                                    15) + comment + "\n"
                                elif endstate == 'B':
                                    line = ai.rjust(6) + aj.rjust(
                                        7) + ak.rjust(6) + funct.rjust(
                                            8) + c1.rjust(15) + k1.rjust(
                                                15) + c2.rjust(15) + fc.rjust(
                                                    10) + comment + "\n"
                            # for dual junction use high force constant and angle of 90 degrees
                            elif len(P2s[bridge_atoms[endstate][i][0]]) == 2:
                                if multiplicity_atoms is not None:
                                    if any(set(atoms).issubset(set(ext)) for ext in multiplicity_atoms):  # at idx 2 P2
                                        continue  # if D1, P1, P2 in angle to keep, do not modify
                                    else:
                                        if endstate == 'A':
                                            line = ai.rjust(6) + aj.rjust(
                                                7) + ak.rjust(6) + funct.rjust(
                                                    8) + c1.rjust(15) + "0.000".rjust(
                                                        10) + c2.rjust(20) + k2.rjust(
                                                            15) + comment + "\n"
                                        elif endstate == 'B':
                                            line = ai.rjust(6) + aj.rjust(
                                                7) + ak.rjust(6) + funct.rjust(
                                                    8) + c1.rjust(15) + k1.rjust(
                                                        15) + c2.rjust(15) + "0.000".rjust(
                                                            10) + comment + "\n"
                                else:
                                    if endstate == 'A':
                                        line = ai.rjust(6) + aj.rjust(7) + ak.rjust(
                                            6) + funct.rjust(8) + dual_angle.rjust(
                                                9) + dual_fc.rjust(16) + c2.rjust(
                                                    20) + k2.rjust(
                                                        15) + comment + "\n"
                                    elif endstate == 'B':
                                        line = ai.rjust(6) + aj.rjust(
                                            7) + ak.rjust(6) + funct.rjust(
                                                8) + c1.rjust(15) + k1.rjust(
                                                    15) + dual_angle.rjust(
                                                        9) + dual_fc.rjust(
                                                            16) + comment + "\n"
            out_file.write(line)

    out_file.write("\n\n")
    return


def process_dhedrals(line, in_f, out_file, P2s=['H8'], multiplicity_atoms=None):
    "remove dihedral angles P3P2P1D1, only keep selected dihedrals of P2P1D1D2"
    multiplicity_atom_sets = []
    if multiplicity_atoms:
        multiplicity_atom_sets = [set(atoms) for atoms in multiplicity_atoms]
    out_file.write(line + "\n")
    while True:
        line = in_f.readline()
        if line.strip() == "":
            out_file.write("\n")
            break
        if line.split()[0] == ";":
            out_file.write(line.strip() + "\n")
        else:
            ai, aj, ak, al, funct, a1, fc1, f1, a2, fc2, f2 = line.strip(
            ).split()[0:11]
            comment = " ".join(line.split()[11:])
            atoms = [ai, aj, ak, al]
            first = line.split()[-1].split('->')[0]
            if first.count("A") > 1 and "D" in first:
                if first.count("D") == 1:
                    if set(atoms) in multiplicity_atom_sets and fc1 != "0": # also apply if exact same atoms also have other dihedral term
                        line = ai.rjust(6) + aj.rjust(7) + ak.rjust(7) + al.rjust(7) + funct.rjust(5) + " " + a1 + " " + "420" + " " + "1" + " " + a2 + " " + fc2 + " " + "1" + " " + comment + "\n"
                    else:
                        line = ai.rjust(6) + aj.rjust(7) + ak.rjust(
                            7
                        ) + al.rjust(7) + funct.rjust(
                            5
                        ) + " " + a1 + " " + "0" + " " + f1 + " " + a2 + " " + fc2 + " " + f2 + " " + comment + "\n"
                elif any(ext in atoms for ext in P2s):
                    if first.count("A") != 2:
                        line = ai.rjust(6) + aj.rjust(7) + ak.rjust(
                            7
                        ) + al.rjust(7) + funct.rjust(
                            5
                        ) + " " + a1 + " " + "0" + " " + f1 + " " + a2 + " " + fc2 + " " + f2 + " " + comment + "\n"
                else:
                    line = ai.rjust(6) + aj.rjust(7) + ak.rjust(7) + al.rjust(
                        7
                    ) + funct.rjust(
                        5
                    ) + " " + a1 + " " + "0" + " " + f1 + " " + a2 + " " + fc2 + " " + f2 + " " + comment + "\n"

            second = line.split()[-1].split('->')[1]
            if second.count("A") > 1 and "D" in second:
                if second.count("D") == 1:
                    if set(atoms) in multiplicity_atom_sets and fc2 != "0": # also apply if exact same atoms also have other dihedral term
                        line = ai.rjust(6) + aj.rjust(7) + ak.rjust(7) + al.rjust(7) + funct.rjust(5) + " " + a1 + " " + fc1 + " " + "1" + " " + a2 + " " + "420" + " " + "1" + " " + comment + "\n"
                    else:
                        line = ai.rjust(6) + aj.rjust(7) + ak.rjust(
                            7
                        ) + al.rjust(7) + funct.rjust(
                            5
                        ) + " " + a1 + " " + fc1 + " " + f1 + " " + a2 + " " + "0" + " " + f2 + " " + comment + "\n"
                elif any(ext in atoms for ext in P2s):
                    if second.count("A") != 2:
                        line = ai.rjust(6) + aj.rjust(7) + ak.rjust(
                            7
                        ) + al.rjust(7) + funct.rjust(
                            5
                        ) + " " + a1 + " " + fc1 + " " + f1 + " " + a2 + " " + "0" + " " + f2 + " " + comment + "\n"
                else:
                    line = ai.rjust(6) + aj.rjust(7) + ak.rjust(7) + al.rjust(
                        7
                    ) + funct.rjust(
                        5
                    ) + " " + a1 + " " + fc1 + " " + f1 + " " + a2 + " " + "0" + " " + f2 + " " + comment + "\n"
            out_file.write(line)
    return


def process_file(in_file, out_file, decouple_params):
    with open(out_file, 'w') as out_f:
        with open(in_file) as in_f:
            for line in iter(in_f.readline, ''):
                if line.strip() == "[ moleculetype ]":
                    print("moleculetype section started!")
                    process_other(line.strip(), in_f, out_f)
                    print("moleculetype section done!")
                elif line.strip() == "[ atoms ]":
                    print("atoms section started!")
                    dummies = process_atoms(line.strip(), in_f, out_f)
                    print("atoms section done!")
                elif line.strip() == "[ bonds ]":
                    print("bonds sectrion started!")
                    bridge, P2s, P3s = process_bonds(line.strip(), in_f, out_f,
                                                     dummies, in_f.tell())
                    print("bonds sectrion done!")
                elif line.strip() == "[ pairs ]":
                    print("pairs section started!")
                    process_other(line.strip(), in_f, out_f)
                    print("pairs section done!")
                elif line.strip() == "[ angles ]":
                    print("angles section started!")
                    # get P2s that do not overlap between attachment point dummies
                    bonds_list = [P2s[key] for key in P2s]
                    freq = Counter(chain.from_iterable(bonds_list))
                    res = {idx for idx in freq if freq[idx] == 1}
                    bonds_list = [[x for x in bond_list if x in res]
                                  for bond_list in bonds_list]

                    multiplicity_atoms = None
                    if not decouple_params['90']:
                        # for dual anchored save dihedrals
                        multiplicity_atoms = []
                        dual_idx = [len(x) == 2 for x in P2s.values()]
                        P1s_dual = list(compress(list(P2s.keys()), dual_idx))
                        for loc in bridge['A']:
                            if loc[0] not in P1s_dual: continue
                            P3 = None
                            D1_A = loc[1]
                            P1 = loc[0]
                            P2_atoms = P2s[P1]
                            for P2_atom in P2_atoms:
                                if P2_atom not in res:
                                    continue
                                if P2_atom in list(P3s.keys()):
                                    P2 = P2_atom
                                    P3 = P3s[P2][0]
                                    multiplicity_atoms.append(
                                        (D1_A, P1, P2, P3))
                                    D1_B = [
                                        x[1] for x in bridge['B']
                                        if x[0] == P1
                                    ][0]
                                    multiplicity_atoms.append(
                                        (D1_B, P1, P2, P3))
                                    break
                            if P3 is None:
                                for P2_atom in P2_atoms:
                                    if P2_atom in list(P3s.keys()):
                                        P2 = P2_atom
                                        P3 = P3s[P2][0]
                                        multiplicity_atoms.append(
                                            (D1_A, P1, P2, P3))
                                        D1_B = [
                                            x[1] for x in bridge['B']
                                            if x[0] == P1
                                        ][0]
                                        multiplicity_atoms.append(
                                            (D1_B, P1, P2, P3))
                                        break
                            if P3 is None:
                                raise KeyError(
                                    "No dihedrals found to anchor dual dummy point, use 90 degree setting instead"
                                )
                    process_angles(line.strip(),
                                   in_f,
                                   out_f,
                                   fc=decouple_params['fc'],
                                   bridge_atoms=bridge,
                                   P2s=P2s,
                                   multiplicity_atoms=multiplicity_atoms)
                    print("angles section done!")
                elif line.strip() == "[ dihedrals ]":
                    print("dihedrals section started!")
                    dihedral_excpt = [suitable[0] for suitable in bonds_list]
                    process_dhedrals(line.strip(),
                                     in_f,
                                     out_f,
                                     P2s=dihedral_excpt,
                                     multiplicity_atoms=multiplicity_atoms)
                    print("dihedrals section done!")
                elif line.strip() == "[ cmap ]":
                    print("cmap section started!")
                    process_other(line.strip(), in_f, out_f)
                    print("cmap section done!")
                elif line.strip() == "":
                    out_f.write(line)
        in_f.close()
    out_f.close()

if __name__ == "__main__":
    to__int = {'fc':'20.92',
               '90': False,
               'name': 'to__int',
              'P2s_dihedral':['H8'],
              'bridge_atoms': {'A': [["C4", "DC10"]],
                               'B': [["C4", "C2"]]},
              'P2s_angle': [['C5', 'H8', 'H9']]}

    ref_int = {'fc':'20.92',
               '90': False,
               'name': 'ref_int',
              'P2s_dihedral':['O3', 'C2'],
              'bridge_atoms': {'A': [["C1", "HV3"], ['C7', 'DN1']],
                               'B': [["C1", "O1"], ['C7', 'N1']]},
              'P2s_angle': [['C2', 'C6'], ['O3', 'C6']]}

    to__ref = {'fc':'20.92',
               '90': False,
               'name': 'to__ref',
              'P2s_dihedral':['H8', 'O1', 'C11'],
              'bridge_atoms': {'A': [["C4", "DC14"], ["C10", "DO1"], ['C12', 'DN1']],
                               'B': [["C4", "C2"], ["C10", "H12"], ['C12', 'N3']]},
              'P2s_angle': [['C5', 'H8', 'H9'], ['C9', 'C11'], ['O1', 'C9']]}

    import json
    for params in [ref_int, to__int, to__ref]: #[ref_int, to__int, to__ref]:
        with open(f'{params["name"]}.json', "w") as outfile:
            json.dump(params, outfile)
        process_file(f'merged_tmp_{params["name"]}.itp', f'merged_tmp_{params["name"]}_decoupled_new_90.itp', decouple_params = params)
