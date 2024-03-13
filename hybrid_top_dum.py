#! /usr/bin/python

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
                   P2s=[['C5', 'H8', 'H9']]):
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
            atoms = comment.split(' ')[1:]
            atoms = [atom.rstrip('x') for atom in atoms]
            for endstate in ['A', 'B']:
                for i in range(len(bridge_atoms[endstate])):
                    if bridge_atoms[endstate][i][0] in atoms and bridge_atoms[
                            endstate][i][1] in atoms:
                        if any(ext in atoms for ext in P2s[i]):
                            # for triple junction set force constant of angle low
                            if len(P2s[i]) == 3:
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
                            elif len(P2s[i]) == 2:
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


def process_dhedrals(line, in_f, out_file, P2s=['H8']):
    "remove dihedral angles P3P2P1D1, only keep selected dihedrals of P2P1D1D2"
    out_file.write(line + "\n")
    while True:
        line = in_f.readline()
        if line.strip() == "":
            out_file.write("\n")
            break
        if line.split()[0] == ";":
            out_file.write(line.strip() + "\n")
        else:
            altered = False
            ai, aj, ak, al, funct, a1, fc1, f1, a2, fc2, f2 = line.strip(
            ).split()[0:11]
            comment = " ".join(line.split()[11:])
            atoms = comment.split(' ')[1:5]
            first = line.split()[-1].split('->')[0]
            if first.count("A") > 1 and "D" in first:
                if any(ext in atoms for ext in P2s):
                    if first.count("A") != 2:
                        line = ai.rjust(6) + aj.rjust(7) + ak.rjust(
                            7
                        ) + al.rjust(7) + funct.rjust(
                            5
                        ) + " " + a1 + " " + "0" + " " + f1 + " " + a2 + " " + fc2 + " " + f2 + " " + comment + "\n"
                        altered = True
                else:
                    line = ai.rjust(6) + aj.rjust(7) + ak.rjust(7) + al.rjust(
                        7
                    ) + funct.rjust(
                        5
                    ) + " " + a1 + " " + "0" + " " + f1 + " " + a2 + " " + fc2 + " " + f2 + " " + comment + "\n"
                    altered = True

            second = line.split()[-1].split('->')[1]
            if second.count("A") > 1 and "D" in second:
                if any(ext in comment for ext in P2s):
                    if second.count("A") != 2:
                        if altered:
                            raise AssertionError(
                                'Dihedral between two dummy atoms is made')
                        line = ai.rjust(6) + aj.rjust(7) + ak.rjust(
                            7
                        ) + al.rjust(7) + funct.rjust(
                            5
                        ) + " " + a1 + " " + fc1 + " " + f1 + " " + a2 + " " + "0" + " " + f2 + " " + comment + "\n"
                else:
                    if altered:
                        raise AssertionError(
                            'Dihedral between two dummy atoms is made')
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
            for line in in_f:
                if line.strip() == "[ moleculetype ]":
                    print("moleculetype section started!")
                    process_other(line.strip(), in_f, out_f)
                    print("moleculetype section done!")
                elif line.strip() == "[ atoms ]":
                    print("atoms section started!")
                    process_other(line.strip(), in_f, out_f)
                    print("atoms section done!")
                elif line.strip() == "[ bonds ]":
                    print("bonds sectrion started!")
                    process_other(line.strip(), in_f, out_f)
                    print("bonds sectrion done!")
                elif line.strip() == "[ pairs ]":
                    print("pairs section started!")
                    process_other(line.strip(), in_f, out_f)
                    print("pairs section done!")
                elif line.strip() == "[ angles ]":
                    print("angles section started!")
                    process_angles(line.strip(), in_f, out_f, fc=decouple_params['fc'], bridge_atoms = decouple_params['bridge_atoms'], P2s = decouple_params['P2s_angle'])
                    print("angles section done!")
                elif line.strip() == "[ dihedrals ]":
                    print("dihedrals section started!")
                    process_dhedrals(line.strip(), in_f, out_f, P2s = decouple_params['P2s_dihedral'])
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
               'name': 'to__int',
              'P2s_dihedral':['H8'],
              'bridge_atoms': {'A': [["C4", "DC10"]],
                               'B': [["C4", "C2"]]},
              'P2s_angle': [['C5', 'H8', 'H9']]}

    ref_int = {'fc':'20.92',
               'name': 'ref_int',
              'P2s_dihedral':['O3', 'C2'],
              'bridge_atoms': {'A': [["C1", "HV3"], ['C7', 'DN1']],
                               'B': [["C1", "O1"], ['C7', 'N1']]},
              'P2s_angle': [['C2', 'C6'], ['O3', 'C6']]}

    to__ref = {'fc':'20.92',
               'name': 'to__ref',
              'P2s_dihedral':['H8', 'O1', 'C11'],
              'bridge_atoms': {'A': [["C4", "DC14"], ["C10", "DO1"], ['C12', 'DN1']],
                               'B': [["C4", "C2"], ["C10", "H12"], ['C12', 'N3']]},
              'P2s_angle': [['C5', 'H8', 'H9'], ['C9', 'C11'], ['O1', 'C9']]}

    import json
    for params in [ref_int, to__int, to__ref]:
        with open(f'{params["name"]}.json', "w") as outfile:
            json.dump(params, outfile)
        process_file(f'merged_tmp_{params["name"]}.itp', f'merged_tmp_{params["name"]}_decoupled_new.itp', decouple_params = params)
