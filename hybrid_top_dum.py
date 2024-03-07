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

def process_angles(line, in_f, out_file):
    out_file.write(line + "\n")
    while True:
        line = in_f.readline()
        if line.strip() == "":
            break
        elif line.split()[0] == ";":
            out_file.write(line)
        else:
            ai, aj, ak, funct, c1, k1, c2, k2 = line.strip().split()[:8]
            comment = " ".join(line.split()[8:])
            if "C2" in comment and "C4" in comment:
                if any(ext in comment for ext in ['C5', 'H8', 'H9']):
                    line = ai.rjust(6) + aj.rjust(7) + ak.rjust(6) + funct.rjust(8) + c1.rjust(15) + "14.85".rjust(10) + c2.rjust(20) + "14.85".rjust(10) + comment + "\n"
            out_file.write(line)
        
    out_file.write("\n\n")
    return

def process_dhedrals(line, in_f, out_file):
    out_file.write(line + "\n")
    i = 0
    while True:
        line = in_f.readline()
        if line.strip() == "":
            out_file.write("\n")
            break
        if line.split()[0] == ";":
            out_file.write(line.strip() + "\n")
        else:
            ai, aj, ak, al, funct, a1, fc1, f1, a2, fc2, f2 = line.strip().split()[0:11]
            comment = " ".join(line.split()[11:])
            second = line.split()[-1].split('->')[1]
            if second.count("A") > 1 and "D" in second:
                if 'H8' not in comment:
                    line = ai.rjust(6) + aj.rjust(7) + ak.rjust(7) + al.rjust(7) + funct.rjust(5) + " " + a1 + " " + "0" + " " + f1 + " " + a2 + " " + "0" + " " + f2 + " " + comment + "\n"
            out_file.write(line)
    return


def process_file(in_file, out_file):
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
                    process_angles(line.strip(), in_f, out_f)
                    print("angles section done!")
                elif line.strip() == "[ dihedrals ]":
                    print("dihedrals section started!")
                    process_dhedrals(line.strip(), in_f, out_f)
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
    process_file()