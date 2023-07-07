import random
import os


rev_alpha = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
rev_bin = {'00': '11', '01': '10', '10': '01', '11': '00'}
alpha2bin = {'A': '00', 'C': '01', 'G': '10', 'T': '11'}
bin2alpha = {'00': 'A', '01': 'C', '10': 'G', '11': 'T'}

def alpha2int(seq: str):
    return int(''.join([alpha2bin[n] for n in seq]), 2)

def int2alpha(num: int, k: int):
    bin_rep = bin(num)[2:].zfill(2*k)
    return ''.join([bin2alpha[bin_rep[2*i:2*(i+1)]] for i in range(k)])

def reverse_alpha(seq: str):
    return ''.join([rev_alpha[n] for n in reversed(seq)])

def reverse_int(num: int, k: int):
    bin_rep = bin(num)[2:].zfill(2*k)
    return int(''.join([rev_bin[bin_rep[2*(i-1):2*i]] for i in range(k, 0, -1)]), 2)

def is_alpha(c):
    return str(c) in rev_alpha

def canonical_alpha2int(seq: str):
    sub = min(seq, reverse_alpha(seq))
    return alpha2int(sub)

def hamming_distance(k1: str, k2: str):
    hamm = 0
    for (i,j) in zip(k1, k2):
        if i != j:
            hamm += 1
    return hamm

def min_pair(a, b):
    return (min(a,b), max(a,b))

def get_assignment_prob(bin_id: float):
    return random.choices(population=[0,1], weights=[1-bin_id, bin_id], k=1)[0]

def get_assignment_round(bin_id: float):
    if bin_id > 0.5:
        return 1
    elif bin_id < 0.5:
        return 0
    else:
        return random.randint(0, 1)

def get_assignment_rand():
    return random.choices(population=[0,1], weights=[0.5, 0.5], k=1)[0]

class Bin_t:
    def __init__(self) -> None:
        self.bin0 = []
        self.bin0len = 0
        self.bin1 = []
        self.bin1len = 0
        self.bin_amb = []
        self.binalen = 0
        return

    def insert_bin_det(self, vid: int, bin_id: int):
        if bin_id == 0:
            self.bin0.append(vid)
            self.bin0len += 1
        elif bin_id == 1:
            self.bin1.append(vid)
            self.bin1len += 1
        else:
            self.bin_amb.append(vid)
            self.binalen += 1
        return
    
    def insert_bin_prob(self, vid: int, bin_id: float):
        self.insert_bin_det(vid, get_assignment_prob(bin_id))
        return
    
    def insert_bin_round(self, vid: int, bin_id: float):
        self.insert_bin_det(vid, get_assignment_round(bin_id))
        return
    
    def insert_all_det(self, vids: list, bin_id: int):
        if bin_id == 0:
            self.bin0.extend(vids)
            self.bin0len += len(vids)
        else:
            self.bin1.extend(vids)
            self.bin1len += len(vids)
        return

    def get_fst_bin(self):
        return self.bin0

    def get_snd_bin(self):
        return self.bin1

    def get_amb_bin(self):
        return self.bin_amb

    def __str__(self) -> str:
        return "Bin 0: {0}\nBin 1: {1}\n".format(self.bin0, self.bin1)

def list_to_string(ids: list, s=""):
    string = s + " - " if s != "" else ""
    for id in ids:
        string += str(id) + ", "
    return string[:-2] if len(string) >= 2 else ""

def gfa_to_fasta(gfa_file, fasta_file):
    os.system(f"""awk '/^S/{{print ">"$2"\\n"$3}}' {gfa_file} > {fasta_file}""")
    return

def get_edges(fasta_file):
    unitigs = {}
    prev_id = ""
    with open(fasta_file, "r") as fd:
        for l in fd:
            if l.startswith(">"):
                prev_id = l[1:-1]
            else:
                unitigs[prev_id] = l[1:-1]
                
        fd.close()
    return unitigs
