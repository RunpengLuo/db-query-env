import os
import sys
from transformer import *

def System(command: str):
    return os.system(command)

def Create(file: str):
    return System("echo "" > " + file)

def Create_Dir(dir: str):
    if not os.path.exists(dir):
        # create new directory
        os.makedirs(dir)
        print("Directory: ", dir, "created.")
        return True
    else:
        print(f"Directory: {dir} exists")
        return False

def line_counter(file: str):
    c = 0
    with open(file, "r") as fd:
        for line in fd:
            if line.startswith(">") or line.startswith("@"):
                c += 1
        fd.close()
    return c

def get_file(file: str, is_fasta: bool):
    if is_fasta:
        return get_fasta(file)
    else:
        return get_fastq(file)

def get_fasta(fasta_file: str):
    res = []
    with open(fasta_file, "r") as fa_fd:
        sid = ""
        seq = ""
        for line in fa_fd:
            if line.startswith(">"):
                # process previous entry
                if sid != "":
                    res.append((sid, seq))
                sid = line.strip()[1:]
                seq = ""
            else:
                seq += line.strip()
        fa_fd.close()
        if sid != "":
            res.append((sid, seq))
    return res

def get_fastq(fastq_file: str):
    reads = []
    counter = 0
    with open(fastq_file, "r") as fd:
        prev_id = ""
        prev_seq = ""
        i = 0
        for l in fd:
            if i == 0:
                assert l.startswith("@")
                counter += 1

                if prev_id != "":
                    reads.append((prev_id, prev_seq))
                
                # read id, get first non-space id
                prev_id = l.strip()
                prev_seq = ""
            elif i == 1:
                prev_seq += l.strip()
            elif i == 2:
                if not l == "+\n":
                    # multiline sed
                    prev_seq += l.strip()
                    continue
            i = (i+1) % 4
        fd.close()
        if prev_id != "":
            reads.append((prev_id, prev_seq))

    print("total reads: ", counter)
    return reads

def reformat_fasta(in_file: str, out_file: str):
    os.system("echo "" > " + out_file)
    with open(in_file, "r") as in_fd:
        with open(out_file, "w") as out_fd:
            sid = ""
            seq = ""
            for line in in_fd:
                if line.startswith(">"):
                    # process previous entry
                    if sid != "":
                        out_fd.write(sid + "\n")
                        out_fd.write(seq + "\n")
                    sid = line.strip()
                    
                    seq = ""
                else:
                    seq += line.strip().upper()
            out_fd.close()
            if sid != "":
                out_fd.write(sid + "\n")
                out_fd.write(seq + "\n")
        in_fd.close()

def store_reads(reads: dict, only_id=True):
    # store basic read informations
    System("echo "" > reads.txt")
    with open("reads.txt", "w") as fd:
        for rid, rseq in reads.items():
            if only_id:
                fd.write(f"{rid}\n")
            else:
                fd.write(f"{rid}:{rseq}\n")
        fd.close()

def get_reads(only_id=True):
    ids = []
    with open("reads.txt", "r") as fd:
        for line in fd:
            if only_id:
                ids.append(line.strip())
            else:
                raise Exception
        fd.close()
    return ids

def read_idx_file(idx_file: str):
    idx_arrs = None
    # read idx_file
    with open(idx_file, "r") as idx_fd:
        num_comp = int(idx_fd.readline().strip())
        idx_arrs = [[] for _ in range(num_comp)]
        glb_idx = 0 # index between components
        for line in idx_fd:
            if line == "\n":
                glb_idx += 1
                continue
            
            [_, sid] = line.strip().split("\t")
            idx_arrs[glb_idx].append(sid)
        idx_fd.close()
    return idx_arrs

def read_db_file(db_file: str):
    glb_kmer_table = {}
    k = -1
    with open(db_file, "r") as db_fd:
        k = int(db_fd.readline().strip())
        for line in db_fd:
            attrs = line.strip().split("\t")
            kmer, idts = attrs[0], attrs[1:]
            kint = alpha2int(kmer)
            assert kint not in glb_kmer_table
            glb_kmer_table[kint] = []
            for idt in idts:
                cidx,sidx = [int(x) for x in idt.split(":")]
                glb_kmer_table[kint].append((cidx, sidx))
        db_fd.close()
    return glb_kmer_table, k


class Node:
    def __init__(self, key="") -> None:
        self.key = key
        self.count = 1
        self.next = None
        self.prev = None

        self.next_none = 0

   
class LinkedList:
    def __init__(self) -> None:
        self.start = Node()
        self.curr = self.start
        self.count = 0
    
    def append(self, n: Node):
        n.prev = self.curr
        self.curr.next = n
        self.count += 1
        self.curr = n
        return True
    
    def source(self):
        return self.start

    def __str__(self) -> str:
        arr = []
        node = self.start.next # exclude root node
        while node != None:
            arr.append(f"{node.key}({node.count})")
            node = node.next
        
        return "->".join(arr)

    def arr_filter(self, fil_key="") -> str:
        arr = []
        node = self.start.next # exclude root node
        while node != None:
            if node.key != fil_key:
                arr.append((node.key, node.count))
            node = node.next
        
        return arr



def merge(dst: Node, src: Node, add_count):
    """
    assume dst.next == src, merge dst and src, preserve curr key
    """
    dst.next = src.next
    if src.next != None:
        src.next.prev = dst
    if add_count:
        dst.count += src.count
    else:
        # next node is del key node
        dst.next_none += src.count
    return dst

    
def iter_merge(ll: LinkedList):
    """
    iterative merge the adjacent node if having same key
    """
    curr = ll.source()

    while curr.next != None:
        next = curr.next
        if curr.key == next.key:
            # print("merging: ", curr.key, next.key, curr.count)
            curr = merge(curr, next, True)
        else:
            curr = next
    # reaching last node, check if first node and last node can merge, soft merge, no link inheritance
    if ll.source().next != None:
        if ll.source().next.key == curr.key:
            # add count to first node
            ll.source().next.count += curr.count
            # remove last node
            curr.prev.next = None
    return ll

def iter_merge_del(ll: LinkedList, del_key: str):
    """
    iterative merge the adjacent node if having same key, or node with del key
    """
    curr = ll.source()

    while curr.next != None:
        next = curr.next
        if curr.key == next.key or next.key == del_key:
            # print("merging: ", curr.key, next.key, curr.count)
            curr = merge(curr, next, next.key != del_key)
        else:
            curr = next
    # reaching last node, check if first node and last node can merge, soft merge, no link inheritance
    if ll.source().next != None:
        if ll.source().next.key == curr.key:
            # add count to first node
            ll.source().next.count += curr.count
            # remove last node
            curr.prev.next = None
    return ll


def iter_del(ll: LinkedList, del_key: str):
    """
    iterative delete the node with `del_key`
    """
    curr = ll.source()
    while curr.next != None:
        next = curr.next
        if next.key == del_key:
            curr = merge(curr, next, False)
        else:
            curr = next
    return ll

def cyc_lis_lds(arr: list, size_arr: int, start_idx: int):
    """
    Dynamic programming to compute LIS and LDS.
    """
    # circular array starting from arr[start_idx]
    dp_lis = [None for _ in range(size_arr)]
    dp_lis[0] = [arr[start_idx]]
    dp_lds = [None for _ in range(size_arr)]
    dp_lds[0] = [arr[start_idx]]
    # arr[start_idx], arr[start_idx + 1], ... arr[(start_idx + size_arr - 1) % size_arr]
    rot_arr = [arr[start_idx]]
    for acc in range(1, size_arr):
        curr_elem = arr[(start_idx + acc) % size_arr]
        rot_arr.append(curr_elem)
        lis_arr = []
        lds_arr = []
        for j in range(0, acc):
            if len(dp_lis[j]) > len(lis_arr):
                if curr_elem > dp_lis[j][-1]:
                    lis_arr = dp_lis[j]
            if len(dp_lds[j]) > len(lds_arr):
                if curr_elem < dp_lds[j][-1]:
                    lds_arr = dp_lds[j]

        dp_lis[acc] = list(lis_arr) # deep copy
        dp_lis[acc].append(curr_elem)
        dp_lds[acc] = list(lds_arr)
        dp_lds[acc].append(curr_elem)
    
    return dp_lis[size_arr - 1], dp_lds[size_arr - 1]
    