import os
import sys
from utils import *
from transformer import *
import time

import pandas as pd


def parse_arg(argv: list):
    argc = len(argv)
    if argc < 5:
        print(f"{argv[0]} build -f<a|q> <k> <out_dir> <ref> <ref> ... <ref>")
        print(f"{argv[0]} query -f<a|q> <db_dir> <out_dir> <query>")
        sys.exit(1)
    _, mode = argv[:2]
    if mode == "build":
        is_fasta = str.endswith(argv[2], "a")
        files = argv[3:] # first file is k-mer size
        return mode, is_fasta, files
    elif mode == "query":
        is_fasta = str.endswith(argv[2], "a")
        files = argv[3:]
        if argc != 6:
            print("abort, either db_dir, out_dir, or query files missing")
            sys.exit(1)
        return mode, is_fasta, files
    else:
        print("abort")
        sys.exit(1)


def build_phase(files: list, is_fasta : bool):
    """
    return: comps_idx.txt, comps_record.txt, comps_db.txt
    """
    print("Build Phase")
    tstart = time.time()
    idx_arrs = []
    comp_records = []
    glb_kmer_table = {}

    [k, db_dir, ref_files] = [int(files[0]), files[1], files[2:]]
    print("k: ", k)

    if db_dir[-1] != "/":
        db_dir += "/"
    if not Create_Dir(db_dir):
        print("Please remove/choose selected/other directory.")
        sys.exit(0)

    out_ref = db_dir + "comps_rname.txt"
    Create(out_ref)
    with open(out_ref, "w") as ref_fd:
        ref_fd.write(";".join(ref_files) + "\n")
        ref_fd.close()


    for ref_file in ref_files:
        print("processing component: " + ref_file)
        comp_idx = len(idx_arrs)
        # 1. get unique kmers within the component
        data = get_file(ref_file, is_fasta)
        kmer_table = {}
        idx_arr = []

        for s_idx, (sid, seq) in enumerate(data):
            idx_arr.append(sid)

            l = len(seq)
            # if idx % 500 == 0:
            #     print(idx)
            for i in range(0, l - k):
                canon_kmer = canonical_alpha2int(seq[i: i+int(k)])
                if canon_kmer not in kmer_table:
                    kmer_table[canon_kmer] = []
                kmer_table[canon_kmer].append(s_idx)
                
        comp_record = [[] for _ in range(len(idx_arr))]
        for kmer, s_idxs in kmer_table.items():
            if len(s_idxs) == 1:
                [s_idx] = s_idxs
                # unique in current column
                if kmer not in glb_kmer_table:
                    glb_kmer_table[kmer] = []
                glb_kmer_table[kmer].append((comp_idx, s_idx))
                comp_record[s_idx].append(kmer)
        
        idx_arrs.append(idx_arr)
        comp_records.append(comp_record)

    # store the results
    out_idx = db_dir + "comps_idx.txt"
    Create(out_idx)
    with open(out_idx, "w") as idx_fd:
        idx_fd.write(f"{len(idx_arrs)}\n") #num_components
        for idx_arr in idx_arrs:
            for idx, sid in enumerate(idx_arr):
                idx_fd.write(f"{idx}\t{sid}\n")
            idx_fd.write("\n")
        idx_fd.close()
    
    out_record = db_dir + "comps_record.txt"
    Create(out_record)
    with open(out_record, "w") as rmap_fd:
        rmap_fd.write(f"{k}\n") #k
        rmap_fd.write(f"{len(comp_records)}\n") #num_components
        for comp_record in comp_records:
            for idx, kmers in enumerate(comp_record):
                kmer_str = "\t".join([int2alpha(s, k) for s in kmers])
                rmap_fd.write(f"{idx}\t{kmer_str}\n")
            rmap_fd.write("\n")
        rmap_fd.close()
    
    out_db = db_dir + "comps_db.txt"
    Create(out_db)
    with open(out_db, "w") as db_fd:
        db_fd.write(f"{k}\n") #k
        for kmer, arr in glb_kmer_table.items():
            ostr = int2alpha(kmer, k)
            for (comp_idx, s_idx) in arr:
                ostr += f"\t{comp_idx}:{s_idx}"
            ostr += "\n"
            db_fd.write(ostr)
        db_fd.close()
    print(f"database init: {time.time() - tstart}s")
    return

def query_phase(files: list, is_fasta : bool):
    print("Query Phase")
    tstart = time.time()
    [db_dir, out_dir, query_file] = files

    if db_dir[-1] != "/":
        db_dir += "/"
    if not os.path.exists(db_dir):
        print("Database directory does not appear: ", db_dir)
        sys.exit(1)

    if out_dir[-1] != "/":
        out_dir += "/"
    if not Create_Dir(out_dir):
        print("Please remove/choose selected/other directory.")
        sys.exit(1)

    idx_file = db_dir + "comps_idx.txt"
    idx_arrs = read_idx_file(idx_file)
    db_file = db_dir + "comps_db.txt"
    glb_kmer_table, k = read_db_file(db_file)
    rname_file = db_dir + "comps_rname.txt"
    rnames = []
    with open(rname_file, "r") as rname_fd:
        rnames = rname_fd.readline().strip().split(";")
        rname_fd.close()
    reads = get_file(query_file, is_fasta)

    temp_rfile = out_dir + "temp.fasta"
    temp_pfile = out_dir + "temp.paf"
    temp_logfile = out_dir + "mmp.log"
    
    total_assembled = [] # contain all parts, and follow the asc/dsc order
    partial_assembled = [] # follow the asc/dsc order, missing some parts
    mis_assembled = [] # missing all parts.
    count_p2t = 0 # number of partial to total conversion

    for ridx, (rid, rseq) in enumerate(reads):
        if ridx % 1001 == 0:
            print("processed reads: ", ridx)
            print(f"time elapsed: {time.time() - tstart}s")

        # process current read
        snode = LinkedList()
        for i in range(0, len(rseq) - k):
            canon_kmer = canonical_alpha2int(rseq[i: i+int(k)])
            insert_key = "N"
            if canon_kmer in glb_kmer_table:
                # found an unique kmer from database, check the source
                if len(glb_kmer_table[canon_kmer]) == 1:
                    kmer_identity = glb_kmer_table[canon_kmer][0]
                    insert_key = f"{kmer_identity[0]},{kmer_identity[1]}"
            snode.append(Node(insert_key))
        iter_merge_del(snode, "N") # modify snode, remove N's and merge adjacent nodes.
        arr = snode.arr_filter()
        size_arr = len(arr)
        lis_seq = []
        lds_seq = []

        # simplify step, compute the LIS and LDS in cyclic manner
        for start_idx in range(size_arr):
            lis, lds = cyc_lis_lds(arr, size_arr, start_idx)
            lis_seq = lis if len(lis) > len(lis_seq) else lis_seq
            lds_seq = lds if len(lds) > len(lds_seq) else lds_seq
        
        if len(lis_seq) == 0 and len(lds_seq) == 0:
            # erroroness read
            mis_assembled.append(rid)
        else:
            lxs_seq = None
            # select longest increasing/decreasing chain
            if len(lis_seq) >= len(lds_seq):
                lxs_seq = lis_seq
            else:
                lds_seq.reverse()
                lxs_seq = lds_seq
            
            num_components = len(idx_arrs)
            assembled_arr = [None for _ in range(num_components)]
            for tp in lxs_seq:
                [j1, j2] = tp[0].split(",")
                assembled_arr[int(j1)] = idx_arrs[int(j1)][int(j2)]
            
            if all([asm_elem != None for asm_elem in assembled_arr]):
                total_assembled.append([rid, assembled_arr])
            else:
                for ri in range(num_components):
                    if assembled_arr[ri] == None:
                        Create(temp_rfile)
                        System(f"echo \">{ridx}\n{rseq}\n\" > {temp_rfile}")
                        System(f"minimap2 --secondary=no {temp_rfile} {rnames[ri]} > {temp_pfile} 2> {temp_logfile}")
                        with open(temp_pfile, "r") as paf_fd:
                            paf_lines = paf_fd.readlines()
                            if len(paf_lines) != 1:
                                # empty file or over alignment
                                assembled_arr[ri] = "*"
                            else:
                                splited = paf_lines[0].strip().split("\t")
                                seg_no = str(splited[0])
                                # ref_no = str(splited[5])
                                assembled_arr[ri] = seg_no
                            paf_fd.close()

                if any([asm_str == "*" for asm_str in assembled_arr]):
                    # belong to partial assembled read
                    partial_assembled.append([rid, assembled_arr])
                else:
                    total_assembled.append([rid, assembled_arr])
                    count_p2t += 1

    System(f"rm {temp_logfile}; rm {temp_pfile}; rm {temp_rfile}")
    out_qry_mis = out_dir + "query_mis.csv"
    Create(out_qry_mis)
    with open(out_qry_mis, "w") as mis_fd:
        for sid in mis_assembled:
            mis_fd.write(f"{sid}\n")
        mis_fd.close()

    out_qry_partial = out_dir +"query_partial.csv"
    Create(out_qry_partial)
    with open(out_qry_partial, "w") as par_fd:
        for [sid, arr] in partial_assembled:
            par_fd.write(f"{sid}," + ",".join(arr) + "\n")
        par_fd.close()
    
    df_partial = pd.read_csv(out_qry_partial, header=None)
    counts = df_partial.groupby(list(df_partial.columns)[1:])[0].count() # group by combination type
    counts.to_csv(out_dir + "query_partial.group.csv")

    out_qry_total = out_dir + "query_total.csv"
    Create(out_qry_total)
    with open(out_qry_total, "w") as tot_fd:
        for [sid, arr] in total_assembled:
            tot_fd.write(f"{sid}," + ",".join(arr) + "\n")
        tot_fd.close()
    
    df_total = pd.read_csv(out_qry_total, header=None)
    counts = df_total.groupby(list(df_total.columns)[1:])[0].count() # group by combination type
    counts.to_csv(out_dir + "query_total.group.csv")

    summary = out_dir + "summary.txt"
    Create(summary)
    with open(summary, "w") as sum_fd:
        sum_fd.write(f"total queried reads: {len(reads)}\n")
        sum_fd.write(f"misassembled: {len(mis_assembled)}\n")
        sum_fd.write(f"partial assembled: {len(partial_assembled)}\n")
        sum_fd.write(f"correctly assembled: {len(total_assembled)}\n")
        sum_fd.write(f"partial to correctly assembled: {count_p2t}\n")
        sum_fd.write(f"time elapsed: {time.time() - tstart}s")
        sum_fd.close()

    
    return


if __name__ == "__main__":
    mode, is_fasta, files = parse_arg(sys.argv)
    run = {"build": build_phase, "query": query_phase}
    run[mode](files, is_fasta)


