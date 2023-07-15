import os
import sys
from utils import *
from transformer import *
import time

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def parse_arg(argv: list):
    argc = len(argv)
    if argc < 3:
        print(f"{argv[0]} build -f<a|q> <k> <out_dir> <ref> <ref> ... <ref>")
        print(f"{argv[0]} query -f<a|q> <db_dir> <out_dir> <query>")
        print(f"{argv[0]} analysis <qry_dir>")
        sys.exit(1)
    _, mode = argv[:2]
    if mode == "build":
        is_fasta = str.endswith(argv[2], "a")
        if argc <= 5:
            print("no reference file(s) provided for build mode")
            sys.exit(1)
        files = argv[3:] # first `file` is k-mer size

        if files[1][-1] != "/":
            files[1] += "/"
        if not Create_Dir(files[1]):
            print("Please remove/choose selected/other directory.")
            sys.exit(0)

        return mode, is_fasta, files
    
    if mode == "query":
        is_fasta = str.endswith(argv[2], "a")
        if argc != 6:
            print("abort, either db_dir, out_dir, or query files missing for query mode")
            sys.exit(1)
        files = argv[3:]
        
        if files[0][-1] != "/":
            files[0] += "/"
        if not os.path.exists(files[0]):
            print("Database directory does not exist.")
            sys.exit(1)

        if files[1][-1] != "/":
            files[1] += "/"
        if not Create_Dir(files[1]):
            print("Please remove/choose selected/other output directory.")
            sys.exit(1)

        return mode, is_fasta, files

    if mode == "analysis":
        if argc != 3:
            print("incorrect parameter for analysis mode")
            sys.exit(1)
        files = argv[2:]
        if files[0][-1] != "/":
            files[0] += "/"
        if not os.path.exists(files[0]):
            print("Query directory does not exist.")
            sys.exit(1)
        return mode, False, files
    print("abort")
    sys.exit(1)


def build_phase(files: list, is_fasta : bool):
    """
    return: comps_idx.txt, comps_record.txt, comps_udb.txt
    """
    print("Build Phase")
    tstart = time.time()
    idx_arrs = []
    comp_records = []
    glb_kmer_table = {}

    [k, db_dir, ref_files] = [int(files[0]), files[1], files[2:]]
    print("k: ", k)

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

            for i in range(0, len(seq) - k):
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
    len_mat = [[0 for _ in range(len(row))] for row in idx_arrs]
    out_udb = db_dir + "comps_udb.txt"
    Create(out_udb)
    # TODO also store the expected # matching kmers to 
    # obtain a upper bound for matches, can be used on line 209
    with open(out_udb, "w") as udb_fd:
        udb_fd.write(f"{k}\n") #k
        for kmer, idy_arr in glb_kmer_table.items():
            # kmer that unique maps to single sub-component
            if len(idy_arr) == 1:
                [(comp_idx, s_idx)] = idy_arr
                len_mat[comp_idx][s_idx] += 1
                udb_fd.write(f"{kmer}\t{comp_idx}:{s_idx}\n")
        udb_fd.close()

    out_idx = db_dir + "comps_idx.txt"
    Create(out_idx)
    with open(out_idx, "w") as idx_fd:
        idx_fd.write(f"{len(idx_arrs)}\n") #num_components
        for comp_idx, idx_arr in enumerate(idx_arrs):
            for s_idx, sid in enumerate(idx_arr):
                idx_fd.write(f"{s_idx}\t{len_mat[comp_idx][s_idx]}\t{sid}\n")
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
    print(f"database init: {time.time() - tstart}s")
    print("Done")
    return

def query_phase(files: list, is_fasta : bool):
    print("Query Phase")
    tstart = time.time()
    [db_dir, out_dir, query_file] = files

    idx_arrs, len_mat = read_idx_file(db_dir + "comps_idx.txt")
    ukmer_table, k = read_db_file(db_dir + "comps_udb.txt")
    rnames = []
    with open(db_dir + "comps_rname.txt", "r") as rname_fd:
        rnames = rname_fd.readline().strip().split(";")
        rname_fd.close()
    reads = get_file(query_file, is_fasta)

    temp_rfile = out_dir + "temp.fasta"
    temp_pfile = out_dir + "temp.paf"
    temp_logfile = out_dir + "mmp.log"
    
    # RESULT
    total_assembled = [] # contain all parts, and follow the asc/dsc order
    partial_assembled = [] # follow the asc/dsc order, missing some parts
    mis_assembled = [] # missing all parts.

    # ANALYSIS
    count_p2t = 0 # number of partial to total conversion
    length_status = [None for _ in range(len(reads))]
    qal_mat = [None for _ in range(len(reads))]

    for ridx, (rid, rseq) in enumerate(reads):
        # if ridx == 100:
        #     break
        if ridx % 1001 == 0:
            print("processed reads: ", ridx)
            print(f"time elapsed: {time.time() - tstart}s")

        # process current read
        snode = LinkedList()
        len_rseq = len(rseq)
        for i in range(0, len_rseq - k):
            canon_kmer = canonical_alpha2int(rseq[i: i+int(k)])
            insert_key = "N"
            if canon_kmer in ukmer_table:
                # found an unique kmer from database, check the source
                kmer_identity = ukmer_table[canon_kmer]
                insert_key = f"{kmer_identity[0]},{kmer_identity[1]}"
            snode.append(Node(insert_key))
        iter_merge_del(snode, "N") # modify snode, remove N's and merge adjacent nodes.

        # at this point, all vertices with key-value pair as (`i:j`, val) is stored in the linkedlist, all
        # None identified kmers have been removed.
        arr = snode.to_arr_noN()
        size_arr = len(arr)
        # convert weight to normalized weight for all vertices
        for a_idx in range(size_arr):
            (ia, ja, wa) = arr[a_idx]
            arr[a_idx][2] = round(wa / len_mat[ia][ja], 3)

        lis_seq = []
        lds_seq = []

        # simplify step, compute the LIS and LDS in cyclic manner
        for start_idx in range(size_arr):
            lis, lds = cyc_lis_lds(arr, size_arr, start_idx)
            lis_seq = lis if len(lis) > len(lis_seq) else lis_seq
            lds_seq = lds if len(lds) > len(lds_seq) else lds_seq
        

        lxs_seq = None
        # select longest increasing/decreasing chain
        if len(lis_seq) >= len(lds_seq):
            lxs_seq = lis_seq
        else:
            lds_seq.reverse()
            lxs_seq = lds_seq


        if len(lxs_seq) == 0:
            # erroroness read
            mis_assembled.append(rid)
            length_status[ridx] = (len_rseq, "MIS")
            qal_mat[ridx] = [rid, ["Nil"]]
            continue

        num_components = len(idx_arrs)
        assembled_arr = [None for _ in range(num_components)]
        qal_arr = ["Nil" for _ in range(num_components)]
        for [j1,j2,kcount] in lxs_seq:
            assembled_arr[j1] = f"{idx_arrs[j1][j2]}"
            qal_arr[j1] = str(kcount)
        qal_mat[ridx] = [rid, qal_arr]
        
        if all([asm_elem != None for asm_elem in assembled_arr]):
            total_assembled.append([rid, assembled_arr])
            length_status[ridx] = (len_rseq, "TOTAL")
            continue

        for ri in range(num_components):
            if assembled_arr[ri] == None:
                Create(temp_rfile)
                System(f"echo \">{ridx}\n{rseq}\n\" > {temp_rfile}")
                System(f"minimap2 --secondary=no {temp_rfile} {rnames[ri]} > {temp_pfile} 2> {temp_logfile}")
                with open(temp_pfile, "r") as paf_fd:
                    paf_lines = paf_fd.readlines()
                    if len(paf_lines) == 1:
                        # perfect alignment
                        splited = paf_lines[0].strip().split("\t")
                        seg_no = str(splited[0])
                        mapq = splited[11]
                        # get full component name
                        for fname in idx_arrs[ri]:
                            if str.startswith(fname, seg_no):
                                assembled_arr[ri] = fname
                                break
                    else:
                        assembled_arr[ri] = "*"
                    paf_fd.close()

        if any([asm_str == "*" for asm_str in assembled_arr]):
            # belong to partial assembled read
            partial_assembled.append([rid, assembled_arr])
            length_status[ridx] = (len_rseq, "PARTIAL")
        else:
            total_assembled.append([rid, assembled_arr])
            length_status[ridx] = (len_rseq, "TOTAL")
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

    out_qry_total = out_dir + "query_total.csv"
    Create(out_qry_total)
    with open(out_qry_total, "w") as tot_fd:
        for [sid, arr] in total_assembled:
            tot_fd.write(f"{sid}," + ",".join(arr) + "\n")
        tot_fd.close()

    out_len_stat = out_dir + "query_len.csv"
    Create(out_len_stat)
    with open(out_len_stat, "w") as len_fd:
        for (rseq_len, rseq_status) in length_status:
            len_fd.write(f"{rseq_len},{rseq_status}\n")
        len_fd.close()

    out_qal_stat = out_dir + "query_qal.csv"
    Create(out_qal_stat)
    with open(out_qal_stat, "w") as qal_fd:
        for [sid, arr] in qal_mat:
            qal_fd.write(f"{sid}," + ",".join(arr) + "\n")
        qal_fd.close()


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
    print("Done")
    return


def analysis_phase(files: list, _):
    print("Analysis Phase")
    tstart = time.time()
    [qry_dir] = files


    out_qry_partial = qry_dir +"query_partial.csv"
    df_partial = pd.read_csv(out_qry_partial, header=None)
    counts = df_partial.groupby(list(df_partial.columns)[1:])[0].count() # group by combination type
    counts.to_csv(qry_dir + "query_partial.group.csv")

    out_qry_total = qry_dir +"query_total.csv"
    df_total = pd.read_csv(out_qry_total, header=None)
    counts = df_total.groupby(list(df_total.columns)[1:])[0].count() # group by combination type
    counts.to_csv(qry_dir + "query_total.group.csv")

    out_len_stat = qry_dir + "query_len.csv"
    df_len = pd.read_csv(out_len_stat, header=None)

    counts = df_len.groupby(list(df_len.columns)[1])
    df_group = pd.DataFrame({"MIS": counts.get_group("MIS")[0], 
                            "PARTIAL": counts.get_group("PARTIAL")[0], 
                            "TOTAL": counts.get_group("TOTAL")[0]})
    _, ax1 = plt.subplots()
    sns.histplot(data=df_group)

    ax1.set_xlabel("Read length (bp)")
    ax1.set_title("Histogram - Count vs Read Length")
    plt.savefig(qry_dir + "hist_len.png")
    print(f"time elapsed: {time.time() - tstart}s")
    print("Done")
    return




if __name__ == "__main__":
    mode, is_fasta, files = parse_arg(sys.argv)
    run = {"build": build_phase, "query": query_phase, "analysis": analysis_phase}
    run[mode](files, is_fasta)


