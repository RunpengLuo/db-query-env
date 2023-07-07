import sys
import os

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print(f"{sys.argv[0]} <fasta> <out>")
        sys.exit(1)
    _, fa, out = sys.argv
    os.system("echo "" > " + out)
    with open(fa, "r") as fa_fd:
        with open(out, "w") as out_fd:
            sid = ""
            seq = ""
            for line in fa_fd:
                if line.startswith(">"):
                    # process previous entry
                    if sid != "":
                        out_fd.write(sid + "\n")
                        out_fd.write(seq + "\n")
                    sid = line.strip()
                    
                    seq = ""
                else:
                    seq += line.strip().upper()
            if sid != "":
                out_fd.write(sid + "\n")
                out_fd.write(seq + "\n")
            out_fd.close()
        fa_fd.close()
