import sys
import time

import seq_align

FILES = ["fastas\HelicoverpaArmigera-cMyc.fasta", "fastas\HomoSapiens-SHH.fasta",
         "fastas\LarimichthysCrocea-cMyc.fasta", "fastas\RatusNorvegicus-SHH.fasta"]

ALIGN_TYPE = ["global", "local", "overlap"]
NUMBER_OF_RUNS = 1

if __name__ == '__main__':
    s = sys.argv
    time_counter = []
    score_counter = []
    for align_type in ALIGN_TYPE:
        for f1 in FILES:
            for f2 in FILES:
                this_s = s.copy() + [f1, f2, "--align_type", align_type]
                sys.argv = this_s
                start = time.perf_counter()
                for i in range(NUMBER_OF_RUNS):
                    score = seq_align.main()
                end = time.perf_counter()
                time_counter.append(str((end - start) / NUMBER_OF_RUNS))
                score_counter.append(str(score))

    with open("timer.txt", "w") as timer:
        timer.write("\n".join(time_counter))

    with open("socre.txt", "w") as score_f:
        score_f.write("\n".join(score_counter))
