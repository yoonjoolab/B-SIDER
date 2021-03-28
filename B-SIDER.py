#!/usr/bin/env python3
# coding: utf-8

import os, sys, sqlite3
import csv, random
import math, statistics
import argparse



class B_SIDER:
    def __init__(self, db="comp_seq_DB.db"):
        self.verbose = True
        if os.path.isfile(db):
            self.vprint("Fragment database: %s\n"%db)
        else:
            self.vprint("The database file %s does not exist"%db)
            sys.exit()
        self.conn = sqlite3.connect(db)
        self.conn.row_factory = sqlite3.Row
        self.cur = self.conn.cursor()
        self.aa_list = list("ACDEFGHIKLMNPQRSTVWY")
        self.scoring_matrix = []
        self.penalty = 100.0 # penalty when no certain amino acid is present

    def comp_seq_search(self, target, parallel, min_frag=3, output=None):
        if parallel:
            self.para_status = "Parallel"
        else:
            self.para_status = "Anti-parallel"
        self.vprint("Searching for the database...\n")
        self.vprint("Target sequence: %s"%target)
        self.vprint("parallelity: %s"%self.para_status)
        self.vprint("Minimum fragment size: %d\n"%min_frag)
        self.target = target
        self.min_frag = min_frag
        window = len(target)
        target_len = len(target)
        self.comp_seq_list = []
        while window >= self.min_frag :
            if window >= self.min_frag:
                for k in range(target_len-window+1):
                    count = 0
                    frag_target = target[k:k+window]
                    front_gap = k
                    back_gap = target_len - len(frag_target) - front_gap
                    new_seq = "-"*front_gap + frag_target + "-"*back_gap
                    self.vprint("Searching fragments: %s"%new_seq, end=" ")
                    self.cur.execute("SELECT * FROM pdb WHERE INSTR(seq1,'%s')>0 AND parallel=%d"%(frag_target, parallel))
                    for l in self.cur.fetchall():
                        pos = l['seq1'].find(frag_target)
                        comp_seq = "-"*front_gap + l['seq2'][pos:pos+window] + "-"*back_gap
                        self.comp_seq_list.append(comp_seq)
                        count += 1
                    self.vprint("(%d matches found)"%count)
            else:
                pass
            window = window - 1
            if window < min_frag :
                break

        self.vprint("\nComplementary sequences:\n")
        for seq in self.comp_seq_list:
            self.vprint(seq)

        if output is not None:
            self.vprint("Complementary sequences are saved in %s"%output)
            output_f = open(output, "w")
            for seq in self.comp_seq_list:
                output_f.write("%s\n"%seq)
            output_f.close()

    def build_score_matrix(self, output=None, background_frequency=None):
        self.vprint("\nBuilding a scoring matrix for %s\n"%self.target)
        if background_frequency is None:
            self.background = {'A': 0.0628, 'C': 0.0193, 'E': 0.0537, 'D': 0.0709, 
                    'G': 0.0513, 'F': 0.0315, 'I': 0.0356, 'H': 0.0224, 
                    'K': 0.0577, 'M': 0.0155, 'L': 0.0637, 'N': 0.0523, 
                    'Q': 0.0343, 'P': 0.0784, 'S': 0.0694, 'R': 0.0427, 
                    'T': 0.0642, 'W': 0.0106, 'V': 0.0499, 'Y': 0.03}
        else:
            if os.path.isfile(background_frequency):
                self.vprint("Background amino acid frequency information: %s"%background_frequency)
            else:
                self.error_print("The background amino acid frequency file (%s) does not exist"%background_frequency)
            background = list(csv.DictReader(open(background_frequency)))
            self.background = dict(map(lambda x: [x["AA"], float(x["freq"])], background))

        self.aa_freq_from_db = []
        self.best_score = [self.penalty] * len(self.target)
        self.best_aa = ["A"] * len(self.target)
        for p in range(len(self.target)):
            pos_aa_freq = dict(map(lambda x: [x, 0], self.aa_list))
            self.aa_freq_from_db.append(pos_aa_freq)
        for seq in self.comp_seq_list:
            for p in range(len(seq)):
                if seq[p] != "-":
                    self.aa_freq_from_db[p][seq[p]] += 1
        for p in range(len(self.target)):
            self.scoring_matrix.append({})
            pos_sum = sum(self.aa_freq_from_db[p].values())
            for aa in self.aa_list:
                if self.aa_freq_from_db[p][aa] == 0:
                    sc = self.penalty
                else:
                    sc = -1*math.log(self.aa_freq_from_db[p][aa] / (pos_sum * self.background[aa]))
                self.scoring_matrix[p][aa] = sc
                if self.best_score[p] > sc:
                    self.best_score[p] = sc
                    self.best_aa[p] = aa
        score_output = [["Pos", "Target"] + self.aa_list]
        score_screen = score_output[:]
        for p in range(len(self.target)):
            score_output.append([p+1, self.target[p]] + list(map(lambda x: self.scoring_matrix[p][x], self.aa_list)))
            score_screen.append([p+1, self.target[p]] + list(map(lambda x: "%.2f"%self.scoring_matrix[p][x], self.aa_list)))
        self.vprint("\nScoring matrix for the target:")
        for l in list(zip(*score_screen)):
            self.vprint("\t".join(map(str, l)))
        if output is not None:
            output_f = open(output, "w")
            score_mtx = csv.writer(output_f)
            #for p in range(len(self.target)):
            for l in list(zip(*score_output)):
                score_mtx.writerow(l)
            output_f.close()
            self.vprint("The scoring matrix was saved in %s"%output)
        print("The best complementarity peptide for %s: %s"%(self.target, "".join(self.best_aa)))
        print("Parallelity: %s"%(self.para_status))
        print("Complementarity score: %f"%(sum(self.best_score)))

    def score_calculator(self, peptide):
        if len(peptide) != len(self.target):
            self.error_print("The sequence lengths between %s (%d) and %s (%d) are difference"%(self.target, len(self.target), peptide, len(peptide)))
        else:
            self.score = 0.0
            for p in range(len(peptide)):
                self.score += self.scoring_matrix[p][peptide[p]]

    def compare_complementarity(self, comp_seq, randnum=10000):
        self.vprint("\nComparing complementarity against %s for %s using %d random sequences\n"%(self.target, comp_seq, randnum))
        self.score_calculator(comp_seq)
        comp_seq_score = self.score
        random_score_list = []
        self.vprint("Random complementarity sequences and scores:\n")
        for n in range(randnum):
            random_peptide = [random.choice(self.aa_list) for i in range(len(comp_seq))]
            random_peptide = "".join(random_peptide)
            self.score_calculator(random_peptide)
            random_score_list.append(self.score)
            self.vprint("%s %f"%(random_peptide, self.score))
        random_score_list.sort()
        for n, s in enumerate(random_score_list):
            if comp_seq_score < s:
                rank = float(n)
                break
        mean = statistics.mean(random_score_list)
        stdev = statistics.stdev(random_score_list)
        z_score = (comp_seq_score - mean) / stdev
        self.vprint("\nSequence: %s"%comp_seq)
        self.vprint("Score: %s"%comp_seq_score)
        self.vprint("Mean (standard deviation) of the random peptides: %f (%f)"%(mean, stdev))
        self.vprint("Rank: %.3f %%"%(n/randnum*100))

    def vprint(self, expr, end="\n"):
        if self.verbose:
            print(expr, end=end)

    def error_print(self, expr):
        sys.stderr.write(expr + "\n")
        sys.exit()

if __name__ == "__main__":
    default_DB = os.path.join(os.getcwd(), "database", "comp_seq_DB.db")
    parser = argparse.ArgumentParser()
    parser.add_argument("-v", "--verbose", dest="verbose", 
            help="Output verbosity (1 for True and 0 for False). Default: True",
            action="store", default=True, type=int)
    parser.add_argument("-t", "--target", dest="target", 
            help="Target sequence",
            action="store", default=None)
    parser.add_argument("-p", "--parallel", dest="parallel", 
            help="Parallelity (1 for parallel sheet and 0 for anti-parallel). Default: False",
            action="store", default=False, type=int)
    parser.add_argument("-d", "--database", dest="db", 
            help="Database (Default: %s)"%default_DB,
            action="store", default=default_DB)
    parser.add_argument("-b", "--background_freq", dest="background_freq", 
            help="Background frequency file (csv). If not given, default values are used",
            action="store", default=None)
    parser.add_argument("-s", "--score_output", dest="score_output", 
            help="Output for the position specific score matrix in CSV (Default: None)",
            action="store", default=None)
    parser.add_argument("-o", "--seq_output", dest="seq_output", 
            help="List of sequences used for the position specific score matrix in text (Default: None)",
            action="store", default=None)
    parser.add_argument("-m", "--min_frag", dest="min_frag", 
            help="The shortest number of residues for search (Default: 3)",
            action="store", default=3, type=int)
    parser.add_argument("-c", "--compare", dest="compare", 
            help="Sequence for comparison. If present, the complementarity of this sequence for the target is compared against random sequences",
            action="store", default=None)
    parser.add_argument("-n", "--randnum", dest="randnum", 
            help="Number of random sequences when a comparing sequence is given (Default: 10,000)",
            action="store", default=10000, type=int)

    args = parser.parse_args()
    if args.target is None:
        sys.stderr.write("A target sequence must be given\n")
        parser.print_help()
        sys.exit()
    verbose = bool(args.verbose)

    a = B_SIDER(db=args.db)
    a.verbose = verbose
    a.comp_seq_search(target=args.target, parallel=args.parallel, min_frag=args.min_frag, output=args.seq_output)
    a.build_score_matrix(output=args.score_output, background_frequency=args.background_freq)

    if args.compare is not None:
        a.compare_complementarity(comp_seq=args.compare, randnum=args.randnum)
