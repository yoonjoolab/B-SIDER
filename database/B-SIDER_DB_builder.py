from pymol import cmd,stored
import sqlite3

one_letter ={'VAL':'V', 'ILE':'I', 'LEU':'L', 'GLU':'E', 'GLN':'Q', 
        'ASP':'D', 'ASN':'N', 'HIS':'H', 'HID':'H', 'HIE':'H', 
        'TRP':'W', 'PHE':'F', 'TYR':'Y', 'ARG':'R', 'LYS':'K', 
        'SER':'S', 'THR':'T', 'MET':'M', 'ALA':'A', 'GLY':'G', 
        'PRO':'P', 'CYS':'C', 'CYX': 'C', 'CYD': 'C'}

def int_res(res):
    try:
        res = int(res)
    except:
        res = int(res[:-1])
    return res

def grouping(residue_list):
    min_ss = 3
    secondary_struct_group = {}
    strand_index = 0
    for chain, res in residue_list:
        if chain not in secondary_struct_group:
            secondary_struct_group[chain] = [[res]]
            strand_index = 0
        else:
            last_res = secondary_struct_group[chain][strand_index][-1]
            res_num_diff = int_res(res) - int_res(last_res)
            if res_num_diff >= 0 and res_num_diff < 2:
                secondary_struct_group[chain][strand_index].append(res)
            else:
                secondary_struct_group[chain].append([res])
                strand_index += 1
    filtered_ss_group = {}
    count = 0
    for chain in secondary_struct_group:
        for res in secondary_struct_group[chain]:
            if len(res) >= min_ss:
                name = "ss%d"%count
                cmd.select(name, "chain %s and resi %s"%(chain, "+".join(res)))
                filtered_ss_group[name] = map(lambda x: [chain, x], res)
                count += 1
    return filtered_ss_group

def pair_match(group):
    ss_name = group.keys()
    pair_list = []
    for i in range(len(ss_name)-1):
        ss1 = ss_name[i]
        for j in range(i+1, len(ss_name)):
            ss2 = ss_name[j]
            d = cmd.distance("%s-%s"%(ss1, ss2), ss1, ss2, mode=2)
            if d == 0.0:
                cmd.delete("%s-%s"%(ss1, ss2))
            else:
                pair_list.append([ss1, ss2])
                #cmd.delete("%s-%s"%(ss1, ss2))
    return pair_list

def parallel_check(group, ss1, ss2):
    hydro_b_cut = 3.5
    parallel = True
    seq1 = False; seq2 = False
    good_pair = False
    res_list1 = group[ss1]
    res_list2 = group[ss2]
    ss1N = "chain %s and resi %s and name ca"%(res_list1[0][0], res_list1[0][1])
    ss2N = "chain %s and resi %s and name ca"%(res_list2[0][0], res_list2[0][1])
    ss2C = "chain %s and resi %s and name ca"%(res_list2[-1][0], res_list2[-1][1])
    dist_ss1N_ss2N = cmd.distance("dist", ss1N, ss2N)
    dist_ss1N_ss2C = cmd.distance("dist", ss1N, ss2C)
    cmd.delete("dist")
    res_pairs = []

    if dist_ss1N_ss2N > dist_ss1N_ss2C:
        parallel = False
    if parallel:
        temp_res_pairs = []
        start_o = False
        start_n = False
        h_bond = True
        new_o_n = None
        for n, res in enumerate(res_list1):
            stored.o_n_pairs = []
            stored.n_o_pairs = []
            atom_neighbor = "(chain %s and resi %s and name o) around %f and (%s and name n)"%(res[0], res[1], hydro_b_cut, ss2)
            cmd.select("atom_neighbor", "%s"%(atom_neighbor))
            cmd.iterate("atom_neighbor","stored.o_n_pairs.append([chain, resi])")
            atom_neighbor = "(chain %s and resi %s and name n) around %f and (%s and name o)"%(res[0], res[1], hydro_b_cut, ss2)
            cmd.select("atom_neighbor", "%s"%(atom_neighbor))
            cmd.iterate("atom_neighbor","stored.n_o_pairs.append([chain, resi])")
            pos_res = res_list1.index(res)
            if len(stored.o_n_pairs) == 1 and len(stored.n_o_pairs) == 0:
                pos_o_n = res_list2.index(stored.o_n_pairs[0])
                if pos_o_n-1 >= 0:
                    temp_res_pairs.append([res, res_list2[pos_o_n-1]])
                    new_o_n = stored.o_n_pairs[0]
                    start_o = True
            elif len(stored.o_n_pairs) == 0 and len(stored.n_o_pairs) == 1:
                pos_n_o = res_list2.index(stored.n_o_pairs[0])
                if pos_n_o+1 < len(res_list2) and pos_n_o+1 > 0:
                    temp_res_pairs.append([res, res_list2[pos_n_o+1]])
            elif len(stored.o_n_pairs) == 0 and len(stored.n_o_pairs) == 0:
                if start_o and res_list2.index(new_o_n) >= n:
                    temp_res_pairs.append([res, new_o_n])
                if start_n:
                    pass
            elif len(stored.o_n_pairs) == 1 and len(stored.n_o_pairs) == 1:
                pos_o_n = res_list2.index(stored.o_n_pairs[0])
                pos_n_o = res_list2.index(stored.n_o_pairs[0])
                if pos_res - 1 >= 0:
                    temp_res_pairs.append([res_list1[pos_res-1], stored.n_o_pairs[0]])
                if pos_n_o+1 < len(res_list2) and pos_n_o+1 > 0:
                    temp_res_pairs.append([res, res_list2[pos_n_o+1]])
                if pos_res+1 < len(res_list1):
                    temp_res_pairs.append([res_list1[pos_res+1], stored.o_n_pairs[0]])
                start_n = True
        for r in temp_res_pairs:
            if r not in res_pairs:
                res_pairs.append(r)
    else:
        start = False
        h_bond = True
        temp_res_pairs = []
        for res in res_list1:
            stored.neighbor = []
            n_neighbor = "(chain %s and resi %s and name n) around %f and (%s and name o)"%(res[0], res[1], hydro_b_cut, ss2)
            cmd.select("n_neighbor", "%s"%(n_neighbor))
            cmd.iterate("n_neighbor","stored.neighbor.append([chain, resi])")
            bonded_res = False
            if len(stored.neighbor) > 0:
                for stored_res in stored.neighbor:
                    o_n_dist = cmd.distance("o_neighbor", "chain %s and resi %s and name o"%(res[0], res[1]), 
                            "chain %s and resi %s and name n"%(stored.neighbor[0][0], stored.neighbor[0][1]))
                    cmd.delete("o_neighbor")
                    if o_n_dist < hydro_b_cut:
                        bonded_res = stored_res
                        start = True
                        break
            if start and bool(bonded_res) == h_bond:
                if bonded_res:
                    temp_res_pairs.append([res, bonded_res])
                else:
                    if res_list2.index(temp_res_pairs[-1][1]) - 1 >= 0:
                        temp_res_pairs.append([res, res_list2[res_list2.index(temp_res_pairs[-1][1])-1]])
                h_bond = not(bool(bonded_res))
        if len(temp_res_pairs) >= 3:
            if not temp_res_pairs[-1][1]:
                temp_res_pairs.pop()
            for r in temp_res_pairs:
                if r not in res_pairs:
                    res_pairs.append(r)
    if len(res_pairs) >= 3:
        print("Parallel?:", parallel, ss1, ss2)
        stored.aa1 = []
        stored.aa2 = []
        for i in res_pairs:
            cmd.select("aa1", "chain %s and resi %s and name ca"%(i[0][0], i[0][1]))
            cmd.iterate("aa1", "stored.aa1.append(resn)")
            cmd.select("aa2", "chain %s and resi %s and name ca"%(i[1][0], i[1][1]))
            cmd.iterate("aa2", "stored.aa2.append(resn)")
            print(i, [stored.aa1[-1], stored.aa2[-1]])
        if not parallel:
            stored.aa2.reverse()
        try:
            aa1 = list(map(lambda x: one_letter[x], stored.aa1))
            aa2 = list(map(lambda x: one_letter[x], stored.aa2))
        except:
            aa1 = []
            aa2 = []
        if len(aa1) == len(stored.aa1) and len(aa1) == len(aa2):
            seq1 = "".join(aa1)
            seq2 = "".join(aa2)
            good_pair = True
        print()
    return good_pair, parallel, int(ss1[2:]), int(ss2[2:]), res_pairs, seq1, seq2

pdb_code = os.path.basename(sys.argv[1])

if len(sys.argv) < 3:
    db_name = "comp_seq_DB.db"
elif len(sys.argv) == 3:
    db_name = sys.argv[2]

con = sqlite3.connect(db_name)
cursor = con.cursor()

cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
if len(cursor.fetchall()) == 0:
    cursor.execute("CREATE TABLE pdb(code text, parallel integer, strand1 integer, strand2 integer, seq1 text, seq2 text)")
    cursor.execute("CREATE TABLE strand(code text, chain text, strand_num integer, range text)")

stored.strand = []
cmd.iterate("(ss S and name ca)","stored.strand.append((chain, resi))")
strand_list = grouping(stored.strand)
strand_pair = pair_match(strand_list)
for i in strand_pair:
    good_pair, parallel, ss1, ss2, res_pairs, seq1, seq2 = parallel_check(strand_list, i[0], i[1])
    if good_pair:
        range1 = ",".join(map(lambda x: x[0][1], res_pairs))
        chain1 = map(lambda x: x[0][0], res_pairs)[0]
        range2 = ",".join(map(lambda x: x[1][1], res_pairs))
        chain2 = map(lambda x: x[1][0], res_pairs)[0]
        if parallel:
            print(parallel, seq1, seq2)
            print(parallel, seq2, seq1)
        else:
            print(parallel, seq1, seq2)
            print(parallel, seq2[::-1], seq1[::-1])
        cursor.execute('SELECT * FROM pdb WHERE (code=? AND strand1=? AND strand2=?)', (pdb_code,ss1, ss2,))
        entry = cursor.fetchone()
        if entry is None:
            cursor.execute('INSERT INTO pdb (code, parallel, strand1, strand2, seq1, seq2) VALUES (?,?,?,?,?,?)', 
                    (pdb_code, int(parallel), ss1, ss2, seq1, seq2,))
            if parallel:
                cursor.execute('INSERT INTO pdb (code, parallel, strand1, strand2, seq1, seq2) VALUES (?,?,?,?,?,?)', 
                        (pdb_code, int(parallel), ss2, ss1, seq2, seq1,))
            else:
                cursor.execute('INSERT INTO pdb (code, parallel, strand1, strand2, seq1, seq2) VALUES (?,?,?,?,?,?)', 
                        (pdb_code, int(parallel), ss2, ss1, seq2[::-1], seq1[::-1],))
            cursor.execute('INSERT INTO strand (code, chain, strand_num, range) VALUES (?,?,?,?)', 
                    (pdb_code, chain1, ss1, range1,))
            cursor.execute('INSERT INTO strand (code, chain, strand_num, range) VALUES (?,?,?,?)', 
                    (pdb_code, chain2, ss2, range2,))
con.commit()
con.close()
