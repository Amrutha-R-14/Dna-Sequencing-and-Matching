from Sequence import Sequence
from print_string_sequence import print_string_sequence

def bm_approx_mat(pat, tx, maxmis):
    patlen = pat.getSize()
    txlen = tx.getSize()

    if patlen == 0 or txlen == 0 or patlen > txlen:
        return []

    pos = []

    i = 0
    while i <= txlen - patlen:
        mismc = 0

        # Separated loop for mismc calculation
        for j in range(patlen):
            if i + j < txlen and pat.getItr(j) != tx.getItr(i + j):
                mismc += 1

        if mismc <= maxmis:
            pos.append(i)

        if i + patlen < txlen:
            lastoc = pat.find(tx.getItr(i + patlen))
            i += max(1, patlen - lastoc)
        else:
            break

    return pos


pattern = Sequence("data/human_seq1.fasta")
seq1 = Sequence("data/human_seq1.fasta")
# Boyer-Moore approximate matching
approx_match_positions = bm_approx_mat(pattern, seq1, maxmis=2)
print("\n\n")
print("Boyer-Moore Approximate Matches:")
print("===============================================================")
if approx_match_positions:
    for pos in approx_match_positions:
        print(f"Approximate match starting at position {pos}")
else:
    print("No Boyer-Moore approximate matches detected.")
print("===============================================================\n")

def create_phylogenic_tree(species):
    def check_valid_sequences(species, keys):
        length = species[keys[0]].getSize()
        for k in keys:
            if length > species[k].getSize():
                length = species[k].getSize()
        return length
        
    def minmatch(seq1, seq2, length):
        return sum(1 for i in range(length) if seq1.getItr(i) != seq2.getItr(i))
    
    def generate_initial_scores(species, keys, num_species, length):
        mismatch_score = {}
        for k in keys:
            mismatch_score[k] = {}
       
        for i in range(num_species):
            for j in range(i+1, num_species):
                s1 = keys[i]
                s2 = keys[j]
                mismatch_score[s1][s2] = minmatch(species[s1], species[s2], length)

        return mismatch_score
    
    def find_min_pos(score, num_species, seq_len):
        min_value = seq_len
        min_pos = (0, 0)

        for i in range(num_species):
            for j in range(i+1, num_species):
                s1 = keys[i]
                s2 = keys[j]
                if score[s1][s2] < min_value:
                    min_value = score[s1][s2]
                    min_pos = (i, j)
        
        return min_pos
    
    def update_score(score, keys, pos1, pos2):
        m = "[" + keys[pos1] + ", " + keys[pos2] + "]";

        mismatch_score = score.copy()
        mismatch_score[m] = {}

        for i in range(num_species):
            if i < pos1:
                mismatch_score[keys[i]][m] = (mismatch_score[keys[i]][keys[pos1]] + mismatch_score[keys[i]][keys[pos2]]) / 2
                del mismatch_score[keys[i]][keys[pos1]]
                del mismatch_score[keys[i]][keys[pos2]]
            elif i > pos1 and i < pos2:
                mismatch_score[keys[i]][m] = (mismatch_score[keys[pos1]][keys[i]] + mismatch_score[keys[i]][keys[pos2]]) / 2
                del mismatch_score[keys[i]][keys[pos2]]
            elif pos2 < i:
                mismatch_score[keys[i]][m] = (mismatch_score[keys[pos1]][keys[i]] + mismatch_score[keys[pos2]][keys[i]]) / 2

        del mismatch_score[keys[pos1]]
        del mismatch_score[keys[pos2]]
        return mismatch_score

    keys = list(species)
    num_species = len(keys)
    seq_length = check_valid_sequences(species, keys)

    if seq_length == 0:
        print("INVALID: All Sequences should be of a length greater than 0.")
        return None

    if num_species == 1:
        return keys[0]
    else:
        score = generate_initial_scores(species, keys, num_species, seq_length)

        while(num_species > 1):
            (i, j) = find_min_pos(score, num_species, seq_length)
            score = update_score(score, keys, i, j)
            keys = list(score)
            num_species = len(keys)
            
        return keys[0]
species = {
    "human1": Sequence("data/human_seq1.fasta"),
    "frog": Sequence("data/blue_poison_arrow_frog.fasta"),
    "elephant": Sequence("data/african_elephant.fasta"),
    "whale": Sequence("data/blue_whale.fasta"),
    "dog": Sequence("data/dog.fasta"),
    "pigeon": Sequence("data/pigeon.fasta"),
    "rabbit": Sequence("data/rabbit.fasta"),
    "human2": Sequence("data/human_seq1.fasta")
}

print("\n\n")
print("Creating Phylogenic Tree:")
print("===============================================================")
print("Input Sizes:")
for s in species:
    print("Length of sequence for species", s, "is:", species[s].getSize())

# Printing Results
tree = create_phylogenic_tree(species)
print("\nTree Representation ([A, B] means A and B are linked genetically)")
print("\nOutput:")
print(tree)
print("===============================================================\n")

"""
    [[elephant, [frog, dog]], [rabbit, [pigeon, [human, whale]]]]
    means 
           __human
        __| 
     __|  |__whale
  __|  |__pigeon
 |  |__rabbit
 |
 |      __frog
 |   __|
 |__|  |__dog
    |__elephant

"""


def sw_algo(seq1, seq2, st):
    r = seq1.getSize()
    c = seq2.getSize()

    # First loop: Creating an empty matrix with zeros
    sm = [[0] * c for _ in range(r)]

    maxs = 0
    maxsp = (0, 0)

    # The rest of the code remains unchanged
    for i in range(1, r):
        for j in range(1, c):
            curr_1 = seq1.getItr(i)
            curr_2 = seq2.getItr(j)

            ms = st[curr_1][curr_2]
            gappen = st["gap"]["gap"]

            diags = sm[i - 1][j - 1] + ms
            ups = sm[i - 1][j] + gappen
            les = sm[i][j - 1] + gappen

            sm[i][j] = max(0, diags, ups, les)

            if sm[i][j] > maxs:
                maxs = sm[i][j]
                maxsp = (i, j)

    alseq1 = ""
    alseq2 = ""
    i, j = maxsp

    while i > 0 and j > 0 and sm[i][j] > 0:
        curr_1 = seq1.getItr(i)
        curr_2 = seq2.getItr(j)

        ms = st[curr_1][curr_2]
        gappen = st["gap"]["gap"]

        if sm[i][j] == sm[i - 1][j - 1] + ms:
            alseq1 = curr_1 + alseq1
            alseq2 = curr_2 + alseq2
            i -= 1
            j -= 1
        elif sm[i][j] == sm[i - 1][j] + gappen:
            alseq1 = curr_1 + alseq1
            alseq2 = '-' + alseq2
            i -= 1
        else:
            alseq1 = '-' + alseq1
            alseq2 = curr_2 + alseq2
            j -= 1

    return alseq1, alseq2

score_table_sw = {
    "A": {"A": 2, "C": -1, "G": -1, "T": -1, "N": -1},
    "C": {"A": -1, "C": 2, "G": -1, "T": -1, "N": -1},
    "G": {"A": -1, "C": -1, "G": 2, "T": -1, "N": -1},
    "T": {"A": -1, "C": -1, "G": -1, "T": 2, "N": -1},
    "N": {"A": -1, "C": -1, "G": -1, "T": -1, "N": -1},
    "gap": {"gap": -1}
}

seq1 = Sequence("data/pigeon.fasta")
seq2 = Sequence("data/dog.fasta")


aligned_seqs = sw_algo(seq1, seq2, score_table_sw)
print("Smith-Waterman Alignment:")
print("Sequence 1:", aligned_seqs[0])
print("\n")
print("Sequence 2:", aligned_seqs[1])




def detmut(seq, mutseq):
    # Function to calculate Hamming distance between two sequences
    def hamdis(seq1, seq2):
        return sum(1 for i in range(seq1.getSize()) if seq1.getItr(i) != seq2.getItr(i))

    # Calculate Hamming distance and sequence length
    hamdist = hamdis(seq, mutseq)
    seqlen = seq.getSize()

    # Print mutation detection using Hamming distance
    print("Hamming Distance between sequences:", hamdist)

    # Calculate and print percentage of mutation
    mutper = (hamdist / seqlen) * 100
    print("Percentage of Mutation:", mutper, "%")

    # Print specific mutations
    mut = []
    for i in range(seqlen):
        if seq.getItr(i) != mutseq.getItr(i):
            mut.append(f"Position {i + 1}: {seq.getItr(i)} -> {mutseq.getItr(i)}")

    if mut:
        print("\nMutations:")
        for mutation in mut:
            print(mutation)
    else:
        print("\nNo mutations detected.")

seq= Sequence("data/pigeon.fasta")
mutseq = Sequence("data/dog.fasta")
detmut(seq, mutseq)


species = {
"h": Sequence("data/human_seq1.fasta"),
"f": Sequence("data/blue_poison_arrow_frog.fasta"),
"e": Sequence("data/african_elephant.fasta"),
"w": Sequence("data/blue_whale.fasta"),
"d": Sequence("data/dog.fasta"),
"p": Sequence("data/pigeon.fasta"),
"r": Sequence("data/rabbit.fasta"),
}
import difflib
def hamming_distance(seq1, seq2):
    # Ensure the sequences are strings
    str_seq1 = str(seq1)
    str_seq2 = str(seq2)
    # Use ndiff to find differences between the sequences
    diff = difflib.ndiff(str_seq1, str_seq2)
    # Count differences that indicate mismatches
    hamdist = sum(1 for item in diff if item.startswith('-'))
    return hamdist
import networkx as nx
import matplotlib.pyplot as plt
def create_weighted_phylogenetic_graph():
    G = nx.DiGraph()
    
    # Define edges and weights
    edges_weights = [
    ("w", "h", 6),
    ("p", "h", 5),
    ("r", "p", 2),
    ("r", "f", 4),
    ("f", "d", 5),
    ("e", "d", 6),
    ("p","w",hamming_distance(species["p"],species["h"])),
    ("r","w",hamming_distance(species["r"],species["w"])),
    ("r","h",hamming_distance(species["r"],species["h"])),
    ("e","f",hamming_distance(species["e"],species["f"]))
    ]
    # Add edges to the graph with weights
    for source, target, weight in edges_weights:
        G.add_edge(source, target, weight=weight)
    return G
# Example usage
weighted_phylogenetic_graph = create_weighted_phylogenetic_graph()
# Draw the graph
pos = nx.spring_layout(weighted_phylogenetic_graph)
labels = nx.get_edge_attributes(weighted_phylogenetic_graph, 'weight')
nx.draw(weighted_phylogenetic_graph, pos, with_labels=True, font_weight='bold')
nx.draw_networkx_edge_labels(weighted_phylogenetic_graph, pos,edge_labels=labels)
plt.show()