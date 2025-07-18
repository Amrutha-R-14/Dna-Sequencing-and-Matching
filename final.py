import sys
import tkinter as tk
from tkinter import filedialog, scrolledtext, Label, Entry
import io
import difflib
import networkx as nx
import matplotlib.pyplot as plt
from print_string_sequence import print_string_sequence
from Sequence import Sequence
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
class DNASequence:
    def __init__(self, data_path):
        self.Data = ""
        with open(data_path, "r") as file1:
            for line in file1.readlines():
                if line.startswith('>'):
                    self.Header = line[:-2]
                else:
                    self.Data = self.Data + line[:-2]
        self.DataSize = len(self.Data)

    def getItr(self, itr):
        return self.Data[itr]

    def getSize(self):
        return self.DataSize

    def find(self, ch):
        return self.Data.rfind(ch)

    def print(self):
        print_string_sequence(self.Data)

def hamming_distance(seq1, seq2):
    diff = difflib.ndiff(str(seq1), str(seq2))
    return sum(1 for item in diff if item.startswith('-'))

def generate_weighted_graph(species):
  G = nx.DiGraph()

  # Define edges and weights based on Hamming distance
  edges_weights = []
  keys = list(species.keys())
  for i in range(len(keys)):
    for j in range(i + 1, len(keys)):
      source = keys[i]
      target = keys[j]
      weight = hamming_distance(species[source], species[target])
      edges_weights.append((source, target, weight))

  # Check if only one sequence exists
  if len(keys) == 1:
     print("Cannot generate a weighted graph with only one sequence as input.")

  # Add edges to the graph with weights
  for source, target, weight in edges_weights:
    G.add_edge(source, target, weight=weight)

  return G


class RedirectText(io.StringIO):
    def __init__(self, widget):
        self.widget = widget

    def write(self, s):
        self.widget.insert(tk.END, s)
        self.widget.see(tk.END)


original_stdout = sys.stdout

# GUI setup
root = tk.Tk()
root.title("DNA SEQUENCING AND MATCHING APPLICATION")

# Function to select a file and display a hint
def select_file_and_display_hint(species_number):
    file_path = filedialog.askopenfilename(title=f"Select Species {species_number+1} File", filetypes=[("FASTA Files", "*.fasta")])
    if file_path:
        output_text.insert(tk.END, f"Selected Species {species_number+1} file: {file_path}\n")
        return file_path
    else:
        output_text.insert(tk.END, f"No file selected for Species {species_number+1}.\n")
        return None
def bm_approx_matching():
    def bm_approx_mat(pat, tx, maxmis):
        patlen = pat.getSize()
        txlen = tx.getSize()

        if patlen == 0 or txlen == 0 or patlen > txlen:
            return []

        pos = []

        i = 0
        while i <= txlen - patlen:
            mismc = 0

            # Separated loop for mismatch calculation
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

    file_path = [select_file_and_display_hint(i) for i in range(2)]
    if file_path:
        pattern = Sequence(file_path[0])
        seq1 = Sequence(file_path[1])
        approx_match_positions = bm_approx_mat(pattern, seq1, maxmis=2)
        print("\n\n")
        print("Boyer-Moore Approximate Matches:")
        print("===============================================================")
        if approx_match_positions:
            for pos in approx_match_positions:
                print(f"Approximate match starting at position {pos}")
        else:
            print("No Boyer-Moore approximate matches detected.")
        print("===============================================================")

# Function to create phylogenetic tree
def create_phylogenetic_tree():
    def create_phylogenic_tree(species, num_files):
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
                for j in range(i + 1, num_species):
                    s1 = keys[i]
                    s2 = keys[j]
                    mismatch_score[s1][s2] = minmatch(species[s1], species[s2], length)

            return mismatch_score

        def find_min_pos(score, num_species, seq_len):
            min_value = seq_len
            min_pos = (0, 0)

            for i in range(num_species):
                for j in range(i + 1, num_species):
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
                    mismatch_score[keys[i]][m] = (mismatch_score[keys[i]][keys[pos1]] + mismatch_score[keys[i]][
                        keys[pos2]]) / 2
                    del mismatch_score[keys[i]][keys[pos1]]
                    del mismatch_score[keys[i]][keys[pos2]]
                elif i > pos1 and i < pos2:
                    mismatch_score[keys[i]][m] = (mismatch_score[keys[pos1]][keys[i]] + mismatch_score[keys[i]][
                        keys[pos2]]) / 2
                    del mismatch_score[keys[i]][keys[pos2]]
                elif pos2 < i:
                    mismatch_score[keys[i]][m] = (mismatch_score[keys[pos1]][keys[i]] + mismatch_score[keys[pos2]][
                        keys[i]]) / 2

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

            while (num_species > 1):
                (i, j) = find_min_pos(score, num_species, seq_length)
                score = update_score(score, keys, i, j)
                keys = list(score)
                num_species = len(keys)

            return keys[0]

    num_files = int(entry_num_files.get())
    file_paths = [select_file_and_display_hint(i) for i in range(num_files)]
    species = {f"Species{i + 1}": Sequence(file_path) for i, file_path in enumerate(file_paths) if file_path}
    if len(species) >= 3:
        print("\n\n")
        print("Creating Phylogenic Tree:")
        print("===============================================================")
        print("Input Sizes:")
        for s, path in zip(species, file_paths):
            print(f"Length of sequence for species {s} is: {Sequence(path).getSize()}")

        # Printing Results
        tree = create_phylogenic_tree(species, num_files)
        print("\nTree Representation ([A, B] means A and B are more likely to be linked genetically)")
        print("\nOutput:")
        print(tree)
        print("===============================================================")

        # Creates a weighted graph
        G = nx.DiGraph()



        # Added nodes
        for node in tree:
            G.add_node(node)

        # Added edges with weights
        for i in range(num_files):
            for j in range(i + 1, num_files):
                species1 = f"Species{i + 1}"
                species2 = f"Species{j + 1}"
                weight = hamming_distance(species[file_paths[i]], species[file_paths[j]])

                G.add_edge(species1, species2, weight=weight)

        # Draw the graph
        pos = nx.spring_layout(G)
        labels = nx.get_edge_attributes(G, 'weight')
        nx.draw(G, pos, with_labels=True)
        nx.draw_networkx_edge_labels(G, pos, edge_labels=labels)
        plt.show()

    else:
        print("Not enough sequences selected for creating a phylogenetic tree.")


# Function to perform Smith-Waterman Alignment
def smith_waterman_alignment():
    def sw_algo(seq1, seq2, st):
        r = seq1.getSize()
        c = seq2.getSize()

       
        sm = [[0] * c for _ in range(r)]

        maxs = 0
        maxsp = (0, 0)

      
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

    file_paths = [select_file_and_display_hint(i) for i in range(2)]  
    seq1 = Sequence(file_paths[0])
    seq2 = Sequence(file_paths[1])
    print("\nSmith-Waterman Alignment:")
    aligned_seqs = sw_algo(seq1, seq2, score_table_sw)
    print("Sequence 1:", aligned_seqs[0])
    print("\nSequence 2:", aligned_seqs[1])

# Function to detect mutations
def detect_mutations():
    def detmut(seq, mutseq):
        def hamdis(seq1, seq2):
            return sum(1 for i in range(min(seq1.getSize(),seq2.getSize())) if seq1.getItr(i) != seq2.getItr(i))

        # Calculate Hamming distance and sequence length
        hamdist = hamdis(seq, mutseq)
        seqlen = min(seq.getSize(),mutseq.getSize())

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

    file_paths = [select_file_and_display_hint(i) for i in range(2)]  
    seq = Sequence(file_paths[0])
    mutseq = Sequence(file_paths[1])
    detmut(seq, mutseq)


def draw_weighted_graph():
    file_paths = [select_file_and_display_hint(i) for i in range(int(entry_num_files.get()))]
    species = {f"Species{i + 1}": DNASequence(file_path) for i, file_path in enumerate(file_paths) if file_path}
    weighted_phylogenetic_graph = generate_weighted_graph(species)

    fig, ax = plt.subplots(figsize=(18, 22))

    pos = nx.spring_layout(weighted_phylogenetic_graph)
    labels = nx.get_edge_attributes(weighted_phylogenetic_graph, 'weight')
    nx.draw(weighted_phylogenetic_graph, pos, with_labels=True, font_weight='bold', ax=ax)
    nx.draw_networkx_edge_labels(weighted_phylogenetic_graph, pos, edge_labels=labels, ax=ax)

    
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.pack()

    
    root.mainloop()



def execute_and_redirect(func):
    # Redirects standard output to the GUI
    output_text.config(state=tk.NORMAL)
    output_text.delete(1.0, tk.END)
    sys.stdout = RedirectText(output_text)

    # Calls the function
    func()

    # Restores the original standard output
    sys.stdout = original_stdout
    output_text.config(state=tk.DISABLED)

# Creates labels and entry widgets for input
label_num_files = Label(root, text="Number of Files for Phylogenetic Tree or Weighted Graph:")
label_num_files.pack(pady=5)
entry_num_files = Entry(root)
entry_num_files.pack(pady=5)


# Creates buttons for each function
bm_button = tk.Button(root, text="Approximate Matching", command=lambda: execute_and_redirect(bm_approx_matching))
bm_button.pack(pady=10)

phylo_button = tk.Button(root, text="Phylogenetic Tree Generation", command=lambda: execute_and_redirect(create_phylogenetic_tree))
phylo_button.pack(pady=10)

sw_button = tk.Button(root, text="Sequence Alignment", command=lambda: execute_and_redirect(smith_waterman_alignment))
sw_button.pack(pady=10)

mut_button = tk.Button(root, text="Mutation Detection", command=lambda: execute_and_redirect(detect_mutations))
mut_button.pack(pady=10)

# Creates buttons for each function
weighted_graph_button = tk.Button(root, text="Weighted Graph Generation", command=lambda: execute_and_redirect(draw_weighted_graph))
weighted_graph_button.pack(pady=10)

# Creates a scrolled text widget for output
output_text = scrolledtext.ScrolledText(root, wrap=tk.WORD, width=80, height=20)
output_text.pack(pady=10)

# Run the GUI
root.mainloop()
