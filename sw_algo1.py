from Sequence import Sequence

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

seq1 = Sequence("data/human_seq1.fasta")
seq2 = Sequence("data/human_seq2.fasta")


aligned_seqs = sw_algo(seq1, seq2, score_table_sw)
print("Smith-Waterman Alignment:")
print("Sequence 1:", aligned_seqs[0])
print("\n")
print("Sequence 2:", aligned_seqs[1])

