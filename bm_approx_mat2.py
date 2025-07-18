from Sequence import Sequence
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