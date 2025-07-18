from Sequence import Sequence
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

seq= Sequence("data/mutated_1.fasta")
mutseq = Sequence("data/mutated_2.fasta")
detmut(seq, mutseq)