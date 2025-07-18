def print_string_sequence(seq):
    PRINT_SIZE = 64
    SEQ_LENGTH = len(seq)
    i = 1

    while(i + PRINT_SIZE < SEQ_LENGTH):
        j = i + PRINT_SIZE - 1
        print(seq[i: j])
        i = j + 1

    j = SEQ_LENGTH
    print(seq[i: j])
