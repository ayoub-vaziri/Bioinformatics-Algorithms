def banded_needleman_wunsch(seq1, seq2, indel, match_score=2, mismatch_penalty=-1, gap_penalty=-1):
    m = len(seq1)
    n = len(seq2)

    # Initialize the score matrix
    score_matrix = [[0] * (n + 1) for _ in range(m + 1)]

    # Initialize the traceback matrix
    traceback_matrix = [[""] * (n + 1) for _ in range(m + 1)]

    # Initialize the scores for the first row and column
    for i in range(1, m + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] + gap_penalty
    for j in range(1, n + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] + gap_penalty

    # Fill in the score and traceback matrices
    for i in range(1, m + 1):
        for j in range(max(1, i - indel), min(n, i + indel) + 1):
            if j > 0 and j <= n:
                match = score_matrix[i - 1][j - 1] + (match_score if seq1[i - 1] == seq2[j - 1] else mismatch_penalty)
                delete = score_matrix[i - 1][j] + gap_penalty
                insert = score_matrix[i][j - 1] + gap_penalty

                # Find the maximum score
                max_score = max(match, delete, insert)

                score_matrix[i][j] = max_score

                # Update the traceback matrix
                if max_score == match:
                    traceback_matrix[i][j] = "M"  # Match
                elif max_score == delete:
                    traceback_matrix[i][j] = "D"  # Delete
                else:
                    traceback_matrix[i][j] = "I"  # Insert

    # Traceback to find the alignment
    align1 = ""
    align2 = ""
    i, j = m, n
    while i > 0 and j > 0:
        if traceback_matrix[i][j] == "M":
            align1 = seq1[i - 1] + align1
            align2 = seq2[j - 1] + align2
            i -= 1
            j -= 1
        elif traceback_matrix[i][j] == "D":
            align1 = seq1[i - 1] + align1
            align2 = "-" + align2
            i -= 1
        else:
            align1 = "-" + align1
            align2 = seq2[j - 1] + align2
            j -= 1

    return align1, align2, score_matrix[m][n]

def read_input(file_path):
    with open(file_path, 'r') as input_file:
        seq1 = input_file.readline().strip()
        seq2 = input_file.readline().strip()
        indel = int(input_file.readline().strip())
    return seq1, seq2, indel
    
# Read input from file
input_file_path = '/home/ayoub/Python/input.txt'
seq1, seq2, indel = read_input(input_file_path)

# Execute the algorithm
alignment_result = banded_needleman_wunsch(seq1, seq2, indel)

print("Sequence Alignment 1:", alignment_result[0])
print("Sequence Alignment 2:", alignment_result[1])
print("Optimal Alignment Score:", alignment_result[2])

# Write the result to a file
output_file_path = '/home/ayoub/Python/output.txt'
with open(output_file_path, 'w') as output_file:
    output_file.write("Sequence Alignment 1: {}\n".format(alignment_result[0]))
    output_file.write("Sequence Alignment 2: {}\n".format(alignment_result[1]))
    output_file.write("Optimal Alignment Score: {}\n".format(alignment_result[2]))
