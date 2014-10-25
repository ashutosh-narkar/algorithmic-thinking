#!/usr/bin/env python
'''
Project 4 - Computing alignments of sequences
'''


def build_scoring_matrix(alphabet, diag_score, off_diag_score, dash_score):
    '''
    Takes as input a set of characters 'alphabet' and three scores
    'diag_score', 'off_diag_score', and 'dash_score'.
    The function returns a dictionary of dictionaries whose entries
    are indexed by pairs of characters in 'alphabet' plus '-'.
    The score for any entry indexed by one or more dashes is 'dash_score'.
    The score for the remaining diagonal entries is 'diag_score'.
    The score for the remaining off-diagonal entries is 'off_diag_score'.

    Although an alignment with two matching dashes is not allowed,
    the scoring matrix should still include an entry for
    two dashes (which will never be used).
    '''
    scoring_matrix = {}

    for char in alphabet:

        # diagonal score
        scoring_matrix[char] = {char: diag_score}

        # dash score
        scoring_matrix[char].update({'-': dash_score})

        # non-diagonal scores
        non_diag_chars = [item for item in alphabet if item not in (char, '-')]
        for non_diag_char in non_diag_chars:
            scoring_matrix[char].update({non_diag_char: off_diag_score})

    # correct dash score
    for item in alphabet:
        if '-' in scoring_matrix:
            scoring_matrix['-'].update({item: dash_score})
        else:
            scoring_matrix['-'] = {item: dash_score}
    scoring_matrix['-'].update({'-': dash_score})

    return scoring_matrix


def compute_alignment_matrix(seq_x, seq_y, scoring_matrix, global_flag):
    '''
    Takes as input two sequences seq_x and seq_y whose elements share a
    common alphabet
    with the scoring matrix 'scoring_matrix'
    Returns the alignment matrix for seq_x and seq_y
    Entries in the alignment matrix will be indexed by their row and column
    with integer indices starting at zero.
    This matrix will be modelled as a list of lists

    If global_flag is True, compute Global Alignment Scores
    If global_flag is False, compute Local Alignment Scores

    To compute a local alignment matrix, whenever a negative value is being
    assigned to alignment_matrix[i][j], replace with 0
    '''

    alignment_matrix = [] # alignment matrix

    alignment_matrix.insert(0, [0]) # s[0][0] = 0

    for ival in range(1, len(seq_x) + 1):
        # accessing the appropriate score from the scoring matrix
        score = scoring_matrix[seq_x[ival - 1]]['-']
        val = alignment_matrix[ival - 1][0] + score
        if val < 0 and not global_flag:
            alignment_matrix.insert(ival, [0])
        else:
            alignment_matrix.insert(ival, [val])

    for jval  in range(1, len(seq_y) + 1):
        # accessing the appropriate score from the scoring matrix
        score = scoring_matrix['-'][seq_y[jval - 1]]
        val = alignment_matrix[0][jval - 1] + score
        if val < 0 and not global_flag:
            alignment_matrix[0].insert(jval, 0)
        else:
            alignment_matrix[0].insert(jval, val)

    for ival in range(1, len(seq_x) + 1):
        for jval  in range(1, len(seq_y) + 1):
            val_1 = alignment_matrix[ival - 1][jval - 1] + \
                    scoring_matrix[seq_x[ival - 1]][seq_y[jval - 1]]
            val_2 = alignment_matrix[ival - 1][jval] + \
                    scoring_matrix[seq_x[ival - 1]]['-']
            val_3 = alignment_matrix[ival][jval - 1] + \
                    scoring_matrix['-'][seq_y[jval - 1]]
            val = max(val_1, val_2, val_3)
            if val < 0 and not global_flag:
                alignment_matrix[ival].insert(jval, 0)
            else:
                alignment_matrix[ival].insert(jval, val)
    return alignment_matrix


def compute_global_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    '''
    Takes as input two sequences seq_x and seq_y whose elements share a common
    alphabet with the scoring matrix 'scoring_matrix'.
    Compute a global alignment of seq_x and seq_y using the global alignment
    matrix 'alignment_matrix'.

    The function returns a tuple of the form (score, align_x, align_y) where
    score is the score of the global alignment align_x and align_y.

    Note that align_x and align_y should have the same length and may include
    the padding character '-'


    **Principle**
    Obtaining the optimal global alignment itself can be done via traceback ,
    which basically means walking back in the matrix from
    entry S[|i||j|] to entry S[0,0] along the path that has score S[|i||j|] ,
    and while doing so, creating the alignment itself.
    where i = len(seq_x) j = len(seq_y)

    '''
    xstr = ''
    ystr = ''

    ival = len(seq_x)
    jval = len(seq_y)

    while ival != 0 and jval != 0:
        if alignment_matrix[ival][jval] == alignment_matrix[ival - 1][jval - 1] + \
           scoring_matrix[seq_x[ival - 1]][seq_y[jval - 1]]:
            xstr = seq_x[ival - 1] + xstr
            ystr = seq_y[jval - 1] + ystr
            ival -= 1
            jval -= 1
        else:
            if alignment_matrix[ival][jval] == alignment_matrix[ival - 1][jval] + \
               scoring_matrix[seq_x[ival - 1]]['-']:
                xstr = seq_x[ival - 1] + xstr
                ystr = '-' + ystr
                ival -= 1
            else:
                xstr = '-' + xstr
                ystr = seq_y[jval - 1] + ystr
                jval -= 1

    while ival != 0:
        xstr = seq_x[ival - 1] + xstr
        ystr = '-' + ystr
        ival -= 1

    while jval != 0:
        xstr = '-' + xstr
        ystr = seq_y[jval - 1] + ystr
        jval -= 1
    return (alignment_matrix[len(seq_x)][len(seq_y)], xstr, ystr)


def compute_local_alignment(seq_x, seq_y, scoring_matrix, alignment_matrix):
    '''
    Takes as input two sequences seq_x and seq_y whose elements share a common
    alphabet with the scoring matrix 'scoring_matrix'.
    Compute a local alignment of seq_x and seq_y using the local alignment matrix
    'alignment_matrix'.

    The function returns a tuple of the form (score, align_x, align_y) where
    score is the score of the local alignment align_x and align_y.
    Note that align_x and align_y should have the same length and may include
    the padding character '-'

    **Principle**
    Start the traceback from the entry in the local alignment matrix that has
    the maximum value over the entire matrix and trace backwards
    similar to compute global alignment
    Stop the traceback when the first entry with value 0 is encountered.
    If the local alignment matrix has more than one entry that has the maximum
    value, any entry with maximum value may be used as the starting entry
    '''

    # find the entry with max value in the local alignment matrix
    val, index = find_max_element_and_index(alignment_matrix)

    ival = index[0]
    jval = index[1]
    xstr = ''
    ystr = ''

    if alignment_matrix[ival][jval] == 0:
        return (0, xstr, ystr)

    while ival != 0 and jval != 0:
        if alignment_matrix[ival - 1][jval - 1] != 0 and \
           alignment_matrix[ival][jval] == alignment_matrix[ival - 1][jval - 1] + \
           scoring_matrix[seq_x[ival - 1]][seq_y[jval - 1]]:
            xstr = seq_x[ival - 1] + xstr
            ystr = seq_y[jval - 1] + ystr
            ival -= 1
            jval -= 1

        elif alignment_matrix[ival - 1][jval - 1] == 0:
            break

        else:
            if alignment_matrix[ival - 1][jval] != 0 and \
               alignment_matrix[ival][jval] == alignment_matrix[ival - 1][jval] + \
               scoring_matrix[seq_x[ival - 1]]['-']:
                xstr = seq_x[ival - 1] + xstr
                ystr = '-' + ystr
                ival -= 1
            elif alignment_matrix[ival - 1][jval] == 0:
                break

            else:
                if alignment_matrix[ival][jval - 1] == 0:
                    break
                xstr = '-' + xstr
                ystr = seq_y[jval - 1] + ystr
                jval -= 1

    # above loop can break without ever evaluating the entire sequence
    # so we need to take the unevaluated part of the sequence and
    # compute a local alignment matrix
    # then use this to compute a global alignment on the unevaluated
    # part of the sequences
    local = compute_alignment_matrix(seq_x[:ival], seq_y[:jval], scoring_matrix, False)
    alignment = compute_global_alignment(seq_x[:ival], seq_y[:jval], scoring_matrix, local)

    xstr = alignment[1] + xstr
    ystr = alignment[2] + ystr

    # if only one of xstr and ystr starts with a '-', we can ignore that matching.
    # Hence we wont get penalized for mismatch
    while xstr.startswith('-') and not ystr.startswith('-') \
          or not xstr.startswith('-') and ystr.startswith('-'):
        xstr = xstr[1:]
        ystr = ystr[1:]

    return (val, xstr, ystr)


def find_max_element_and_index(matrix):
    '''
    Given a matrix modelled as a list of lists, return the
    max element in the matrix and its indices
    '''
    max_val = 0
    index = [0, 0]
    for row in range(len(matrix)):
        if max(matrix[row]) > max_val:
            max_val = max(matrix[row])
            index[0] = row
            index[1] = matrix[row].index(max_val)
    return (max_val, index)


def main():
    '''
    tests
    '''
    alphabet = set(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l',
                    'm', 'n', 'o', 'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x',
                    'y', 'z', '-'])
    scoring_matrix = build_scoring_matrix(alphabet, 2, -1, -1)
    #print scoring_matrix
    #xstrtest = 'happypedestrianwalker'
    #ystrtest = 'sadpedesxtriandriver'
    xstrtest = 'actact'
    ystrtest = 'agcta'

    # global alignment
    #global_alignment_matrix = compute_alignment_matrix(xstrtest, ystrtest,
    #                                                   scoring_matrix, True)
    #global_alignment = compute_global_alignment(xstrtest, ystrtest,
    #                                  scoring_matrix, global_alignment_matrix)
    #print global_alignment

    # local alignment
    local_alignment_matrix = compute_alignment_matrix(xstrtest, ystrtest,
                                                      scoring_matrix, False)
    #print local_alignment_matrix
    local_alignment = compute_local_alignment(xstrtest, ystrtest, scoring_matrix,
                                              local_alignment_matrix)
    print local_alignment


if __name__ == '__main__':
    main()
