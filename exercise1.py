import argparse
from itertools import product


def nucl_to_number(nucl):
    a = 1
    t = 2
    g = 3
    c = 0
    if nucl == 'A':
        number = a
    elif nucl == 'T':
        number = t
    elif nucl == 'G':
        number = g
    elif nucl == 'C':
        number = c
    else:
        raise Exception('Wrong character in sequence: {}'.format(nucl))
    return number

def number_to_nucl(number):
    a = 1
    t = 2
    g = 3
    c = 0
    if number == a:
        nucl = 'A'
    elif number == t:
        nucl = 'T'
    elif number == g:
        nucl = 'G'
    elif number == c:
        nucl = 'C'
    else:
        raise Exception('Wrong character in sequence: {}'.format(number))
    return nucl
def kmer_hashing(kmer):
    hash = 0
    multiplier = 1
    for index, nucleotide in enumerate(kmer[::-1]):
        number = str(nucl_to_number(nucleotide))
        if index > 0:
            multiplier = multiplier * 4
        hash += int(number) * multiplier
    return hash

def hash_to_kmer(hash, length, restored_kmer = None, divider = 4):
    if restored_kmer is None:
        restored_kmer = []
    residue = hash % divider
    main = hash // divider
    restored_kmer.append(number_to_nucl(residue))
    if main < divider:
        restored_kmer.append(number_to_nucl(main))
        if len(restored_kmer) < length:
            restored_kmer.append(number_to_nucl(0))
    else:
        hash_to_kmer(main, length, restored_kmer)
    return restored_kmer[::-1]

def hash_table_preparation(sequence, length, divider = 4):
    size = len(list(product(['A', 'T', 'G', 'C'], repeat=length)))
    table = [[None, None, None, None] for i in range(size)]
    current_hash = kmer_hashing(sequence[0:length])
    table[current_hash][0] = 1
    table[current_hash][1] = 0
    table[current_hash][2] = 0
    table[current_hash][3] = 0
    multiplier = divider**(length-1)
    previous_start_nucl = sequence[0]
    for nucl_i in range(length, len(sequence)):
        previous_start_number = nucl_to_number(previous_start_nucl) * multiplier
        current_hash = (current_hash - previous_start_number) * divider + nucl_to_number(sequence[nucl_i])
        if table[current_hash][0]:
            table[current_hash][0] += 1
        else:
            table[current_hash][0] = 1
        if table[current_hash][1] is None:
            table[current_hash][1] = nucl_i-length+1
        if table[current_hash][2] is None or table[current_hash][2] < nucl_i:
            table[current_hash][2] = nucl_i-length+1
            if table[current_hash][1] == table[current_hash][2]:
                table[current_hash][3] = table[current_hash][1]
            else:
                table[current_hash][3] = (table[current_hash][2] - table[current_hash][1])/table[current_hash][0]
        previous_start_nucl = sequence[nucl_i - length + 1]
    return table


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='gather kmer stats over given data')
    parser.add_argument('-k', '--kmer_length', type=int, default=11)
    args = parser.parse_args()

    length = args.kmer_length
    sequence = str(input('Input sequence:\n'))
    if len(sequence) < length:
        raise Exception('k-mer size cannot be greater then sequence length!')
    table = hash_table_preparation(sequence, length)
    for index in range(len(table)):
        if table[index][0]:
            print('kmer: ',''.join(hash_to_kmer(index, length)),' - ', '\t'.join([str(i) for i in table[index]]))