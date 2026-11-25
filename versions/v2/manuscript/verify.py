import hashlib
import itertools
import json
from functools import lru_cache
from math import comb, gcd
from pathlib import Path


def literal_count(occurrences):
    size = 3 * occurrences + 2
    total = 0
    for permutation in itertools.permutations(range(size)):
        valid = True
        for start in range(0, 3 * occurrences, 3):
            if not (permutation[start] < permutation[start + 3]
                    < permutation[start + 2] < permutation[start + 1]
                    < permutation[start + 4]):
                valid = False
                break
        total += valid
    return total


def poset_count(occurrences):
    return prescribed_count([2] * (occurrences - 1))


def prescribed_count(overlaps):
    starts = [0]
    for overlap in overlaps:
        starts.append(starts[-1] + 5 - overlap)
    size = starts[-1] + 5
    predecessors = [0] * size
    for start in starts:
        chain = [start, start + 3, start + 2, start + 1, start + 4]
        for lower, upper in zip(chain, chain[1:]):
            predecessors[upper] |= 1 << lower
    full = (1 << size) - 1

    @lru_cache(maxsize=None)
    def extensions(chosen):
        if chosen == full:
            return 1
        total = 0
        for vertex, required in enumerate(predecessors):
            if not chosen & (1 << vertex) and chosen & required == required:
                total += extensions(chosen | (1 << vertex))
        return total

    count = extensions(0)
    states = extensions.cache_info().currsize
    extensions.cache_clear()
    return count, states


def transfer_counts(limit):
    row = [1] * (2 * limit + 2)
    counts = [row[3]]
    for step in range(1, limit):
        next_row = [0] * (len(row) - 2)
        running = 0
        for argument in range(1, len(next_row)):
            running += argument * row[argument + 2]
            next_row[argument] = running
        row = next_row
        counts.append(row[3])
    return counts


def convolution(left, right, degree):
    return [sum(left[index] * right[power - index]
                for index in range(power + 1)) % 4
            for power in range(degree + 1)]


def series_checks(degree):
    catalan = [comb(2 * index, index) // (index + 1)
               for index in range(degree + 1)]
    squared = convolution(catalan, catalan, degree)
    series_v = [(3 * entry - (2 if index == 0 else 0)) % 4
                for index, entry in enumerate(squared)]
    series_u = [1] + series_v[:-1]
    square_v = convolution(series_v, series_v, degree)
    for power in range(degree + 1):
        value = series_v[power]
        if power >= 1:
            value += 2 * series_v[power - 1]
        if power >= 2:
            value += square_v[power - 2]
        assert value % 4 == (1 if power == 0 else 0)
    powers = [[1] + [0] * degree]
    for exponent in range(1, 11):
        powers.append(convolution(powers[-1], series_v, degree))
    family = []
    for argument in range(20):
        family.append(convolution(series_u, powers[argument // 2], degree)
                      if argument % 2 == 0 else powers[argument // 2 + 1])
    for argument in range(17):
        assert family[argument][0] == 1
        for power in range(1, degree + 1):
            adjacent = (family[argument + 1][power]
                        + (2 * argument + 2) * family[argument + 1][power - 1]
                        - (2 * argument + 3) * family[argument + 2][power - 1])
            assert family[argument][power] == adjacent % 4
            recurrence = ((2 * argument + 1) * family[argument + 1][power - 1]
                          - sum(family[index][power - 1]
                                for index in range(1, argument + 1)))
            assert family[argument][power] == recurrence % 4


def modular_model_checks(counts):
    records = []
    for modulus in range(2, 101):
        if gcd(modulus, 6) != 1:
            continue
        bound = (modulus - 1) // 6 if modulus % 6 == 1 else (modulus + 1) // 6
        cutoff = next(height for height in range(1, bound + 1)
                      if (36 * height ** 2 - 1) % modulus == 0)
        vector = [1] + [0] * (cutoff - 1)
        for count in counts:
            assert vector[0] == count % modulus
            following = [0] * cutoff
            for height, weight in enumerate(vector):
                following[height] += (12 * height + 6) * weight
                if height + 1 < cutoff:
                    following[height + 1] += weight
                if height:
                    following[height - 1] += (36 * height ** 2 - 1) * weight
            vector = [entry % modulus for entry in following]
        records.append({'modulus': modulus, 'cutoff': cutoff, 'bound': bound})
    denominator = [1, 4, 20, 20, 15]
    numerator = [1, 10, 15, 15]

    def continuant(first, last):
        previous, current = [0], [1]
        for height in range(first, last + 1):
            following = current[:] + [0]
            for power, coefficient in enumerate(current):
                following[power + 1] -= (12 * height + 6) * coefficient
            if height > first:
                for power, coefficient in enumerate(previous):
                    following[power + 2] -= (36 * height ** 2 - 1) * coefficient
            previous, current = current, [coefficient % 25 for coefficient in following]
        return current

    assert continuant(0, 3) == denominator
    assert continuant(1, 3) == numerator
    for degree in range(len(counts)):
        coefficient = sum(denominator[index] * counts[degree - index]
                          for index in range(min(degree, 4) + 1)) % 25
        assert coefficient == (numerator[degree] if degree < 4 else 0)
    return {'moduli': records, 'first_j': 1, 'last_j': len(counts),
            'modulo_25_rational_series_degree_inclusive': len(counts) - 1,
            'modulo_25_determinants': {'numerator': numerator, 'denominator': denominator},
            'matches': True}


def separator_vertices(overlaps):
    starts = [0]
    for overlap in overlaps:
        starts.append(starts[-1] + 5 - overlap)
    size = starts[-1] + 5
    predecessors = [0] * size
    for start in starts:
        chain = [start, start + 3, start + 2, start + 1, start + 4]
        for lower, upper in zip(chain, chain[1:]):
            predecessors[upper] |= 1 << lower
    for pivot in range(size):
        for vertex in range(size):
            if predecessors[vertex] & (1 << pivot):
                predecessors[vertex] |= predecessors[pivot]

    def is_chain(subset):
        cardinality = subset.bit_count()
        comparable_pairs = sum((predecessors[vertex] & subset).bit_count()
                               for vertex in range(size) if subset & (1 << vertex))
        return comparable_pairs == cardinality * (cardinality - 1) // 2

    separators = []
    for vertex in range(size):
        lower = predecessors[vertex]
        upper = sum(1 << other for other in range(size)
                    if predecessors[other] & (1 << vertex))
        if (lower | upper).bit_count() == size - 1 and not is_chain(lower) and not is_chain(upper):
            separators.append(vertex)
    return separators


def separator_checks():
    cases = []
    for length in range(7):
        for overlaps in itertools.product([1, 2], repeat=length):
            blocks = [1]
            for overlap in overlaps:
                if overlap == 1:
                    blocks.append(1)
                else:
                    blocks[-1] += 1
            vertices = separator_vertices(overlaps)
            expected = sum(block >= 2 for block in blocks) >= 2
            assert bool(vertices) == expected
            cases.append({'overlaps': list(overlaps), 'blocks': blocks,
                          'separator_vertices_zero_based': vertices})
    for windows in range(1, 17):
        assert not separator_vertices([2] * (windows - 1))
    return {'maximum_windows': 7, 'words_checked': len(cases),
            'pure_overlap_two_maximum_windows': 16, 'cases': cases, 'matches': True}


def generic_pattern_count(pattern, overlap, windows):
    relative_positions = sorted(range(len(pattern)), key=lambda position: pattern[position])
    stride = len(pattern) - overlap
    size = len(pattern) + (windows - 1) * stride
    predecessors = [0] * size
    for occurrence in range(windows):
        chain = [occurrence * stride + position for position in relative_positions]
        for lower, upper in zip(chain, chain[1:]):
            predecessors[upper] |= 1 << lower
    full = (1 << size) - 1

    @lru_cache(maxsize=None)
    def extensions(chosen):
        if chosen == full:
            return 1
        return sum(extensions(chosen | (1 << vertex))
                   for vertex, required in enumerate(predecessors)
                   if not chosen & (1 << vertex) and chosen & required == required)

    actual = extensions(0)
    extensions.cache_clear()
    return actual


def catalan_sibling_checks():
    records = []
    for pattern, overlap in [('1324', 2), ('14253', 3), ('136245', 3),
                             ('142536', 4), ('1376245', 3), ('369247815', 5),
                             ('162738495', 7)]:
        values = []
        for windows in range(1, 9):
            actual = generic_pattern_count(pattern, overlap, windows)
            assert actual == comb(2 * windows, windows) // (windows + 1)
            values.append(actual)
        records.append({'pattern': pattern, 'overlap': overlap, 'counts_j_1_to_8': values})
    return {'cases': records, 'matches': True}


def equivalent_pattern_checks(counts):
    records = []
    for pattern in ['14325', '52341', '25314', '41352']:
        values = [generic_pattern_count(pattern, 2, windows) for windows in range(1, 9)]
        assert values == counts[:8]
        records.append({'pattern': pattern, 'overlap': 2, 'counts_j_1_to_8': values})
    return {'cases': records, 'matches': True}


def main():
    limit = 256
    counts = transfer_counts(limit)
    literal = {str(occurrences): literal_count(occurrences)
               for occurrences in [1, 2]}
    independent = []
    for occurrences in range(1, 9):
        count, ideals = poset_count(occurrences)
        assert count == counts[occurrences - 1]
        if str(occurrences) in literal:
            assert count == literal[str(occurrences)]
        independent.append({'j': occurrences, 'count': count, 'order_ideals': ideals})
    for occurrences, count in enumerate(counts, start=1):
        binary = occurrences + 1
        dyadic = binary & (binary - 1) == 0
        assert count % 5 == 1
        assert count % 35 == (1 if occurrences % 2 else 6)
        assert count % 10 == (1 if dyadic else 6)
        assert (count % 4 != 0) == (binary.bit_count() <= 2)
        if occurrences >= 2:
            catalan = comb(2 * occurrences, occurrences) // (occurrences + 1)
            assert count % 4 == 3 * catalan % 4
            digit_count = binary.bit_count()
            assert catalan % 4 == (1 if digit_count == 1 else 2 if digit_count == 2 else 0)
            residue = {0: 16, 1: 11, 2: 6}[catalan % 4]
            assert count % 20 == residue
            assert count % 140 == (36 * pow(6, occurrences - 1, 140) + 35 * catalan) % 140
    for argument in range(100):
        first = [0, 1, 3, 1, 0]
        stable = [0, 1, 1, 1, 0]
        assert (sum(range(1, argument + 1)) % 5) == first[argument % 5]
        for profile in [first, stable]:
            value = sum(weight * profile[(weight + 2) % 5]
                        for weight in range(1, argument + 1))
            assert value % 5 == stable[argument % 5]
    series_checks(64)
    degree = 64
    fraction = [1] + [0] * degree
    for level in range(33, -1, -1):
        denominator = [1, -(12 * level + 6)] + [
            -(36 * (level + 1) ** 2 - 1) * fraction[power - 2]
            for power in range(2, degree + 1)]
        inverse = [1]
        for power in range(1, degree + 1):
            inverse.append(-sum(denominator[index] * inverse[power - index]
                                for index in range(1, power + 1)) % 20)
        fraction = inverse
    assert fraction == [count % 20 for count in counts[:degree + 1]]
    mixed = []
    for length in range(4):
        for overlaps in itertools.product([1, 2], repeat=length):
            blocks = [1]
            for overlap in overlaps:
                if overlap == 1:
                    blocks.append(1)
                else:
                    blocks[-1] += 1
            expected = 1
            for block in blocks:
                expected *= counts[block - 1]
            actual, ideals = prescribed_count(overlaps)
            assert actual == expected
            assert actual % 35 == pow(6, sum(block - 1 for block in blocks), 35)
            assert actual % 10 == (1 if all((block + 1).bit_count() == 1
                                           for block in blocks) else 6)
            assert (actual % 4 != 0) == (sum((block + 1).bit_count() - 1
                                            for block in blocks) <= 1)
            mixed.append({'overlaps': list(overlaps), 'count': actual})
    known = [1, 6, 71, 1266, 30206, 902796, 32420011, 1359292626,
             65164480466, 3515569641156, 210779736073446, 13903319821066836,
             1000559812125494076, 78012524487061315416, 6550837823204594551731,
             589404446176366002280146]
    assert counts[:16] == known
    result = {
        'arithmetic': 'exact integers and coefficient arithmetic modulo four',
        'literal_permutations': literal,
        'comparison_poset_counts': independent,
        'transfer_and_congruence_range': {'first_j': 1, 'last_j': limit},
        'source_numerical_table': {'first_j': 1, 'last_j': 16, 'matches': True},
        'period_five_profiles': {'first_k': 0, 'last_k': 99, 'matches': True},
        'series': {'degree_inclusive': 64, 'first_k': 0, 'last_k': 16, 'matches': True},
        'mixed_overlaps': {'maximum_windows': 4, 'cases': mixed, 'matches': True},
        'jacobi_fraction': {'depth': 34, 'degree_inclusive': degree, 'modulus': 20, 'matches': True},
        'extended_congruences': {'moduli': [35, 140], 'first_j': 1,
                                'last_j': limit, 'modulo_140_first_j': 2, 'matches': True},
        'modular_models': modular_model_checks(counts),
        'separators': separator_checks(),
        'catalan_siblings': catalan_sibling_checks(),
        'equivalent_patterns': equivalent_pattern_checks(counts),
        'catalan_residue_classification': {'first_j': 2, 'last_j': limit,
                                          'possible_residues_modulo_4': [0, 1, 2], 'matches': True},
        'script_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'status': 'pass',
        'scope': 'Finite checks supplement the written proof for all positive j.'
    }
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
