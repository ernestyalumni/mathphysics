"""Check Z from the links against all transfer directions at exp(-2 beta) = 1/2.

This checks conventions, not the strong-coupling bound. Use integers until the
final normalization. Each frustrated plaquette supplies a factor 1/2.
"""
from fractions import Fraction
from itertools import combinations, product


def geometry(sides):
    sites = list(product(*(range(n) for n in sides)))
    index = {x: i for i, x in enumerate(sites)}
    steps = [[index[x[:j] + ((x[j] + 1) % n,) + x[j+1:]]
              for j, n in enumerate(sides)] for x in sites]
    return len(sites), steps


def face_masks(volume, steps, dim):
    # XOR, not union: a link traversed twice in Z2 makes no contribution.
    def link(s, i):
        return 1 << (s * dim + i)
    return [link(s, i) ^ link(steps[s][i], j) ^ link(steps[s][j], i) ^ link(s, j)
            for s in range(volume) for i, j in combinations(range(dim), 2)]


def frustrated(u, masks):
    return sum((u & m).bit_count() % 2 for m in masks)


def direct(sides):
    v, steps = geometry(sides)
    faces = face_masks(v, steps, 4)
    # Scaled by 2^(6v) to clear all plaquette denominators; average over 2^(4v) fields.
    numerator = sum(1 << (6*v - frustrated(u, faces)) for u in range(1 << (4*v)))
    return Fraction(numerator, 1 << (10*v))


def sliced(sides, axis):
    t = sides[axis]
    if t not in (1, 2):
        raise ValueError('this check expands only traces of the first and second powers')
    v, steps = geometry(sides[:axis] + sides[axis+1:])
    states = range(1 << (3*v))
    faces = face_masks(v, steps, 3)
    spatial = [frustrated(a, faces) for a in states]
    # Every temporal-link field supplies an endpoint difference on each spatial link.
    differences = [sum((((u >> s) ^ (u >> steps[s][j])) & 1) << (3*s + j)
                       for s in range(v) for j in range(3)) for u in range(1 << v)]

    def weight(a, b):
        # Charge both spatial half-weights on the incoming slice instead: this
        # preserves the product around a closed walk and avoids square roots.
        return sum(1 << (6*v - spatial[a] - (a ^ b ^ g).bit_count()) for g in differences)

    if t == 1:
        trace = sum(weight(a, a) for a in states)
    else:
        trace = sum(weight(a, b) * weight(b, a) for a, b in product(states, repeat=2))
    return Fraction(trace, 1 << (10*v*t))


if __name__ == '__main__':
    for sides in [(1, 1, 1, 1), (2, 1, 1, 1), (2, 2, 1, 1)]:
        z = direct(sides)
        for axis in range(4):
            trace = sliced(sides, axis)
            if trace != z:
                raise AssertionError((sides, axis, z, trace))
        print(sides, z, 'all four axes match')
