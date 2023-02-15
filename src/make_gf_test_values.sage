from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
import itertools as it

F = GF(2)
R.<X> = F[]
p = X^128 + X^7 + X^2 + X + 1
L = GF(2^128, 'z', p)
L.inject_variables()

q = X^8 + X^4 + X^3 + X + 1
K = GF(2^8, 'y', q)
K.inject_variables()


H = Hom(K, L)
f = FiniteFieldHomomorphism_generic(H)

K_values =  [
       K.from_integer(0x4b),
       K.from_integer(0x27),
       K.from_integer(0xb1),
       K.from_integer(0xfc),
   ]

L_values =  [
       L.from_integer(0x616ce329b8aee6b6c752890eaec5fdca),
       L.from_integer(0xc9d404e824cf265e31f38f85087a788f),
       L.from_integer(0x7eceb30086f94c398a6881ea95c275b2),
       L.from_integer(0x41a376dbe3eebf96a8f49dff52f12e1c),
   ]

L_sums = [x + y for x, y in it.product(values, repeat=2)]
L_prods = [x * y for x, y in it.product(values, repeat=2)]
L_invs = [x^-1 for x in values]
K_in_L_embeddings = [f(v) for v in K_values]

code_embedding_constants = '''
    pub const GF2P8_EMBEDDING_POX: [Self; 8] = [
        Self::ONE,
        {}
    ];
'''.format(',\n        '.join(f'Self::from_u128(0x{(f(y^i)).to_integer():032x})' for i in range(1, 8)))

print(code_embedding_constants)

L_code_values = '''
    const GF2P128_VALUES: [u128; 4] = [
        {}
    ];
    const GF2P128_SUMS: [u128; 16] = [
        {}
    ];
    const GF2P128_PRODS: [u128; 16] = [
        {}
    ];
    const GF2P128_INVS: [u128; 4] = [
        {}
    ];
    const GF2P8_IN_GF2P128_EMBEDDINGS: [u128; 4] = [
        {}
    ];
'''.format(
        ',\n        '.join(f'0x{v.to_integer():032x}' for v in L_values),
        ',\n        '.join(f'0x{v.to_integer():032x}' for v in L_sums),
        ',\n        '.join(f'0x{v.to_integer():032x}' for v in L_prods),
        ',\n        '.join(f'0x{v.to_integer():032x}' for v in L_invs),
        ',\n        '.join(f'0x{v.to_integer():032x}' for v in K_in_L_embeddings),
           )

print(L_code_values)

L_code_test_sums = '''
        let values = GF2P128_VALUES.map(F::from_u128);
        let sums = GF2P128_SUMS.map(F::from_u128);
        {}
'''.format('\n        '.join(f'assert_eq!(values[{i}] + values[{j}], sums[{4 * i + j}]);'
                             for i, j in it.product(range(4), repeat=2)))

print(L_code_test_sums)

L_code_test_prods = '''
        let values = GF2P128_VALUES.map(F::from_u128);
        let prods = GF2P128_PRODS.map(F::from_u128);
        {}
'''.format('\n        '.join(f'assert_eq!(values[{i}] * values[{j}], prods[{4 * i + j}]);'
                             for i, j in it.product(range(4), repeat=2)))

print(L_code_test_prods)

L_code_test_invs = '''
        let values = GF2P128_VALUES.map(F::from_u128);
        let invs = GF2P128_INVS.map(F::from_u128);
        {}
'''.format('\n        '.join(f'assert_eq!(values[{i}].invert().unwrap(), invs[{i}]);'
                             for i in range(4)))

print(L_code_test_invs)
