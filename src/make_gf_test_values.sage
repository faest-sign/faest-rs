from sage.rings.finite_rings.hom_finite_field import FiniteFieldHomomorphism_generic
import itertools as it

F = GF(2)
R.<X> = F[]
p = X^128 + X^7 + X^2 + X + 1
L = GF(2^128, 'z', p)
L.inject_variables()

q8 = X^8 + X^4 + X^3 + X + 1
K8 = GF(2^8, 'y8', q8)
K8.inject_variables()
g8 = y8 + 1

q7 = X^7 + X + 1
K7 = GF(2^7, 'y7', q7)
K7.inject_variables()
g7 = y7

q9 = X^9 + X^4 + 1
K9 = GF(2^9, 'y9', q9)
K9.inject_variables()
g9 = y9

q10 = X^10 + X^3 + 1
K10 = GF(2^10, 'y10', q10)
K10.inject_variables()
g10 = y10

q11 = X^11 + X^2 + 1
K11 = GF(2^11, 'y11', q11)
K11.inject_variables()
g11 = y11


H = Hom(K8, L)
f = FiniteFieldHomomorphism_generic(H)

K7_values =  [
    K7.from_integer(0x7f), K7.from_integer(0x65), K7.from_integer(0x59), K7.from_integer(0x71)
]
K8_values =  [
    K8.from_integer(0x4b), K8.from_integer(0x27), K8.from_integer(0xb1), K8.from_integer(0xfc),
]
K9_values =  [
    K9.from_integer(0x12a), K9.from_integer(0x1f4), K9.from_integer(0x1a7), K9.from_integer(0x02b)
]
K10_values =  [
    K10.from_integer(0x0fd), K10.from_integer(0x235), K10.from_integer(0x3cb), K10.from_integer(0x0c7)
]
K11_values =  [
    K11.from_integer(0x1fb), K11.from_integer(0x113), K11.from_integer(0x112), K11.from_integer(0x544)
]

def make_tables(name, K, g, values):
    if len(K) <= 256:
        bitlen = 8
        hexlen = 2
        basic_type = 'u8'
    else:
        bitlen = 16
        hexlen = 3
        basic_type = 'u16'

    name_caps = name.upper()

    log_table = [0] + [
        log(K.from_integer(i), g) for i in range(1, len(K))
    ]
    anti_log_table = [g^i for i in range(len(K) - 1)]

    code_tables = f'''
        const LOG_TABLE: [{basic_type}; {len(K)}] = [
            {', '.join(f'0x{v:02x}' for v in log_table)}
        ];
        const ANTI_LOG_TABLE: [Self; {len(K)-1}] = [
            {', '.join(f'Self(0x{v.to_integer():02x})' for v in anti_log_table)}
        ];
    '''

    print(code_tables)

    sums = [x + y for x, y in it.product(values, repeat=2)]
    prods = [x * y for x, y in it.product(values, repeat=2)]
    invs = [x^-1 for x in values]

    code_values = f'''
        const {name_caps}_VALUES: [{name}; 4] = [
            {(',' + chr(0xa) + '            ').join(f'{name}(0x{v.to_integer():0{hexlen}x})' for v in values)}
        ];
        const {name_caps}_SUMS: [{name}; 16] = [
            {(',' + chr(0xa) + '            ').join(f'{name}(0x{v.to_integer():0{hexlen}x})' for v in sums)}
        ];
        const {name_caps}_PRODS: [{name}; 16] = [
            {(',' + chr(0xa) + '            ').join(f'{name}(0x{v.to_integer():0{hexlen}x})' for v in prods)}
        ];
        const {name_caps}_INVS: [{name}; 4] = [
            {(',' + chr(0xa) + '            ').join(f'{name}(0x{v.to_integer():0{hexlen}x})' for v in invs)}
        ];
    '''

    print(code_values)


def make_tests(name, K):
    if len(K) <= 256:
        bitlen = 8
        hexlen = 2
        basic_type = 'u8'
    else:
        bitlen = 16
        hexlen = 3
        basic_type = 'u16'

    name_caps = name.upper()
    table_var_prefix = '$tab'

    code_test_sums = f'''
            let values = {table_var_prefix}_values;
            let sums = {table_var_prefix}_sums;
            {(chr(0xa) + '            ').join(f'assert_eq!(values[{i}] + values[{j}], sums[{4 * i + j}]);'
                                 for i, j in it.product(range(4), repeat=2))}
    '''

    print(code_test_sums)

    code_test_prods = f'''
            let values = {table_var_prefix}_values;
            let prods = {table_var_prefix}_prods;
            {(chr(0xa) + '            ').join(f'assert_eq!(values[{i}] * values[{j}], prods[{4 * i + j}]);'
                                 for i, j in it.product(range(4), repeat=2))}
    '''

    print(code_test_prods)

    code_test_invs = f'''
            let values = {table_var_prefix}_values;
            let invs = {table_var_prefix}_invs;
            {(chr(0xa) + '            ').join(f'assert_eq!(values[{i}].invert().unwrap(), invs[{i}]);'
                                 for i in range(4))}
    '''

    print(code_test_invs)

#  make_tables('GF2p7', K7, g7, K7_values)
#  make_tables('GF2p8', K8, g8, K8_values)
make_tables('GF2p9', K9, g9, K9_values)
make_tables('GF2p10', K10, g10, K10_values)
make_tables('GF2p11', K11, g11, K11_values)

#  make_tests('GF2p8', K8)

#  L_values =  [
#         L.from_integer(0x616ce329b8aee6b6c752890eaec5fdca),
#         L.from_integer(0xc9d404e824cf265e31f38f85087a788f),
#         L.from_integer(0x7eceb30086f94c398a6881ea95c275b2),
#         L.from_integer(0x41a376dbe3eebf96a8f49dff52f12e1c),
#     ]
#
#  L_sums = [x + y for x, y in it.product(values, repeat=2)]
#  L_prods = [x * y for x, y in it.product(values, repeat=2)]
#  L_invs = [x^-1 for x in values]
#  K_in_L_embeddings = [f(v) for v in K_values]
#
#  code_embedding_constants = '''
#      pub const GF2P8_EMBEDDING_POX: [Self; 8] = [
#          Self::ONE,
#          {}
#      ];
#  '''.format(',\n        '.join(f'Self::from_u128(0x{(f(y^i)).to_integer():032x})' for i in range(1, 8)))
#
#  print(code_embedding_constants)
#
#  L_code_values = '''
#      const GF2P128_VALUES: [u128; 4] = [
#          {}
#      ];
#      const GF2P128_SUMS: [u128; 16] = [
#          {}
#      ];
#      const GF2P128_PRODS: [u128; 16] = [
#          {}
#      ];
#      const GF2P128_INVS: [u128; 4] = [
#          {}
#      ];
#      const GF2P8_IN_GF2P128_EMBEDDINGS: [u128; 4] = [
#          {}
#      ];
#  '''.format(
#          ',\n        '.join(f'0x{v.to_integer():032x}' for v in L_values),
#          ',\n        '.join(f'0x{v.to_integer():032x}' for v in L_sums),
#          ',\n        '.join(f'0x{v.to_integer():032x}' for v in L_prods),
#          ',\n        '.join(f'0x{v.to_integer():032x}' for v in L_invs),
#          ',\n        '.join(f'0x{v.to_integer():032x}' for v in K_in_L_embeddings),
#             )
#
#  print(L_code_values)
#
#  L_code_test_sums = '''
#          let values = GF2P128_VALUES.map(F::from_u128);
#          let sums = GF2P128_SUMS.map(F::from_u128);
#          {}
#  '''.format('\n        '.join(f'assert_eq!(values[{i}] + values[{j}], sums[{4 * i + j}]);'
#                               for i, j in it.product(range(4), repeat=2)))
#
#  print(L_code_test_sums)
#
#  L_code_test_prods = '''
#          let values = GF2P128_VALUES.map(F::from_u128);
#          let prods = GF2P128_PRODS.map(F::from_u128);
#          {}
#  '''.format('\n        '.join(f'assert_eq!(values[{i}] * values[{j}], prods[{4 * i + j}]);'
#                               for i, j in it.product(range(4), repeat=2)))
#
#  print(L_code_test_prods)
#
#  L_code_test_invs = '''
#          let values = GF2P128_VALUES.map(F::from_u128);
#          let invs = GF2P128_INVS.map(F::from_u128);
#          {}
#  '''.format('\n        '.join(f'assert_eq!(values[{i}].invert().unwrap(), invs[{i}]);'
#                               for i in range(4)))
#
#  print(L_code_test_invs)
