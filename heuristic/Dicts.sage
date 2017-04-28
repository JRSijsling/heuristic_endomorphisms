"""
 *  Dictionary conversions
 *
 *  Copyright (C) 2016-2017
 *            Edgar Costa      (edgarcosta@math.dartmouth.edu)
 *            Davide Lombardo  (davide.lombardo@math.u-psud.fr)
 *            Jeroen Sijsling  (jeroen.sijsling@uni-ulm.de)
 *
 *  See LICENSE.txt for license details.
"""

_index_dict_ = dict()
# Magma indices for Lattice
_index_dict_['field'] = 1
_index_dict_['structure'] = 2
# Magma indices for field key
_index_dict_['seq'] = 1
_index_dict_['magma'] = 2
# Magma indices for OverField
_index_dict_['representation'] = 1
_index_dict_['algebra'] = 2
_index_dict_['description'] = 3
# Magma indices for representation key
_index_dict_['tangent'] = 1
_index_dict_['homology'] = 2
_index_dict_['approx'] = 3
# Magma indices for algebra key
_index_dict_['alg_QQ'] = 1
_index_dict_['alg_ZZ'] = 2
_index_dict_['alg_RR'] = 3
_index_dict_['alg_ST'] = 4
# Sage indices for description key
_index_dict_['factors_QQ'] = 0
_index_dict_['desc_ZZ'] = 1
_index_dict_['desc_RR'] = 2
_index_dict_['sato_tate'] = 3
# Sage indices for a factor_QQ
_index_dict_['albert_type'] = 0
_index_dict_['base_field'] = 1
_index_dict_['dim_sqrt'] = 2
_index_dict_['disc'] = 3
# Sage indices for desc_ZZ key
_index_dict_['index'] = 0
_index_dict_['is_eichler'] = 1
# Decomposition. Note the compatility of the first field with what went before.
_index_dict_['field'] = 1
_index_dict_['idem'] = 2
_index_dict_['factor'] = 3
_index_dict_['proj'] = 4

def sagify_description(desc_list):
    return eval(repr(magma.SagifyDescription(desc_list)))

def dict_lattice(lattice):
    dicts = [ ]
    for tup in lattice:
        dikt = dict()
        dikt['field'] = dict_field(tup[_index_dict_['field']])
        dikt['structure'] = dict_structure(tup[_index_dict_['structure']])
        dicts.append(dikt)
    return dicts

def desc_lattice(lattice):
    descs = [ ]
    for tup in lattice:
        desc = [ ]
        desc.append(sagify_description(tup[_index_dict_['field']][_index_dict_['seq']]))
        desc.append(sagify_description(tup[_index_dict_['structure']][_index_dict_['description']]))
        descs.append(desc)
    return descs

def dict_structure(structure):
    dikt = dict()
    dikt['representation'] = dict_rep(structure[_index_dict_['representation']])
    dikt['algebra'] = dict_alg(structure[_index_dict_['algebra']])
    desc = sagify_description(structure[_index_dict_['description']])
    dikt['description'] = dict_description(desc)
    return dikt

def desc_structure(structure):
    return sagify_description(structure[_index_dict_['description']])

def dict_field(field):
    dikt = dict()
    dikt['seq'] = field[_index_dict_['seq']]
    dikt['magma'] = field[_index_dict_['magma']]
    return dikt

def dict_rep(rep):
    dicts = [ ]
    for tup in rep:
        dikt = dict()
        dikt['tangent'] = tup[_index_dict_['tangent']]
        dikt['homology'] = tup[_index_dict_['homology']]
        dikt['approx'] = tup[_index_dict_['approx']]
        dicts.append(dikt)
    return dicts

def dict_alg(rep):
    dikt = dict()
    dikt['alg_QQ'] = rep[_index_dict_['alg_QQ']]
    dikt['alg_ZZ'] = rep[_index_dict_['alg_ZZ']]
    dikt['alg_RR'] = rep[_index_dict_['alg_RR']]
    dikt['alg_ST'] = rep[_index_dict_['alg_ST']]
    return dikt

def dict_description(desc):
    dikt = dict()
    dikt['factors_QQ'] = [ dict_factor_QQ(factor_QQ) for factor_QQ in desc[_index_dict_['factors_QQ']] ]
    dikt['desc_ZZ'] = dict_desc_ZZ(desc[_index_dict_['desc_ZZ']])
    dikt['desc_RR'] = dict_desc_RR(desc[_index_dict_['desc_RR']])
    dikt['sato_tate'] = desc[_index_dict_['sato_tate']]
    return dikt

def dict_factor_QQ(factor_QQ):
    dikt = dict()
    dikt['albert_type'] = factor_QQ[_index_dict_['albert_type']]
    dikt['base_field'] = factor_QQ[_index_dict_['base_field']]
    dikt['dim_sqrt'] = factor_QQ[_index_dict_['dim_sqrt']]
    dikt['disc'] = factor_QQ[_index_dict_['disc']]
    return dikt

def dict_desc_ZZ(desc_ZZ):
    dikt = dict()
    dikt['index'] = desc_ZZ[_index_dict_['index']]
    dikt['is_eichler'] = desc_ZZ[_index_dict_['is_eichler']]
    return dikt

def dict_desc_RR(desc_RR):
    dikt = dict()
    dikt['factors'] = desc_RR
    return dikt
