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
_index_dict_['gens_ZZ'] = 2
_index_dict_['desc_RR'] = 3
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

def sagify_description(desc_list):
    return eval(repr(magma.SagifyDescription(desc_list)))

def dict_lattice(lattice):
    dicts = [ ]
    for tup in lattice:
        dict_to_fill = dict()
        dict_to_fill['field'] = dict_field(tup[_index_dict_['field']])
        dict_to_fill['structure'] = dict_structure(tup[_index_dict_['structure']])
        dicts.append(dict_to_fill)
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
    dict_to_fill = dict()
    dict_to_fill['representation'] = dict_rep(structure[_index_dict_['representation']])
    dict_to_fill['algebra'] = dict_alg(structure[_index_dict_['algebra']])
    desc = sagify_description(structure[_index_dict_['description']])
    dict_to_fill['description'] = dict_description(desc)
    return dict_to_fill

def desc_structure(structure):
    return sagify_description(structure[_index_dict_['description']])

def dict_field(field):
    dict_to_fill = dict()
    dict_to_fill['seq'] = field[_index_dict_['seq']]
    dict_to_fill['magma'] = field[_index_dict_['magma']]
    return dict_to_fill

def dict_rep(rep):
    dict_to_fill = dict()
    dict_to_fill['tangent'] = rep[_index_dict_['tangent']]
    dict_to_fill['homology'] = rep[_index_dict_['homology']]
    dict_to_fill['approx'] = rep[_index_dict_['approx']]
    return dict_to_fill

def dict_alg(rep):
    dict_to_fill = dict()
    dict_to_fill['alg_QQ'] = rep[_index_dict_['alg_QQ']]
    dict_to_fill['gens_ZZ'] = rep[_index_dict_['gens_ZZ']]
    dict_to_fill['desc_RR'] = rep[_index_dict_['desc_RR']]
    return dict_to_fill

def dict_description(desc):
    dict_to_fill = dict()
    dict_to_fill['factors_QQ'] = [ dict_factor_QQ(factor_QQ) for factor_QQ in desc[_index_dict_['factors_QQ']] ]
    dict_to_fill['desc_ZZ'] = dict_desc_ZZ(desc[_index_dict_['desc_ZZ']])
    dict_to_fill['desc_RR'] = dict_desc_RR(desc[_index_dict_['desc_RR']])
    dict_to_fill['sato_tate'] = desc[_index_dict_['sato_tate']]
    return dict_to_fill

def dict_factor_QQ(factor_QQ):
    dict_to_fill = dict()
    dict_to_fill['albert_type'] = factor_QQ[_index_dict_['albert_type']]
    dict_to_fill['base_field'] = factor_QQ[_index_dict_['base_field']]
    dict_to_fill['dim_sqrt'] = factor_QQ[_index_dict_['dim_sqrt']]
    dict_to_fill['disc'] = factor_QQ[_index_dict_['disc']]
    return dict_to_fill

def dict_desc_ZZ(desc_ZZ):
    dict_to_fill = dict()
    dict_to_fill['index'] = desc_ZZ[_index_dict_['index']]
    dict_to_fill['is_eichler'] = desc_ZZ[_index_dict_['is_eichler']]
    return dict_to_fill

def dict_desc_RR(desc_RR):
    dict_to_fill = dict()
    dict_to_fill['factors'] = desc_RR
    return dict_to_fill
