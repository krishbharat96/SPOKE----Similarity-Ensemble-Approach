#!/usr/bin/python2
from __future__ import print_function

import csv
import math
from decimal import Decimal as D

Func_types = ["AC50", "GI50", "LD50", "ED50", "ID50", "pD'2", "pD2", "pA2", "Log AC50", "Log GI50", "Log LD50", "-Log AC50", "-Log GI50", "-Log LD50"]
Func_data = ["identity", "identity", "identity", "identity", "identity", "minus_log", "minus_log", "minus_log", "plus_log", "plus_log", "plus_log", "minus_log", "minus_log", "minus_log"]
FUNCTIONAL_TYPES = zip(Func_types, Func_data)
FUNCTIONAL_TYPES_ACT = dict(FUNCTIONAL_TYPES)
BOUNDING_TYPES = {'IC60':'identity', 'IC70':'identity',
                  'IC80':'identity', 'IC90':'identity',
                  'IC95':'identity', 'IC99':'identity',
                  'Log IC60':'plus_log', 'Log IC70':'plus_log',
                  'Log IC80':'plus_log', 'Log IC90':'plus_log',
                  'Log IC95':'plus_log', 'Log IC99':'plus_log',
                  '-Log IC60':'minus_log', '-Log IC70':'minus_log',
                  '-Log IC80':'minus_log', '-Log IC90':'minus_log',
                  '-Log IC95':'minus_log', '-Log IC99':'minus_log'}
STD_TYPES = {'Ki':'identity', 'Kd':'identity',
             'IC50':'identity', 'pKi':'minus_log',
             'pKd':'minus_log', 'pIC50':'minus_log',
             '-Log Ki':'minus_log', '-Log Kd': 'minus_log',
             '-Log IC50': 'minus_log', '-Log KD': 'minus_log',
             'Log 1/Ki': 'minus_log', 'Log 1/Kd': 'minus_log',
             'Log 1/IC50': 'minus_log', 'log(1/Ki)': 'minus_log',
             'log(1/Kd)': 'minus_log', 'log(1/IC50)': 'minus_log',
             'Log Ki': 'plus_log', 'Log Kd': 'plus_log', 'Log IC50':
             'plus_log', 'logKi': 'plus_log', 'logKd': 'plus_log',
             'logIC50': 'plus_log', 'EC50': 'identity', 'pEC50':'minus_log',
             '-Log EC50': 'minus_log', 'Log 1/EC50': 'minus_log',
             'log(1/EC50)': 'minus_log', 'Log EC50': 'plus_log',
             'logEC50': 'plus_log'}

UNITS2NM = {'M': 1e9, 'mM': 1e6, 'uM': 1e3, 'nM':1, 'pM':1e-3, 'fM':1e-6, 'fmol/ml':1e-3, 'pmol/ml':1, 'nmol/ml': 1e3, 'umol/ml':1e9, 'mmol/ml':1e12, 'M-1':1e-9, 'NULL':1}
units_arr = UNITS2NM.keys()

def minus_log(x):
    return 10**(9-x)

def plus_log(x):
    if x > 2:
	x = -x
    return 10**(9+x)

def identity(x):
    return x

COMBINED_TYPE = {}
COMBINED_TYPE.update(FUNCTIONAL_TYPES_ACT)
COMBINED_TYPE.update(STD_TYPES)
COMBINED_TYPE.update(BOUNDING_TYPES)
test_arr = COMBINED_TYPE.keys()

sqlite_query = """
    SELECT td.chembl_id,
           md.chembl_id,
           act.standard_units,
           act.standard_value,
           act.standard_relation,
           act.standard_type,
           act.doc_id,
           assays.src_id,
           assays.confidence_score,
           assays.curated_by,
           td.target_type,
           td.organism
      FROM activities as act
      JOIN assays
            ON act.assay_id = assays.assay_id
      JOIN molecule_dictionary as md
            ON act.molregno = md.molregno
      JOIN target_dictionary as td
            ON assays.tid = td.tid
      WHERE not assays.curated_by = 'Autocuration'
            AND assays.confidence_score >= 4; 
"""

def get_assays_sqlite(ChEMBL_DB_Path):
    import sqlite3
    cdb = sqlite3.connect(ChEMBL_DB_Path)
    c = cdb.cursor()
    records = c.execute(ChEMBL_DB_Path)
    
    count = 0
    fin_dict = dict()
    for line in records:
        count = count + 1
        if (count%5000) == 0:
            print (str(count) + "/727442")
        if line[2].strip() in units_arr and line[11].strip() == 'Homo sapiens' and line[5].strip() in test_arr and (line[10].strip() == 'SINGLE PROTEIN' or line[10].strip() == 'PROTEIN COMPLEX'):
            if not line[6].strip() in ['7', '15']:
                targ_id = line[0].strip()
                cmpd_id = line[1].strip()
                concat_id = cmpd_id + "|" + targ_id

                test_type = line[5].strip()
                unit = line[2].strip()
                val = float(line[3])
                units_converted_val = val*float(UNITS2NM[unit])
                conversion_type = COMBINED_TYPE[test_type]
                
                final_val = 0
                if conversion_type == 'minus_log':
                    final_val = minus_log(units_converted_val)
                elif conversion_type == 'plus_log':
                    final_val = plus_log(units_converted_val)
                else:
                    final_val = identity(units_converted_val)

                value_in_M = final_val * 1e-9
                value_in_neglog = -(D(value_in_M).log10())
                
                bound_rel = line[4].strip()
                confidence = line[8].strip()
                curated_by = line[9].strip()

                fin_arr_elem = [str(value_in_neglog), str(bound_rel),
                                str(confidence), str(curated_by)]
                fin_arr_string = str(value_in_neglog) + "|" + str(bound_rel) + "|" + str(confidence) + "|" + str(curated_by) + "|" + str(test_type)

                if not concat_id in fin_dict.keys():
                    fin_dict.update({concat_id:[fin_arr_string]})
                else:
                    fin_dict[concat_id].append(fin_arr_string)

return fin_dict               

    
    

                            
