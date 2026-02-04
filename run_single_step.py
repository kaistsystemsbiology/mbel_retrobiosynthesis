import os
import argparse
import subprocess
from glob import glob

import numpy as np
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem

rulebook_dir = '/your_directory/mbel_retrobiosynthesis/reaction_rules/RuleBook.txt'
rule_info_dir = '/your_directory/mbel_retrobiosynthesis/reaction_rules/Rule_Information.txt'


def argument_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_inchi', required=True,  help='InChI of target product for retrobiosynthesis')
    parser.add_argument('-o', '--output_dir', required=True, help='Output directory')
    parser.add_argument('-n', '--step_num', required=False, type=int, default=1, help='Output directory')
    return parser


def run_reaction_inchi(target_inchi, rules, inchi=True):
    if inchi:
        target_mol = Chem.MolFromInchi(target_inchi)
    else:
        target_mol = Chem.MolFromSmiles(target_inchi)
    target_mol = Chem.AddHs(target_mol)
    
    results = {}
    for idx, rule in rules.items():
        try:
            ps = rule.RunReactants((target_mol,))
        except:
            1
        else:
            if len(ps) == 0:
                continue
            else:
                results[idx] = [item[0] for item in ps]
    return results


def run_cycle(target_inchi, rules, inchi):
    
    results = run_reaction_inchi(target_inchi, rules=rules, inchi=inchi)
    sub2idx = {}
    for idx, lst in results.items():
        for subs in lst:
            tmp_inchi = []
            smiles_long = Chem.MolToSmiles(subs)
            smiles_long = smiles_long.split('.')
            smiles_long.sort()
            tmp_inchi = [Chem.MolFromSmiles(item) for item in smiles_long]
            tmp_inchi = [Chem.MolToInchi(item) for item in tmp_inchi]
            tmp_inchi.sort()
            tmp_inchi = ';'.join(tmp_inchi)
            if tmp_inchi not in sub2idx:
                sub2idx[tmp_inchi] = set()
            sub2idx[tmp_inchi].add(idx)
    
    # DeepRFC
    
    intermediates = dict()
    for key, val in sub2idx.items():
        for inchi in key.split(';'):
            if inchi not in intermediates:
                intermediates[inchi] = set()
            intermediates[inchi].update(val)
    
    return results, sub2idx, intermediates



if __name__ == '__main__':
    '''
    ! python run_single_step.py -i 'InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6?/m1/s1' -o './output' -n 1
    '''
    
    parser = argument_parser()
    options = parser.parse_args()
    target_inchi = options.input_inchi
    output_dir = options.output_dir
    step_num = options.step_num
    
    if output_dir[-1] == '/':
        output_dir = output_dir[:-1]
        
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
        
    rule_info = pd.read_csv(rule_info_dir, sep='\t', index_col='IDX')
    df_rules = pd.read_csv(rulebook_dir, sep='\t', index_col='IDX')
    rules = {}
    for i, row in df_rules.iterrows():
        rules[i] = AllChem.ReactionFromSmarts(row['SMARTS'])
        
    print(f'Number of rules: {len(rules)}')
    
    inchi = False
    if 'InChI=' in target_inchi:
        target_smi = Chem.MolToSmiles(Chem.MolFromInchi(target_inchi))
        inchi = True
    else:
        target_smi = target_inchi
    
    reaction_results, sub2idx, intermediates = run_cycle(target_inchi, rules=rules, inchi=inchi)
    
    print(f'Number of matched rules: {len(reaction_results)}')
    print(f'Number of predicted reactions: {len(sub2idx)}')
    print(f'Number of predicted substrates: {len(intermediates)}')
    
    
    idx2inchi = {}
    for i, inchi in enumerate(intermediates):
        idx2inchi[i] = inchi
    
    with open(f'{output_dir}/tmp_deeprfc_input.txt', 'w') as fp:
        fp.write('ID	Reactant	Product\n')
        cnt = 0 
        for cnt, inchi in idx2inchi.items():
            try:
                tmp_mol = Chem.MolFromInchi(inchi)
            except:
                print(inchi)
                Chem.MolFromInchi(inchi)
            else:
                smi = Chem.MolToSmiles(tmp_mol)
            fp.write(f'{cnt}\t{smi}\t{target_smi}\n')
            
            
    subprocess.call(
        f"python deeprfc.py -i '{output_dir}/tmp_deeprfc_input.txt' -o {output_dir}/rfc_result",
        shell=True,
        stderr=subprocess.STDOUT
    )
    
    output_rfc = pd.read_csv(f'{output_dir}/rfc_result/result.txt', sep='\t')
    
    intermediates_feasible = {}
    feasible_inchis = set()
    for idx in output_rfc[output_rfc['Feasibility']=='feasible'][:]['ID']:
        inchi = idx2inchi[idx]
        intermediates_feasible[inchi] = intermediates[inchi]
        feasible_inchis.add(inchi)
        
        
    print(f'Number of feasible, predicted substrates: {len(feasible_inchis)}')
        

    sub2idx_processed = {}
    for key, val in sub2idx.items():
        substrates = key.split(';')
        if len(set(substrates) - feasible_inchis) > 0:
            continue
        sub2idx_processed[key] = val
        
    print(f'Number of feasible, predicted reactions: {len(sub2idx_processed)}')

    tmp_save = []
    for subs, rule_idxs in sub2idx_processed.items():
        rule_idxs = list(rule_idxs)
        rule_idxs.sort()
        rule_idxs = ';'.join([str(i) for i in rule_idxs])
        tmp_save.append([subs, rule_idxs])
        
    tmp_save = pd.DataFrame(tmp_save)
    
    
    result_equation = []
    for i, row in tmp_save.iterrows():
        subs, idxs = row
        smart_idxs = idxs.split(';')
        smart_idxs = [int(i) for i in smart_idxs]

        dict_coeff = {}
        for sub in subs.split(';'):
            dict_coeff[sub] = subs.count(sub)
        eqn_sub = []
        for sub, coeff in dict_coeff.items():
            eqn_sub.append(f'{float(coeff)} {sub}')
        eqn_sub = ' + '.join(eqn_sub)

        for idx in smart_idxs:
            df_tmp = rule_info[rule_info['SMARTS_IDX']==idx]
            for j, sub_row in df_tmp.iterrows():
                rule_id, legacy, diameter, smarts, score, ec, eqn, _ = sub_row
                for eqn in df_tmp['Equation']:
                    new_eqn = eqn + f' + 1.0 {target_inchi}'
                    new_eqn = eqn_sub + ' + ' + new_eqn

                    result_equation.append([subs, idx, new_eqn, rule_id, legacy, diameter, smarts, score, ec])

    df_equation = pd.DataFrame(result_equation)
    df_equation.columns = ['Substrates', 'SMARTS_IDX', 'Equation', 'Rule_ID', 'Legacy', 'Diameter', 'SMARTS', 'Score' ,'EC_number']
    df_equation.to_csv(f'{output_dir}/Results.txt', sep='\t')