"""Identify optimal order of XFST replace rules (mal-rules) to avoid reversions.
"""
import argparse
from collections import defaultdict
from itertools import combinations
from pprint import pprint
import re
import subprocess
from typing import DefaultDict

from matplotlib import pyplot as plt
import networkx as nx

# INTERESTING CYCLES TO INVESTIGATE
# ['prijti', 'j2i', 'SRc', 'i2j']

parser = argparse.ArgumentParser(description='Identify optimal ordering of '
                                             'XFST replace rules (mal-rules) '
                                             'to avoid reversions.')
parser.add_argument('files', type=str, nargs='*', default=[],
                    help='XFST regex file(s) containing replace rules.')
# parser.add_argument('--sum', dest='accumulate', action='store_const',
#                     const=sum, default=max,
#                     help='sum the integers (default: find the max)')

args = parser.parse_args()

ex_changes = {'V': [('a', 'b')],
              'W': [('b', 'a')],
              'X': [('b', 'c')],
              'Y': [('c', 'd')],
              'Z': [('d', 'a')]
              }

def affix_indices(fnames):
    shared_fname_prefix = 0
    while True:
        if len({fname[shared_fname_prefix] for fname in fnames}) == 1:
            shared_fname_prefix += 1
        else:
            break
    shared_fname_suffix = 0
    while True:
        if len({fname[shared_fname_suffix-1] for fname in fnames}) == 1:
            shared_fname_suffix -= 1
        else:
            if shared_fname_suffix == 0:
                shared_fname_suffix = None
            break
    return shared_fname_prefix, shared_fname_suffix

def foma_dot_graph_from_files(fnames) -> nx.DiGraph:
    """Use foma to generate a dot file from the composition
    of all the rules. Return networkx DiGraph of dot file.
    """
    shared_fname_prefix, shared_fname_suffix = affix_indices(fnames)
    foma_src = ''
    for i, fname in enumerate(fnames):
        display_name = fname[shared_fname_prefix:shared_fname_suffix]
        # foma_src += f'# {display_name}\n'
        if i == 0:
            foma_src += 'regex '
        cleaned_lines = []
        with open(fname) as f:
            for line in f:
                if line.startswith('#') or line.strip() == '':
                    continue
                else:
                    cleaned = line.strip().rstrip(' \n\t;')
                    if cleaned:
                        cleaned_lines.append(cleaned)
        cleaned_file = ' '.join(cleaned_lines)
        if not cleaned_file.startswith('[') or not cleaned_file.endswith(']'):
            cleaned_file = f'[ {cleaned_file} ]'
        cleaned_file = f'{cleaned_file} .o.\n'
        foma_src += cleaned_file
    assert ';' not in foma_src
    foma_src = foma_src.rstrip('.o\n') + ';\nprint dot\nview net\n'
    # foma_src = foma_src.replace('.o. ', '.o.\n\t')
    print(foma_src)
    escaped_quote = '\\"'
    dot = subprocess.run('foma', input=foma_src,
                         text=True, capture_output=True, shell=True)
    dot = re.search(r'(digraph.*)', dot.stdout, re.S).group(1).strip()
    print(dot)


def changes_from_files(fnames) -> dict:
    shared_fname_prefix, shared_fname_suffix = affix_indices(fnames)
    change_dict = {}
    for fname in fnames:
        display_name = fname[shared_fname_prefix:shared_fname_suffix]
        with open(fname) as f:
            change_dict[display_name] = re.findall(r'(\w+)\s*\(<-\)\s*(\w+)', f.read())
    return change_dict


def src_and_tgt_symbols(change_dict):
    src_symbols = DefaultDict(list)
    tgt_symbols = DefaultDict(list)
    for rule_name, changes in change_dict.items():
        for tgt, src in changes:
            src_symbols[src].append(rule_name)
            tgt_symbols[tgt].append(rule_name)
    return src_symbols, tgt_symbols


def parallelize_rules(rule_sets, u, v):
    print('input rule_sets', rule_sets)
    print('input nodes:', u, v)
    output_sets = []
    new_set = set()
    for s in rule_sets:
        print('current_set:', s)
        print('new_set:', new_set)
        if u in s or v in s:
            new_set.add(u)
            new_set.add(v)
            for w in s:
                new_set.add(w)
        else:
            output_sets.append(s)
    if new_set:
        output_sets.append(new_set)
    else:
        output_sets.append({u, v})
    print('=' * 79)
    return output_sets


def get_pos_from_dot(dot_str):
    """Get y-coordinates from from dot-format string."""
    dot_str = re.sub(r'digraph {(.*)}', r'\1', dot_str, re.S)
    entries = [e for e in dot_str.split(';')
               if 'change="' not in e]
    y_coords = []
    for e in entries:
        match = re.search(r'(\w+)\s+\[.*?pos="[0-9.]+,([0-9.]+)', e, re.S)
        if match:
            y_coords.append(match.groups())
        else:
            print('no match:', repr(e))
    return y_coords
    

def rule_ordering(rule_sets, node_pos, maximal_feeding=True):
    y_coords = sorted(node_pos, key=lambda x: x[1], reverse=maximal_feeding)
    already_added = set()
    ordering = []
    for node, y in node_pos:
        if node in already_added:
            continue
        i = None
        for j, s in enumerate(rule_sets):
            if node in s:
                i = j
                break
        if i is not None:  # if node is in a group
            group = rule_sets.pop(i)
            for n in group:
                already_added.add(n)
            ordering.append(group)
        else:
            already_added.add(node)
            ordering.append(node)
    return ordering


if __name__ == '__main__':
    rule_sets = []
    if args.files:
        change_dict = changes_from_files(args.files)
        src_symbols, tgt_symbols = src_and_tgt_symbols(change_dict)
    else:
        change_dict = ex_changes
        src_symbols, tgt_symbols = src_and_tgt_symbols(change_dict)
        
    g = nx.MultiDiGraph()
    removed_g = nx.MultiDiGraph()
    for rule_name, changes in change_dict.items():
        for tgt, src in changes:
            for bleed_name in src_symbols.get(tgt, []):
                g.add_edge(rule_name, bleed_name, change=f'{src}->{tgt}')

    while remaining_cycles := sorted((c for c in nx.simple_cycles(g) if len(c) > 1), key=len):
        pprint(nx.to_dict_of_dicts(g))
        cycles_len = len(remaining_cycles)
        if not cycles_len:
            break
        else:
            print(g, ':', cycles_len, 'cycles.')
            if cycles_len < 10:
                print('Less than 10 cycles remaining:', remaining_cycles)
            min_cycle_len = min(len(c) for c in remaining_cycles)
            min_len_cycles = [c for c in remaining_cycles
                              if len(c) == min_cycle_len]
            print('Cycles of length', min_cycle_len, ':', min_len_cycles)
            response = input('Remove edges (Y/n)? ')
            if response not in ('', 'Y', 'y', 'yes', 'Yes', 'YES'):
                break
            else:
                for c in min_len_cycles:
                    if min_cycle_len <= 2:
                        edges = tuple(e for e in zip(c, c[1:]))
                    else:
                        edges = tuple(e for e in zip(c, c[1:] + [c[0]]))
                    for i, (u, v) in enumerate(edges):
                        print(f'({i}) {u}->{v}')
                    print('Current parallel rules:', rule_sets)
                    response = input('Which edges should be removed (space-separated) or (Q)uit ([enter] for all)? ')
                    if response in ('Q', 'q', 'Quit', 'quit', 'QUIT'):
                        break
                    elif response in ('', 'y', 'yes', 'Yes', 'Y', 'YES'):
                        response = [i for i in range(len(edges))]
                    else:
                        response = [int(i) for i in response.split(' ')]
                    for i, (u, v) in enumerate(edges):
                        if i in response:
                            print(f'removing the edge ({u}, {v})...')
                            while g.has_edge(u, v):
                                g.remove_edge(u, v)
                            while g.has_edge(v, u):
                                g.remove_edge(v, u)
                            rule_sets = parallelize_rules(rule_sets, u, v)
    nx.draw_networkx(g, pos=nx.spring_layout(g))
    # plt.show()

    # nx.drawing.nx_pydot.write_dot(g, 'network.dot')
    g_pydot = nx.drawing.nx_pydot.to_pydot(g)
    g_pydot.write_png('network.png')
    gv_dot_str = g_pydot.create_dot().decode('utf-8')
    print('dot str:\n', gv_dot_str)
    node_pos = get_pos_from_dot(gv_dot_str)
    print('node_pos:\n', node_pos)
    print('rule_sets:\n', rule_sets)
    ordering = rule_ordering(rule_sets, node_pos)
    print('ordering:\n', ordering)




    # print(gv_dot_bytes)
    # with open('/tmp/gv_dot_string.dot', 'wb') as f:
    #     f.write(gv_dot_bytes)
    # g_from_gv = nx.drawing.nx_pydot.read_dot('/tmp/gv_dot_string.dot')
    # print(sorted(g_from_gv.nodes(), key=lambda x: x['pos'][]))
    # print(nx.to_dict_of_dicts(g_from_gv))
