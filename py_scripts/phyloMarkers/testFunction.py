import argparse
import subprocess,os
import numpy as np
import regex as re
import pandas as pd
import os.path

def createSubtree(tree,paralogsFamily):
    with open(tree,"r") as f:
        tree = f.readline()
    f.close()
    treeBeginning = paralogsFamily[0]
    treeEnding = paralogsFamily[-1]
    treeEnding = treeEnding.replace("|", "\\|")
    subTree = ""
    regex_pattern = r'\((.*?)' + re.escape(treeBeginning)
    subTree = re.split(regex_pattern, tree)
    if len(subTree) > 1:
        subTree =  subTree[-1].strip()
    match = re.search(treeEnding, subTree)
    if match:
        substring_after_regex = subTree[match.end():]
        first_comma_position = substring_after_regex.find(',')
        first_semicolon_position = substring_after_regex.find(';')
        first_separator_position = min(first_comma_position, first_semicolon_position)
        if first_separator_position != -1:
            text_before_separator = subTree[:match.end() + first_separator_position]
            subTree = text_before_separator.strip()
    subTree = treeBeginning + subTree
    subTree = subTree + ')'
    subTreeCount1 = subTree.count('(')
    subTreeCount2 = subTree.count(')')
    diff = max(subTreeCount2,subTreeCount1) - min(subTreeCount2,subTreeCount1)
    if subTreeCount1 > subTreeCount2:
        subTree = subTree + ')' * diff
    else:
        subTree = '(' * diff + subTree
    subTree = subTree + ';'
    """
    with open(pathOutput+"/subTree.tree","w") as f:
        f.write(subTree)
    f.close()
    """
    print(subTree)
    return


tree = "output/synaptotagminTest/HMM/C2:C2/MSA/Tree.tree"
paralogsFamily = ['c_hem|XLOC_007395', 'm_vir|scaffold204.g30.t1', 'ast|STRG.26091'] 
createSubtree(tree,paralogsFamily)
