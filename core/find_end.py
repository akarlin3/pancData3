import re
import sys

def check_file(filename):
    with open(filename, 'r', encoding='utf-8') as f:
        lines = f.readlines()
        
    stack = []
    
    # MATLAB keywords that start a block
    # Note: we need to handle cases where they are at the start of a line or after spaces
    # and not part of another word.
    block_starts = re.compile(r'^\s*(if|for|while|switch|try|parfor|function)\b')
    block_ends = re.compile(r'^\s*end\b')
    
    for i, line in enumerate(lines):
        line = line.split('%')[0] # remove comments
        if not line.strip():
            continue
            
        # check starts
        match = block_starts.search(line)
        if match:
            keyword = match.group(1)
            stack.append((keyword, i+1))
            
        # check ends
        if block_ends.search(line):
            if stack:
                stack.pop()
            else:
                print(f"Extra END at line {i+1}")
                return
                
    if stack:
        for keyword, line_num in stack:
            print(f"Unmatched {keyword} at line {line_num}")
    else:
        print("All matched!")

check_file('load_dwi_data.m')
