import hashlib
import math
import os
import csv
from q1 import solve_sudoku
from typing import List
from tqdm import tqdm
import random


def is_valid_sudoku(original:List[List[int]], grid: List[List[int]]) -> bool:
    
    n = len(grid)
    sqrt_n = int(math.sqrt(n))
    if sqrt_n * sqrt_n != n:
        raise ValueError("Grid size is not square with square subgrids.")

    def check_unique(values):
        nums = [v for v in values if v != 0]
        return len(set(nums)) == len(nums)

    for i in range(n):
        if not check_unique(grid[i]) or not check_unique([grid[j][i] for j in range(n)]):
            return False

    for br in range(0, n, sqrt_n):
        for bc in range(0, n, sqrt_n):
            block = [grid[r][c] for r in range(br, br + sqrt_n) for c in range(bc, bc + sqrt_n)]
            if not check_unique(block):
                return False
    return True


puzzles = []

with open('testcases', 'r') as f:
    lines = f.readlines()

# Filter valid lines of length 81
valid_lines = [line.strip() for line in lines if len(line.strip()) == 81]

# Randomly sample 500 lines without replacement
sampled_lines = random.sample(valid_lines, min(500, len(valid_lines)))

# Convert each line to a 9x9 grid
for line in sampled_lines:
    grid = [
        [int(c) if c.isdigit() else 0 for c in line[i*9:(i+1)*9]]
        for i in range(9)
    ]
    puzzles.append(grid)



passed = 0

for i, puzzle in enumerate(tqdm(puzzles, desc="Solving puzzles")):
    solved = solve_sudoku(puzzle)
    if is_valid_sudoku(puzzle, solved):
        passed += 1
    else:
        print(f"❌ Test case {i + 1} failed!")

print(f"\n✅ {passed}/{len(puzzles)} test cases passed.")

