"""
sudoku_solver.py

Implement the function `solve_sudoku(grid: List[List[int]]) -> List[List[int]]` using a SAT solver from PySAT.
"""
import math
from pysat.formula import CNF
from pysat.solvers import Solver
from typing import List
def assign(r,c,d,N):
      return N*N*(r-1)+(c-1)*N+(d-1)+1 # to assign a unique number to each varible in the  CNF

def solve_sudoku(grid: List[List[int]]) -> List[List[int]]:
    """Solves a Sudoku puzzle using a SAT solver. Input is a 2D grid with 0s for blanks."""
    cnf=CNF()
    N=len(grid)
    #Each cell has at least one diigit
    for r in range(1,N+1):
         for c in range(1,N+1):
              l=[]
              for d in range(1,N+1):
                   l.append(assign(r,c,d,N))
              cnf.append(l)
    #each cell has atmost one digit
    for r in range(1,N+1):
       for c in range(1,N+1):
            for d in range(1,N+1):
                 for d1 in range(1,N+1):
                      if(d!=d1):
                          cnf.append([-assign(r,c,d,N),-assign(r,c,d1,N)])
    #each digit appears once per a row
    for r in range(1,N+1):
         for d in range(1,N+1):
              l=[]
              for c in range(1,N+1):
                   l.append(assign(r,c,d,N))
              cnf.append(l)              
    #each digit appears once per a column
    for c in range(1,N+1):
         for d in range(1,N+1):
              l=[]
              for r in range(1,N+1):
                   l.append(assign(r,c,d,N))
              cnf.append(l)
    #each digit appears once per a sub grid
    for d in range(1,N+1):
         l=[]
         for r in range(1,int(math.sqrt(N))+1):
              for c in range(1,int(math.sqrt(N))+1):
                   l.append(assign(r,c,d,N))
         cnf.append(l)   

    #fixed clues from the input
    for r in range(1,N+1):
         for c in range(1,N+1):
              if grid[r-1][c-1]!=0:
                   cnf.append([assign(r,c,grid[r-1][c-1],N)])
    # TODO: implement encoding and solving using PySAT
    
    return grid