"""
sudoku_solver.py

Implement the function `solve_sudoku(grid: List[List[int]]) -> List[List[int]]` using a SAT solver from PySAT.
"""
import math
from pysat.formula import CNF
from pysat.solvers import Solver
from typing import List
def assign(r,c,d,N):
      return N*N*(r-1)+(c-1)*N+d # to assign a unique number to each varible in the  CNF
def decode(n,N):
     n=n-1
     d=(n%N)
     n=n//N
     c=(n%N)
     n=n//N
     r=(n%N)
     return r+1,c+1,d+1

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
     #each digit appears atmost once per row
    for r in range(1,N+1):
         for d in range(1,N+1):
              l=[]
              for c in range(1,N+1):
                   l.append(assign(r,c,d,N))
              for i in range(1,len(l)):
                   for j in range(i+1,len(l)):
                         cnf.append([-l[i],-l[j]])     

    #each digit appears atmost once in the column
    for c in range(1,N+1):
      for d in range(1,N+1):
           l=[]
           for r in range(1,N+1):
                l.append(assign(r,c,d,N))
           for i in range(1,len(l)):
                for j in range(i+1,len(l)):
                     cnf.append([-l[i],-l[j]])              
                                     
    #each digit appears once per a sub grid
    sq_rt = int(math.sqrt(N))
    for d in range(1,N+1):
         for r in range(0,N):
              for c in range(0,N):
                   i=(r//sq_rt)*sq_rt
                   j=(c//sq_rt)*sq_rt
                   l=[]
                   for r1 in range(i+1,i+sq_rt+1):
                        for c1 in range(j+1,j+sq_rt+1):
                             l.append(assign(r1,c1,d,N))
                   cnf.append(l)
                   for i in range(len(l)):
                        for j in range(i+1,len(l)):
                             cnf.append([-l[i],-l[j]])          
    for r in range(1,N+1):
         for c in range(1,N+1):
              if grid[r-1][c-1]!=0:
                   cnf.append([assign(r,c,grid[r-1][c-1],N)])
                   
    # TODO: implement encoding and solving using PySAT
    # Initialize a SAT solver 
    with Solver(name='glucose3') as solver:
        solver.append_formula(cnf.clauses)
        if solver.solve():
            model = solver.get_model()
          #   print("SAT solution:", model)
            for item in model:
               if item>0:
                  r,c,d=decode(item,N)
                  grid[r-1][c-1]=d
        else:
            print("UNSAT")
    return grid