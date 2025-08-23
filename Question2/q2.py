"""
Sokoban Solver using SAT (Boilerplate)
--------------------------------------
Instructions:
- Implement encoding of Sokoban into CNF.
- Use PySAT to solve the CNF and extract moves.
- Ensure constraints for player movement, box pushes, and goal conditions.

Grid Encoding:
- 'P' = Player
- 'B' = Box
- 'G' = Goal
- '#' = Wall
- '.' = Empty space
"""

from pysat.formula import CNF
from pysat.solvers import Solver

# Directions for movement
DIRS = {'U': (-1, 0), 'D': (1, 0), 'L': (0, -1), 'R': (0, 1)}


class SokobanEncoder:
    def __init__(self, grid, T):
        """
        Initialize encoder with grid and time limit.

        Args:
            grid (list[list[str]]): Sokoban grid.
            T (int): Max number of steps allowed.
        """
        self.grid = grid
        self.T = T
        self.N = len(grid)
        self.M = len(grid[0])

        self.goals = []
        self.boxes = []
        self.player_start = None

        # TODO: Parse grid to fill self.goals, self.boxes, self.player_start
        self._parse_grid()

        self.num_boxes = len(self.boxes)
        self.cnf = CNF()

    def _parse_grid(self):
        """Parse grid to find player, boxes, and goals."""
        # TODO: Implement parsing logic
        for i in range(len(self.grid)):
            for j in range(len(self.grid[0])):
                if self.grid[i][j]=='P':
                     self.player_start=(i,j)
                if self.grid[i][j]=='B':
                     self.boxes.append((i,j))
                if self.grid[i][j]=='G':
                     self.goals.append((i,j))           
        
    # ---------------- Variable Encoding ----------------
    def var_player(self, x, y, t):
        """
        Variable ID for player at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        return (0 << 20) | (x << 8) | (y << 4) | t+1


    def var_box(self, b, x, y, t):
        """
        Variable ID for box b at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        
        return (1 << 20) | (b << 12) | (x << 8) | (y << 4) | t+1

    # ---------------- Encoding Logic ----------------
    def encode(self):
        """
        Build CNF constraints for Sokoban:
        - Initial state
        - Valid moves (player + box pushes)
        - Non-overlapping boxes
        - Goal condition at final timestep
        """
        # TODO: Add constraints for:
        # 1. Initial conditions  
        self.cnf.append([self.var_player(self.player_start[0],self.player_start[1],0)])
        no_of_boxes=len(self.boxes)
        for i in range(no_of_boxes):
            self.cnf.append([self.var_box(i,self.boxes[i][0],self.boxes[i][1],0)])
        no_of_goals=len(self.goals)

        # 2. Player movement
        for i in range(0,self.N):
            for j in range(0,self.M):
                for t in range(0,self.T):
                    l=[]
                    if i==0 and j!=0 and j!=self.M-1:
                      p1=self.var_player(i+1,j,t+1)
                      l.append(p1)
                      p2=self.var_player(i,j+1,t+1)
                      l.append(p2)
                      p3=self.var_player(i,j-1,t+1)
                      l.append(p3)
                      l.append(-self.var_player(i, j, t))
                    elif i==0 and j==0:
                       p1=self.var_player(i+1,j,t+1)
                       p2=self.var_player(i,j+1,t+1)
                       l.append(p1)
                       l.append(p2)
                       l.append(-self.var_player(i, j, t))
                    elif i==0 and j==self.M-1:
                       p2=self.var_player(i+1,j,t+1)
                       p3=self.var_player(i,j-1,t+1)
                       l.append(p2)
                       l.append(p3)
                       l.append(-self.var_player(i, j, t))
                    elif i!=0 and j==0:
                       p1=self.var_player(i+1,j,t+1)
                       p2=self.var_player(i-1,j,t+1)
                       p3=self.var_player(i,j+1,t+1)
                       l.append(p1)
                       l.append(p2)
                       l.append(p3)
                       l.append(-self.var_player(i, j, t))
                    elif i!=0 and j!=0:
                       p1=self.var_player(i+1,j,t+1)
                       p2=self.var_player(i,j+1,t+1)
                       p3=self.var_player(i,j-1,t+1)
                       p4=self.var_player(i-1,j,t+1)     
                       l.append(p1)
                       l.append(p2)
                       l.append(p3)
                       l.append(p4)  
                       l.append(-self.var_player(i, j, t))
                    self.cnf.append(l)   

        # 3. Box movement (push rules)
        for b in range(no_of_boxes):
          for i in range(self.N):
              for j in range(self.M):
                  for t in range(self.T):
                      box_next = self.var_box(b, i, j, t+1)
                      box_curr = self.var_box(b, i, j, t)
                      
                      # Box stays in same position (no push)
                      no_push_clause = [-box_next, box_curr]
                      
                      # Add conditions that no push happened from any direction
                      # Push from left (player at (i,j-2) pushes box at (i,j-1) to (i,j))
                      if j > 0:
                          no_push_clause.append(-self.var_player(i, j-2, t))
                          no_push_clause.append(-self.var_box(b, i, j-1, t))
                      
                      # Push from right (player at (i,j+2) pushes box at (i,j+1) to (i,j))
                      if j < self.M-1:
                          no_push_clause.append(-self.var_player(i, j+2, t))
                          no_push_clause.append(-self.var_box(b, i, j+1, t))
                      
                      # Push from top (player at (i-2,j) pushes box at (i-1,j) to (i,j))
                      if i > 0:
                          no_push_clause.append(-self.var_player(i-2, j, t))
                          no_push_clause.append(-self.var_box(b, i-1, j, t))
                      
                      # Push from bottom (player at (i+2,j) pushes box at (i+1,j) to (i,j))
                      if i < self.N-1:
                          no_push_clause.append(-self.var_player(i+2, j, t))
                          no_push_clause.append(-self.var_box(b, i+1, j, t))
                      
                      self.cnf.append(no_push_clause)
                      
                      # Case 2: Box pushed from LEFT (player at (i,j-2) pushes box at (i,j-1) to (i,j))
                      if j > 0 and j < self.M-1:  # Can push from left and target is valid
                          push_left_clause = [
                              -box_next,
                              self.var_player(i, j-2, t),    # Player behind box
                              self.var_box(b, i, j-1, t),    # Box was to left
                              -self.var_box(b, i, j, t)      # Box wasn't already here
                          ]
                          # Ensure target position is not blocked by other boxes
                          for other_b in range(no_of_boxes):
                              if other_b != b:
                                  push_left_clause.append(-self.var_box(other_b, i, j, t))
                          self.cnf.append(push_left_clause)
                      
                      # Case 3: Box pushed from RIGHT (player at (i,j+2) pushes box at (i,j+1) to (i,j))
                      if j < self.M-1 and j > 0:  # Can push from right and target is valid
                          push_right_clause = [
                              -box_next,
                              self.var_player(i, j+2, t),    # Player behind box
                              self.var_box(b, i, j+1, t),    # Box was to right
                              -self.var_box(b, i, j, t)      # Box wasn't already here
                          ]
                          for other_b in range(no_of_boxes):
                              if other_b != b:
                                  push_right_clause.append(-self.var_box(other_b, i, j, t))
                          self.cnf.append(push_right_clause)
                      
                      # Case 4: Box pushed from TOP (player at (i-2,j) pushes box at (i-1,j) to (i,j))
                      if i > 0 and i < self.N-1:  # Can push from top and target is valid
                          push_top_clause = [
                              -box_next,
                              self.var_player(i-2, j, t),    # Player behind box
                              self.var_box(b, i-1, j, t),    # Box was above
                              -self.var_box(b, i, j, t)      # Box wasn't already here
                          ]
                          for other_b in range(no_of_boxes):
                              if other_b != b:
                                  push_top_clause.append(-self.var_box(other_b, i, j, t))
                          self.cnf.append(push_top_clause)
                      
                      # Case 5: Box pushed from BOTTOM (player at (i+2,j) pushes box at (i+1,j) to (i,j))
                      if i < self.N-1 and i > 0:  # Can push from bottom and target is valid
                          push_bottom_clause = [
                              -box_next,
                              self.var_player(i+2, j, t),    # Player behind box
                              self.var_box(b, i+1, j, t),    # Box was below
                              -self.var_box(b, i, j, t)      # Box wasn't already here
                          ]
                          for other_b in range(no_of_boxes):
                              if other_b != b:
                                  push_bottom_clause.append(-self.var_box(other_b, i, j, t))
                          self.cnf.append(push_bottom_clause)
                          
        # 4. Non-overlap constraints
        for i in range(0,self.N):
            for j in range(0,self.M):
                for t in range(0,self.T):
                    p=self.var_player(i,j,t)
                    for b in range(no_of_boxes):
                        b=self.var_box(b,i,j,t)
                        self.cnf.append([-p,-b])
                    for a in range(no_of_boxes):
                        for b in range(a+1,no_of_boxes):
                            a=self.var_box(a,i,j,t)
                            b=self.var_box(b,i,j,t)
                            self.cnf.append([-a,-b])   
        # 5. Goal conditions
        for (g, h) in self.goals:
            l = []
            for b in range(no_of_boxes):
                l.append(self.var_box(b,g,h, self.T - 1))
            self.cnf.append(l)   

        # 6. Other conditions
        for t in range(0,self.T):
            l=[]
            for i in range(0,self.N):
                for j in range(0,self.M):
                   l.append(self.var_player(i,j,t))
            self.cnf.append(l)

        for t in range(0,self.T):
            for i1 in range(0,self.N):
                for i2 in range(0,self.N):
                    for j1 in range(0,self.M):
                        for j2  in range(0,self.M):
                            if not (i1==i2 and j1==j2):
                               self.cnf.append([-self.var_player(i1,j1,t),-self.var_player(i2,j2,t)])
        for t in range(0, self.T):
           for b in range(no_of_boxes):   # loop over boxes
               l = []
               for i in range(0, self.N):
                   for j in range(0, self.M):
                       l.append(self.var_box(b, i, j, t))
               self.cnf.append(l)

        for t in range(0, self.T):
          for b in range(no_of_boxes):
              for i1 in range(0, self.N):
                  for j1 in range(0, self.M):
                      for i2 in range(0, self.N):
                          for j2 in range(0, self.M):
                              if not (i1 == i2 and j1 == j2):
                                  self.cnf.append([
                                      -self.var_box(b, i1, j1, t),
                                      -self.var_box(b, i2, j2, t)
                                  ])    
                                                 

def decode(model, encoder):
    """
    Decode SAT model into list of moves ('U', 'D', 'L', 'R').

    Args:
        model (list[int]): Satisfying assignment from SAT solver.
        encoder (SokobanEncoder): Encoder object with grid info.

    Returns:
        list[str]: Sequence of moves.
    """
    N, M, T = encoder.N, encoder.M, encoder.T
    # TODO: Map player positions at each timestep to movement directions
    l=set(v for v in model if v >0)
    pos={}
    for t in range(T+1):
      for i in range(N):
          for j in range(M):
              var = encoder.var_player(i, j, t)
              if var in l:
                  pos[t] = (i, j)
   
    moves=[]
    for t in range(1, T+1):
        if t not in pos or (t-1) not in pos:
            continue
        i1, j1 = pos[t-1]
        i2, j2 = pos[t]
        if i2 == i1-1 and j2 == j1:  # up
            moves.append("U")
        elif i2 == i1+1 and j2 == j1:  # down
            moves.append("D")
        elif i2 == i1 and j2 == j1-1:  # left
            moves.append("L")
        elif i2 == i1 and j2 == j1+1:  # right
            moves.append("R")
    return moves

def solve_sokoban(grid, T):
    """
    DO NOT MODIFY THIS FUNCTION.

    Solve Sokoban using SAT encoding.

    Args:
        grid (list[list[str]]): Sokoban grid.
        T (int): Max number of steps allowed.

    Returns:
        list[str] or "unsat": Move sequence or unsatisfiable.
    """
    encoder = SokobanEncoder(grid, T)
    cnf = encoder.encode()

    with Solver(name='g3') as solver:
        solver.append_formula(encoder.cnf)
        if not solver.solve():
            return []

        model = solver.get_model()
        if not model:
            return []

        return decode(model, encoder)


