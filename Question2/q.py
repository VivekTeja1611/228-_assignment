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

        self.var_map = {} #why this is there
        self.var_count = 0  #why this is there
        self.no_of_boxes = len(self.boxes)
        self.cnf = CNF()
    def _var(self, kind, *args):
        """Creates or retrieves a unique integer variable for a given propositional atom."""
        key = (kind,) + args
        if key not in self.var_map:
            self.var_count += 1
            self.var_map[key] = self.var_count
        return self.var_map[key]
    
    def _is_free(self, i, j):
        """Checks if a cell is within bounds and not a wall."""
        return 0 <= i < self.N and 0 <= j < self.M and self.grid[i][j] != '#'    
    def _parse_grid(self):
        """Parse grid to find player, boxes, and goals."""
        # TODO: Implement parsing logic
        for r in range(self.N):
            for c in range(self.M):
                if self.grid[r][c] == 'P':
                    self.player_start = (r, c)
                elif self.grid[r][c] == 'B':
                    self.boxes.append((r, c))
                elif self.grid[r][c] == 'G':
                    self.goals.append((r, c))
        #note that in above one u didn't take the walls location  
    def var_player(self,x,y,t):
       """
       Variable ID for player at (x, y) at time t.
       """
       if not self._is_free(x, y):
           return None
       return self._var('P', x, y, t)
    
    def var_box(self,b,x,y,t):
     """
     Variable ID for box b at (x, y) at time t.
     """
     if not self._is_free(x, y):
        return None
     return self._var('B', b, x, y, t)

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
        #player start position
        for i in range(self.N):
            for j in range(self.M):
                if self.var_player(i,j,0):
                    if (i,j)==self.player_start:
                        self.cnf.append([self.var_player(i,j,0)])
                    else:    
                        self.cnf.append([-self.var_player(i,j,0)])
        #same goes for boxes at t=0                
        for b in range(self.no_of_boxes):
            for i in range(self.N):
                for j in range(self.M):
                    if self.var_box(b,i,j,0):
                        if(i,j)==self.boxes[b]:
                            self.cnf.append([self.var_box(b,i,j,0)])
                        else:
                            self.cnf.append([-self.var_box(b,i,j,0)])   
        
        # 2. Player movement      # 3. Box movement (push rules)
        for t in range(self.T):
            #at each state the player mush move to one of the adjacent cells
            for i in range(self.N):
                for j in range(self.M):
                    l=[]
                    if self._is_free(i, j):
                       for delta_i,delta_j in DIRS.values():
                           if self._is_free(i+delta_i, j+delta_j):
                               l.append(self.var_player(i+delta_i,j+delta_j,t+1))
                       self.cnf.append([-self.var_player(i,j,t)]+l)
            # a move can only happend if there are no wall and the pushed bpx has a empty cell to occupy
            for i in range(self.N):
                for j in range(self.M):
                    if self._is_free(i, j):
                        for di,dj in DIRS.values():
                            new_player_i,new_player_j=i+di,j+dj
                            if not self._is_free(new_player_i, new_player_j):
                                continue
                            box_new_i,box_new_j=new_player_i+di,new_player_j+dj
                            #there are towawyas to get blockd..either awall or another box
                            for b in range(self.no_of_boxes):
                                if not self._is_free(box_new_i, box_new_j):
                                    #blocked by wall
                                    self.cnf.append([-self.var_player(i,j,t),-self.var_box(b,new_player_i,new_player_j,t),-self.var_player(new_player_i,new_player_j,t+1)])
                                else:
                                    #blocked by box
                                    for b1 in range(self.no_of_boxes):
                                        self.cnf.append([-self.var_player(i,j,t),-self.var_box(b,new_player_i,new_player_j,t),-self.var_box(b1,box_new_i,box_new_j,t),-self.var_player(new_player_i,new_player_j,t+1)])
            for i in range(self.N):
               for j in range(self.M):
                   if not self._is_free(i, j): 
                       continue
                   for di,dj in DIRS.values():
                       new_i,new_j=i+di,j+dj
                       if not self._is_free(new_i, new_j):
                           continue
                       new_box_i,new_box_j=new_i+di,new_j+dj
                       if not self._is_free(new_box_i, new_box_j):
                           continue
                       for b in range(self.no_of_boxes):
                           p_t=self.var_player(i,j,t)
                           p_t1=self.var_player(new_i,new_j,t+1)
                           b_t=self.var_box(b,new_i,new_j,t)
                           b_t1=self.var_box(b,new_box_i,new_box_j,t+1)
                           if p_t and p_t1 and b_t and b_t1: 
                               self.cnf.append([-p_t,-p_t1,-b_t,b_t1])

            for b in range(self.no_of_boxes):
                for i in range(self.N):
                    for j in range(self.M):
                         if  self._is_free(i, j):
                             l=[]
                             for di,dj in DIRS.values():
                                 pi,pj=i-di,j-dj
                                 if self._is_free(pi, pj):
                                       v = self._var('push', b, i, j, t, pi, pj)
                                       p_prev=self.var_player(pi,pj,t)
                                       p_nex=self.var_player(i,j,t+1)
                                       self.cnf.append([-v,p_prev])
                                       self.cnf.append([-v,p_nex])
                                       self.cnf.append([v,-p_prev,-p_nex])
                                       l.append(v)
                             self.cnf.append([-self.var_box(b,i,j,t),self.var_box(b,i,j,t+1)]+l)      

        # 4. Non-overlap constraints
        for t in range(self.T+1):
            l=[]
            #atleast one postion
            for i in range(self.N):
                for j in range(self.M):
                    if not self._is_free(i, j):
                       continue
                    else:
                       if self.var_player(i,j,t):
                           l.append(self.var_player(i,j,t))
            self.cnf.append(l)               
            #atmost one position condition               
            for i in range(len(l)):
                for j in range(i+1,len(l)):
                    self.cnf.append([-l[i],-l[j]])
            #now similarly for boxes            
            for b in range(self.no_of_boxes):
                l=[]
                #atleast one postion
                for i in range(self.N):
                    for j in range(self.M):
                        if not self._is_free(i, j):
                           continue
                        else:
                           if self.var_box(b,i,j,t):
                               l.append(self.var_box(b,i,j,t))
                self.cnf.append(l)
                #atmost one position condition               
                for i in range(len(l)):
                    for j in range(i+1,len(l)):
                        self.cnf.append([-l[i],-l[j]])    
            for i in range(self.N):
                for j in range(self.M):
                     if not self._is_free(i, j):
                          continue
                     #overlap constraint between box and player
                     for b in range(self.no_of_boxes):
                          self.cnf.append([-self.var_player(i,j,t),-self.var_box(b,i,j,t)])
                     #overlap constraint between the boxes
                     for b1 in range(self.no_of_boxes):
                         for b2 in range(b1+1,self.no_of_boxes):
                                self.cnf.append([-self.var_box(b1,i,j,t),-self.var_box(b2,i,j,t)])
            
                                        
        # 5. Goal conditions
        for b in range(self.no_of_boxes):
            l=[]
            for (i,j) in self.goals:
                if  self._is_free(i, j):
                      l.append(self.var_box(b,i,j,self.T))
            self.cnf.append(l)

        # 6. Other conditions
        return self.cnf  
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
    if not model:
        return -1

    rev_map = {v: k for k, v in encoder.var_map.items()}
    player_pos = [None] * (encoder.T + 1)

    for var in model:
        if var > 0 and var in rev_map:
            kind, *args = rev_map[var]
            if kind == 'P':
                r, c, t = args
                player_pos[t] = (r, c)

    moves = []
    for t in range(encoder.T):
        if player_pos[t] is None or player_pos[t+1] is None:
            # Should not happen if the model is valid
            return -1 
        r1, c1 = player_pos[t]
        r2, c2 = player_pos[t+1]
        dr, dc = r2 - r1, c2 - c1
        
        move_found = False
        for move_char, (mr, mc) in DIRS.items():
            if (dr, dc) == (mr, mc):
                moves.append(move_char)
                move_found = True
                break
        if not move_found:
            # Player did not move to an adjacent square, invalid model
             return -1
    
    return "".join(moves)

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
        solver.append_formula(cnf)
        if not solver.solve():
            return -1

        model = solver.get_model()
        if not model:
            return -1

        return decode(model, encoder)