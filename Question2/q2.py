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

        self.var_map = {}
        self.var_count = 0
        self.num_boxes = len(self.boxes)
        self.cnf = CNF()

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
                    
    
    def _var(self, kind, *args):
        """Creates or retrieves a unique integer variable for a given propositional atom."""
        key = (kind,) + args
        if key not in self.var_map:
            self.var_count += 1
            self.var_map[key] = self.var_count
        return self.var_map[key]
    
    def _is_free(self, r, c):
        """Checks if a cell is within bounds and not a wall."""
        return 0 <= r < self.N and 0 <= c < self.M and self.grid[r][c] != '#'
    # ---------------- Variable Encoding ----------------
    def var_player(self, r, c, t):
        """
         Variable ID for player P at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        if not self._is_free(r, c):
            return None
        return self._var('P', r, c, t)

    def var_box(self, b, r, c, t):
        """
         Variable ID for box b at (x, y) at time t.
        """
        # TODO: Implement encoding scheme
        if not self._is_free(r, c):
            return None
        return self._var('B', b, r, c, t)

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
        for r in range(self.N):
            for c in range(self.M):
                v = self.var_player(r, c, 0)
                if v:
                    self.cnf.append([v] if (r, c) == self.player_start else [-v])
        
        for b in range(self.num_boxes):
            for r in range(self.N):
                for c in range(self.M):
                    v = self.var_box(b, r, c, 0)
                    if v:
                        self.cnf.append([v] if (r, c) == self.boxes[b] else [-v])

        # 2. player movement and box movement (push rules)
        for t in range(self.T):
            # Player must move to an adjacent cell.
            for r, c in [(r, c) for r in range(self.N) for c in range(self.M) if self._is_free(r,c)]:
                p_t = self.var_player(r, c, t)
                possible_moves = [self.var_player(r + dr, c + dc, t + 1) for dr, dc in DIRS.values() if self._is_free(r + dr, c + dc)]
                self.cnf.append([-p_t] + possible_moves)

            # A move is only possible if there are no walls and inside the grid.
            for r, c in [(r, c) for r in range(self.N) for c in range(self.M) if self._is_free(r,c)]:
                p_t = self.var_player(r, c, t)
                for dr, dc in DIRS.values():
                    nr, nc = r + dr, c + dc #player position after push
                    if not self._is_free(nr, nc): continue
                    p_t1 = self.var_player(nr, nc, t + 1)
                    
                    # box postion after a push
                    br, bc = nr + dr, nc + dc
                    
                    # A move is blocked if it's a push, but the destination is a wall or another box.
                    is_blocked_by_wall = not self._is_free(br, bc)
                    for b1_idx in range(self.num_boxes):
                        b_in_way = self.var_box(b1_idx, nr, nc, t)
                        if is_blocked_by_wall:
                            # P(r,c,t) ^ B(b1,nr,nc,t) => ~P(nr,nc,t+1)  (blocked by wall)
                            self.cnf.append([-p_t, -b_in_way, -p_t1])
                        else:
                            for b2_idx in range(self.num_boxes):
                                b_at_dest = self.var_box(b2_idx, br, bc, t)
                                # P(r,c,t) ^ B(b1,nr,nc,t) ^ B(b2,br,bc,t) => ~P(nr,nc,t+1) (blocked by another box)
                                self.cnf.append([-p_t, -b_in_way, -b_at_dest, -p_t1])

            # Action Consequences: Player's move causes boxes to move.
            for r, c in [(r, c) for r in range(self.N) for c in range(self.M) if self._is_free(r,c)]:
                p_t = self.var_player(r, c, t)
                for dr, dc in DIRS.values():
                    nr, nc = r + dr, c + dc
                    if not self._is_free(nr, nc): continue
                    p_t1 = self.var_player(nr, nc, t + 1)
                    
                    br, bc = nr + dr, nc + dc
                    if not self._is_free(br, bc): continue
                    
                    for b_idx in range(self.num_boxes):
                        b_t = self.var_box(b_idx, nr, nc, t)
                        b_t1_new = self.var_box(b_idx, br, bc, t + 1)
                        # P(r,c,t) ^ P(nr,nc,t+1) ^ B(b,nr,nc,t) => B(b,br,bc,t+1)
                        self.cnf.append([-p_t, -p_t1, -b_t, b_t1_new])

            # A box stays put unless pushed.
            for b_idx in range(self.num_boxes):
                for r, c in [(r, c) for r in range(self.N) for c in range(self.M) if self._is_free(r,c)]:
                    b_t = self.var_box(b_idx, r, c, t)
                    b_t1 = self.var_box(b_idx, r, c, t + 1)
                    
                    # A box at (r,c) is pushed if player comes from behind and moves onto it.
                    push_conditions = []
                    for dr, dc in DIRS.values():
                        pr, pc = r - dr, c - dc
                        if self._is_free(pr, pc):
                            p_behind = self.var_player(pr, pc, t)
                            p_on = self.var_player(r, c, t + 1)
                            # aux variable for (p_behind AND p_on)
                            aux = self._var('push', b_idx, r, c, t, pr, pc)
                            self.cnf.append([-aux, p_behind])
                            self.cnf.append([-aux, p_on])
                            self.cnf.append([aux, -p_behind, -p_on])
                            push_conditions.append(aux)

                    # B(t) AND NOT(was_pushed) => B(t+1)
                    # ~B(t) OR was_pushed OR B(t+1)
                    self.cnf.append([-b_t, b_t1] + push_conditions)
        # 3. Non-overlap constraints
        for t in range(self.T + 1):
            # Player is in exactly one position.
            player_vars = [self.var_player(r, c, t) for r in range(self.N) for c in range(self.M) if self._is_free(r, c)]
            self.cnf.append(player_vars) # At least one
            for i in range(len(player_vars)):
                for j in range(i + 1, len(player_vars)):
                    self.cnf.append([-player_vars[i], -player_vars[j]]) # At most one

            # Each box is in exactly one position.
            for b_idx in range(self.num_boxes):
                box_vars = [self.var_box(b_idx, r, c, t) for r in range(self.N) for c in range(self.M) if self._is_free(r,c)]
                self.cnf.append(box_vars) # At least one
                for i in range(len(box_vars)):
                    for j in range(i + 1, len(box_vars)):
                        self.cnf.append([-box_vars[i], -box_vars[j]]) # At most one
            
            # No two objects can be in the same cell.
            for r in range(self.N):
                for c in range(self.M):
                    if not self._is_free(r,c): continue
                    
                    # Player and any box cannot overlap.
                    p_var = self.var_player(r, c, t)
                    for b_idx in range(self.num_boxes):
                        b_var = self.var_box(b_idx, r, c, t)
                        self.cnf.append([-p_var, -b_var])
                    
                    # No two boxes can overlap.
                    for b1 in range(self.num_boxes):
                        for b2 in range(b1 + 1, self.num_boxes):
                            b1_var = self.var_box(b1, r, c, t)
                            b2_var = self.var_box(b2, r, c, t)
                            self.cnf.append([-b1_var, -b2_var])            

        # 4. Goal State Constraint
        # At time T, every box must be on a goal cell.
        for b_idx in range(self.num_boxes):
            goal_vars = [self.var_box(b_idx, gr, gc, self.T) for gr, gc in self.goals if self._is_free(gr, gc)]
            self.cnf.append(goal_vars)
            
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