"""
Sokoban Solver using traditional search algorithms
--------------------------------------
This implementation uses basic Python concepts to solve the Sokoban puzzle
without using SAT solver. It uses BFS to find the solution.

Grid Encoding:
- 'P' = Player
- 'B' = Box
- 'G' = Goal
- '#' = Wall
- '.' = Empty space
"""

from collections import deque, defaultdict
from copy import deepcopy

class SokobanStater:
    def __init__(self, grid, player_pos, boxes, move_sequence=""):
        self.grid = grid
        self.player_pos = player_pos
        self.boxes = set(boxes)  # Convert to set for faster lookup
        self.move_sequence = move_sequence
        self.hash_val = None

    def __hash__(self):
        if self.hash_val is None:
            # Hash based on player position and frozen set of box positions
            self.hash_val = hash((self.player_pos, frozenset(self.boxes)))
        return self.hash_val

    def __eq__(self, other):
        return (self.player_pos == other.player_pos and 
                self.boxes == other.boxes)

class SokobanSolver:
    def __init__(self, grid, T):
        self.grid = [list(row) for row in grid]
        self.T = T
        self.N = len(grid)
        self.M = len(grid[0])
        self.goals = set()
        self.boxes = set()
        self.player_start = None
        self.directions = {
            'U': (-1, 0),
            'D': (1, 0),
            'L': (0, -1),
            'R': (0, 1)
        }
        self._parse_grid()

    def _parse_grid(self):
        """Find initial positions of player, boxes, and goals"""
        for r in range(self.N):
            for c in range(self.M):
                cell = self.grid[r][c]
                if cell == 'P':
                    self.player_start = (r, c)
                    self.grid[r][c] = '.'  # Convert to empty space for easier checking
                elif cell == 'B':
                    self.boxes.add((r, c))
                    self.grid[r][c] = '.'
                elif cell == 'G':
                    self.goals.add((r, c))

    def _is_valid_move(self, pos):
        """Check if position is within bounds and not a wall"""
        r, c = pos
        return (0 <= r < self.N and 
                0 <= c < self.M and 
                self.grid[r][c] != '#')

    def _get_next_pos(self, pos, direction):
        """Get next position after moving in given direction"""
        return (pos[0] + direction[0], pos[1] + direction[1])

    def _is_deadlock(self, boxes):
        """Check if boxes are in a deadlock position"""
        for box in boxes:
            if box in self.goals:
                continue

            r, c = box
            # Check if box is in corner
            horizontal_blocked = False
            vertical_blocked = False

            # Check horizontal movement
            if (not self._is_valid_move((r, c-1)) or 
                not self._is_valid_move((r, c+1))):
                horizontal_blocked = True

            # Check vertical movement
            if (not self._is_valid_move((r-1, c)) or 
                not self._is_valid_move((r+1, c))):
                vertical_blocked = True

            # If both directions are blocked and not on goal, it's a deadlock
            if horizontal_blocked and vertical_blocked:
                return True

        return False

    def _get_next_states(self, state):
        """Generate all possible next states from current state"""
        next_states = []
        
        for move, direction in self.directions.items():
            # New player position
            new_player_pos = self._get_next_pos(state.player_pos, direction)
            
            # Check if move is valid
            if not self._is_valid_move(new_player_pos):
                continue

            # If there's a box in new position
            if new_player_pos in state.boxes:
                # Calculate where box would move
                new_box_pos = self._get_next_pos(new_player_pos, direction)
                
                # Check if box can be moved
                if (self._is_valid_move(new_box_pos) and 
                    new_box_pos not in state.boxes):
                    # Move box
                    new_boxes = set(state.boxes)
                    new_boxes.remove(new_player_pos)
                    new_boxes.add(new_box_pos)
                    
                    # Check for deadlock in new position
                    if not self._is_deadlock(new_boxes):
                        new_state = SokobanState(
                            self.grid,
                            new_player_pos,
                            new_boxes,
                            state.move_sequence + move
                        )
                        next_states.append(new_state)
            else:
                # Just move player
                new_state = SokobanState(
                    self.grid,
                    new_player_pos,
                    state.boxes,
                    state.move_sequence + move
                )
                next_states.append(new_state)

        return next_states

    def _is_goal_state(self, boxes):
        """Check if all boxes are on goals"""
        return all(box in self.goals for box in boxes)

    def solve(self):
        """
        Solve the Sokoban puzzle using BFS.
        Returns move sequence or -1 if no solution found.
        """
        if len(self.boxes) > len(self.goals):
            return -1  # Impossible case: more boxes than goals

        # Create initial state
        initial_state = SokobanState(self.grid, self.player_start, self.boxes)
        
        # Use BFS to find solution
        queue = deque([initial_state])
        visited = set([initial_state])
        
        while queue:
            current_state = queue.popleft()
            
            # Check if we've exceeded move limit
            if len(current_state.move_sequence) > self.T:
                continue
            
            # Check if we've reached goal state
            if self._is_goal_state(current_state.boxes):
                return current_state.move_sequence
            
            # Generate and process next states
            for next_state in self._get_next_states(current_state):
                if next_state not in visited:
                    visited.add(next_state)
                    queue.append(next_state)
        
        return -1  # No solution found

def solve_sokoban(grid, T):
    """
    Main function to solve Sokoban puzzle.
    Returns move sequence or -1 if unsolvable.
    """
    solver = SokobanSolver(grid, T)
    return solver.solve()
