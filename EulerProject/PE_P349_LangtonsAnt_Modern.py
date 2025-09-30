"""
Description:
  Project Euler Problem #349
  Langton's ant
  https://projecteuler.net/problem=349
  Method: The brute force method is improved by:
    1. Encoding the move of the ant based on the direction and color
    2. Detecting the repeating pattern in the highway by an algorithm which combines KMP with bit manipulation

Assumptions:
  See the docstring of PE_P349_LangtonsAnt.py.

Method:
  Similar to the brute force solution, I will assume that P2 number of moves (e.g. 12000)
  would be enough to detect the highway pattern.
  Hence, in the first part I will let the ant move P2 times
  and store the data required to detect the highway pattern and to determine the black cell count.

  Then, I will perform the pattern detection and determine the blacck cell count
  using the prepared data.
  
  The move of the ant is encoded by:
    1. 4 directions are simulated by:
      * 0: east (x+1, y)
      * 1: south (x, y+1)
      * 2: west (x-1, y)
      * 3: north (x, y-1)
    2. The color of the cell after the move:
      * 0: white
      * 1: black
    
    Hence, we have 8 values that requires 3 bits.

  I will store the move codes in a bitarray: move_codes.
  This bitarray will be used to detect the pattern by using the bitwise equality operation.

  Another data to be stored is the move ids for one of the 8 move codes (e.g. 0).
  This array would help to activate the KMP algorithm.
  Lets assume that the move ids for code 0 are: [80, 90, 100, 110, 124, 228, 332, ...].
  This means that code 0 is achieved at the end of this move ids.
  when we analyze the occurances, it looks like there exists a pattern with the length of 10.
  But, 124 does not satisfy the pattern condition.
  10 satisfies the pattern condition 3 times but the 4th one is not achieved.
  Similarly, 104 satisfies the pattern condition 3 times as well.
  If it satisfies the conditin n times (e.g. 10) there may exist a pattern starting from the 110th move.
  We already know that the length of the pattern is 104.
  So it will satisfy the condition.
  Up to here, I use a modified form of KMP algorithm.
  In summary, from this algorithm we achieved two numbers: 10 and 104.
  Lets assume that both satisfies the pattern condition m times.
  I will perform the pattern detection
  for each number achieved in this algorithm (i.e. 10 and 104).

  Remember that I stored the move codes in a bitarray.
  The above algorithm is an initial step to speed up the pattern detection.
  In order to detect the pattern, we need to inspect
  the remaining codes in between the two occurances of code 0 (e.g. between moves 80 and 90).
  For this purpose, I will use bit operations instead of the KMP algorithm.
  I will create slices from the bitarray and inspect for bitwise equality.
  The pattern is detected if the equality holds for m times:
    move_codes[80,90] == move_codes[91,101] == move_codes[102,112] == move_codes[113,123] == ...
  
  




Space Complexity:
  The problem requests the black cell count after N = 1E18 moves.
  The travel in this approach contains
  P1 moves in the best case and P2 moves in the worst case.
  The row and column-wise size of the travel is limited to Z (ARRAY_SIZE_GRID).
  See Nomenclature section for P1, P2, and Z.
  
  The arrays listed in the Method section are the sources used.
  They are static and allocated by P2 and Z.
  Hence, the space complexity is S(aP2) + S(bZ^2) where
  a and b are constants and P2, Z are independent.
  Hence, the space complexity is linear for P2 but quadratic for Z,
  which is not memory critical as P2 and Z are too small compared to N.

Time Complexity:
  Let's assume that P2 > T (i.e., the algorithm succeeds).
  
  We can split the algorithm into two parts:
    1. Travel of the ant with P1 moves without pattern detection.
    2. Travel of the ant with [T - P1] moves with pattern detection.
  
  For the 1st part, the time complexity is O(k1P1) where
  k1 is a constant for the actions performed for the ant's move.
  
  Pattern requirement states that a sequence of move IDs must repeat
  the requested amount of times (n).
  n is initialized in the main routine.
  
  In the best case scenario, P1 >= T.
  Hence, it's expected that after P1 moves, the highway is already started
  and the pattern would be detected after the 1st move.
  In this case, the time complexity of the 2nd part is related
  to the length of the real pattern, which is too small compared to T.
  Hence, for this case (i.e., the best case),
  the time complexity of the 2nd part can be neglected.
  
  However, in the worst case (P1 < T), the travel would continue up to T moves.
  The pattern detection is processed once in w moves
  in order to reduce the runtime for the worst case.
  
  Let's assume that the current move is the p-th occurrence of a move ID.
  Let's also assume that, on average, each move ID is repeated once in r moves.
  Hence, p = [i / r], where
    i is the current move index, which is between P1 and T.
  
  The detection algorithm loops through the occurrences of the current move ID
  from the p-th till the n-th in the backward direction.
  Hence, the 1st loop runs for [(p - n) ~ p] times.
  
  The 2nd loop performs the comparison of the move IDs
  to check whether they are repeating n times or not.
  As we assumed that each move ID is repeated once in r moves,
  the number of moves between the two occurrences is r
  for the 1st cycle of the 1st loop.
  The move count becomes [2 * r] for the 2nd cycle and [j * r] for the j-th cycle.
  
  Hence, the ranges for the two loops are:
    1. p = [i / r]
    2. njr, where j is [1, p]
  
  Hence, the time complexity for pattern detection is:
    Sum(njr), where j is [1, p]
  which yields:
    [O(T^3) - O(P1^3)] / w ~ O(T^3) / w.
  
  Now, the time complexity is cubic.
  
  Increasing the value of w will decrease the runtime as expected.
  The time complexity evaluates to O(T^2) when w approaches T, 
  which yields a quadratic time complexity for the worst case.
  
  Average time complexity calculation is not possible, as a mathematical relation cannot be defined between P1, T, n, and w.
  
  In summary:
    - Best case: Omega(P1) where P1 >= T
    - Worst case: O(T^3) where P1 < T
    - Average: NA
  
  In order to see the effect of having a cubic time complexity, the following configuration is executed:
    - P1 = 1
    - w = 1
    - n = 10
  
  The runtime is 1680 seconds in a workstation.
  Nevertheless, as stated in the CAUTION section, the results are wrong!
  
  The runtime is less than a second for the best case.
  
  Note that, the algorithm contains only the basic operations like bit equality, increment, decrement, list/array random access, etc.
"""

import numpy as np

# The size limit for the containers indexed by the grid IDs
ARRAY_SIZE_GRID = np.uint16(512)

# The size limit for the containers indexed by the reduced move ID (see get_move_ids)
ARRAY_SIZE_MOVE_ID_REDUCED = np.uint8(8)

# The size limit for the containers indexed by the move indices
ARRAY_SIZE_MOVE_INDEX = np.uint16(30000)

# Stores the colors of the cells of the grid
# GRID_CELL_COLORS[i_row][i_clm]: np.bool_:
# True if the cell is black
GRID_CELL_COLORS = np.zeros(
  shape=(ARRAY_SIZE_GRID, ARRAY_SIZE_GRID),
  dtype=np.bool_)

# Stores the reduced move IDs indexed by move indices:
# MOVE_INDEX_TO_ID_FULL[move_index]: uint8:
# See get_move_ids for the definition of the reduced move ID
MOVE_INDEX_TO_ID_REDUCED = np.zeros(
  shape=(ARRAY_SIZE_MOVE_INDEX),
  dtype=np.uint8)

# Stores the full move IDs indexed by move indices:
# MOVE_INDEX_TO_ID_FULL[move_index]: uint32:
# See get_move_ids for the definition of the full move id
MOVE_INDEX_TO_ID_FULL = np.zeros(
  shape=(ARRAY_SIZE_MOVE_INDEX),
  dtype=np.uint32)

# Stores the occurrence/count of the reduced move ids indexed by the reduced move id:
# MOVE_ID_TO_REDUCED_OCCURRENCE[move_id_reduced]: uint16:
# The number of occurrences of move_id_reduced during the travel
# One of the values in the range: [0, ARRAY_SIZE_MOVE_INDEX]:
# 0: If move_id_reduced has no occurrence
# ARRAY_SIZE_MOVE_INDEX: If all moves have the same id (move_id_reduced)
MOVE_ID_TO_REDUCED_OCCURRENCE = np.zeros(
  shape=(ARRAY_SIZE_MOVE_ID_REDUCED),
  dtype=np.uint16)

# Stores the move indices indexed by the reduced move id and the occurrence/count of it:
# MOVE_ID_TO_REDUCED_INDEX[move_id_reduced, i_occurrence]: uint16:
# The index of the move corresponding to the ith occurrence of move_id_reduced.
# A move id may be repeated during the travel of the ant.
# Hence, the array is two-dimensional
# where the 2nd index is for the occurrence of the move id.
MOVE_ID_TO_REDUCED_INDEX = np.zeros(
  shape=(ARRAY_SIZE_MOVE_ID_REDUCED, ARRAY_SIZE_MOVE_INDEX),
  dtype=np.uint16)

# Stores the black cell count indexed by the move index:
# MOVE_INDEX_TO_BLACK_COUNT[move_index]: uint16:
# The number of black cells at the ith move
MOVE_INDEX_TO_BLACK_COUNT = np.zeros(
  shape=(ARRAY_SIZE_MOVE_INDEX),
  dtype=np.uint16)

# Stores the rotation relations based on the reduced move id.
# ROTATE_ORIENTATIONS_REDUCED[move_id_reduced][0]: np.int8: X-orientation after the rotation
# ROTATE_ORIENTATIONS_REDUCED[move_id_reduced][1]: np.int8: Y-orientation after the rotation
# See get_move_ids for the definition of the reduced move id
ROTATE_ORIENTATIONS_REDUCED = [
  [np.int8(1), np.int8(0)],  # [WHITE] North -> East
  [np.int8(0), np.int8(-1)], # [WHITE] East -> South
  [np.int8(0), np.int8(1)],  # [WHITE] West -> North
  [np.int8(-1), np.int8(0)], # [WHITE] South -> West
  [np.int8(-1), np.int8(0)], # [BLACK] North -> West
  [np.int8(0), np.int8(1)],  # [BLACK] East -> North
  [np.int8(0), np.int8(-1)], # [BLACK] West -> South
  [np.int8(1), np.int8(0)]]  # [BLACK] South -> East

def get_move_ids(
    current_cell_color,
    current_dir_x,
    current_dir_y,
    current_row,
    current_clm):
  """
  Description:
    Each move of the travel is assigned to two ids
    one excluding the grid information while the other including.
    Hence, the 1st move id contains the following information:
      a. The current colour of the cell before flipping,
      b. The orientation of the ant in x-direction before rotation
      c. The orientation of the ant in y-direction before rotation
    The three input parameters yield 8 possible values for the 1st move id
    as the 2nd and 3rd are constrained by NEWS (North, East, ...) directions.
  
    The 2nd move id, additionally, contains the grid locations:
      a. Current grid location in x-direction
      b. Current grid location in y-direction
    The five input parameters yield [8 * ARRAY_SIZE_GRID ^ 2] possible values
    for the 2nd move id.
  
  Parameters:
    current_cell_color: np.bool_
      The color of the cell before flipping
      True if black
    current_dir_x: np.int8
      X-direction of the ant before the rotation
      One of the following: 0, 1, -1
    current_dir_y: np.int8
      Y-direction of the ant before the rotation
      One of the following: 0, 1, -1
    current_row: np.uint16
      Current grid location in x-direction
    current_clm: np.uint16
      Current grid location in y-direction
  
  Returns:
    move_id_reduced: np.uint8
      The reduced id of the input move excluding the grid information:
      One of the values in the range: [0, 8]
    move_id_full: np.uint32
      The full id of the input move including the grid information:
      One of the values in the range: [0, 8 * ARRAY_SIZE_GRID ^ 2]
  
  Modifies:
    None
  """
  move_id_dir = 0  # u == 0 and v == 1 by definition (North)
  if current_dir_x == 1:  # v == 0 by definition (East)
    move_id_dir = 1
  elif current_dir_x == -1:  # v == 0 by definition (West)
    move_id_dir = 2
  elif current_dir_y == -1:  # u == 0 by definition (South)
    move_id_dir = 3
  
  coeff_dir = np.uint8(4)
  move_id_reduced = np.uint8(coeff_dir * current_cell_color + move_id_dir)
  
  coeff_clm = np.uint32(8)
  coeff_row = np.uint32(coeff_clm * np.uint32(ARRAY_SIZE_GRID))
  move_id_full_clm = np.uint32(
    coeff_clm * np.uint32(current_clm - 1) +
    np.uint32(move_id_reduced))
  move_id_full = np.uint32(
    coeff_row * np.uint32(current_row - 1) +
    move_id_full_clm)
  return move_id_reduced, move_id_full

def detect_pattern(current_move_index, pattern_repeat_count_req):
  """
  Description:
    Let's assume that the current move is the kth occurrence of the reduced move id, rmi.
    The method inspects the previous occurrences of rmi
    to determine whether one of them forms a pattern repeated n times.
    See Nomenclature section of the module docstring for the definition of n.
  
  Parameters:
    current_move_index: int
      The index of the current move
    pattern_repeat_count_req: np.uint8
      The required number of repeats for a pattern to be accepted
  
  Returns:
    move_index_pattern_start: uint16
      The move index where the pattern starts formation
    move_index_pattern_end: uint16
      The move index where the pattern ends
  
  Method:
    The method is based on the occurrences of a reduced move id.
    Consider the reduced move ids as letters (a, b, c, ...).
    Consider the sequence of reduced move ids as a string: ...abacadabacadabacada
    Consider "a" is the reduced move id corresponding to the input move.
    As you can see, abacad is the pattern we are looking for,
    but in that pattern, "a" is repeated many times.
    Consider, the last "a" is the kth occurrence of "a".
    The method inspects the (k-1)th, (k-2)th, ... nth occurrences of "a"
    if any satisfies the pattern rules:
      (k-1)th: ad: Fail
      (k-2)th: acad: Fail
      (k-3)th: abacad: Success
    See the Nomenclature section in the module docstring for the definition of "n".
    
    1. Get the reduced move id for the input move index.
    2. Loop through the previous occurrences of the reduced move id,
       starting from the last one till the nth one.
    
    For the ith cycle of the 1st loop:
    3. Assume the ith previous occurrence of the reduced move id is the 1st move of the pattern.
    4. Determine the length of the pattern:
       The distance from the ith previous occurrence to the input move minus 1.
    5. Run another loop with the range of the input pattern repeat requirement.
    
    For the jth cycle of the 2nd loop:
    6. Run another loop with a range of pattern length.
    
    For the kth cycle of the 3rd loop:
    7. Inspect the reduced and the full move ids to see whether a pattern exists.
    
    8. Return the pattern move indices if found, otherwise return None.
  
  Modifies:
    None
  """
  # Get the current move id
  move_id_reduced = MOVE_INDEX_TO_ID_REDUCED[current_move_index]
  
  # Loop through the previous occurrences of the input move id,
  # starting from the last one till the 1st one: range(last, 1st, -1)
  _range = range(
    MOVE_ID_TO_REDUCED_OCCURRENCE[move_id_reduced] - 1,
    pattern_repeat_count_req,
    -1)
  for pattern_move_id_count in _range:
    # Assume the ith previous occurrence is the 1st move of the pattern
    pattern_start_move_index_ith = MOVE_ID_TO_REDUCED_INDEX[
      move_id_reduced,
      pattern_move_id_count]
  
    # Determine the length of the pattern:
    # The distance from the ith previous occurrence to the input move minus 1
    pattern_length = current_move_index - pattern_start_move_index_ith
  
    # Inspect if there exist enough moves for pattern inspection
    pattern_start_move_index_1st = (
      pattern_start_move_index_ith -
      pattern_length * (pattern_repeat_count_req + 1))
    if pattern_start_move_index_1st < 1:
      return None, None
  
    # Inspect if a pattern satisfying the pattern requirements exists
    check_pattern = inspect_pattern_repeated_req(
      pattern_repeat_count_req,
      pattern_length,
      pattern_start_move_index_ith)
  
    # Return the pattern move indices if the pattern requirements are satisfied
    if check_pattern:
      move_index_pattern_start = pattern_start_move_index_ith
      move_index_pattern_end = move_index_pattern_start + pattern_length - 1
      return move_index_pattern_start, move_index_pattern_end
  
  # The current move index does not satisfy the pattern requirements
  return None, None

def inspect_pattern_repeated_req(
    pattern_repeat_count_req,
    pattern_length,
    pattern_start_move_index_ith):
  """
  Description:
    Inspects the reduced and the full move ids to check whether a sequence of the 
    pattern_length number of moves starting at the ith move index repeats n times.

  Parameters:
    pattern_repeat_count_req: np.uint8
      The required number of repeats for a pattern to be accepted
    pattern_length: np.uint16
      The length of the inspected pattern
    pattern_start_move_index_ith: np.uint16
      The move index where the inspected pattern starts

  Returns:
    bool: True if a pattern repeats for pattern_repeat_count_req

  Modifies:
    None
  """
  # Run a loop with the range of the input pattern repeat requirement
  for i_pattern_repeat in range(pattern_repeat_count_req):
    # The pattern starting move index for the jth pattern repeat
    pattern_start_move_index_jth = (
      pattern_start_move_index_ith -
      pattern_length * (i_pattern_repeat + 1))

    # Inspect the reduced and the full move ids for the jth pattern repeat
    if not inspect_pattern_once(
        pattern_length,
        pattern_start_move_index_ith,
        pattern_start_move_index_jth):
      return False

  # The pattern detection not failed -> A pattern repeating n times found
  return True

def inspect_pattern_once(
    pattern_length,
    pattern_start_move_index_ith,
    pattern_start_move_index_jth):
  """
  Description:
    Consider two sequences of moves with the same length (pattern_length):
      1. Starting from pattern_start_move_index_ith
      2. Starting from pattern_start_move_index_jth

    Inspects all corresponding moves from the two sequences
    whether they satisfy the following pattern requirements:
      1. The cell colour must be the same
      2. The row-wise shift must be the same
      3. The clm-wise shift must be the same
    The above checks can be performed using MOVE_INDEX_TO_ID_FULL
    as it contains all the three information.

  Parameters:
    pattern_length: np.uint16
      The length of the inspected pattern
    pattern_start_move_index_ith: np.uint16
      The move index where the inspected pattern starts
    pattern_start_move_index_jth: np.uint16
      The move index where the inspected pattern repeat starts

  Returns:
    bool: True if a pattern repeats once

  Modifies:
    None
  """
  # Get the difference between the corresponding 1st pattern move ids
  diff_full_1st = np.int32(
    np.int32(MOVE_INDEX_TO_ID_FULL[pattern_start_move_index_ith]) -
    np.int32(MOVE_INDEX_TO_ID_FULL[pattern_start_move_index_jth]))

  # Run a loop with the range of the input pattern length
  for i in range(1, pattern_length):
    # The difference between the corresponding pattern move ids must be the same
    diff_full_current = np.int32(
      np.int32(MOVE_INDEX_TO_ID_FULL[pattern_start_move_index_ith + i]) -
      np.int32(MOVE_INDEX_TO_ID_FULL[pattern_start_move_index_jth + i]))
    if diff_full_current != diff_full_1st:
      return False

  # None of the corresponding moves failed -> A pattern found
  return True

def perform_limited_travel(
    initials,
    travel_move_count_limit,
    pattern_detection_start_move_index,
    pattern_detection_range,
    pattern_repeat_count_req):
  """
  Description:
    Performs a travel of the ant in order to detect the highway pattern.
    The travel is limited as the global arrays are static.
  
    Starts executing the pattern detection after (P1)th move.
  
    P1 is used to delay the pattern detection
    to skip the arbitrary travel region.
  
    See Nomenclature section of the module docstring for the definition of P1.
  
    Fills the global arrays during the travel.

  Parameters:
    initials: list[]
      initials[0]: np.uint16: Initial row id
      initials[1]: np.uint16: Initial column id
      initials[2]: np.int8: Initial X-direction
      initials[3]: np.int8: Initial Y-direction
    travel_move_count_limit: np.uint16
      A limit value for the move count in order to prevent an infinite loop
      in case of a failure in the pattern detection procedure.
    pattern_detection_start_move_index: np.uint16
      This variable is used to delay the pattern detection.
    pattern_detection_range: np.uint16
      Perform pattern detection after pattern_detection_start_move_index
      in every pattern_detection_range moves.
      See module docstring Time Complexity section
    pattern_repeat_count_req: np.uint8
      The required number of repeats for a pattern to be accepted

  Returns:
    move_index_pattern_start: uint16
      The move index where the pattern starts formation
    move_index_pattern_end: uint16
      The move index where the pattern ends

  Modifies:
    The global variables are initialized right after the module docstring:
      GRID_CELL_COLORS
      MOVE_INDEX_TO_ID_REDUCED
      MOVE_INDEX_TO_ID_FULL
      MOVE_ID_TO_REDUCED_OCCURRENCE
      MOVE_ID_TO_REDUCED_INDEX
      MOVE_INDEX_TO_BLACK_COUNT
  """
  # Initialize the travel
  row = initials[0]
  clm = initials[1]
  dir_x = initials[2]
  dir_y = initials[3]
  black_count = 0

  # Run the ant until the pattern is detected
  pattern_detection_counter = 0
  for move_index in range(travel_move_count_limit):
    # Flip the cell colour
    GRID_CELL_COLORS[row][clm] = not GRID_CELL_COLORS[row][clm]
    if GRID_CELL_COLORS[row][clm]:
      black_count += 1
    else:
      black_count -= 1
  
    # Get the move ids before moving the ant
    move_id_reduced, move_id_full = get_move_ids(
      GRID_CELL_COLORS[row][clm], dir_x, dir_y, row, clm)
  
    # Fill the global arrays
    MOVE_INDEX_TO_ID_REDUCED[move_index] = move_id_reduced
    MOVE_INDEX_TO_ID_FULL[move_index] = move_id_full
    MOVE_INDEX_TO_BLACK_COUNT[move_index] = black_count
    MOVE_ID_TO_REDUCED_OCCURRENCE[move_id_reduced] += 1
    MOVE_ID_TO_REDUCED_INDEX[
      move_id_reduced,
      MOVE_ID_TO_REDUCED_OCCURRENCE[move_id_reduced]] = move_index
  
    # Inspect if the pattern with the required repeat count is detected
    if (
        move_index >= pattern_detection_start_move_index and
        pattern_detection_counter == 0):
      move_index_pattern_start, move_index_pattern_end = detect_pattern(
        move_index, pattern_repeat_count_req)
      if move_index_pattern_start is not None:
        return move_index_pattern_start, move_index_pattern_end
  
    # Move the ant
    [dir_x, dir_y] = ROTATE_ORIENTATIONS_REDUCED[int(move_id_reduced)]
    row += dir_x
    clm += dir_y
  
    # Increment the pattern detection counter
    pattern_detection_counter += 1
    if pattern_detection_counter == pattern_detection_range:
      pattern_detection_counter = 0
  
  # Pattern detection has failed
  return None, None

def determine_black_count(
    move_count_req,
    move_index_pattern_start,
    move_index_pattern_end):
  """
  Description:
    Determines the black cell count for the whole travel of the ant.

  Method:
    The total number of black cells for the whole travel of the ant
    is the sum of the black cell counts in the arbitrary and highway regions.
  
    The black cells in the highway region can be calculated by...
    determining the black cells in the repeating pattern.
    
    However, the total move count of the highway may not be divided to the
    repeating pattern move count evenly.
    
    If so, the count of the black cells must be calculated
    for the remaining moves as well.

  Parameters:
    move_count_req: np.uint64
      The move count requirement by Project Euler Problem #349
    move_index_pattern_start: uint16
      The move index where the pattern starts formation
    move_index_pattern_end: uint16
      The move index where the pattern ends for a single pattern repeat
  
  Returns:
    np.uint64: The total number of the black cells for the whole travel of the ant
  
  Modifies:
    None
  """
  # Determine the move indices
  move_index_arb_start = np.uint64(0)
  move_index_arb_end = np.uint64(move_index_pattern_start - 1)
  
  # Determine the move counts
  move_count_pattern = np.uint64(
    move_index_pattern_end - move_index_pattern_start + 1)
  move_count_arb = np.uint64(move_index_arb_end - move_index_arb_start + 1)
  move_count_highway = np.uint64(move_count_req - move_count_arb)
  move_count_remaining = np.uint64(move_count_highway % move_count_pattern)
  
  # Determine the black cell counts for the whole travel of the ant
  black_count_pattern = np.uint64(
    MOVE_INDEX_TO_BLACK_COUNT[move_index_pattern_end] -
    MOVE_INDEX_TO_BLACK_COUNT[move_index_pattern_start - 1])
  black_count_arb = np.uint64(MOVE_INDEX_TO_BLACK_COUNT[move_index_arb_end])
  black_count_highway = np.uint64(
    np.uint64((move_count_highway - move_count_remaining) / move_count_pattern) *
    black_count_pattern)
  if move_count_remaining == 0:
    return np.uint64(black_count_arb + black_count_highway)
  
  # Determine the move indices for the remaining moves
  move_index_remaining_start = int(move_index_pattern_start)
  move_index_remaining_end = (
    move_index_remaining_start + int(move_count_remaining) - 1)
  
  # Determine the black cell count for the remaining moves
  black_count_remaining = np.uint64(
    MOVE_INDEX_TO_BLACK_COUNT[move_index_remaining_end] -
    MOVE_INDEX_TO_BLACK_COUNT[move_index_remaining_start - 1])
  
  # Return the total number of the black cells for the whole travel of the ant
  return np.uint64(black_count_arb + black_count_highway + black_count_remaining)

import time

def main():
  """
  Description:
    The main function
  
  Method:
    See the module docstring
  
  Assumptions:
    See the module docstring
  
  Parameters:
    None
  
  Returns:
    np.uint64: The total number of the black cells for the whole travel of the ant
  """
  t0 = time.time()
  
  # Initialize the ant.
  initial_row = np.uint16(ARRAY_SIZE_GRID / 2)
  initial_clm = np.uint16(ARRAY_SIZE_GRID / 2)
  initial_dir_x = np.int8(0)
  initial_dir_y = np.int8(-1)
  
  # Set the required move count
  move_count_req = np.uint64(10 ** 18)
  
  # Set a limit value for the move count in order to prevent an infinite loop
  # in case of a failure in the pattern detection procedure.
  travel_move_count_limit = np.uint16(30000)
  
  # This variable is used to delay the start of the pattern detection procedure.
  pattern_detection_start_move_index = np.uint16(10000)
  
  # The required number of repeats for a pattern to be accepted
  pattern_repeat_count_req = np.uint8(10)
  
  # Perform pattern detection after pattern_detection_start_move_index
  # in every pattern_detection_range moves.
  pattern_detection_range = 100
  
  # Perform the limited travel and get the pattern move indices
  move_index_pattern_start, move_index_pattern_end = perform_limited_travel(
    [initial_row, initial_clm, initial_dir_x, initial_dir_y],
    travel_move_count_limit,
    pattern_detection_start_move_index,
    pattern_detection_range,
    pattern_repeat_count_req)
  
  t1 = time.time()
  print(t1 - t0)
  
  # The pattern detection has failed
  if move_index_pattern_start is None:
    return None
  
  # Determine the total number of the black cells for the whole travel of the ant
  return determine_black_count(
    move_count_req,
    move_index_pattern_start,
    move_index_pattern_end)

if __name__ == '__main__':
  print ("Langton's ant problem:")
  
  black_count = main()
  if black_count is None:
    print ("ERROR: Pattern detection has failed")
  else:
    print ("The number of the black cells for 10^18 moves: " + str(black_count))
