# -*- coding: utf-8 -*-
"""
A simple solution for the ProjectEuler problem:
    ID: 14
    Name: Longest Collatz sequence

Problem definition: https://projecteuler.net/problem=14
    The following iterative sequence is defined for the set of positive integers:

    n -> n/2 (n is even)
    n -> 3n + 1 (n is odd)

    Using the rule above and starting with 13, we generate the following sequence:

    13 -> 40 -> 20 -> 10 -> 5 -> 16 -> 8 -> 4 -> 2 -> 1
    It can be seen that this sequence (starting at 13 and finishing at 1) contains 10 terms.
    Although it has not been proved yet (Collatz Problem),
    it is thought that all starting numbers finish at 1.

    Which starting number, under one million, produces the longest chain?

    NOTE: Once the chain starts the terms are allowed to go above one million.

Solutions:
    CAUTION:
        BOTH SOLUTIONS ASSUME THAT THE COLLATZ CONJECTURE IS PROVEN.
        HENCE, THE PROBLEM IS BOUNDED (THE LOOP ENDS UP WITH 1).
        BUT THEORITICALLY, THERE IS NO PROOF FOR THE CONJECTURE.
        AS ITS NOT PROVEN YET, THE CODE MUST INSPECT CONVERGENCE.
        THE CONVERGENCE CRITERIA IS NOT IMPPLEMENTED IN THIS SHEET FOR SIMPLICITY.
        HENCE, THE PROGRAM MAY FALL INTO AN INFINITE LOOP,
        IF THERE EXIST EXCEPTIONAL NUMBERS FOR THE CONJECTURE (A DISPROVE).

    Solution 1 (Direct solution):
        Description:
            Create a loop from 3 to the input max number,
            For each number in the loop, determine the sequence ending with 1
            Determine the max length of the sequences

            The sequence is determined with a recursive function

        Time complexity:
            Sequence determination will run for each cycle of the loop.
            Hence, the time complexity is: O(N) * T
            where T stands for the worst case for the determination of the sequence.
            However, T can not be determined as its equivalent
            to prove (or disprove) the Collatz conjecture.
            If the conjecture exists, the upper bound of the sequence can somehow be determined.
            Then, T can be calculated.

        Space complexity:
            No storage needed: O(1)

    Solution 2 (Single visit solution):
        Description:
            During the calculations, store two values for each number:
                The number of the values visitted till the current number,
                The number of the values remaining to 1.

            For each new value, inspect if any value stored for the number of the visitted numbers:
                If a value stored:
                    If the new value is higher than the stored,
                    set the new value and continue for the remaining numbers to 1.
                    Otherwise, continue with the new cycle.

            With this method, each number is visitted only once.

            The sequence is determined with a recursive function

        Time complexity:
            Cannot be determined (see time complexity for Solution1).

        Space complexity:
            No limitation for the space complexity as there is no upper bound for the problem.
            Hence, the code may fail due to the memory problems.

@author: baris.albayrak.ieee@gmail.com
"""

import numpy as np
import time

def find_next_val(val):
    '''
    Description:
        Find the next value in the sequence

    Parameters:
        val: int:
            The number for which the sequence is being determined

    Outputs:
        int: The next number in the sequence
    '''
    if val % 2:
        return 3 * val + 1
    return int(val / 2)








'''
*******************************************
*******************************************
SOLUTION 1
*******************************************
*******************************************
'''

def find_sequence_1(current_val, counter):
    '''
    Description:
        The sequence length determination for Solution 1 in the module docstring.

    Parameters:
        current_val : int
            DESCRIPTION:
                The number for which the sequence is determined
        counter : int
            DESCRIPTION:
                The number of the visitted values up to the input number.
                Required due to the recursion

    Outputs:
        int:
            DESCRIPTION: The length of the sequence for the input number
    '''
    counter += 1
    next_val = find_next_val(current_val)
    if next_val == 1:
        return counter

    return counter + find_sequence_2(next_val, 0)

def solution_1(limit_val):
    '''
    Description:
        Solution 1 in the module docstring.

    Parameters:
        limit_val : int:
            The upper bound of the problem

    Outputs:
        c_max: int:
            The number for which the length of the sequence is max
        n: int:
            The length of the longest sequence
    '''
    c_max = 2
    for i in range(3, limit_val):
        c_count = find_sequence_1(i, 0)
        if c_count > c_max:
            c_max = c_count
            n = i

    return n, c_max














'''
*******************************************
*******************************************
SOLUTION 2
*******************************************
*******************************************
'''

ARRAY_BOUND = int(1e7)
SEQUENCE1 = np.zeros([ARRAY_BOUND, 2], dtype=int) # Array store
SEQUENCE2 = {} # Dictionary store. For numbers higher than the array upper bound

def get_visitted_remaining(val):
    '''
    Description:
        Get the number of visitted points and the number of the remaining points
        from the stored data if exist.

    Parameters:
        val : int:
            The number for which the sequence is being determined

    Outputs:
        [0]: int:
            The number of visitted points:
                Zero if the input number is not visitted yet
        [1]: int:
            The number of the remaining points:
                Zero if the input number is not visitted yet
    '''
    if val < ARRAY_BOUND:
        return SEQUENCE1[val]
    if val in SEQUENCE2.keys():
        return SEQUENCE2[val]
    return [0, 0]

def set_visitted_remaining(val, visitteds, remainings):
    '''
    Description:
        Set the number of visitted points and the number of the remaining points

    Parameters:
        val : int:
            The number for which the sequence is being determined
        visitteds : int:
            The number of visitted points
        remainings : int:
            The number of the remaining points

    Outputs:
        void
    '''
    if val < ARRAY_BOUND:
        SEQUENCE1[val][0] = visitteds
        SEQUENCE1[val][1] = remainings

    SEQUENCE2[val] = [visitteds, remainings]

def find_sequence_2(current_val, counter):
    '''
    Description:
        The sequence length determination for Solution 2 in the module docstring.

    Parameters:
        current_val : int
            DESCRIPTION:
                The number for which the sequence is determined
        counter : int
            DESCRIPTION:
                The number of the visitted values up to the input number.
                Required due to the recursion

    Outputs:
        int:
            DESCRIPTION: The length of the sequence for the input number
    '''
    counter += 1
    next_val = find_next_val(current_val)
    if next_val == 1:
        return counter

    totals = get_visitted_remaining(next_val)
    if counter < totals[0]:
        return 0

    totals[0] = counter
    if totals[1] > 0:
        return counter + totals[1]

    remainings = find_sequence_2(next_val, 0)
    set_visitted_remaining(next_val, counter, remainings)
    return counter + remainings

def solution_2(limit_val):
    '''
    Description:
        Solution 2 in the module docstring.

    Parameters:
        limit_val : int:
            The upper bound of the problem

    Outputs:
        c_max: int:
            The number for which the length of the sequence is max
        n: int:
            The length of the longest sequence
    '''
    c_max = 2
    for i in range(3, limit_val):
        if get_visitted_remaining(i)[0] > 0:
            continue

        c_count = find_sequence_2(i, 0)
        if c_count > c_max:
            c_max = c_count
            n = i

    return n, c_max












if __name__ == '__main__':
    MAX = int(1e6)

    # Solution 1
    t0 = time.time()
    n, c_max = solution_1(MAX)
    t1 = time.time()

    print('Upper Bound: ' + str(MAX))

    print('\nSolution 1:')
    print('Sequence length: ' + str(n))
    print('Number with longest sequence: ' + str(c_max))
    print('Runtime: ' + str(t1 - t0))

    # Solution 2
    t0 = time.time()
    n, c_max = solution_2(MAX)
    t1 = time.time()

    print('\nSolution 2:')
    print('Sequence length: ' + str(n))
    print('Number with longest sequence: ' + str(c_max))
    print('Runtime: ' + str(t1 - t0))
