#!/usr/bin/env python
#
# Copyright 2024 John Kendrick & Andrew Burnett
#
# This file is part of PDielec
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the MIT License
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#
# You should have received a copy of the MIT License along with this program, if not see https://opensource.org/licenses/MIT
#
"""The checkcsv command read in two csv files and check that the float contents are the same.

It is used by the :mod:`~PDielec.pdmake` command to check the validity of the reference calculations in the test suite.
"""
import sys
from os.path import splitext

from termcolor import colored

global separator
global threshold

def compare_lines( line1, line2, swapped):
  """Compare two lines of text to find and report differences in numerical values.

  Parameters
  ----------
  line1 : str
      The first line of text to be compared.
  line2 : str
      The second line of text to be compared.
  swapped : bool
      A flag indicating if the comparison allows for the first numerical value to be 
      considered as swapped, i.e., if an absolute difference of 1 between the first
      numerical values in each line should be ignored.

  Returns
  -------
  tuple : A tuple containing the following elements:
      - store_error : bool (A flag indicating if an error was stored during comparison.)
      - nerrors : int (The total number of errors detected in the comparison.)
      - max_percentage_error : float (The maximum percentage error found during comparisons.)
      - keep_word1 : str (The word from line1 that resulted in the maximum percentage error.)
      - keep_word2 : str (The word from line2 that resulted in the maximum percentage error.)

  Notes
  -----
  - The function makes use of two globals: `separator`, a string defining how to split
    the lines into words, and `threshold`, a float defining the percentage error
    threshold above which a discrepancy is considered an error.
  - The function compares corresponding numerical values in `line1` and `line2`. Non-numeric 
    words are ignored in the error calculation, but they contribute to the discrepancy count
    if one word is numeric and the corresponding word is not.
  - Special consideration is given to the word "processors" in the fourth position; 
    if present and matching in both lines, it's an immediate pass, regardless of subsequent content.
  - The function also allows for a tolerance in the comparison of the first numeric values 
    of the lines if `swapped` is True, ignoring a discrepancy of exactly 1.

  """    
  global separator
  global threshold
  store_error = False
  nerrors = 0
  percentage_error = 0.0
  max_percentage_error = 0.0
  keep_word1 = ""
  keep_word2 = ""
  len1 = len(line1.split(separator))
  len2 = len(line2.split(separator))
  if len1 == len2 and len1 > 4:
      word1 = line1.split(separator)[3]
      word2 = line2.split(separator)[3]
      if word1 == word2  and word1 == "processors":
          # Ignore this line 'cos the number of processors used depends on the machine
          return (store_error, nerrors, max_percentage_error, keep_word1, keep_word2)
      # end if word1 == word2 and len1 > 4
  #end if len1 == len2
  for index,(word1,word2) in enumerate(zip(line1.split(separator),line2.split(separator))):
    if isfloat(word1):
        if not isfloat(word2):
            nerrors +=1
            store_error = True
        else:
            float1 = float(word1)
            float2 = float(word2)
            error = abs(float1 - float2)
            error2 = abs(abs(float1 - float2) -1.0)
            # This complicated test is there so that if the lines have been swapped and the first number
            # is an index number we ignore the index in the comparison
            # EG:  File1                                File2
            #      5 0.99999                            5 0.88888
            #      6 0.88888                            6 0.99999
            if not swapped or index > 1 or error2 > 1.0e-10:
                percentage_error = 100.0*error/(abs(float1)+abs(float2)+2.0)/2.0
                if percentage_error > threshold:
                    nerrors += 1
                    store_error = True
                    if max_percentage_error < percentage_error:
                        keep_word1 = word1
                        keep_word2 = word2
                        max_percentage_error = percentage_error
                # end if
             #end if error_swapped
        #end if not
    #end if isfloat
  # end for word1, word2
  return (store_error, nerrors, max_percentage_error, keep_word1, keep_word2)


def isfloat(value):
  """Check if the given value can be converted to float.

  Parameters
  ----------
  value : any
      The value to be checked.

  Returns
  -------
  bool
      Returns True if value can be converted to float, False otherwise.

  """    
  try:
    float(value)
    return True
  except ValueError:
    return False

def main():
    """Compare two files line by line to check for significant changes based on a provided threshold.

    This script takes at least two arguments that are file names along with optional arguments for
    separator (`-sep`) and threshold value (`-thresh`). It compares numbers from corresponding lines
    in the two files. A significant change is flagged if the specified comparison formula exceeds
    a given threshold. Defaults to a threshold of 3.0E-2 and treats spaces as the default separator
    for `.txt` and `.out` files, and commas for `.csv` files.

    Parameters
    ----------
    None explicitly declared. The function operates on command-line arguments passed to the script.

    - `-sep` : str, optional
        Defines the separator to be used for splitting entries in a line.
        Default is None, which then defaults to spaces for `.txt` and `.out` files,
        and commas for `.csv` files.
    - `-thresh` : float, optional
        The threshold for flagging a difference between the corresponding numbers
        from the two files as significant. Default is 3.0E-2.

    Returns
    -------
    tuple
        A tuple containing the number of errors found (`nerrors`), the line number with
        the largest percentage difference (`keep_line_number`), the two values with the
        largest difference (`keep_word1` and `keep_word2`), and the maximum percentage
        difference found (`max_percentage_error`).

    Raises
    ------
    FileNotFoundError
        If one of the input files does not exist.

    Notes
    -----
    The command-line arguments need to be parsed using `sys.argv`.
    The comparison formula used is `2*abs(f1-f2)/(f1+f2+2)`.

    """    
    global separator
    global threshold
    # Start processing the directories
    if len(sys.argv) <= 1 :
        print("checkcsv file file2 [-sep blanks] [-thresh 3.0E-2]", file=sys.stderr)
        print("         Compare the two files for any significant changes", file=sys.stderr)
        print('         The default separator is "blanks" for txt/out files', file=sys.stderr)
        print('         and "," for csv files', file=sys.stderr)
        print("         Numbers f1 and f2 in each file are compared", file=sys.stderr)
        print("         an error is flagged if 2*abs(f1-f2)/(f1+f2+2) > threshold", file=sys.stderr)
        print("         The default value for threshold is 3.0E-2 ", file=sys.stderr)
        sys.exit()

    threshold = 3.0E-2
    separator = None
    tokens = sys.argv[1:]
    ntokens = len(tokens)-1
    itoken = -1
    files = []
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        if token == "-sep":
            itoken +=1
            separator = tokens[itoken]
            if separator == "blanks":
                separator = None
        elif token == "-thresh":
            itoken +=1
            threshold = float(tokens[itoken])
        else:
            files.append(tokens[itoken])
        # end if
    # end while
    file1 = files[0]
    file2 = files[1]
    file_name,extension = splitext(file1)
    if extension == ".csv":
        separator = ","
    elif extension in (".txt", ".out"):
        separator = None
    # Read the command line again
    files = []
    while itoken < ntokens:
        itoken += 1
        token = tokens[itoken]
        if token == "-sep":
            itoken +=1
            separator = tokens[itoken]
            if separator == "blanks":
                separator = None
        elif token == "-thresh":
            itoken +=1
            threshold = float(tokens[itoken])
        else:
            files.append(tokens[itoken])
        # end if
    # end while
    keep_line_number = 0
    keep_word1 = 0
    keep_word2 = 0
    nerrors = 0
    store_error = False
    max_percentage_error = 0.0
    compare_next_line = False
    with open(file1) as fd1, open(file2) as fd2:
        lines1 = fd1.readlines()
        lines2 = fd2.readlines()
        l1 = len(lines1)
        l2 = len(lines2)
        if l1 != l2:
            max_percentage_error = 100.0
            nerrors = nerrors +1
        for line_number,(line1,line2) in enumerate(zip(lines1,lines2)):
            if not compare_next_line:
                store_error, nerror, percentage_error, keep_word1c, keep_word2c = compare_lines(line1,line2,False)
                max_percentage_error = max(max_percentage_error,percentage_error)
                if store_error:
                    # if there was an error let's see if two lines are swapped
                    line1a = line1
                    line2a = line2
                    compare_next_line = True
            else:
                # There was an error on the previous line, see if lines are swapped
                store_error1, nerror1, percentage_error1, keep_word1a, keep_word2a = compare_lines(line1,line2a,True)
                store_error2, nerror2, percentage_error2, keep_word1b, keep_word2b = compare_lines(line2,line1a,True)
                compare_next_line = False
                keep_word1 = keep_word1a
                keep_word2 = keep_word2a
                keep_line_number = line_number
                if store_error1 or store_error2:
                    # There is still an error on the swapped lines so we have an error on line1 and line1a
                    nerrors += nerror
                    if percentage_error > max_percentage_error:
                        max_percentage_error = percentage_error
                        keep_word1 = keep_word1c
                        keep_word2 = keep_word2c
                        keep_line_number = line_number - 1
                    store_error, nerror, percentage_error, keep_word1c, keep_word2c = compare_lines(line1,line2,False)
                    nerrors += nerror
                    if percentage_error > max_percentage_error:
                        max_percentage_error = percentage_error
                        keep_word1 = keep_word1c
                        keep_word2 = keep_word2c
                        keep_line_number = line_number
                else:
                    # The swapped lines are the same so there is no error
                    store_error = False
                # endif store_error1...
            # end if not
        #end for line1, line2
    # Close open files
    if nerrors > 0:
      print("  "+colored("ERRORS:","red")+f"({nerrors}) LARGEST ON LINE {keep_line_number} OF {file1}({keep_word1}) and {file2}({keep_word2}) -- max %error={max_percentage_error}")
    else:
      print("  "+colored("OK:","blue")+f" {file1} = {file2} -- max %error={max_percentage_error}")
    # end
    return (nerrors,keep_line_number,keep_word1,keep_word2,max_percentage_error)

if __name__=="__main__":
    main()
