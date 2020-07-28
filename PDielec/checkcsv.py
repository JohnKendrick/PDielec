#!/usr/bin/env python
"""Read in two csv files and check that the float contents are the same"""
from __future__ import print_function
import sys
from os.path import splitext
from termcolor import colored
global separator
global threshold

def compare_lines( line1, line2, swapped):
  global separator
  global threshold
  store_error = False
  nerrors = 0
  percentage_error = 0.0
  max_percentage_error = 0.0
  keep_word1 = ''
  keep_word2 = ''
  index = 0
  len1 = len(line1.split(separator))
  len2 = len(line2.split(separator))
  if len1 == len2:
      if len1 > 4:
          word1 = line1.split(separator)[3]
          word2 = line2.split(separator)[3]
          if word1 == word2  and word1 == "processors":
              # Ignore this line 'cos the number of processors used depends on the machine
              return (store_error, nerrors, max_percentage_error, keep_word1, keep_word2)
          # end if word1 == word2
      # end if len1 > 4
  #end if len1 == len2
  for word1,word2 in zip(line1.split(separator),line2.split(separator)):
    index += 1
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
  try:
    float(value)
    return True
  except ValueError:
    return False


# Start processing the directories
if len(sys.argv) <= 1 :
    print('checkcsv file file2 [-sep blanks] [-thresh 3.0E-2]', file=sys.stderr)
    print('         Compare the two files for any significant changes', file=sys.stderr)
    print('         The default separator is "blanks" for txt/out files', file=sys.stderr)
    print('         and "," for csv files', file=sys.stderr)
    print('         Numbers f1 and f2 in each file are compared', file=sys.stderr)
    print('         an error is flagged if 2*abs(f1-f2)/(f1+f2+2) > threshold', file=sys.stderr)
    print('         The default value for threshold is 3.0E-2 ', file=sys.stderr)
    exit()

threshold = 3.0E-2
separator = None
tokens = sys.argv[1:]
ntokens = len(tokens)-1
itoken = -1
files = []
while itoken < ntokens:
    itoken += 1
    token = tokens[itoken]
    if token == '-sep':
        itoken +=1
        separator = tokens[itoken]
        if separator == 'blanks':
            separator = None
    elif token == '-thresh':
        itoken +=1
        threshold = float(tokens[itoken])
    else:
        files.append(tokens[itoken])
    # end if
# end while
file1 = files[0]
file2 = files[1]
file_name,extension = splitext(file1)
if extension == '.csv':
    separator = ','
elif extension == '.txt':
    separator = None
elif extension == '.out':
    separator = None
# Read the command line again
files = []
while itoken < ntokens:
    itoken += 1
    token = tokens[itoken]
    if token == '-sep':
        itoken +=1
        separator = tokens[itoken]
        if separator == 'blanks':
            separator = None
    elif token == '-thresh':
        itoken +=1
        threshold = float(tokens[itoken])
    else:
        files.append(tokens[itoken])
    # end if
# end while
fd1 = open(file1)
fd2 = open(file2)
#
line1_errors = []
line2_errors = []
line_numbers_with_error = []
nerrors = 0
store_error = False
max_percentage_error = 0.0
max_percentage_error2= 0.0
line_number = 0
compare_next_line = False
for line1,line2 in zip(fd1,fd2):
  line_number += 1
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
      # endif
   # end if not
#end for line1, line2
if nerrors > 0:
  print('  '+colored('ERRORS:','red')+"({}) LARGEST ON LINE {} OF {}({}) and {}({}) -- max %error={}".format(nerrors, keep_line_number, file1,keep_word1,file2,keep_word2,max_percentage_error))
else:
  print('  '+colored('OK:','blue')+" {} = {} -- max %error={}" .format(file1,file2,max_percentage_error))
