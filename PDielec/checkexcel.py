#!/usr/bin/env python
"""Read in two csv files and check that the float contents are the same"""
from __future__ import print_function
import sys
from termcolor import colored
from openpyxl import load_workbook

global threshold

# Start processing the directories
if len(sys.argv) <= 1 :
    print('checkexcel file file2 [-thresh 1.0E-3] [-f]', file=sys.stderr)
    print('         Compare the two excel files for any significant changes', file=sys.stderr)
    print('         Numbers f1 and f2 in each file are compared', file=sys.stderr)
    print('         an error is flagged if 2*abs(f1-f2)/(f1+f2+2) > threshold', file=sys.stderr)
    print('         The default value for threshold is 1.0E-3 ', file=sys.stderr)
    print('         -f forces a full comparison of the sheets ', file=sys.stderr)
    print('            by default Settings is not included    ', file=sys.stderr)
    exit()

threshold = 1.0E-3
separator = None
tokens = sys.argv[1:]
ntokens = len(tokens)-1
itoken = -1
files = []
full = False
while itoken < ntokens:
    itoken += 1
    token = tokens[itoken]
    if token == '-thresh':
        itoken +=1
        threshold = float(tokens[itoken])
    if token == '-f':
        full = True
    else:
        files.append(tokens[itoken])
    # end if
# end while

file1 = files[0]
file2 = files[1]

#
# Open the work books
#
wb1 = load_workbook(filename=file1, read_only=True)
wb2 = load_workbook(filename=file2, read_only=True)
#
# Initialise variables
#
error = None
nerrors = 0
max_percentage_error = 0.0
#
# Loop over sheets
#
sheets = ['Main','Scenarios','Powder Molar Absorption (cells)','Powder Absorption','Powder Real Permittivity','Powder Imaginary Permittivity', 'Powder ATR Reflectance', 'Analysis','Crystal R_p','Crystal R_s','Crystal T_p','Crystal T_s','Real Crystal Permittivity','Imag Crystal Permittivity']
if full:
    sheets.append('Settings')
for sheet in sheets:
    if not sheet in  wb1 :
        continue
    if not sheet in  wb2 :
        continue
    print('Checking sheet ',sheet)
    ws1 = wb1[sheet]
    ws2 = wb2[sheet]
    max_rows1 = ws1.max_row
    max_rows2 = ws2.max_row
    max_columns1 = ws1.max_column
    max_columns2 = ws2.max_column
    if max_rows1 != max_rows2:
        print('Error - the number rows in the sheet are not the same',sheet,max_rows1,max_rows2)
        nerrors += 1
        continue
    if max_columns1 != max_columns2:
        print('Error - the number columns in the sheet are not the same',sheet,max_columns1,max_columns2)
        nerrors += 1
        continue
    #
    # Loop over rows
    #
    row_index = 0
    for row1,row2 in zip(ws1.rows, ws2.rows):
        row_index += 1
        col_index = 0
        #
        # Loop over cells
        #
        for cell1,cell2, in zip(row1,row2):
            col_index += 1
            value1 = cell1.value
            value2 = cell2.value
            if cell1 is None and cell2 is None:
                pass
            elif cell1 is None or cell2 is None:
                nerrors += 1
                percentage_error = 0.0
                error = (sheet, row_index, col_index, value1, value2, percentage_error)
            else:
                if value1 != value2:
                    if cell1.data_type == 'n' and cell2.data_type == 'n':
                        #
                        # Flag an error which is numeric
                        #
                        percentage_error = 100.0*abs(2.0*(value1 - value2)/(abs(value1) + abs(value2)+2))
                        if percentage_error > 100.*threshold:
                            nerrors += 1
                            if percentage_error > max_percentage_error:
                                max_percentage_error = percentage_error
                                error = (sheet, row_index, col_index, value1, value2, percentage_error)
                            # if percentage error
                        else:
                            if percentage_error > max_percentage_error:
                                max_percentage_error = percentage_error
                        # if percentage error > threshold
                    else:
                        #
                        # This is a non-numeric error
                        #
                        nerrors += 1
                        percentage_error = 0.0
                        error = (sheet, row_index, col_index, value1, value2, percentage_error)
                    # if cell1.data_type
                # if value1 != value2
            # if cell1 is none
        # for cell1, cell2
    # for rowq, row
# for sheet
if error is not None:
    sheet,row,col,value1,value2,max_percentage_error = error
    print('  '+colored('ERRORS:','red')+'({}) LARGEST ON ROW,COL {},{} OF SHEET {}, {}({}) and {}({}) -- max %error={}'.format(nerrors, row,col,sheet,file1,value1,file2,value2,max_percentage_error))
elif nerrors > 0:
    print('  '+colored('ERRORS:','red')+'({}) Dimensions of spreadsheet were wrong                                    '.format(nerrors))
else:
    print('  '+colored('OK:','blue')+" {} = {} -- max %error={}" .format(file1,file2,max_percentage_error))
