import xlsxwriter as xlsx
class SpreadSheetManager():
    def __init__(self, filename):
        # Manage a spread sheet for PDielec / PDGui
        self.workbook = xlsx.Workbook(filename)
        self.tab_names = ['Main', 'Settings', 'Scenarios', 'Molar Absorption','Absorption', 'Real Permittivity', 'Imaginary Permittivity', 'ATR Reflectance','Analysis']
        self.worksheets = {}
        # Positions points to where we write to next
        self.positions  = {}
        self.max_col    = {}
        self.max_row    = {}
        for tab in self.tab_names:
            self.worksheets[tab] = self.workbook.add_worksheet(tab)
            self.positions[tab] = (0,0)
            self.max_col[tab] = 0
            self.max_row[tab] = 0
        self.name = 'Main'

    def selectWorkSheet(self,name):
        self.name = name

    def writeNextRow(self,items, row=None, col=None, check=''):
        oldRow,oldCol = self.positions[self.name]
        if col is None:
            col = oldCol
        if row is None:
            row = oldRow
        if col == 0:
            print('We have a problem, col is 0')
        self.write(row,0,check)
        for item in items:
            #print('writing ', row, col, item)
            self.write(row, col, item)
            col += 1
        row += 1

    def write(self,row,col,item):
        self.worksheets[self.name].write(row,col,item)
        self.max_col[self.name] = max(self.max_col[self.name], col)
        self.max_row[self.name] = max(self.max_row[self.name], row)
        self.positions[self.name] = (row+1,col+1)

    def delete(self):
        for row in range(0,self.max_row[self.name]):
            for col in range(0, self.max_col[self.name]):
                self.worksheets[self.name].write(row,col,'')
        self.positions[self.name] = (-1,0)
        self.max_col[self.name] = 0
        self.max_row[self.name] = 0

    def close(self):
        self.workbook.close()
