import xlsxwriter as xlsx
class SpreadSheetManager():
    def __init__(self, filename):
        # Manage a spread sheet for PDielec / PDGui
        self.workbook = xlsx.Workbook(filename)
        self.closed = False
        self.tab_names = ['Main', 'Settings', 'Scenarios',
                'Powder Molar Absorption (cells)',
                'Powder Molar Absorption (atoms)',
                'Powder Molar Absorption (mols)',
                'Powder Absorption',
                'Powder Real Permittivity',
                'Powder Imaginary Permittivity',
                'Powder ATR Reflectance',
                'Analysis',
                'Crystal R_p',
                'Crystal R_s',
                'Crystal T_p',
                'Crystal T_s',
                'Real Crystal Permittivity',
                'Imag Crystal Permittivity' ]
        self.worksheets = {}
        # Positions points to where we write to next
        self.positions  = {}
        self.max_col    = {}
        self.max_row    = {}
        self.opened     = {}
        for tab in self.tab_names:
            self.opened[tab] = False
        self.name = 'Main'
        self.openWorkSheet(self.name)

    def openWorkSheet(self,tab):
        if self.opened[tab]:
            return
        self.worksheets[tab] = self.workbook.add_worksheet(tab)
        self.positions[tab] = (0,0)
        self.max_col[tab] = 0
        self.max_row[tab] = 0
        self.opened[tab] = True

    def selectWorkSheet(self,name):
        self.name = name
        if not self.opened[name]:
            self.openWorkSheet(self.name)

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
        if not self.closed:
            self.workbook.close()
        self.closed = True
