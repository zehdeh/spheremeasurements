from openpyxl import Workbook
from openpyxl.formatting.rule import ColorScaleRule
from openpyxl.utils import (_get_column_letter)

def writeReport(filePath, measures, results):
	wb = Workbook()
	ws = wb.active

	rule = ColorScaleRule(start_type='percentile', start_value=10, start_color='00FF00', mid_type='percentile', mid_value=50, mid_color='FFFF00', end_type='percentile', end_value=90, end_color='FF0000')

	for i,measure in enumerate(measures):
		c = ws.cell(row = 1, column = i+1)
		c.value = measure.name
		if measure.color:
			ws.conditional_formatting.add(_get_column_letter(i+1) + '2:' + _get_column_letter(i+1) + str(len(results)+1), rule)
	
	for i, result in enumerate(results):
		for j, val in enumerate(result):
			c = ws.cell(row = 2+i, column = j+1)
			try:
				c.value = val
			except ValueError:
				c.value = float(val)

	wb.save(filePath)
