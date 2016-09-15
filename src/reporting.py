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
			ws.conditional_formatting.add(_get_column_letter(j) + '2:' + _get_column_letter(j) + str(results.shape[0]+1), rule)

	wb.save(filePath)
