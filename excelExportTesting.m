function excelExportTesting

exc = actxserver('Excel.Application');

exc.Visible = 1;

eWorkbook = exc.Workbooks.Add(-4167);
eWorkbook.WorkSheets.Add;

% eSheets = exc.ActiveWorkbook.Sheets;
eSheets = eWorkbook.Sheets;
eSheets.Item(1).Name = 'ephys';
eSheets.Item(2).Name = 'imaging';
% eSheet1 = eSheets.get('Item',1);
eSheet1 = eSheets.get('Item','ephys');
eSheet1.Activate

eWorkbook.Worksheets.Count

A = rand(10,10);
% range = exc.Application.ConvertFormula('=SUM(R1C1:R10C10)',-4150,1)

% eActivesheetRange = get(exc.Activesheet,'Range','A1:B2');
% eActivesheetRange.Value = A;
cellsRange1 = get(eSheet1,'Cells',1,1);
cellsRange2 = get(eSheet1,'Cells',10,10);
range = get(eSheet1,'Range',cellsRange1,cellsRange2);
range.Value = A;

fig1 = figure;
ax1=axes(fig1);
plot(ax1,1:10,exp(1:10))

% https://www.mathworks.com/matlabcentral/fileexchange/24424-xlswritefig
% ezt használtam alapul, akár át is lehetne venni
hgexport(fig1,'-clipboard')
Paste(exc.Activesheet,get(exc.Activesheet,'Range','L1','O4'))



currPath = cd
% SaveAs(eWorkbook,[currPath,'\myfile.xlsx'])
eWorkbook.SaveAs([currPath,'\myfile2.xlsx'])

eWorkbook.Saved = 1;
Close(eWorkbook)

Quit(exc)
delete(exc)