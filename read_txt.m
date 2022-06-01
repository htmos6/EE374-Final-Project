opts = detectImportOptions("C:\Users\Legion\Documents\MATLAB\EE374_project\datas\library.csv");
opts.VariableTypes = "string";
R = readtable("datas\library.csv", opts)
