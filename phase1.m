close all; clear all; clc;

text_path = 'C:\Users\Legion\Documents\MATLAB\EE374_project\datas\Input_file_example1.txt';
text_pathh = 'C:\Users\Legion\Documents\MATLAB\EE374_project\datas\Input_file_example2.txt';
library_path  = 'C:\Users\Legion\Documents\MATLAB\EE374_project\datas\library.csv';
[S_base, V_base, number_of_circuit, number_of_bundle, bundle_distance, length, conductor_name, outside_diameter, R_ac, GMR_conductor] = e237441_p1(text_path, library_path)


function [S_base, V_base, number_of_circuit, number_of_bundle, bundle_distance, length, conductor_name, outside_diameter, R_ac, GMR_conductor] = e237441_p1(text_path, library_path)
lib = readtable(library_path);                                                                  % read excel file as a table

lib.Properties.VariableNames([2 3 4 5 6 7 8 ]) = {'Aluminum_Area_m2' 'Strand' 'Layers_Of_Aluminum' 'Outside_Diameter_m' 'DC_Resistance_20C_ohm_per_m' 'AC_50Hz_Resistance_20C_ohm_per_m' 'GMR_m'};

lib.Aluminum_Area_m2 = (5.06707479*1e-10).*(lib.Aluminum_Area_m2);                              % convert cmil to m^2
lib.Outside_Diameter_m = (lib.Outside_Diameter_m).*(0.0254);                                    % convert inch to m
lib.DC_Resistance_20C_ohm_per_m = (lib.DC_Resistance_20C_ohm_per_m).*(0.0032808399);            % convert inch to ohm/(1000*ft) to ohm/m
lib.AC_50Hz_Resistance_20C_ohm_per_m = (lib.AC_50Hz_Resistance_20C_ohm_per_m).*(39370.0787);    % convert inch to ohm/(1000*mil) to ohm/m
lib.GMR_m = lib.GMR_m.*(0.3048);                                                                % convert ft to m

lines = strsplit(fileread(text_path), {'\r', '\n'});                                            % Read your input from ".txt" files

S_base = str2double((lines{2}.'))*1e6;                                                          % From ".txt" file, extract variable "Complex Power"
V_base = str2double((lines{4}.'))*1e3;                                                          % From ".txt" file, extract variable "Voltage"
number_of_circuit = str2double((lines{6}.'));                                                   % From ".txt" file, extract variable "Number of circuits"
number_of_bundle = str2double((lines{8}.'));                                                    % From ".txt" file, extract variable "Number of bundles"
bundle_distance = str2double((lines{10}.'));                                                    % From ".txt" file, extract variable "Bundle Distance"
length = str2double((lines{12}.'))*1e3;                                                         % From ".txt" file, extract variable "Length"
conductor_name = lines{14};

%%%%%%% Define your new row names %%%%%%%
code_words = ["Waxwing";"Ostrich";"Linnet";"Ibis";"Hawk";"Dove";"Rook";"Drake";"Rail";"Cardinal";"Bluejay";"Pheasant";"Plover";"Falcon";"Bluebird"];
                       
lib.Properties.RowNames = code_words;                                                           % Assign your new row names to the table
lib.CodeWord=[];                                                                                % delete codeword column to find necessary row easily

writetable(lib,'lib_new.csv','WriteRowNames',true,'Delimiter',' ');  

outside_diameter = lib{conductor_name,4};                                                       % For a given conductor name, find its corresponding "outside-diameter"
R_ac = lib{[conductor_name],6};                                                                 % For a given conductor name, find its corresponding "AC-resistance"
GMR_conductor = lib{[conductor_name],7};                                                        % For a given conductor name, find its corresponding "GMR-conductor-length"

strand = char(lib{[conductor_name],2});                                                                                           
strand = myDivider(strand);                                                                     % Send "strand" variable to function, function converts it to double by parsing
end


function y = myDivider(x)                                                                       % If we need to use "strand" variable, it converts string data type to double i.e. "21/3" --> 21/3 = 7
    x = strsplit(x,'/');
    y = str2double(x{1})/str2double(x{2});
end
