% To do list:
%     Under if 1:
%         Calculate
%             R_pu 
%             B_pu 
%     Under if 2:
%         Calculate
%             R_pu 
%             B_pu 
%             X_pu

close all; clear all; clc;

text_path = 'C:\Users\Legion\Documents\MATLAB\EE374_project\datas\Input_file_example1.txt';
text_path_2 = 'C:\Users\Legion\Documents\MATLAB\EE374_project\datas\Input_file_example2.txt';
library_path  = 'C:\Users\Legion\Documents\MATLAB\EE374_project\datas\library.csv';

[S_base, V_base, number_of_circuit, number_of_bundle, bundle_distance, length, conductor_name, outside_diameter, R_ac, GMR_conductor] = e237441_p1(text_path, library_path); % take necessary parameters from phase 1 function in order not to repeat same code again 

[R_pu, X_pu, B_pu] = e237441_p2(text_path_2, library_path);


function [R_pu, X_pu, B_pu] = e237441_p2(text_path, library_path)

    [S_base, V_base, number_of_circuit, number_of_bundle, bundle_distance, length, conductor_name, outside_diameter, R_ac, GMR_conductor] = e237441_p1(text_path, library_path);
    frequency = 50; % frequency is given as 50Hz
    conductor_radius = outside_diameter/2.0000; % Divide outside diameter by 2 which gives the radius
    [GMR_bundle, r_eq] = bundle_GMR_calculator(conductor_radius, GMR_conductor, bundle_distance, number_of_bundle); % find your GMR of the bundle in another function to repeat many number of possibilities

    lines = strsplit(fileread(text_path), {'\r', '\n'}); % read your text file to collect necessary position informations
    if (number_of_circuit == 1)
        c1_c = [str2double((lines{16}.')),str2double((lines{17}.'))]; % extract position of circuit 1 phase c
        c1_a = [str2double((lines{19}.')),str2double((lines{20}.'))]; % extract position of circuit 1 phase a
        c1_b = [str2double((lines{22}.')),str2double((lines{23}.'))]; % extract position of circuit 1 phase b
       
        [d_ab, d_bc, d_ac] = distance_finder(c1_a, c1_b, c1_c); % calculate distance between phases according to phase positions
        L = L_calculator(GMR_bundle, d_ab, d_bc, d_ac)*length; % calculate total inductance by multiplying length
        X_L = 2*pi*frequency*L; % calculate total reactance
        X_pu = X_pu_calculator(S_base, V_base, X_L); % calculate total reactance in per unit system
        
    elseif (number_of_circuit == 2)
        c1_c = [str2double((lines{16}.')),str2double((lines{17}.'))]; % extract position of circuit 1 phase c
        c1_a = [str2double((lines{19}.')),str2double((lines{20}.'))]; % extract position of circuit 1 phase a
        c1_b = [str2double((lines{22}.')),str2double((lines{23}.'))]; % extract position of circuit 1 phase b
        c2_c = [str2double((lines{25}.')),str2double((lines{26}.'))]; % extract position of circuit 2 phase c
        c2_a = [str2double((lines{28}.')),str2double((lines{29}.'))]; % extract position of circuit 2 phase a
        c2_b = [str2double((lines{31}.')),str2double((lines{32}.'))]; % extract position of circuit 2 phase b
    end
end

%%% PER UNIT REACTANCE CALCULATOR %%%
function X_pu = X_pu_calculator(S_base, V_base, X_L)
    S_base = S_base/1000; % in terms of KVA 
    V_base = V_base/1000; % in terms of KVA 
    X_pu = (X_L*S_base)/(((V_base)^2)*1000); % calculate X_pu
end

%%% INDUCTANCE CALCULATOR %%%
function L = L_calculator(GMR, d_ab, d_bc, d_ac)
    GMD = (d_ab*d_bc*d_ac)^(1/3); % calculate GMD
    L = (2*1e-7)*(log(GMD/GMR)); % calculate inductance
end

%%% DISTANCE CALCULATOR BETWEEN PHASES %%%
function [d_ab, d_bc, d_ac] = distance_finder(c1_a, c1_b, c1_c)
    d_ab = ((c1_a(1)-c1_b(1))^2 + (c1_a(2)-c1_b(2))^2 )^(0.5);
    d_bc = ((c1_b(1)-c1_c(1))^2 + (c1_b(2)-c1_c(2))^2 )^(0.5);
    d_ac = ((c1_a(1)-c1_c(1))^2 + (c1_a(2)-c1_c(2))^2 )^(0.5);
end

%%% BUNDLE GMR CALCULATOR %%%
function [GMR_bundle, r_eq] = bundle_GMR_calculator(conductor_radius, GMR_conductor, bundle_distance, number_of_bundle)
    switch number_of_bundle
        case 1
            r_eq = conductor_radius;
            GMR_bundle = GMR_conductor;
        case 2 
            r_eq = (bundle_distance * conductor_radius)^(1/2);  
            GMR_bundle = (bundle_distance * GMR_conductor)^(1/2);
        case 3 
            r_eq = (((bundle_distance)^2) * conductor_radius)^(1/3);
            GMR_bundle = (((bundle_distance)^2) * GMR_conductor)^(1/3);
        case 4
            r_eq = ((((bundle_distance)^3) * conductor_radius)^(1/4))*1.09;
            GMR_bundle = ((((bundle_distance)^3) * GMR_conductor)^(1/4))*1.09;
        case 5
            r_eq = ((((bundle_distance)^4) * conductor_radius)^(1/5))*1.46956;
            GMR_bundle = ((((bundle_distance)^4) * GMR_conductor)^(1/5))*1.46956;
        case 6
            r_eq = ((((bundle_distance)^5) * conductor_radius)^(1/6))*1.34800;
            GMR_bundle = ((((bundle_distance)^5) * GMR_conductor)^(1/6))*1.34800;
        case 8
            r_eq = ((((bundle_distance)^7) * conductor_radius)^(1/8))*1.638703;
            GMR_bundle = ((((bundle_distance)^7) * GMR_conductor)^(1/8))*1.638703;
    end
end

%%% PHASE 1 FUNCTION %%%
function [S_base, V_base, number_of_circuit, number_of_bundle, bundle_distance, length, conductor_name, outside_diameter, R_ac, GMR_conductor] = e237441_p1(text_path, library_path)
    lib = readtable(library_path); % read excel file as a table

    lib.Properties.VariableNames([2 3 4 5 6 7 8 ]) = {'Aluminum_Area_m2' 'Strand' 'Layers_Of_Aluminum' 'Outside_Diameter_m' 'DC_Resistance_20C_ohm_per_m' 'AC_50Hz_Resistance_20C_ohm_per_m' 'GMR_m'};

    lib.Aluminum_Area_m2 = (5.06707479*1e-10).*(lib.Aluminum_Area_m2); % convert cmil to m^2
    lib.Outside_Diameter_m = (lib.Outside_Diameter_m).*(0.0254); % convert inch to m
    lib.DC_Resistance_20C_ohm_per_m = (lib.DC_Resistance_20C_ohm_per_m).*(0.0032808399); % convert inch to ohm/(1000*ft) to ohm/m
    lib.AC_50Hz_Resistance_20C_ohm_per_m = (lib.AC_50Hz_Resistance_20C_ohm_per_m).*(0.000621371192);% convert inch to ohm/(1000*mil) to ohm/m
    lib.GMR_m = lib.GMR_m.*(0.3048); % convert ft to m

    lines = strsplit(fileread(text_path), {'\r', '\n'}); % Read your input from ".txt" files

    S_base = str2double((lines{2}.'))*1e6; % From ".txt" file, extract variable "Complex Power"
    V_base = str2double((lines{4}.'))*1e3; % From ".txt" file, extract variable "Voltage"
    number_of_circuit = str2double((lines{6}.')); % From ".txt" file, extract variable "Number of circuits"
    number_of_bundle = str2double((lines{8}.')); % From ".txt" file, extract variable "Number of bundles"
    bundle_distance = str2double((lines{10}.')); % From ".txt" file, extract variable "Bundle Distance"
    length = str2double((lines{12}.'))*1e3; % From ".txt" file, extract variable "Length"
    conductor_name = lines{14};
    
    code_words = ["Waxwing";"Ostrich";"Linnet";"Ibis";"Hawk";"Dove";"Rook";"Drake";"Rail";"Cardinal";"Bluejay";"Pheasant";"Plover";"Falcon";"Bluebird"];  %%%%%%% Define your new row names %%%%%%%

    lib.Properties.RowNames = code_words; % Assign your new row names to the table
    lib.CodeWord=[]; % delete codeword column to find necessary row easily

    writetable(lib,'lib_new.csv','WriteRowNames',true,'Delimiter',' ')  

    outside_diameter = lib{conductor_name,4}; % For a given conductor name, find its corresponding "outside-diameter"
    R_ac = lib{[conductor_name],6}; % For a given conductor name, find its corresponding "AC-resistance"
    GMR_conductor = lib{[conductor_name],7}; % For a given conductor name, find its corresponding "GMR-conductor-length"

    strand = char(lib{[conductor_name],2});                                                                                           
    strand = myDivider(strand); % Send "strand" variable to function, function converts it to double by parsing
end

%%% STRING DIVISION FUNCTION FOR PHASE 1 %%% 
function y = myDivider(x) % If we need to use "strand" variable, it converts string data type to double i.e. "21/3" --> 21/3 = 7
    x = strsplit(x,'/');
    y = str2double(x{1})/str2double(x{2});
end
