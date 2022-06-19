close all; clear all; clc;

text_path = 'C:\Users\Legion\Documents\MATLAB\EE374_project\datas\Input_file_example1.txt';
text_path_2 = 'C:\Users\Legion\Documents\MATLAB\EE374_project\datas\Input_file_example2.txt';
library_path  = 'C:\Users\Legion\Documents\MATLAB\EE374_project\datas\library.csv';

[R_pu, X_pu, B_pu] = exxxx_p2(text_path, library_path)

function [R_pu, X_pu, B_pu] = exxxx_p2(text_path, library_path)
    format long;
    %%% DATA EXTRACTION FUNCTION %%%
    function [S_base, V_base, number_of_circuit, number_of_bundle, bundle_distance, length, conductor_name, outside_diameter, R_ac, GMR_conductor] = data_extraction(text_path, library_path)
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
        R_dc = lib{[conductor_name],5}; % For a given conductor name, find its corresponding "DC-resistance"
        R_ac = lib{[conductor_name],6}; % For a given conductor name, find its corresponding "AC-resistance"
        GMR_conductor = lib{[conductor_name],7}; % For a given conductor name, find its corresponding "GMR-conductor-length"

        strand = char(lib{[conductor_name],2});                                                                                           
        strand = myDivider(strand); % Send "strand" variable to function, function converts it to double by parsing
    end
    
    %%% STRING DIVISION FUNCTION FOR DATA EXTRACTION FUNCTION %%% 
    function y = myDivider(x) % If we need to use "strand" variable, it converts string data type to double i.e. "21/3" --> 21/3 = 7
        x = strsplit(x,'/');
        y = str2double(x{1})/str2double(x{2});
    end
    
    %%% PER UNIT REACTANCE AND SUSCEPTANCE CALCULATOR FUNCTION %%%
    function [R_pu, X_pu, B_pu] = RXB_pu_calculator(S_base, V_base, X_L, B_C, R_total)
        S_base = S_base/1000; % in terms of KVA 
        V_base = V_base/1000; % in terms of KVA 
        R_pu = (R_total*S_base)/(((V_base)^2)*1000); % ??????????????????????????
        X_pu = (X_L*S_base)/(((V_base)^2)*1000); % calculate X_pu
        B_pu = (B_C*(((V_base)^2)*1000))/S_base; % calculate B_pu
    end
    
    %%% DISTANCE/GMR/GMD CALCULATOR BETWEEN PHASES FOR SINGLE CIRCUIT FUNCTION %%%
    function [L, C] = distance_finder_1_cct(c1_a, c1_b, c1_c, GMR, r_eq, length)
        %%% FIND DISTANCES BETWEEN POINTS %%%
        d_ab = ((c1_a(1)-c1_b(1))^2 + (c1_a(2)-c1_b(2))^2 )^(0.5);
        d_bc = ((c1_b(1)-c1_c(1))^2 + (c1_b(2)-c1_c(2))^2 )^(0.5);
        d_ac = ((c1_a(1)-c1_c(1))^2 + (c1_a(2)-c1_c(2))^2 )^(0.5);

        %%% CALCULATE GMD %%%
        GMD = (d_ab*d_bc*d_ac)^(1/3);

        %%% CALCULATE INDUCTANCE & CAPACITANCE %%%
        L = (2*1e-7)*(log(GMD/GMR))*length; % calculate inductance
        C = ((2*pi*8.85418782*1e-12)/log(GMD/r_eq))*length; % calculate capacitance
    end
    
    %%% DISTANCE/GMR/GMD CALCULATOR BETWEEN PHASES FOR DOUBLE CIRCUIT FUNCTION %%%
    function [L, C] = distance_finder_2_cct(c1_a, c1_b, c1_c, c2_a, c2_b, c2_c, GMR_bundle, r_eq, length)
        %%% FIND DISTANCES BETWEEN POINTS %%%
        d_c1a_c2a = ((c1_a(1)-c2_a(1))^2 + (c1_a(2)-c2_a(2))^2 )^(0.5);
        d_c1b_c2b = ((c1_b(1)-c2_b(1))^2 + (c1_b(2)-c2_b(2))^2 )^(0.5);
        d_c1c_c2c = ((c1_c(1)-c2_c(1))^2 + (c1_c(2)-c2_c(2))^2 )^(0.5);

        d_c1a_c1b = ((c1_a(1)-c1_b(1))^2 + (c1_a(2)-c1_b(2))^2 )^(0.5); % RED LINES
        d_c1a_c2b = ((c1_a(1)-c2_b(1))^2 + (c1_a(2)-c2_b(2))^2 )^(0.5); % RED LINES
        d_c1a_c1c = ((c1_a(1)-c1_c(1))^2 + (c1_a(2)-c1_c(2))^2 )^(0.5); % YELLOW LINES
        d_c1a_c2c = ((c1_a(1)-c2_c(1))^2 + (c1_a(2)-c2_c(2))^2 )^(0.5); % YELLOW LINES

        d_c1b_c1c = ((c1_b(1)-c1_c(1))^2 + (c1_b(2)-c1_c(2))^2 )^(0.5); % PURPLE LINES
        d_c1b_c2c = ((c1_b(1)-c2_c(1))^2 + (c1_b(2)-c2_c(2))^2 )^(0.5); % PURPLE LINES
        d_c1b_c2a = ((c1_b(1)-c2_a(1))^2 + (c1_b(2)-c2_a(2))^2 )^(0.5); % RED LINES

        d_c1c_c2b = ((c1_c(1)-c2_b(1))^2 + (c1_c(2)-c2_b(2))^2 )^(0.5); % PURPLE LINES
        d_c1c_c2a = ((c1_c(1)-c2_a(1))^2 + (c1_c(2)-c2_a(2))^2 )^(0.5); % YELLOW LINES

        d_c2a_c2b = ((c2_a(1)-c2_b(1))^2 + (c2_a(2)-c2_b(2))^2 )^(0.5); % RED LINES
        d_c2a_c2c = ((c2_a(1)-c2_c(1))^2 + (c2_a(2)-c2_c(2))^2 )^(0.5); % YELLOW LINES

        d_c2b_c2c = ((c2_b(1)-c2_c(1))^2 + (c2_b(2)-c2_c(2))^2 )^(0.5); % PURPLE LINES

        %%%%%%%% FIND DISTANCES TO CALCULATE EARTH EFFECT %%%%%%%% 
        %%% _G  means symmetric of the point to ground %%%
        c1_a_G = [c1_a(1), -c1_a(2)]; 
        c1_b_G = [c1_b(1), -c1_b(2)]; 
        c1_c_G = [c1_c(1), -c1_c(2)]; 
        c2_a_G = [c2_a(1), -c2_a(2)];
        c2_b_G = [c2_b(1), -c2_b(2)];
        c2_c_G = [c2_c(1), -c2_c(2)]; 

        %%% CALCULATE TURQUOISE COLOR LINES %%%
        turquoise_a = (((c2_c(1)-c1_c_G(1))^2 + (c2_c(2)-c1_c_G(2))^2 )^(0.5));
        turquoise_b = (((c1_c(1)-c1_c_G(1))^2 + (c1_c(2)-c1_c_G(2))^2 )^(0.5));
        turquoise_c = (((c1_c(1)-c2_c_G(1))^2 + (c1_c(2)-c2_c_G(2))^2 )^(0.5));
        turquoise_d = (((c2_c(1)-c2_c_G(1))^2 + (c2_c(2)-c2_c_G(2))^2 )^(0.5));
        turquoise = (turquoise_a * turquoise_b * turquoise_c * turquoise_d)^(1/4);

        %%% CALCULATE DARK GREEN COLOR LINES %%%
        dark_green_a = (((c1_a(1)-c1_a_G(1))^2 + (c1_a(2)-c1_a_G(2))^2 )^(0.5));
        dark_green_b = (((c1_a(1)-c2_a_G(1))^2 + (c1_a(2)-c2_a_G(2))^2 )^(0.5));
        dark_green_c = (((c2_a(1)-c2_a_G(1))^2 + (c2_a(2)-c2_a_G(2))^2 )^(0.5));
        dark_green_d = (((c2_a(1)-c1_a_G(1))^2 + (c2_a(2)-c1_a_G(2))^2 )^(0.5));
        dark_green = (dark_green_a * dark_green_b * dark_green_c * dark_green_d)^(1/4);

        %%% CALCULATE CREAM COLOR LINES %%%
        cream_a = (((c1_b(1)-c1_b_G(1))^2 + (c1_b(2)-c1_b_G(2))^2 )^(0.5));
        cream_b = (((c1_b(1)-c2_b_G(1))^2 + (c1_b(2)-c2_b_G(2))^2 )^(0.5));
        cream_c = (((c2_b(1)-c1_b_G(1))^2 + (c2_b(2)-c1_b_G(2))^2 )^(0.5));
        cream_d = (((c2_b(1)-c2_b_G(1))^2 + (c2_b(2)-c2_b_G(2))^2 )^(0.5));
        cream = (cream_a * cream_b * cream_c * cream_d)^(1/4);

        %%% CALCULATE YELLOW COLOR LINES %%%
        yellow_a = (((c1_a(1)-c1_c_G(1))^2 + (c1_a(2)-c1_c_G(2))^2 )^(0.5));
        yellow_b = (((c1_a(1)-c2_c_G(1))^2 + (c1_a(2)-c2_c_G(2))^2 )^(0.5));
        yellow_c = (((c2_a(1)-c1_c_G(1))^2 + (c2_a(2)-c1_c_G(2))^2 )^(0.5));
        yellow_d = (((c2_a(1)-c2_c_G(1))^2 + (c2_a(2)-c2_c_G(2))^2 )^(0.5));
        yellow = (yellow_a * yellow_b * yellow_c * yellow_d)^(1/4);

        %%% CALCULATE RED COLOR LINES %%%
        red_a = (((c1_a(1)-c1_b_G(1))^2 + (c1_a(2)-c1_b_G(2))^2 )^(0.5));
        red_b = (((c1_a(1)-c2_b_G(1))^2 + (c1_a(2)-c2_b_G(2))^2 )^(0.5));
        red_c = (((c2_a(1)-c1_b_G(1))^2 + (c2_a(2)-c1_b_G(2))^2 )^(0.5));
        red_d = (((c2_a(1)-c2_b_G(1))^2 + (c2_a(2)-c2_b_G(2))^2 )^(0.5));
        red = (red_a * red_b * red_c * red_d)^(1/4);

        %%% CALCULATE BLUE COLOR LINES %%%
        blue_a = (((c1_b(1)-c1_c_G(1))^2 + (c1_b(2)-c1_c_G(2))^2 )^(0.5));
        blue_b = (((c1_b(1)-c2_c_G(1))^2 + (c1_b(2)-c2_c_G(2))^2 )^(0.5));
        blue_c = (((c2_b(1)-c1_c_G(1))^2 + (c2_b(2)-c1_c_G(2))^2 )^(0.5));
        blue_d = (((c2_b(1)-c2_c_G(1))^2 + (c2_b(2)-c2_c_G(2))^2 )^(0.5));
        blue = (blue_a * blue_b * blue_c * blue_d)^(1/4);

        %%%%%%%% GMR/GMD/R_EQ/CAPACITANCE WITH EARTH EFFECT CALCULATIONS %%%%%%%%
        %%% CALCULATE GMR %%%
        GMR_AA = (d_c1a_c2a * GMR_bundle)^0.5;
        GMR_BB = (d_c1b_c2b * GMR_bundle)^0.5;
        GMR_CC = (d_c1c_c2c * GMR_bundle)^0.5;
        GMR_3_PHASE = (GMR_AA * GMR_BB * GMR_CC)^(1/3);

        %%% CALCULATE R_EQ %%%
        R_EQV_AA = (d_c1a_c2a * r_eq)^0.5; % RED LINE MULTIPLICATIONS
        R_EQV_BB = (d_c1b_c2b * r_eq)^0.5; % YELLOW LINE MULTIPLICATIONS
        R_EQV_CC = (d_c1c_c2c * r_eq)^0.5; % PURPLE LINE MULTIPLICATIONS
        R_EQV_3_PHASE = (R_EQV_AA * R_EQV_BB * R_EQV_CC)^(1/3);

        %%% CALCULATE GMD %%%
        GMD_AB = (d_c1a_c2b * d_c1a_c1b * d_c1b_c2a * d_c2a_c2b)^(1/4); % RED LINE MULTIPLICATIONS
        GMD_AC = (d_c1a_c1c * d_c1a_c2c * d_c1c_c2a * d_c2a_c2c)^(1/4); % YELLOW LINE MULTIPLICATIONS
        GMD_BC = (d_c1b_c1c * d_c1b_c2c * d_c1c_c2b * d_c2b_c2c)^(1/4); % PURPLE LINE MULTIPLICATIONS
        GMD_3_PHASE = (GMD_AB * GMD_AC * GMD_BC)^(1/3);

        %%% CAPACITANCE DUE TO EARTH EFFECT %%%
        num = (yellow * red * blue)^(1/3);
        denum = (turquoise * dark_green * cream)^(1/3);

        %%% CALCULATE INDUCTANCE & CAPACITANCE %%%
        L = (2*1e-7)*(log(GMD_3_PHASE/GMR_3_PHASE))*length; % calculate inductance
        C = ((2*pi*8.85418782*1e-12)/(log(GMD_3_PHASE/R_EQV_3_PHASE)-log(num/denum)))*length;
    end

    %%% BUNDLE GMR & R_EQ CALCULATOR FUNCTION %%%
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
                r_eq = ((((bundle_distance)^3) * conductor_radius)^(1/4))*1.0905;
                GMR_bundle = ((((bundle_distance)^3) * GMR_conductor)^(1/4))*1.0905;
            case 5
                r_eq = ((((bundle_distance)^4) * conductor_radius)^(1/5))*1.21225; 
                GMR_bundle = ((((bundle_distance)^4) * GMR_conductor)^(1/5))*1.21225;
            case 6
                r_eq = ((((bundle_distance)^5) * conductor_radius)^(1/6))*1.34800;
                GMR_bundle = ((((bundle_distance)^5) * GMR_conductor)^(1/6))*1.34800;
            case 7
                r_eq = ((((bundle_distance)^6) * conductor_radius)^(1/7))*1.49116;
                GMR_bundle = ((((bundle_distance)^6) * GMR_conductor)^(1/7))*1.49116;    
            case 8
                r_eq = ((((bundle_distance)^7) * conductor_radius)^(1/8))*1.638703;
                GMR_bundle = ((((bundle_distance)^7) * GMR_conductor)^(1/8))*1.638703;
        end
    end
                                   %%% MAIN PART OF PHASE 2 %%%
    %%% FUNCTIONS DEFINED AS NESTED NOT TO REPEAT SAME PROCESSES AGAIN AND AGAIN %%%
    
    [S_base, V_base, number_of_circuit, number_of_bundle, bundle_distance, length, conductor_name, outside_diameter, R_ac, GMR_conductor] = data_extraction(text_path, library_path);
    
    frequency = 50; % frequency is given as 50Hz
    conductor_radius = outside_diameter/2.0000; % Divide outside diameter by 2 which gives the radius
    [GMR_bundle, r_eq] = bundle_GMR_calculator(conductor_radius, GMR_conductor, bundle_distance, number_of_bundle); % find your GMR of the bundle in another function to repeat many number of possibilities
    lines_of_txt_file = strsplit(fileread(text_path), {'\r', '\n'}); % does not depend on number of circuit, hence extract outside the if blocks.
    
    if (number_of_circuit == 1)
        c1_c = [str2double((lines_of_txt_file{16}.')),str2double((lines_of_txt_file{17}.'))]; % extract position of circuit 1 phase c
        c1_a = [str2double((lines_of_txt_file{19}.')),str2double((lines_of_txt_file{20}.'))]; % extract position of circuit 1 phase a
        c1_b = [str2double((lines_of_txt_file{22}.')),str2double((lines_of_txt_file{23}.'))]; % extract position of circuit 1 phase b
        
        [L,C] = distance_finder_1_cct(c1_a, c1_b, c1_c, GMR_bundle, r_eq, length); % calculate distance between phases according to phase positions
        
        X_L = 2*pi*frequency*L; % calculate total reactance
        B_C = 2*pi*frequency*C; % calculate total susceptance
        R_total = (R_ac)*length/number_of_bundle; % calculate total resistance
        
        [R_pu, X_pu, B_pu] = RXB_pu_calculator(S_base, V_base, X_L, B_C, R_total); % calculate necessary parameters in per unit system

    elseif (number_of_circuit == 2)
        c1_c = [str2double((lines_of_txt_file{16}.')),str2double((lines_of_txt_file{17}.'))]; % extract position of circuit 1 phase c
        c1_a = [str2double((lines_of_txt_file{19}.')),str2double((lines_of_txt_file{20}.'))]; % extract position of circuit 1 phase a
        c1_b = [str2double((lines_of_txt_file{22}.')),str2double((lines_of_txt_file{23}.'))]; % extract position of circuit 1 phase b
        c2_c = [str2double((lines_of_txt_file{25}.')),str2double((lines_of_txt_file{26}.'))]; % extract position of circuit 2 phase c
        c2_a = [str2double((lines_of_txt_file{28}.')),str2double((lines_of_txt_file{29}.'))]; % extract position of circuit 2 phase a
        c2_b = [str2double((lines_of_txt_file{31}.')),str2double((lines_of_txt_file{32}.'))]; % extract position of circuit 2 phase b
        
        [L, C] = distance_finder_2_cct(c1_a, c1_b, c1_c, c2_a, c2_b, c2_c, GMR_bundle, r_eq, length);
    
        X_L = 2*pi*frequency*L; % calculate total reactance
        B_C = 2*pi*frequency*C; % calculate total susceptance
        R_total = (R_ac)*length/(number_of_bundle*number_of_circuit); % calculate total resistance
        
        [R_pu, X_pu, B_pu] = RXB_pu_calculator(S_base, V_base, X_L, B_C, R_total); % calculate necessary parameters in per unit system
    end
end
