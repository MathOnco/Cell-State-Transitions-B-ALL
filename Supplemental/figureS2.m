clc;clear;close all;

four_colors = [123,78,155;
          162,99,169;
          219,156,185;
          242,221,219]/255;

colors = [four_colors;four_colors;four_colors;four_colors];

%% read in flow data
filename = '../Data/markov.csv';

AllIDs = getStringRow(filename,'PatientID');
DP = getStringRow(filename,'DiseaseProgression');
S = getStringRow(filename,'Specimen');
B = getStringRow(filename,'BCR::ABL1/like')


matrix = getMarkov(filename);
Xlabels={'M_{11}','M_{12}','M_{13}','M_{14}','M_{21}','M_{22}','M_{23}','M_{24}','M_{31}','M_{32}','M_{33}','M_{34}','M_{41}','M_{42}','M_{43}','M_{44}'};


BCRpositive_Matrixs = matrix(:, find(strcmp(DP, 'Diagnosis') & (strcmp(B, 'BCR::ABL1-like') | strcmp(B, 'BCR::ABL1') ) ));
BCRnegative_Matrixs = matrix(:, find(strcmp(DP, 'Diagnosis') & (strcmp(B, 'Negative')  ) ));
remission_Matrixs = matrix(:, find(strcmp(DP, 'Remission')    ));

figure(1);
    violinplot(BCRpositive_Matrixs',Xlabels,'ViolinColor',colors,...
        'DataStyle', 'scatter',... % scatter, none
        'ShowNotches', false,...
        'ViolinAlpha',0.75,...
        'HalfViolin','full',...% left, full
        'Bandwidth', 0.25, ...
        'Width', 0.4, ...
        'QuartileStyle','boxplot'); 


    ylim([0 1.15]); xlim([0,17]); clean(); resize(3.5,1); hold on; setfontsize(28);
    ylabel('cell fate transition probability')

figure(2);
    violinplot(BCRnegative_Matrixs',Xlabels,'ViolinColor',colors,...
        'DataStyle', 'scatter',... % scatter, none
        'ShowNotches', false,...
        'ViolinAlpha',0.75,...
        'HalfViolin','full',...% left, full
        'Bandwidth', 0.25, ...
        'Width', 0.4, ...
        'QuartileStyle','boxplot'); 


    ylim([0 1.15]); xlim([0,17]); clean(); resize(3.5,1); hold on; setfontsize(28);
    ylabel('cell fate transition probability')

figure(3);
    violinplot(remission_Matrixs',Xlabels,'ViolinColor',colors,...
        'DataStyle', 'scatter',... % scatter, none
        'ShowNotches', false,...
        'ViolinAlpha',0.75,...
        'HalfViolin','full',...% left, full
        'Bandwidth', 0.25, ...
        'Width', 0.4, ...
        'QuartileStyle','boxplot'); 


    ylim([0 1.15]); xlim([0,17]); clean(); resize(3.5,1); hold on; setfontsize(28);
    ylabel('cell fate transition probability')


printPNG(figure(1),'../plot/S2A.png')
printPNG(figure(2),'../plot/S2B.png')
printPNG(figure(3),'../plot/S2C.png')






function idx = getRowIndex(filename, target_string)
    % Open file and read first column
    fid = fopen(filename, 'r');
    idx = 1;
    while ~feof(fid)
        line = fgetl(fid);
        if ~isempty(line)
            parts = strsplit(line, ',');
            if strcmp(parts{1}, target_string)
                fclose(fid);
                return;
            end
        end
        idx = idx + 1;
    end
    fclose(fid);
    idx = -1; % Return -1 if string not found
end

function row = getStringRow(filename,Variable)
    idx = getRowIndex(filename,Variable);
    % Read the idx row (excluding first value) as a string array
    fid = fopen(filename, 'r');
    for i = 1:idx
        tline = fgetl(fid);
    end
    row = strsplit(tline, ',');
    row = row(2:end); % Exclude first value
    fclose(fid);
end

function matrix = getMarkov(filename)
    % Read rows 21-24 (excluding first column) as numeric matrix
    fid = fopen(filename, 'r');
    for i = 1:19
        tline = fgetl(fid); % Skip first 19 lines
    end
    
    matrix = [];
    for i = 1:16
        tline = fgetl(fid);
        row = str2double(strsplit(tline, ','));
        matrix(i,:) = row(2:end); % Exclude first column
    end
    fclose(fid);
end
