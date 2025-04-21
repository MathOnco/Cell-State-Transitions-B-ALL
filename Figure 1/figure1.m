clc;clear;close all;

colors = [123,78,155;
          162,99,169;
          219,156,185;
          242,221,219]/255;


Xlabels={'CD34+ CD38-','CD34+ CD38+','CD34- CD38+','CD34- CD38-'};


filename = '../Data/flow.csv';


DP = getStringRow(filename,'DiseaseProgression');
S = getStringRow(filename,'Specimen');


matrix = getFlow(filename);


% Diagnosis PB
figure(1);
sorted = matrix(:, find(strcmp(DP, 'Diagnosis') & strcmp(S, 'PB')));
violinplot(sorted',Xlabels,'bandwidth',0.1,'ViolinColor',colors, 'DataStyle', 'scatter', 'ShowNotches', false, 'ViolinAlpha',0.5, 'QuartileStyle','boxplot'); ylim([0 1]); 


% Diagnosis BM
figure(2);
sorted = matrix(:, find(strcmp(DP, 'Diagnosis') & strcmp(S, 'BM')));
violinplot(sorted',Xlabels,'bandwidth',0.1,'ViolinColor',colors, 'DataStyle', 'scatter', 'ShowNotches', false, 'ViolinAlpha',0.5, 'QuartileStyle','boxplot'); ylim([0 1]); 

% Relapse BM
figure(3);
sorted = matrix(:, find(strcmp(DP, 'Relapse') & strcmp(S, 'BM')));
violinplot(sorted',Xlabels,'bandwidth',0.1,'ViolinColor',colors, 'DataStyle', 'scatter', 'ShowNotches', false, 'ViolinAlpha',0.5, 'QuartileStyle','boxplot'); ylim([0 1]); 

% Remission BM
figure(4);
sorted = matrix(:, find(strcmp(DP, 'Remission') & strcmp(S, 'BM')));
violinplot(sorted',Xlabels,'bandwidth',0.1,'ViolinColor',colors, 'DataStyle', 'scatter', 'ShowNotches', false, 'ViolinAlpha',0.5, 'QuartileStyle','boxplot'); ylim([0 1]); 

printPNG(figure(1),'../plot/1E.png')
printPNG(figure(2),'../plot/1F.png')
printPNG(figure(3),'../plot/1G.png')
printPNG(figure(4),'../plot/1H.png')




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

function matrix = getFlow(filename)
    % Read rows 21-24 (excluding first column) as numeric matrix
    fid = fopen(filename, 'r');
    for i = 1:19
        tline = fgetl(fid); % Skip first 20 lines
    end
    
    matrix = [];
    for i = 1:4
        tline = fgetl(fid);
        row = str2double(strsplit(tline, ','));
        matrix(i,:) = row(2:end); % Exclude first column
    end
    fclose(fid);
end
