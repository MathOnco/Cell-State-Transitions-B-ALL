clc;clear;close all;

PBcolor = [195,253,196]/255;
BMcolor = [239,183,182]/255;




%% read in flow data
filename = '../Data/markov.csv';

AllIDs = getStringRow(filename,'PatientID');
DP = getStringRow(filename,'DiseaseProgression');
S = getStringRow(filename,'Specimen');


matrix = getMarkov(filename);
Xlabels={'M_{11}','M_{12}','M_{13}','M_{14}','M_{21}','M_{22}','M_{23}','M_{24}','M_{31}','M_{32}','M_{33}','M_{34}','M_{41}','M_{42}','M_{43}','M_{44}'};


BM_Ids = AllIDs(find(strcmp(DP, 'Diagnosis') & strcmp(S, 'BM')));
BM_Matrixs = matrix(:, find(strcmp(DP, 'Diagnosis') & strcmp(S, 'BM') ));

PB_Ids = AllIDs(find(strcmp(DP, 'Diagnosis') & strcmp(S, 'PB')));
PB_Matrixs = matrix(:, find(strcmp(DP, 'Diagnosis') & strcmp(S, 'PB') ));

% Find overlapping patient IDs between diagnosis and relapse
overlapping_ids = intersect(BM_Ids, PB_Ids);

% Find indices of overlapping IDs in Diagnosis and Relapse arrays
BM_indices = [];
PB_indices = [];

j = 1;
for i = 1:length(overlapping_ids)

    array = find(strcmp(BM_Ids, overlapping_ids{i}));
    BM_indices = [BM_indices,array(1)];
    
    array = find(strcmp(PB_Ids, overlapping_ids{i}));
    PB_indices = [PB_indices,array(1)];    
end

MatricesBM = BM_Matrixs(:,BM_indices);
MatricesPB = PB_Matrixs(:,PB_indices);



    violinplot(MatricesPB',Xlabels,'ViolinColor',PBcolor,...
        'DataStyle', 'scatter',... % scatter, none
        'ShowNotches', false,...
        'ViolinAlpha',0.75,...
        'HalfViolin','left',...% left, full
        'Bandwidth', 0.25, ...
        'Width', 0.4, ...
        'QuartileStyle','boxplot'); 

    violinplot(MatricesBM',Xlabels,'ViolinColor',BMcolor,...
        'DataStyle', 'scatter',... % scatter, none
        'ShowNotches', false,...
        'ViolinAlpha',0.5,...
        'HalfViolin','right',...% left, full
        'Bandwidth', 0.25, ...
        'Width', 0.4, ...
        'QuartileStyle','boxplot'); 

    ylim([0 1.15]); 


    hold on;
    clean();

    ylabel('cell fate transition probability')

    h = ttest2(MatricesPB',MatricesBM');

    for i = 1:1:length(h)
        if (h(i) == 0)
            text(i-0.2,1.05,"N.S.",'FontSize',18);
        else
            text(i-0.2,1.05,"*",'FontSize',18);
        end
    end

xlim([0,17])
resize(3.5,1)
setfontsize(28);

printPNG(figure(1),'../plot/S1.png')






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
