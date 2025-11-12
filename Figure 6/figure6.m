clc;clear;close all;
MS = 45;
lw = 2;
FS = 24;

% DIAGNOSIS_COLOR = [70,153,144]/255;
% RELAPSE_COLOR = [61,118,175]/255;
% REMISSION_COLOR = [165,35,211]/255;

DIAGNOSIS_COLOR = [71, 154, 145]/255;
REMISSION_COLOR = [180, 6, 217]/255;
RELAPSE_COLOR = [32, 121, 180]/255;


colors = [RELAPSE_COLOR;REMISSION_COLOR];

DiseaseProgression2 = {'Relapse','Remission'};

f = 1;
for DP2 = DiseaseProgression2
    


    
    SECOND_COLOR = REMISSION_COLOR;
    
    if (strcmp(DP2,'Relapse'))
        SECOND_COLOR = RELAPSE_COLOR;
    end


    
    filename = '../Data/markov.csv';
    
    AllIDs = getStringRow(filename,'PatientID');
    DP = getStringRow(filename,'DiseaseProgression');
    S = getStringRow(filename,'Specimen');
    
    matrix = getMatrix(filename);
    
    Diagnosis_Ids = AllIDs(find(strcmp(DP, 'Diagnosis') & strcmp(S, 'BM')));
    Diagnosis_Matrixs = matrix(:, find(strcmp(DP, 'Diagnosis') & strcmp(S, 'BM') ));
    
    DP2_Ids = AllIDs(find(strcmp(DP, char(DP2))));
    DP2_Matrixs = matrix(:, find(strcmp(DP, char(DP2)) ));
    
    % Find overlapping patient IDs between diagnosis and relapse
    overlapping_ids = intersect(Diagnosis_Ids, DP2_Ids);
    
    % Find indices of overlapping IDs in Diagnosis and Relapse arrays
    diagnosis_indices = [];
    relapse_indices = [];
    
    j = 1;
    for i = 1:length(overlapping_ids)
    
        array = find(strcmp(Diagnosis_Ids, overlapping_ids{i}));
        diagnosis_indices = [diagnosis_indices,array(1)];
        
        array = find(strcmp(DP2_Ids, overlapping_ids{i}));
        relapse_indices = [relapse_indices,array(1)];    
    end
    
    diagnosis = Diagnosis_Matrixs(:,diagnosis_indices);
    relapse = DP2_Matrixs(:,relapse_indices);
    
    
    figure(f); hold on;
    xvals = 1:1:16;
    delta = 0.45;
    
    for i = 1:1:length(diagnosis_indices)
            
        for j = 1:1:length(xvals)
            plot([xvals(j),xvals(j)+delta],[diagnosis(j,i),relapse(j,i)],'-','Color',black(),'LineWidth',lw);
        end
    end
    for i = 1:1:length(diagnosis_indices)


        MS = 15;

        plot(xvals,diagnosis(:,i),'o','MarkerFaceColor',DIAGNOSIS_COLOR,'MarkerEdgeColor','black','LineWidth',1.8,'MarkerSize',MS);
        plot(xvals+delta,relapse(:,i),'d','MarkerFaceColor',SECOND_COLOR,'MarkerEdgeColor','black','LineWidth',1.8,'MarkerSize',MS);
    


        
    end
    
    Xlabels={'M_{11}','M_{12}','M_{13}','M_{14}','M_{21}','M_{22}','M_{23}','M_{24}','M_{31}','M_{32}','M_{33}','M_{34}','M_{41}','M_{42}','M_{43}','M_{44}'};
    
    
    %% clean up figure
    clean();
    resize(3.5,1);
    setfontsize(FS);
    xlim([0.5,17])
    xticks(xvals);
    xticklabels(Xlabels);
    ylim([0 1.15])
    addSignificance(diagnosis',relapse')
    ylabel('cell fate transition probability');

    f = f + 1;
end





printPNG(figure(1),'../plot/6B.png')
printPNG(figure(2),'../plot/6D.png')


function [] = addSignificance(diagnosis,relapse)    

    h = zeros(1,size(relapse,2));
    for i = 1:1:size(relapse,2)
        [h1,p1] = ttest2(diagnosis(:,i),relapse(:,i)) % ranksum, signrank, ttest2
        h(i) = h1;
    end
    
    for i = 1:1:length(h)
        if (h(i) == 0)
            text(i-0.2,1.05,"N.S.",'FontSize',26);
        else
            text(i-0.2,1.05,"*",'FontSize',32);
        end
    end
end


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

function matrix = getMatrix(filename)
    % Read rows 21-24 (excluding first column) as numeric matrix
    fid = fopen(filename, 'r');
    for i = 1:19
        tline = fgetl(fid); % Skip first 20 lines
    end
    
    matrix = [];
    for i = 1:16
        tline = fgetl(fid);
        row = str2double(strsplit(tline, ','));
        matrix(i,:) = row(2:end); % Exclude first column
    end
    fclose(fid);
end
