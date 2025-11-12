clc;clear;close all;
gray = black()+(220/255);
RESIZE = 2.8;
FS = 30;
BANDWIDTH = 0.1;
MRDpositivecolor = [227,206,169]/255;
MRDnegativecolor = [154,185,186]/255;
RELAPSE_POSITIVE = [58,104,66]/255;
RELAPSE_NEGATIVE = [191,114,112]/255;

filename = '../Data/markov.csv';
DP = getStringRow(filename,'DiseaseProgression');
B = getStringRow(filename,'BCR::ABL1/like');
MRDflow = getStringRow(filename,'MRD(flow)');
MRDmolec = getStringRow(filename,'MRD(molecular)');
Relapse3 = getStringRow(filename,'3YrRelapse');


matrix = getMarkov(filename);

M11 = matrix(1,:);

[M11,indices] = sort(M11,'descend');

B = B(indices);
DP = DP(indices);
MRDflow = MRDflow(indices);
MRDmolec = MRDmolec(indices);
Relapse3 = Relapse3(indices);

%% create waterfall for BCR
figure(1);
i = 1;
violin_BCR = [];
for row = 1:1:length(B)
    nan_row = zeros(1,3) + NaN; % labels={'-','+','RM'};
    j = -1;
    if strcmp(DP{row},'Remission')
        color = gray;
        j = 3;
    elseif strcmp(B{row},'BCR::ABL1')
        color = blue(); % blue is positive
        j = 2;
    elseif strcmp(B{row},'BCR::ABL1-like')
        color = [166, 209, 229]/255; % light blue
        j = 2;
    else
        color = red(); % red is negative
        j = 1;
    end

    bar(i,M11(row),'FaceColor',color,'EdgeColor',color); hold on;

    % used for violins later:
    nan_row(j) = M11(row);
    violin_BCR = [violin_BCR; nan_row];

        i = i+1;
end
clean(); xticks([]); resize(RESIZE,1); ylim([0 1.2]);

%% create waterfall for MRD
figure(2);
i = 1;
violin_MRD = [];
for row = 1:1:length(B)

    if strcmp(MRDflow{row},'NaN')

    else
        nan_row = zeros(1,3) + NaN; % labels={'-','+','RM'};
        j = -1;
        if strcmp(DP{row},'Remission')
            color = gray;
            j = 3;
        elseif strcmp(MRDflow{row},'Positive (> 10^-3)')
            color = MRDpositivecolor;
            j = 2;
        else
            color = MRDnegativecolor;
            j = 1;
        end
    
        bar(i,M11(row),'FaceColor',color,'EdgeColor',color); hold on;
    
        % used for violins later:
        nan_row(j) = M11(row);
        violin_MRD = [violin_MRD; nan_row];

        i = i+1;

    end

        
end
clean(); xticks([]); resize(RESIZE,1); ylim([0 1.2]);


%% create waterfall for RELAPSE
figure(3);
i = 1;
violin_RELAPSE = [];
for row = 1:1:length(B)

    if strcmp(Relapse3{row},'NaN')

    else
        nan_row = zeros(1,3) + NaN; % labels={'-','+','RM'};
        j = -1;
        if strcmp(DP{row},'Remission')
            color = gray;
            j = 3;
        elseif strcmp(Relapse3{row},'Yes')
            color = RELAPSE_POSITIVE;
            j = 2;
        else
            color = RELAPSE_NEGATIVE;
            j = 1;
        end

        bar(i,M11(row),'FaceColor',color,'EdgeColor',color); hold on;

        % used for violins later:
        nan_row(j) = M11(row);
        violin_RELAPSE = [violin_RELAPSE; nan_row];

        i = i+1;

    end


end
clean(); xticks([]); resize(RESIZE,1); ylim([0 1.2]);





    
% 
    %% VIOLIN
    colors = [red();blue();gray];
    labels={'-','+','RM'};

    figure(11); hold on;
    violinplot(violin_BCR,labels,'ViolinColor',colors,'DataStyle', 'scatter','ShowNotches', false,'ViolinAlpha',0.75,'HalfViolin','full','Bandwidth', BANDWIDTH,'Width', 0.2,'QuartileStyle','boxplot'); 
    xlim([0,4]);
    setfontsize(FS);
    ylim([0 1.2]); 

    [h,p] = addSignificance(violin_BCR(:,1),violin_BCR(:,2), 1, 2);
    [h,p] = addSignificance(violin_BCR(:,1),violin_BCR(:,3), 1, 3);
    [h,p] = addSignificance(violin_BCR(:,2),violin_BCR(:,3), 2, 3);


    %% VIOLIN
    colors = [MRDnegativecolor;MRDpositivecolor;gray];
    labels={'-','+','RM'};

    figure(12); hold on;
    violinplot(violin_MRD,labels,'ViolinColor',colors,'DataStyle', 'scatter','ShowNotches', false,'ViolinAlpha',0.75,'HalfViolin','full','Bandwidth', BANDWIDTH,'Width', 0.2,'QuartileStyle','boxplot'); 
    xlim([0,4]);
    setfontsize(FS);
    ylim([0 1.2]); 

    [h,p] = addSignificance(violin_MRD(:,1),violin_MRD(:,2), 1, 2);
    [h,p] = addSignificance(violin_MRD(:,1),violin_MRD(:,3), 1, 3);
    [h,p] = addSignificance(violin_MRD(:,2),violin_MRD(:,3), 2, 3);


    %% VIOLIN
    colors = [RELAPSE_NEGATIVE;RELAPSE_POSITIVE;gray];
    labels={'-','+','RM'};

    figure(13); hold on;
    violinplot(violin_RELAPSE,labels,'ViolinColor',colors,'DataStyle', 'scatter','ShowNotches', false,'ViolinAlpha',0.75,'HalfViolin','full','Bandwidth', BANDWIDTH,'Width', 0.2,'QuartileStyle','boxplot'); 
    xlim([0,4]);
    setfontsize(FS);
    ylim([0 1.2]); 

    [h,p] = addSignificance(violin_RELAPSE(:,1),violin_RELAPSE(:,2), 1, 2);
    [h,p] = addSignificance(violin_RELAPSE(:,1),violin_RELAPSE(:,3), 1, 3);
    [h,p] = addSignificance(violin_RELAPSE(:,2),violin_RELAPSE(:,3), 2, 3);



printPNG(figure(1),'../plot/5B.png')
printPNG(figure(2),'../plot/5E.png')
printPNG(figure(3),'../plot/5H.png')

printPNG(figure(11),'../plot/5C.png')
printPNG(figure(12),'../plot/5F.png')
printPNG(figure(13),'../plot/5I.png')







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%











function [h,p] = addSignificance(col1,col2, x1, x2)    

    % col1 = col1(bool_vec',:)
    % col2 = col2(bool_vec',:)

    [h,p] = ttest2(col1,col2);
    
    y0 = 1.05;
    dy = (x2-x1)*.07 - 0.08;
    dx = 0.02;
    plot([x1+dx,x2-dx],[y0+dy,y0+dy],'-k','LineWidth',2.5); hold on;

    if (h == 0)
        text(x1+(x2-x1)/2-.15,y0+0.04+dy,"N.S.",'FontSize',18);
    else
        text(x1+(x2-x1)/2,y0+0.02+dy,"*",'FontSize',24);
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
