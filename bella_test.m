% choose anaylsis files 

homedir = [getenv('USERPROFILE')];
if ~ispc
    homedir = getuserdir
end

exptDir = fullfile(homedir,'Dropbox (University of Oregon)','UO-Sylwestrak Lab','Opioids','TestScript','data');
%exptDir = '~/Dropbox (University of Oregon)/UO-Sylwestrak Lab/Opioids/Morphine_FosGAD_BH/data/';

%% Get Analysis Directory
dataDates = { '20220318', '20220329', '20220331'};

analysisFiles = [];
for d=1:numel(dataDates)
imageDate = dataDates{d};
dataDir =  fullfile(exptDir, imageDate, 'processed');
aFiles = dir(fullfile(dataDir, '**/*_coexpr.*'));  %get list of files and folders in any subfolder
analysisFiles = [analysisFiles; aFiles];
end

%% Extract info from filenames
%Extract Dates
dates = cellfun(@(n) split(n,'_'),{analysisFiles.name},'UniformOutput',false);
dates = cell2mat(cellfun(@(n) str2double(n{1}),dates,'UniformOutput',false)');

%Extract Section Number
section = cellfun(@(n) split(n,'_'),{analysisFiles.name},'UniformOutput',false);
section = cellfun(@(n) n{4},section,'UniformOutput',false);
section = cellfun(@(n) split(n,'-'),section,'UniformOutput',false);
section = cellfun(@(n) str2num(n{2}),section,'UniformOutput',false)';

%Extract Drug Treatment 
drug = cellfun(@(n) split(n,'_'),{analysisFiles.name},'UniformOutput',false);
drug = cellfun(@(n) n{4},drug,'UniformOutput',false);
drug = cellfun(@(n) split(n,'-'),drug,'UniformOutput',false);
drug = cellfun(@(n) n{1},drug,'UniformOutput',false)';
drug1=categorical(); drug2=categorical();
treatment=categorical();

%Recombine Treatment to make labels for graph
for i=1:numel(drug)
    idx = find(isstrprop(drug{i},'upper'));
    drug1(i) = drug{i}(idx(1):idx(2)-1);
    drug2(i) = drug{i}(idx(2):end);
    treatment(i) = [drug{i}(idx(1):idx(2)-1) '+' drug{i}(idx(2):end)];
end

%% Load in Cell Count Data
nFiles = numel(analysisFiles);

%Set up empty 
G = table(dates,section);
CoExpr = cell(nFiles,1);
Coord = cell(nFiles,1);
M_F = cell(nFiles,1);
L_F = cell(nFiles,1);
M_G = cell(nFiles,1);
L_G = cell(nFiles,1);
GADCoords = cell(nFiles,1);
coexpr = cell(nFiles,1);
H_F= cell(nFiles,1);
H_G= cell(nFiles,1);

%Iterate through experimental dates
for d = 1:nFiles
    %Load original fos coordinates
    cellFile = dir([analysisFiles(d).folder '/*_cells.mat'] );
    load([cellFile.folder '/' cellFile.name]);
    %Clean up coordinate array to remove NaNs
    idxNAN = find(~isnan(nCoord(:,1)),1,'last');
    Coord{d} = nCoord(1:idxNAN);
    
    %Load fos positional analysis
    cellFile = dir([analysisFiles(d).folder '/*_analysis.mat'] );
    load([cellFile.folder '/' cellFile.name]);
    M_F{d} = inMHb;
    L_F{d} = inLHb;
    H_F{d} = inMHb+inLHb;
    
    %Load GAD coordinates
    cellFile = dir([analysisFiles(d).folder '/*_GADcells.mat'] );
    load([cellFile.folder '/' cellFile.name]);
    %Clean up coordinate array to remove NaNs
    idxNAN = find(~isnan(GADCoord(:,1)),1,'last');
    GADCoords{d} = GADCoord(1:idxNAN);
    
    %Load fos positional analysis
    cellFile = dir([analysisFiles(d).folder '/*_GADanalysis.mat'] );
    load([cellFile.folder '/' cellFile.name]);
    M_G{d} = GADinMHb;
    L_G{d} = GADinLHb;
    H_G{d} = GADinMHb+GADinLHb;
    
    %load Coexpression Data
    cellFile = dir([analysisFiles(d).folder '/*coexpr.mat'] );
    load([cellFile.folder '/' cellFile.name]);
    %Clean up coordinate array to remove NaNs
    idxNAN = find(~isnan(nCoord(:,1)),1,'last');
    coexpr(isnan(coexpr))=0;
    CoExpr{d} = logical(coexpr(1:idxNAN));
  
end

%Add Data to table 'G'
G.Drug1 = drug1';
G.Drug2 = drug2';
G.Treatment = treatment';
G.Coord_Fos = Coord;
G.Coord_Gad = GADCoords;
G.inMHb_Fos = cellfun(@logical, M_F, 'UniformOutput', false);
G.inLHb_Fos = cellfun(@logical, L_F, 'UniformOutput', false);
G.inHb_Fos = cellfun(@logical, H_F, 'UniformOutput', false);
G.inMHb_GAD = cellfun(@logical, M_G, 'UniformOutput', false);
G.inLHb_GAD = cellfun(@logical, L_G, 'UniformOutput', false);
G.inHb_GAD = cellfun(@logical, H_G, 'UniformOutput', false);
G.CoExpr = CoExpr;

%% Plot Fos+ neurons vs Treatment Condition
%Calculate Stats

Morph_rows = G(G.Drug1== 'Morph',:)
Sal_rows = G(G.Drug1== 'Sal',:)

Morph_MHbFos = cellfun(@sum, Morph_rows.inMHb_Fos)
Morph_LHbFos = cellfun(@sum, Morph_rows.inLHb_Fos)
Sal_MHbFos = cellfun(@sum, Sal_rows.inMHb_Fos)
Sal_LHbFos = cellfun(@sum, Sal_rows.inLHb_Fos)


%Loop through each anatomical group and run two-way anova to determine if
%there is an effect of morphine

T = table(); 
MorphTable = table();
meanMorph_MHbFos = mean(Morph_MHbFos)
meanMorph_LHbFos = mean(Morph_LHbFos)
meanSal_MHbFos = mean(Sal_MHbFos)
meanSal_LHbFos = mean(Sal_LHbFos)

boxmorphdata = horzcat(Morph_MHbFos, Morph_LHbFos);
boxsaldata = horzcat(Sal_MHbFos, Sal_LHbFos)

data = horzcat(meanMorph_MHbFos, meanMorph_LHbFos,meanSal_MHbFos, meanSal_LHbFos)
sdev = horzcat(std(Morph_MHbFos), std(Morph_LHbFos), std(Sal_MHbFos), std(Sal_LHbFos))

[h, p, ci] = ttest2(boxmorphdata,boxsaldata)
p1 = p 

%Rearrange order of conditions in data set and labels
keySet = {'m_MHb','m_LHb','s_MHb','s_LHb'};
valueSet = [1 2 3 4];
t = containers.Map(keySet,valueSet);
[~, idx] = sort(cellfun(@(n) t(n),keySet));
data = data(:,idx);
err = sdev(:,idx);


% Make a Figure
subplot(3,3,1)
b = bar(data,'k');
hold on
errorbar(b.XData,data,err,'k','LineStyle','none','LineWidth',1.5)

ylim([0 max(data(:))+5])
xticklabels(keySet)
xtickangle(90)
uniformFigureProps()
title('Number of Fos+ cells In Hb', 'FontSize', 9)
ylabel('Number of Fos+ cells In Hb', 'FontSize', 9)



% make a box and whisker plot 
boxmorphdata = horzcat(Morph_MHbFos, Morph_LHbFos);
boxsaldata = horzcat(Sal_MHbFos, Sal_LHbFos)

% Left axes
subplot(3,3,2)
boxplot(boxmorphdata, 'Labels', {'Mhb', 'LHb'})
title('Morphine Treated animals')
ylabel('Number of Fos+ cells')

% Right axes
subplot(3,3,3)
boxplot(boxsaldata, 'Labels', {'Mhb', 'LHb'})
title('Saline Treated animals')
ylabel('Number of Fos+ cells')

%% Determine # of GAD+ Fos+ cells
%Calculate Stats
clear stat stats
Morph_rows = G(G.Drug1== 'Morph',:)
Sal_rows = G(G.Drug1== 'Sal',:)

Morph_MHbFos = cell2mat(cellfun(@(x,y) sum(x*y), Morph_rows.CoExpr,Morph_rows.inMHb_Fos, 'UniformOutput',false));
Morph_LHbFos = cell2mat(cellfun(@(x,y) sum(x*y), Morph_rows.CoExpr,Morph_rows.inLHb_Fos, 'UniformOutput',false));
Sal_MHbFos = cell2mat(cellfun(@(x,y) sum(x*y), Sal_rows.CoExpr,Sal_rows.inMHb_Fos, 'UniformOutput',false));
Sal_LHbFos = cell2mat(cellfun(@(x,y) sum(x*y), Sal_rows.CoExpr,Sal_rows.inLHb_Fos, 'UniformOutput',false));


%Loop through each anatomical group and run two-way anova to determine if
%there is an effect of morphine

T = table(); 
MorphTable = table();
meanMorph_MHbFos = mean(Morph_MHbFos);
meanMorph_LHbFos = mean(Morph_LHbFos);
meanSal_MHbFos = mean(Sal_MHbFos);
meanSal_LHbFos = mean(Sal_LHbFos);

data = horzcat(meanMorph_MHbFos, meanMorph_LHbFos,meanSal_MHbFos, meanSal_LHbFos);
sdev = horzcat(std(Morph_MHbFos), std(Morph_LHbFos), std(Sal_MHbFos), std(Sal_LHbFos));

[h, p, ci] = ttest2(boxmorphdata,boxsaldata)
p2 = p
% Make a bar plot 
%Rearrange order of conditions in data set and labels
keySet = {'m_MHb','m_LHb','s_MHb','s_LHb'};
valueSet = [1 2 3 4];
t = containers.Map(keySet,valueSet);
[~, idx] = sort(cellfun(@(n) t(n),keySet));
data = data(:,idx);
err = sdev(:,idx);


% Make a Figure
subplot(3,3,4)
b = bar(data,'k');
hold on
errorbar(b.XData,data,err,'k','LineStyle','none','LineWidth',1.5)

ylim([0 max(data(:))+5])
xticklabels(keySet)
xtickangle(90)
uniformFigureProps()
title('Number of Fos+ GAD+ cells In Hb', 'FontSize', 9)
ylabel('Number of Fos+ GAD+ cells In Hb', 'FontSize', 9)




% make a box and whisker plot 
boxmorphdata = horzcat(Morph_MHbFos, Morph_LHbFos);
boxsaldata = horzcat(Sal_MHbFos, Sal_LHbFos)

% Left axes
subplot(3,3,5)
boxplot(boxmorphdata, 'Labels', {'Mhb', 'LHb'})
title('Morphine Treated animals')
ylabel('Number of Fos+ Gad+ cells')

% Right axes
subplot(3,3,6)
boxplot(boxsaldata, 'Labels', {'Mhb', 'LHb'})
title('Saline Treated animals')
ylabel('Number of Fos+ Gad+ cells')

%%  Determine % of GAD+ that are Fos+
%Calculate Stats
clear stat stats
%G.CoExpr is a logical array if the Fos+ coexpresses GAD
%G.inHb_Fos is a logical array if the Fos+ is in the Hb
%G.inHb_GAD is a logical array of all GAD+ neurons if they are in the Hb
%G.CoExpr * G.inHb_Fos limit to coexpressing cells in the Hb
%Dividing by sum(G.inHb_GAD) normalizes for all GAD+ cells in HB
%therefore giving the %of GAD+ neurons expression Fos

Morph_rows = G(G.Drug1== 'Morph',:);
Sal_rows = G(G.Drug1== 'Sal',:);

Morph.GADwithFosHb= cell2mat(cellfun(@(x,y,z) sum(x*y)./sum(z)*100, Morph_rows.CoExpr,Morph_rows.inHb_Fos,Morph_rows.inHb_GAD, 'UniformOutput',false));
Morph_cellGroupings = {Morph.GADwithFosHb};

Sal.GADwithFosHb= cell2mat(cellfun(@(x,y,z) sum(x*y)./sum(z)*100, Sal_rows.CoExpr,Sal_rows.inHb_Fos,Sal_rows.inHb_GAD, 'UniformOutput',false));
Sal_cellGroupings = {Sal.GADwithFosHb};

%Summary Data
data = horzcat(mean(Morph.GADwithFosHb), mean(Sal.GADwithFosHb));
err = horzcat(std(Morph.GADwithFosHb), std(Sal.GADwithFosHb)) 
[h, p, ci] = ttest2(Morph.GADwithFosHb, Sal.GADwithFosHb)
p3 = p 
%Rearrange order of conditions in data set and labels
keySet = {'Morph+Sal', 'Sal+Sal'};
valueSet = [1 2];
t = containers.Map(keySet,valueSet);
[~, idx] = sort(cellfun(@(n) t(n),keySet));
data = data(:,idx);
err = err(:,idx);

% Make a Figure
subplot(3,3,7)
b = bar(data,'k');
hold on
errorbar(b.XData,data,err,'k','LineStyle','none','LineWidth',1.5)
  

% ylim([0 max(data(:))+5])
xticklabels(keySet)
xtickangle(90)
% uniformFigureProps()
title('% of GAD+ that are Fos+', 'FontSize', 9)
ylabel('(Fos+GAD)/GAD*100 In Hb', 'FontSize', 9)