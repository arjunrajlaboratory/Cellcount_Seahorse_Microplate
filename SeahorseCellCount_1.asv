% counting cells using matlab's built-in functions
% This comment added 26-oct-2021 8:25pm
% 
% based on colonycounting_v2:
%   count_cells_all_scans's algorithm2:
%   d2utils.makeDAPImask
%   nucleiTable.findNuclei
%
% but with a gaussian filter added, and a junk mask
%% To read multiple TIF files
%addpath(genpath('~/Documents/Github/')); don't need to set this anymore
input_dir ='/Users/lucianncuenca/Dropbox (RajLab)/luciann/AllData/Seahorse assay/Seahorse_DAPI/E19_20211214/';
output_dir = '/Users/lucianncuenca/Dropbox (RajLab)/luciann/AllData/Seahorse assay/Seahorse_DAPI/E19_20211214/Cell count/';
filePattern = fullfile(input_dir, '*.tif');
all_files = dir(filePattern);
%find(contains(vertcat({all_files.name})','E9')) % find a text string
%Cell number
cell_count = cell(length(all_files), 2);

for k = 1: length(all_files)
%for k = 40
    baseFileName = all_files(k).name;
    fullFileName = fullfile(all_files(k) .folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName);
    numNuclei=CountCells(fullFileName, true,'rescaleMinMax',[0 12000],'maxNucleusArea',400,'gaussianBlurSigma',0.75);
    cell_count{k,2}=numNuclei;
    cell_count{k,1}=baseFileName;
end
export_table = cell2table(cell_count);
writetable(export_table, strcat(output_dir, 'cell_count.csv'));
