% Performs analysis on the results of
% script_parameterAnalysis_kuramotoSheet.m
%
% Currently performs:
%   Spectral Clustering using spect_clust.m
%   Separates useful information like dw/dwabs

% Note: This script is full of hardcoded crap, adjust values/code to whatever is in the
% inputFile (e.g. amount of iterations)
% Load the file and check its contents, then try a single iteration before
% committing to the whole thing.
%% Load File
clear
[filename, outputFolder] = uigetfile('*.mat','Select the file containing the parameterAnalysis data');

inputFile    = [outputFolder filename];
outputFile   = [outputFolder 'analysis.mat'];

file         = matfile(inputFile, 'Writable', false);

%% Rest of the program
outFile      = matfile(outputFile, 'Writable', true);
code         = [mfilename('fullpath') '.m'];
copyfile(code,outputFolder);

N = 10; M = 10;

init_struct = struct('dw',{},'dwabs',{},'parameters',{},'cc',{},'c',{},'bic',{},'aic',{}); %Make sure the output struct and init_struct have the same order, it really likes to bitch about this.
init_struct(N,M).cc = [];
outFile.out1 = init_struct;
outFile.out2 = init_struct;
outFile.out3 = init_struct;

for ss=1:N
    tic
    data1 = file.out1(ss,:);
    data2 = file.out2(ss,:);
    data3 = file.out3(ss,:);
    
    output = cell(3,M);
    
    for TT=1:M
        states    = cell(1,3);
        states{1} = data1(1,TT).state;
        states{2} = data2(1,TT).state;
        states{3} = data3(1,TT).state;
        
        output{1,TT}.dw     = data1(1,TT).dw;
        output{1,TT}.dwabs  = data1(1,TT).dwabs;
        output{1,TT}.parameters = data1(1,TT).parameters;
        output{2,TT}.dw     = data2(1,TT).dw;
        output{2,TT}.dwabs  = data2(1,TT).dwabs;
        output{2,TT}.parameters = data2(1,TT).parameters;
        output{3,TT}.dw     = data3(1,TT).dw;
        output{3,TT}.dwabs  = data3(1,TT).dwabs;
        output{3,TT}.parameters = data3(1,TT).parameters;
        
        for ii=1:3
            if any(any(isnan(states{ii})))
                output{ii,TT}.cc = NaN;
                output{ii,TT}.c = NaN;
                output{ii,TT}.bic = NaN;
                output{ii,TT}.aic = NaN;
                continue
            end           
            
            x = states{ii}(:,9001:end);
            cc = calcSimilarity(x);
            c = spect_clust(cc, 20);
            [bic, aic] = baic(x, c);

            output{ii,TT}.cc = cc;
            output{ii,TT}.c = c;
            output{ii,TT}.bic = bic;
            output{ii,TT}.aic = aic;
            
            clear x cc rho c bic aic
        end
        
        clear states
    end
    
    outFile.out1(ss,:) = output{1,:};
    outFile.out2(ss,:) = output{2,:};
    outFile.out3(ss,:) = output{3,:};
    
    clear output
    iterTime = toc;
    disp(['Iteration ' mat2str(ss) ' took ' mat2str(iterTime) ' seconds.'])
end
