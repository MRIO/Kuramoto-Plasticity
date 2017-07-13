% This script is full of hardcoded crap, adjust it to whatever is in the
% inputFile (e.g. amount of iterations)

clear

outputFolder = 'C:\Users\Jeroen\Documents\GitHub\OliveTree\data\script_parameterAnalysis\[2017_7_12_11_56_31.602]\';
inputFile    = 'C:\Users\Jeroen\Documents\GitHub\OliveTree\data\script_parameterAnalysis\[2017_7_12_11_56_31.602]\[2017_7_12_11_56_31.602].mat';
outputFile   = 'C:\Users\Jeroen\Documents\GitHub\OliveTree\data\script_parameterAnalysis\[2017_7_12_11_56_31.602]\analysis.mat';

file         = matfile(inputFile, 'Writable', false);
outFile      = matfile(outputFile, 'Writable', true);
code         = [mfilename('fullpath') '.m'];
copyfile(code,outputFolder);

init_struct = struct('cc',{},'c',{},'bic',{},'aic',{});
init_struct(20,20).cc = [];
outFile.out1 = init_struct;
outFile.out2 = init_struct;
outFile.out3 = init_struct;

for ss=1:20
    tic
    data1 = file.out1(ss,:);
    data2 = file.out2(ss,:);
    data3 = file.out3(ss,:);
    
    output = cell(20,3);
    
    for TT=1:20
        states = cell(1,3);
        states{1} = data1(1,TT).state;
        states{2} = data2(1,TT).state;
        states{3} = data3(1,TT).state;
        
        for ii=1:3
            x = sin(states{ii});
            cc = partialcorr(x', mean(x)'); % partial corherelation w.r.t the mean
            rho = (cc + 1)/2; % normalize rho to get an affinity matrix 0<=rho<=1
            rho = (rho+rho')/2; % rho had better be symmetric
            c = spect_clust(rho, 20);
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
    disp(['Iteration took ' mat2str(iterTime) ' seconds. This was iteration ' mat2str(ss)]) 
end
