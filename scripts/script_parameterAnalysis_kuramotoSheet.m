clear

%%% File/Folder Management %%%
outputFolder = 'C:\Users\Jeroen\Documents\GitHub\OliveTree\data\script_parameterAnalysis\';
name = strrep(mat2str(clock),' ','_');
outputFolder2 = [outputFolder name];
mkdir(outputFolder2);
code = [mfilename('fullpath') '.m'];
copyfile(code,outputFolder2);
outputFile = [outputFolder2 '\' name '.mat'];

file = matfile(outputFile,'Writable', true);

disp(' ')
disp(['Starting time ' name])
disp(['Output will be saved to: ' outputFile '.mat'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% STATIC PARAMS
N = 10; M = 10;
scaling = 4*pi;

connectivity = {'euclidean' 'all to all' 'null'};
radius = 2;
init_scaling = 0.6;

sigmoid = [10 0.5];
stepval = 0.501;

time = 10;
steps = 10000;
decay = 0;

training = [0 10*2*pi];
training_time = 0.5;
training_signal = [];

record_adjacency = 0;

% DYNAMIC PARAMS - [initial (minus 1 step), final, steps]
sigP  = [0 10 20];
TauP  = [0.001 0.01 20];

totalIter = sigP(3)*TauP(3);
disp(['Simulation will run for ' mat2str(totalIter) ' iterations.'])

for ss = 1:sigP(3)
    for TT = 1:TauP(3)
        tic;
        curIter = (ss-1)*TauP(3)+TT;
        disp(' ')
        disp(['Starting iteration ' mat2str(curIter) ' of ' mat2str(totalIter) ' iterations.'])
      
        
        sigVal = sigP(1) + (sigP(2) - sigP(1))./sigP(3)*ss;
        TauVal = TauP(1) + (TauP(2) - TauP(1))./TauP(3)*TT;
        
        disp(['TauVal = ' mat2str(TauVal) ', SigVal = ' mat2str(sigVal) '.'])
        
        % DEPENDENT VALUES
        median = 10;
        interval = [median-sigVal median+sigVal];
        oscillators = interval(1) + (interval(2)-interval(1)).*rand(N*M,1);
        oscillators = oscillators*2*pi;
        
        plasticity = {'STDP' 0.3 [1 1] [TauVal TauVal]};
        
        out1 = kuramotoSheet([N M],scaling,'plotme', 0, 'record_adjacency', record_adjacency', 'connectivity', connectivity{1}, 'radius', radius, 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity,'sigmoid',sigmoid, 'init_scaling',init_scaling, 'decay', decay,'training', training, 'training_signal', training_signal,'training_time', training_time, 'stepval',stepval);
        out2 = kuramotoSheet([N M],scaling,'plotme', 0, 'record_adjacency', record_adjacency', 'connectivity', connectivity{2}, 'radius', radius, 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity,'sigmoid',sigmoid, 'init_scaling',init_scaling, 'decay', decay,'training', training, 'training_signal', training_signal,'training_time', training_time, 'stepval',stepval);
        out3 = kuramotoSheet([N M],scaling,'plotme', 0, 'record_adjacency', record_adjacency', 'connectivity', connectivity{3}, 'radius', radius, 'oscillators',oscillators, 'time', time, 'dt', time/steps, 'plasticity', plasticity,'sigmoid',sigmoid, 'init_scaling',init_scaling, 'decay', decay,'training', training, 'training_signal', training_signal,'training_time', training_time, 'stepval',stepval);
        
        fields = {'phase', 'orderparameter'};
        out1 = rmfield(out1,fields);
        out2 = rmfield(out2,fields);
        out3 = rmfield(out3,fields);
        
        file.out1(ss,TT) = out1;
        file.out2(ss,TT) = out2;
        file.out3(ss,TT) = out3;
        
        clear out1 out2 out3 sigVal TauVal interval oscillators plasticity
        timeIter = toc;
        disp(['Last iteration took ' mat2str(timeIter) ' seconds. Estimated time remaining is ' mat2str(floor((totalIter-curIter)*timeIter/60)) ' seconds.'])
    end
end