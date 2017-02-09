% This script generates a periodically stimulated inferior olive                                                                             
% 
% preface script call with parameters
% 
% example: HPGCGPU_scripts

rng(0,'twister')

% [=================================================================]
%  script parameters
% [=================================================================]

savehist = true;
saveappliednoise = true;
testing = 0;

% [=================================================================]
%  simulation parameters
% [=================================================================]

dt = 0.02;
simtime = 2000;

if ~exist('simtype')    	; simtype  = 'gallop'	  ; end
if ~exist('conntype')   	; conntype = 'iso'	  ; end
if ~exist('numruns')    	; numruns  = 1   	  ; end
if ~exist('sametoall')  	; sametoall = 0.1 	  ; end
if ~exist('tau')	    	; tau = 20  		  ; end
if ~exist('noisemu')	   	; noisemu = 0		  ; end
if ~exist('noisesig')    	; noisesig = 0		  ; end
if ~exist('gaps')	    	; gaps = [0 0.04]     ; end
if ~exist('gapcomp')    	; gapcomp = 0 		  ; end
if ~exist('moreoscillation'); moreoscillation = 0 ; end
if ~exist('nameprefix')  	; nameprefix = 'test' ; end
if ~exist('randampa')  		;  randampa = 0	      ; end
if ~exist('seed')  			;  seed = 0		      ; end
	
 


displaytext = [simtype '_' conntype '_' num2str(numruns) '_' num2str(sametoall)];

% [=================================================================]
%  % create network
% [=================================================================]

netsize = [50 1 1];
% netsize = [3 30 30];
	noneurons = prod(netsize);

plotthis  = 0;
rd = 3;
meannoconn = 8;

normleak  = 1;
randomize = 1;
scaling   = 1;
maxiter	  = 1;
somatapositions = [];
randomize = 1;
symmetrize = 1;



if ~exist(conntype);conntype = 'iso';end
switch conntype
	case 'iso'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [0 0 0 0], normleak);
	case 'cluster'
		W  = createW('3d_chebychev', netsize, rd, scaling, randomize, plotthis, maxiter, meannoconn, somatapositions, symmetrize, [1 100 1 0], normleak);
end


disp('[=================================================================]')
disp( ['using:' conntype])
disp( ['netsize:' num2str(netsize)])
disp( ['gap conditions:' num2str(gaps   )])
disp( ['sametoall:' num2str(sametoall)])
disp( ['simtime:' num2str(simtime) 'ms'])
disp('[=================================================================]')



% [=================================================================]
%  % noise levels
% [=================================================================]
noise_level_transients = [0 0 0 0];

noise_level = [1/tau noisesig noisemu 0];
noise_level = [];

% [=================================================================]
%  create perturbations AMPA
% [=================================================================]% 

pulses = [1 5 9 13];
inputrad  = 4;
switch simtype
	case 'gallop'
		triggersS = [5250:650:simtime];
		triggersL = [5650:650:simtime];
		triggers = union(triggersS, triggersL);
		% maybe repair, perchance remove: triggers = sort(unique(bsxfun(@plus, pulses,triggers')))
		spont = 0;
	case '1Hz'
		triggers = [5000:1000:simtime];
		triggers = [5000:1000:simtime];
		% triggers = sort(unique(bsxfun(@plus, pulses,triggers')))
		spont = 0;
	case 'spont'
		triggers = [];
		spont = 1;
	otherwise
		disp('simtype not recognized')
		disp('check string')
end


if not(spont)
	pert.mask{1}  	  = create_input_mask(netsize, 'dist_to_point', 'radius', inputrad,'cell_coordinates', W.coords,'projection_center', netsize/2,'synapseprobability',1,'plotme',plotthis);
	pert.amplitude{1} = 2;
	pert.type{1}	  = 'ampa';
	pert.duration{1}  = 1;
else
	pert = [];
end

% [=================================================================]
%  create perturbations OU noise
% [=================================================================]

% pert.mask{2}  	   = create_input_mask(netsize, 'random', 'synapse_probability', .5);;
% pert.amplitude{2} = 2;
% pert.type{2}	  = 'ou_noise';
% pert.duration{2}  = 1;

% th =	 1/5 ; % decay time parameter
% mu = 	 0 ; % pA
% sig = 	 0 ; % pA
% mix =    1;

% pert.param{nm}(1)  = 1/5 ;
% pert.param{nm}(2)  = 0  ;
% pert.param{nm}(3)  = sig ;
% pert.param{nm}(4)  = 0.1;

%     ____  __  ___   __   _____ ______  ___
%    / __ \/ / / / | / /  / ___//  _/  |/  /
%   / /_/ / / / /  |/ /   \__ \ / // /|_/ / 
%  / _, _/ /_/ / /|  /   ___/ // // /  / /  
% /_/ |_|\____/_/ |_/   /____/___/_/  /_/   
                                          
if testing
	ttime1 = 3; ttime2 =3 ; pks = 5;  pko = 5;
else	 
	ttime1 = 100;
	ttime2 = 300;
end

if ~exist('transients')
	
	neurons = createDefaultNeurons(noneurons,'celltypes','randomized','gapcompensation',0);
	
	[transients] = IOnet( 'networksize', netsize ,'time',ttime1,'delta',dt,'cell_parameters', neurons ,'W',W.W*gaps(1)/meannoconn,'ou_noise', noise_level_transients, 'sametoall',sametoall);
	[continuation] = IOnet( 'networksize', netsize ,'time',ttime2,'delta',dt,'cell_parameters', neurons ,'W',W.W*gaps(1)/meannoconn,'ou_noise', noise_level_transients, 'sametoall',sametoall, 'tempState',transients.lastState);


end

if ~spont
		[v pks] = findpeaks(mean(continuation.networkHistory.V_soma),'minpeakdistance',80);
		pko = pks(1);
		pert.triggers{1}  = pko + sort(unique(bsxfun(@plus, pulses,triggers')));
end


s = 0; 
for g = gaps


	% [=================================================================]
	%  % create neurons
	% [=================================================================]

	rng(0,'twister')

	neurons = createDefaultNeurons(noneurons,'celltypes','randomized','gapcompensation',gapcomp);

	if randampa
		neurons.gbar_ampa_soma = .08*ones(noneurons,1) + .02*rand(noneurons,1);
	end

	if moreoscillation
		neurons.g_CaL = neurons.g_CaL+.05;
	end

			for n = 1:numruns
				s = s+1;

				seed = seed+1;

				displaytext = [simtype '_' conntype '_' num2str(n) '_' num2str(sametoall)];

				noise_level(4) = seed;
				simresults{s} = IOnet('networksize', netsize,'time',simtime,'delta',dt,'cell_parameters',neurons,'tempState',transients.lastState,'W',W.W*g/meannoconn ,'ou_noise', noise_level , 'perturbation', pert,'sametoall',sametoall,'saveappliednoise',saveappliednoise, 'displaytext',displaytext);

				simresults{s}.spikes  = spikedetect(simresults{s});	
				simresults{s}.W = W;

				if not(savehist) & s>1
			   		simresults{s}.networkHistory = [];
			   		
			    end

			end

			seed = 0;
			
end

eval(['save periodic_ampa_' nameprefix num2str(s) '_' conntype '_' num2str(gaps) '_' simtype '_' num2str(simtime) '_' num2str(numruns) '_' date ' -v7.3'])

	