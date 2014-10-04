global rte;
%rte = 10;
if ~exist('loadDirectory'), loadDirectory = '/research-projects/tantra/tiwarip/imrt/'; end
%loadDirectory='~/imrt/';
disp('Loading data...');
eval(['load ' loadDirectory 'case',int2str(patID), '.mat']); %_9b_s1_th01
global planC; % (input parameter now)
indexS = planC{end};
caseParams = getCaseParams(patID);
eval(['init_pat_xx_', int2str(patID)]);
getBoundaryVoxel
voxelC

numTargets = length(targets);
