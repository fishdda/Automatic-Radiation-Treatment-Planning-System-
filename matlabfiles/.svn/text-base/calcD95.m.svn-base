%patID=30;
if ~exist('loadDirectory'), loadDirectory = '/research-projects/tantra/tiwarip/imrt/'; end
disp('Loading data...');
%loadDirectory='~/imrt/';
eval(['load ' loadDirectory 'case',int2str(patID), '.mat']); %_9b_s1_th01

global planC; % (input parameter now)
indexS = planC{end};
choppedIM = size(planC{planC{end}.IM},2);
initFileName = sprintf('init_case_%d',patID);
initwdExt = strcat(initFileName,'.m');
if(exist(initwdExt,'file')==0)
    getParametersHN(planC,patID);
end
eval(initFileName)

getAllVoxel
numTargets = length(targets);


influenceM = getGlobalInfluenceM2(planC{indexS.IM}(choppedIM).IMDosimetry, unique([targets step2OARs step3OARs]));
numPBs = size(influenceM,2);
% wVI = ones(numPBs,1);
% wVII = ones(numPBs,1);
% wVIII = ones(numPBs,1);
% wVIV = ones(numPBs,1);
maxstep=4;
i=1;
IMNumber = size(planC{indexS.IM},2);
IM = planC{planC{end}.IM}(choppedIM).IMDosimetry;
while i<=maxstep
    d=[];
    if(i==1)
        wVI(numPBs+1:end) = [];
        d = influenceM*wVI;
	planC{indexS.IM}(IMNumber).IMDosimetry.solutions(i).beamletWeights = wVI;
    elseif(i==2)
        wVII(numPBs+1:end) = [];
        d = influenceM*wVII;
	planC{indexS.IM}(IMNumber).IMDosimetry.solutions(i).beamletWeights = wVII;
    elseif(i==3)
        wVIII(numPBs+1:end) = [];
        d = influenceM*wVIII;
	planC{indexS.IM}(IMNumber).IMDosimetry.solutions(i).beamletWeights=wVIII;
    elseif(i==4)
        wVIV(numPBs+1:end) = [];
        d = influenceM*wVIV;
	planC{indexS.IM}(IMNumber).IMDosimetry.solutions(i).beamletWeights = wVIV;
    end
    [dose3DM] = inflateSampledDoseM(d,Skin,max(sampleRates));
    doseName = sprintf('Step %i',i);
    dose2CERR(dose3DM,[],doseName,'CERR test','Test PB distribution.','UniformCT',[],'no', IM.assocScanUID); %assocScanUID added for CERR3
    planC{indexS.IM}(choppedIM).IMDosimetry.solution(i).doseArray = sparse(dose3DM(:));
    i=i+1;
end

disp('Saving data...');
pfname = sprintf('%scase%d_ipopt.mat',loadDirectory,patID);
save(pfname,'planC');

doseNum = size(planC{planC{end}.dose},2);
l = length(targets)+length(step2OARs);
metrics=zeros(l,1);
k=1;
for tar = [targets ]
    
    metrics(k) = Dx(planC, tar, doseNum, 95);
    k=k+1;
end

for i=1:length(step2OARs)
    metrics(k) = MOHx(planC, step2OARs(i), doseNum, 5);
    k=k+1;
end
% for i=1:length(step3OARs)
%     d = influenceM(allVoxelC{step3OARs(i)},:) * wV / doseScale;
%     metrics(k) = min(d);
%     k=k+1;
%     metrics(k) = mean(d);
%     k=k+1;
%     metrics(k) = max(d);
%     k=k+1;
% end
