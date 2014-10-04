initialize
%influenceM = getGlobalInfluenceM2(planC{indexS.IM}(choppedIM).IMDosimetry, unique([targets step2OARs step3OARs]));
%dlmwrite('inlog.txt',influenceM);

numPBs = size(influenceM,2);
%w = zeros(numPBs+numTargets,1);
D=[];
for tar = 1:numTargets
    D(end+1:end+numVoxels(targets(tar)),1) = (1/numVoxels(targets(tar))) *((influenceM(voxelC{targets(tar)},:)* w(1:numPBs)) - dosePrescribedTarget(tar)); 
end
for i=1:length(D) 
    D(i) = D(i)*D(i);
end
dlmwrite('log.txt',D);
t = w(numPBs+1:end);
for i=1:length(t) 
    t(i) = t(i)*t(i); 
end
val = sum(D) +  sum(t);
