%% This code is to generate the domain-filling and evenly-spaced streamlines aligning with layering orientations
%% This code is adapted from the publicaly available stress visualizer "3D-TSV" (https://github.com/Junpeng-Wang-TUM/3D-TSV)
%% Author: Junpeng Wang (junpeng.wang@tum.de)
%% Version: 2021-08-05
%% Update:  2023-11-10
function ShowRank3byStreamlines(alphaList, thetaList, meshInfo, PSLsDensityCtrl)
	global topRankNparas_;
	global meshStruct_;
	global followingLayerTabgent_;
	global snappingOpt_;
	global interceptionThreshold_;
	
	meshStruct_ = meshInfo;
	topRankNparas_ = zeros(meshStruct_.numElements, 9);
	topRankNparas_(:, [1 4 7]) = alphaList;
	topRankNparas_(:, [2 3 5 6 8 9]) = [cos(thetaList(:,1)) sin(thetaList(:,1)) ...
		cos(thetaList(:,2)) sin(thetaList(:,2)) cos(thetaList(:,3)) sin(thetaList(:,3))];		
	followingLayerTabgent_ = 1;
	interceptionThreshold_ = 6;
	snappingOpt_ = 0;
	
	Setup();
	rankLayerOpt = [1 1 1];
	CreatePSLsViaTSV(PSLsDensityCtrl, rankLayerOpt);
	ShowStreamlines();
end

function Setup()
	global meshStruct_;
	global nodeCoords_;
	global eleSize_;
	global tracingStepWidth_;

	global layer_1_StreamlinePool_;
	global layer_2_StreamlinePool_;
	global layer_3_StreamlinePool_;
	
	nodeCoords_ = meshStruct_.nodeCoords(meshStruct_.nodMapBack,:);
	eleSize_ = meshStruct_.eleSize(1);
	layer_1_StreamlinePool_ = StreamlineStruct();
	layer_2_StreamlinePool_ = StreamlineStruct();
	layer_3_StreamlinePool_ = StreamlineStruct();
    tracingStepWidth_ = eleSize_/2;
end

function CreatePSLsViaTSV(resCtrl, rankLayerOpt)
	global eleSize_;
	global nodeCoords_;
	global meshStruct_;
	
	global minimumEpsilon_;
	global mergeTrigger_;
	global seedPointsHistory_;
	global seedPoints_;
	global seedPointsValence_;

	global startCoord_;
	global relaxedFactor_;
	
	global layer_1_StreamlinePool_;
	global layer_2_StreamlinePool_;
	global layer_3_StreamlinePool_;
	global layer_1_CoordList_; layer_1_CoordList_ = [];
    global layer_2_CoordList_; layer_2_CoordList_ = [];
    global layer_3_CoordList_; layer_3_CoordList_ = [];
	global seedPoints_layer_1_;
	global seedPoints_layer_2_;
	global seedPoints_layer_3_;

	%%1. Create Seed Points
	minimumEpsilon_ = min(meshStruct_.boundingBox(2,:)-meshStruct_.boundingBox(1,:))/resCtrl;
	validNodes = zeros((meshStruct_.resX+1)*(meshStruct_.resY+1),1);
	validNodes(meshStruct_.nodMapBack) = (1:length(meshStruct_.nodMapBack))';
	validNodes = reshape(validNodes, meshStruct_.resY+1, meshStruct_.resX+1);
    
    if min([meshStruct_.resX meshStruct_.resX])<100
        seedDensCtrl = 1;
    else
        seedDensCtrl = max(ceil(minimumEpsilon_/eleSize_/10), 2);
    end	
	sampledNodes = validNodes(seedDensCtrl+1:seedDensCtrl:meshStruct_.resY+1-seedDensCtrl, ...
		seedDensCtrl+1:seedDensCtrl:meshStruct_.resX+1-seedDensCtrl);	
	sampledNodes = reshape(sampledNodes, numel(sampledNodes), 1);
	sampledNodes(0==sampledNodes) = [];
	sampledNodes = setdiff(sampledNodes, meshStruct_.nodesOnBoundary);
	seedPointsHistory_ = nodeCoords_(sampledNodes,:);

	%%2. Initialize PSL pool
	relaxedFactor_ = 1;
	startCoord_ = sum(meshStruct_.boundingBox,1)/2;
	mergeTrigger_ = minimumEpsilon_;
	seedPoints_ = seedPointsHistory_;
	numSeedPoints = size(seedPoints_,1);	
    seedPointsValence_ = zeros(numSeedPoints, 1);

	%%3. Streamline sampling via TSV
	disp('Performing TSV to Generate Streamlines');
	%%3.1 Layer 1
	if rankLayerOpt(1)
		seedPoints_layer_1_ = seedPointsHistory_;
		seedValence_layer_1 = seedPointsValence_;
		its = 0;
		looper = sum(seedValence_layer_1);	
		while looper<1*numSeedPoints
			its = its + 1;
			unFinishedSpps = find(seedValence_layer_1<1);
			spp = unFinishedSpps(1);
			if 0==looper
				[~, spp] = min(vecnorm(startCoord_-seedPoints_layer_1_,2,2));
			else				
				[~, tarPos] = min(vecnorm(startCoord_-seedPoints_layer_1_(unFinishedSpps,:),2,2));					
				spp = unFinishedSpps(tarPos);		
            end
			seed = seedPoints_layer_1_(spp,:);	
			seedValence_layer_1(spp,1) = 1;
			iStreamline = Have1moreStreamline(seed, 'ONE');		
			if 0==iStreamline.length
				looper = sum(seedValence_layer_1); 
				disp([' Iteration: ' sprintf('%4i',its) ' Progress: ' sprintf('%6i',looper) ' Total: ' sprintf('%6i',1*numSeedPoints)]);
				continue; 
			end			
			layer_1_StreamlinePool_(end+1,1) = iStreamline;				
			layer_1_CoordList_(end+1:end+iStreamline.length,:) = iStreamline.phyCoordList;
			sppsEmptyLayer1Valence = find(0==seedValence_layer_1);
			if ~isempty(sppsEmptyLayer1Valence)
				[potentialDisList1, potentialPosList1] = GetDisListOfPointList2Curve(seedPoints_layer_1_(...
						sppsEmptyLayer1Valence,:), iStreamline.phyCoordList);					
				potentialSolidSpps1 = find(potentialDisList1<relaxedFactor_);
				if ~isempty(potentialSolidSpps1)
					spps2BeMerged = sppsEmptyLayer1Valence(potentialSolidSpps1);
					seedPoints_layer_1_(spps2BeMerged,:) = potentialPosList1(potentialSolidSpps1,:);								
					seedValence_layer_1(spps2BeMerged,1) = 1;			
				end
			end	
			looper = sum(seedValence_layer_1);
			disp([' Iteration: ' sprintf('%4i',its) ' Progress: ' sprintf('%6i',looper) ' Total: ' sprintf('%6i',1*numSeedPoints)]);					
		end
	end
	
	%%3.2 Layer 2
	if rankLayerOpt(2)
		seedPoints_layer_2_ = seedPointsHistory_;
		seedValence_layer_2 = seedPointsValence_;
		its = 0;
		looper = sum(seedValence_layer_2);	
		while looper<1*numSeedPoints
			its = its + 1;
			unFinishedSpps = find(seedValence_layer_2<1);
			spp = unFinishedSpps(1);
			if 0==looper
				[~, spp] = min(vecnorm(startCoord_-seedPoints_layer_2_,2,2));
			else				
				[~, tarPos] = min(vecnorm(startCoord_-seedPoints_layer_2_(unFinishedSpps,:),2,2));					
				spp = unFinishedSpps(tarPos);		
            end
			seed = seedPoints_layer_2_(spp,:);	
			seedValence_layer_2(spp,1) = 1;
			iStreamline = Have1moreStreamline(seed, 'TWO');		
			if 0==iStreamline.length
				looper = sum(seedValence_layer_2); 
				disp([' Iteration: ' sprintf('%4i',its) ' Progress: ' sprintf('%6i',looper) ' Total: ' sprintf('%6i',1*numSeedPoints)]);
				continue; 
			end			
			layer_2_StreamlinePool_(end+1,1) = iStreamline;				
			layer_2_CoordList_(end+1:end+iStreamline.length,:) = iStreamline.phyCoordList;
			sppsEmptyLayer1Valence = find(0==seedValence_layer_2);
			if ~isempty(sppsEmptyLayer1Valence)
				[potentialDisList2, potentialPosList2] = GetDisListOfPointList2Curve(seedPoints_layer_2_(...
						sppsEmptyLayer1Valence,:), iStreamline.phyCoordList);					
				potentialSolidSpps2 = find(potentialDisList2<relaxedFactor_);
				if ~isempty(potentialSolidSpps2)
					spps2BeMerged = sppsEmptyLayer1Valence(potentialSolidSpps2);
					seedPoints_layer_2_(spps2BeMerged,:) = potentialPosList2(potentialSolidSpps2,:);								
					seedValence_layer_2(spps2BeMerged,1) = 1;			
				end
			end	
			looper = sum(seedValence_layer_2);
			disp([' Iteration: ' sprintf('%4i',its) ' Progress: ' sprintf('%6i',looper) ' Total: ' sprintf('%6i',1*numSeedPoints)]);					
		end
	end	

	%%3.2 Layer 3
	if rankLayerOpt(3)
		seedPoints_layer_3_ = seedPointsHistory_;
		seedValence_layer_3 = seedPointsValence_;
		its = 0;
		looper = sum(seedValence_layer_3);	
		while looper<1*numSeedPoints
			its = its + 1;
			unFinishedSpps = find(seedValence_layer_3<1);
			spp = unFinishedSpps(1);
			if 0==looper
				[~, spp] = min(vecnorm(startCoord_-seedPoints_layer_3_,2,2));
			else				
				[~, tarPos] = min(vecnorm(startCoord_-seedPoints_layer_3_(unFinishedSpps,:),2,2));					
				spp = unFinishedSpps(tarPos);		
            end
			seed = seedPoints_layer_3_(spp,:);	
			seedValence_layer_3(spp,1) = 1;
			iStreamline = Have1moreStreamline(seed, 'THREE');		
			if 0==iStreamline.length
				looper = sum(seedValence_layer_3); 
				disp([' Iteration: ' sprintf('%4i',its) ' Progress: ' sprintf('%6i',looper) ' Total: ' sprintf('%6i',1*numSeedPoints)]);
				continue; 
			end			
			layer_3_StreamlinePool_(end+1,1) = iStreamline;				
			layer_3_CoordList_(end+1:end+iStreamline.length,:) = iStreamline.phyCoordList;
			sppsEmptyLayer1Valence = find(0==seedValence_layer_3);
			if ~isempty(sppsEmptyLayer1Valence)
				[potentialDisList3, potentialPosList3] = GetDisListOfPointList2Curve(seedPoints_layer_3_(...
						sppsEmptyLayer1Valence,:), iStreamline.phyCoordList);					
				potentialSolidSpps3 = find(potentialDisList3<relaxedFactor_);
				if ~isempty(potentialSolidSpps3)
					spps2BeMerged = sppsEmptyLayer1Valence(potentialSolidSpps3);
					seedPoints_layer_3_(spps2BeMerged,:) = potentialPosList3(potentialSolidSpps3,:);								
					seedValence_layer_3(spps2BeMerged,1) = 1;			
				end
			end	
			looper = sum(seedValence_layer_3);
			disp([' Iteration: ' sprintf('%4i',its) ' Progress: ' sprintf('%6i',looper) ' Total: ' sprintf('%6i',1*numSeedPoints)]);					
		end
	end	
	
	%%Compact PSLs sets
	minLength = 2;
	PSLs2Bdiscarded = []; index = 0;
	for ii=1:length(layer_1_StreamlinePool_)
		if layer_1_StreamlinePool_(ii).length<minLength
			index = index + 1;
			PSLs2Bdiscarded(index) = ii;
		end
	end
	layer_1_StreamlinePool_(PSLs2Bdiscarded) = [];
	
	PSLs2Bdiscarded = []; index = 0;
	for ii=1:length(layer_2_StreamlinePool_)
		if layer_2_StreamlinePool_(ii).length<minLength
			index = index + 1;
			PSLs2Bdiscarded(index) = ii;
		end
	end
	layer_2_StreamlinePool_(PSLs2Bdiscarded) = [];
	
	PSLs2Bdiscarded = []; index = 0;
	for ii=1:length(layer_3_StreamlinePool_)
		if layer_3_StreamlinePool_(ii).length<minLength
			index = index + 1;
			PSLs2Bdiscarded(index) = ii;
		end
	end
	layer_3_StreamlinePool_(PSLs2Bdiscarded) = [];
end

function iPSL = Have1moreStreamline(seed, psDir)
	global snappingOpt_;
	switch psDir
		case 'ONE'
			[iPSL, ~, ~] = GeneratePrincipalStressLines_elementWise(seed, psDir);
		case 'TWO'
			[~, iPSL, ~] = GeneratePrincipalStressLines_elementWise(seed, psDir);
		case 'THREE'
			[~, ~, iPSL] = GeneratePrincipalStressLines_elementWise(seed, psDir);			
	end
	if snappingOpt_, iPSL = CroppingPSLifNeeded(iPSL, psDir); end	
end

function [potentialDisList, potentialPosList] = GetDisListOfPointList2Curve(pointList, curveLine)
	global mergeTrigger_;
	disT = (curveLine(:,1) - pointList(:,1)').^2;
	disT = disT + (curveLine(:,2) - pointList(:,2)').^2;
	disT = sqrt(disT);	
	[minVal, minValPos] = min(disT,[],1);
	potentialDisList = minVal';
	potentialDisList = potentialDisList/mergeTrigger_;
	potentialPosList = curveLine(minValPos,:);	
end

function tarPSL = CroppingPSLifNeeded(srcPSL, psDir)
	global mergeTrigger_;
	global layer_1_CoordList_; global layer_2_CoordList_; global layer_3_CoordList_;
	tarPSL = srcPSL;
	if 5>=srcPSL.length, return; end
	disThreshold = 1.5;
	relaxedThreshold = 0.2;
	switch psDir
		case 'ONE'
			if isempty(layer_1_CoordList_), return; end
			srcCoordList = layer_1_CoordList_;		
		case 'TWO'
			if isempty(layer_2_CoordList_), return; end
			srcCoordList = layer_2_CoordList_;					
		case 'THREE'
			if isempty(layer_3_CoordList_), return; end
			srcCoordList = layer_3_CoordList_;		
	end
	
	if srcPSL.midPointPosition == srcPSL.length || srcPSL.midPointPosition == 1
		if 1==srcPSL.midPointPosition
			tarCoordList = tarPSL.phyCoordList;
			disT = (srcCoordList(:,1) - tarCoordList(:,1)').^2;
			disT = disT + (srcCoordList(:,2) - tarCoordList(:,2)').^2;	
			disT = sqrt(disT);
			miniDisList2SrcPSL = min(disT);
			tarPositions = find(miniDisList2SrcPSL<mergeTrigger_/disThreshold);
			%if isempty(tarPositions), return; end
			if length(tarPositions)/size(tarCoordList,1)<relaxedThreshold, return; end
			if length(tarPositions) == size(tarCoordList,1), tarPSL = StreamlineStruct(); return; end
			startPos = 1; endPos = min(tarPositions);
		else
			tarCoordList = flip(tarPSL.phyCoordList,1);
			disT = (srcCoordList(:,1) - tarCoordList(:,1)').^2;
			disT = disT + (srcCoordList(:,2) - tarCoordList(:,2)').^2;	
			disT = sqrt(disT);
			miniDisList2SrcPSL = min(disT);
			tarPositions = find(miniDisList2SrcPSL<mergeTrigger_/disThreshold);
			%if isempty(tarPositions), return; end
			if length(tarPositions)/size(tarCoordList,1)<relaxedThreshold, return; end
			if length(tarPositions) == size(tarCoordList,1), tarPSL = StreamlineStruct(); return; end
			startPos = srcPSL.length - min(tarPositions) + 1; endPos = srcPSL.length;
		end	
	else
		tarCoordList = tarPSL.phyCoordList(srcPSL.midPointPosition:srcPSL.length,:);
		disT = (srcCoordList(:,1) - tarCoordList(:,1)').^2;
		disT = disT + (srcCoordList(:,2) - tarCoordList(:,2)').^2;	
		disT = sqrt(disT);
		miniDisList2SrcPSL = min(disT);
		tarPositions = find(miniDisList2SrcPSL<mergeTrigger_/disThreshold);		
		%if isempty(tarPositions)
		if length(tarPositions)/size(tarCoordList,1)<relaxedThreshold
			endPos = srcPSL.length;
		elseif length(tarPositions) == size(tarCoordList,1)
			endPos = srcPSL.midPointPosition;
		else
			endPos = srcPSL.midPointPosition+min(tarPositions)-1;
		end
		
		tarCoordList = flip(tarPSL.phyCoordList(1:srcPSL.midPointPosition,:),1);
		disT = (srcCoordList(:,1) - tarCoordList(:,1)').^2;
		disT = disT + (srcCoordList(:,2) - tarCoordList(:,2)').^2;	
		disT = sqrt(disT);
		miniDisList2SrcPSL = min(disT);
		tarPositions = find(miniDisList2SrcPSL<mergeTrigger_/disThreshold);		
		%if isempty(tarPositions)
		if length(tarPositions)/size(tarCoordList,1)<relaxedThreshold
			startPos = 1;
		elseif length(tarPositions) == size(tarCoordList,1)
			startPos = srcPSL.midPointPosition;
		else		
			startPos = size(tarCoordList,1) - min(tarPositions) + 1;
		end
	end
	tarPSL.eleIndexList = srcPSL.eleIndexList(startPos:endPos,:);
	tarPSL.phyCoordList = srcPSL.phyCoordList(startPos:endPos,:);
	tarPSL.alphaList = srcPSL.alphaList(startPos:endPos,:);
	tarPSL.length = length(tarPSL.eleIndexList);
	tarPSL.midPointPosition = 1;		
end

function [layer_1_Streamline, layer_2_Streamline, layer_3_Streamline] = GeneratePrincipalStressLines_elementWise(initialSeed, typePSL)
	global tracingStepWidth_;
    global meshStruct_;
	limiSteps = ceil(1.5*norm(meshStruct_.boundingBox(2,:)-meshStruct_.boundingBox(1,:))/tracingStepWidth_);
	layer_1_Streamline = StreamlineStruct();
	layer_2_Streamline = StreamlineStruct();
	layer_3_Streamline = StreamlineStruct();
	
	%%1. Spot the Starting Point	
	[eleIndex, phyCoord, layerPara, opt] = PreparingForTracing_elementWise(initialSeed);
	if ~opt, return; end
	
	%%2. Compute PSL(s)
	switch typePSL
		case 'ONE'
			psDir = [2 3];
			layer_1_Streamline = ComputePSL_elementWise(eleIndex, phyCoord, layerPara, psDir, layerPara(psDir), limiSteps);			
		case 'TWO'
			psDir = [5 6];			
			layer_2_Streamline = ComputePSL_elementWise(eleIndex, phyCoord, layerPara, psDir, layerPara(psDir), limiSteps);	
		case 'THREE'
			psDir = [8 9];			
			layer_3_Streamline = ComputePSL_elementWise(eleIndex, phyCoord, layerPara, psDir, layerPara(psDir), limiSteps);					
		case 'ALL'
			psDir = [2 3];
			layer_1_Streamline = ComputePSL_elementWise(eleIndex, phyCoord, layerPara, psDir, layerPara(psDir), limiSteps);	
			psDir = [5 6];
			layer_2_Streamline = ComputePSL_elementWise(eleIndex, phyCoord, layerPara, psDir, layerPara(psDir), limiSteps);
			psDir = [8 9];			
			layer_3_Streamline = ComputePSL_elementWise(eleIndex, phyCoord, layerPara, psDir, layerPara(psDir), limiSteps);					
	end
end

function iStreamline = ComputePSL_elementWise(eleIndex, phyCoord, layerPara, psDir, iniDir, limiSteps)
	global tracingStepWidth_;
	global followingLayerTabgent_;
	if followingLayerTabgent_
		iniDir = [-iniDir(2) iniDir(1)];
	end
	iStreamline = StreamlineStruct();
	PSLphyCoordList = phyCoord;
	PSLeleIndexList = eleIndex;
	PSLprincipalStressList = layerPara;
	%% tracing along first direction (v1)
	nextPoint = phyCoord + tracingStepWidth_*iniDir;
	[phyCoordList, eleIndexList, alphaList] = TracingPSL_elementWise(nextPoint, iniDir, eleIndex, psDir, limiSteps);
	PSLphyCoordList = [PSLphyCoordList; phyCoordList];
	PSLeleIndexList = [PSLeleIndexList; eleIndexList];
	PSLprincipalStressList = [PSLprincipalStressList; alphaList];
	%% tracing along second direction (-v1)			
	nextPoint = phyCoord - tracingStepWidth_*iniDir;
	[phyCoordList, eleIndexList, alphaList] = TracingPSL_elementWise(nextPoint, -iniDir, eleIndex, psDir, limiSteps);
	if size(phyCoordList,1) > 1
		phyCoordList = flip(phyCoordList);
		eleIndexList = flip(eleIndexList);
		alphaList = flip(alphaList);				
	end
	PSLphyCoordList = [phyCoordList; PSLphyCoordList];
	PSLeleIndexList = [eleIndexList; PSLeleIndexList];
	PSLprincipalStressList = [alphaList; PSLprincipalStressList];
	iStreamline.midPointPosition = size(phyCoordList,1)+1;	
	%%2.3 finish Tracing the current major PSL			
	iStreamline.length = size(PSLphyCoordList,1);
	iStreamline.eleIndexList = PSLeleIndexList;
	iStreamline.phyCoordList = PSLphyCoordList;
	iStreamline.alphaList = PSLprincipalStressList;
	if 2==psDir(1)
		iStreamline.type = 'ONE';
	elseif 5==psDir(1)
		iStreamline.type = 'TWO';
	elseif 8==psDir(1)
		iStreamline.type = 'THREE';	
	else
		error('Wrong Input!');
	end
end

function [eleIndex, phyCoord, layerPara, opt] = PreparingForTracing_elementWise(initialSeed)
    global meshStruct_;
    global nodeCoords_; 
	global topRankNparas_;
	eleIndex = 0;
	phyCoord = 0; 
	layerPara = 0;
	if 3==size(initialSeed,2)
		opt = 1;
		formatedSeed = initialSeed;
		eleIndex = formatedSeed(1,1);
	elseif 2==size(initialSeed,2)
		[eleIndex, paraCoordinates, opt] = FindAdjacentElement(initialSeed);		
		if opt
			formatedSeed = [eleIndex paraCoordinates];
		else
			return;
		end
	else
		error('Wrong Input!');
	end
	NIdx = meshStruct_.eNodMat(eleIndex,:)';
	eleNodeCoords = nodeCoords_(NIdx,:);
	paraCoord = formatedSeed(1, 2:3);
	shapeFuncs = ShapeFunction(paraCoord(1), paraCoord(2));	
	phyCoord = shapeFuncs*eleNodeCoords;	
	layerPara = topRankNparas_(eleIndex,:);	
end

function [phyCoordList, eleIndexList, alphaList] = TracingPSL_elementWise(nextPoint, iniDir, elementIndex, typePSL, limiSteps)
	%% Tracing the PSL by 2-nd order Runge-Kutta Scheme 
	global tracingStepWidth_; 
	global topRankNparas_;
	global followingLayerTabgent_;
	phyCoordList = zeros(limiSteps,2);
	eleIndexList = zeros(limiSteps,1);
	alphaList = zeros(limiSteps,9);
	index = 0;	
	dirDetectOpt = 1;
	[elementIndex, ~, bool1] = FindAdjacentElement(nextPoint);		
	while 1==bool1
		index = index + 1; if index > limiSteps, index = index-1; break; end
		layerPara = topRankNparas_(elementIndex,:);
		if followingLayerTabgent_
			tmp = layerPara;
			layerPara([2 5 8]) = -tmp([3 6 9]);
			layerPara([3 6 9]) = tmp([2 5 8]);
		end
		if dirDetectOpt
			nextDir = DirectionSelecting_elementWise(iniDir, layerPara(typePSL), -layerPara(typePSL));
		else
			dirList = [layerPara([2 3]); -layerPara([2 3]); layerPara([5 6]); -layerPara([5 6]); layerPara([8 9]); -layerPara([8 9])];
			nextDir = DirectionSelecting_elementWise_new(iniDir, dirList);
		end
		if 0 == AngleTerminationCondition_elementWise(iniDir, nextDir), index = index-1; break; end	
		iniDir = nextDir;
		phyCoordList(index,:) = nextPoint;
		eleIndexList(index,:) = elementIndex;
		alphaList(index,:) = layerPara;					
		nextPoint = nextPoint + tracingStepWidth_*iniDir;
		[elementIndex, ~, bool1] = FindAdjacentElement(nextPoint);
	end	
	phyCoordList = phyCoordList(1:index,:);
	eleIndexList = eleIndexList(1:index,:);
	alphaList = alphaList(1:index,:);				
end

function val = AngleTerminationCondition_elementWise(dirct1, dirct2)
    global interceptionThreshold_;
	angle = acos((dirct1*dirct2') / (norm(dirct1)*norm(dirct2)));
	if angle > pi/interceptionThreshold_
		val = 0;
	else
		val = 1;
	end
end

function targetDirection = DirectionSelecting_elementWise(originalVec, Vec1, Vec2)
	angle1 = acos(originalVec*Vec1');
	angle2 = acos(originalVec*Vec2');
	if angle1 < angle2
		targetDirection = Vec1;
	else
		targetDirection = Vec2;
	end
end

function targetDirection = DirectionSelecting_elementWise_new(originalVec, dirList)
	angList = acos(dirList*originalVec');
	[~, minAngPos] = min(angList);
	targetDirection = dirList(minAngPos,:);
end

function val = StreamlineStruct()
	val = struct(...
		'ith',						0, 	...
		'type',						[],	...
		'length',					0,	...
		'midPointPosition',			0,	...		
		'phyCoordList',				[], ...
		'eleIndexList',				[], ...
		'alphaList',				[] ...
	);	
end

function [nextElementIndex, paraCoordinates, opt] = FindAdjacentElement(physicalCoordinates)
	global eleSize_;
	global nodeCoords_; 
	global meshStruct_;
	Lbound = meshStruct_.boundingBox(1,:);
	nextElementIndex = 0; paraCoordinates = []; opt = 0;
	
	physicalCoordinates = physicalCoordinates - Lbound;
	if 0==physicalCoordinates(1)
		eleX = 1;				
	else
		eleX = ceil(physicalCoordinates(1)/eleSize_);
		if eleX<1 || eleX>meshStruct_.resX, return; end
	end
	if 0==physicalCoordinates(2)
		eleY = 1;
	else
		eleY = ceil(physicalCoordinates(2)/eleSize_);
		if eleY<1 || eleY>meshStruct_.resY, return; end
	end				
	
	tarEle = meshStruct_.resY*(eleX-1)+(meshStruct_.resY-eleY+1);
	nextElementIndex = meshStruct_.eleMapForward(tarEle);
	if nextElementIndex	
		opt = 1;
		relatedNodes = meshStruct_.eNodMat(nextElementIndex,:);
		relatedNodeCoords = nodeCoords_(relatedNodes',:)-Lbound;
		paraCoordinates = 2*(physicalCoordinates - relatedNodeCoords(1,:)) / eleSize_ - 1;
	end
end

function N = ShapeFunction(s, t)
	%				   	   __s (parametric coordinate system)
	%				  	 /-t
	%				*4			*3
	%			*1			*2
	%
	%				nodes
	s = s(:);
	t = t(:);
	N = zeros(size(s,1), 4);
	N(:,1) = 0.25*(1-s).*(1-t);
	N(:,2) = 0.25*(1+s).*(1-t);
	N(:,3) = 0.25*(1+s).*(1+t);
	N(:,4) = 0.25*(1-s).*(1+t);	
end

function ShowStreamlines()
	global meshStruct_;
	global nodeCoords_;
	global topRankNparas_;
	global layer_1_StreamlinePool_;
	global layer_2_StreamlinePool_;
	global layer_3_StreamlinePool_;
	
	lw = 2;
	figure;	axHandle_ = gca;
	
	%%Show Density Layout
	xPatchs = nodeCoords_(:,1); xPatchs = xPatchs(meshStruct_.eNodMat');
	yPatchs = nodeCoords_(:,2); yPatchs = yPatchs(meshStruct_.eNodMat');
	layerThickness = topRankNparas_(:, [1 4 7]);
	cPatchs = 1 - prod(1-layerThickness,2);
    [~,pos] = min(cPatchs); cPatchs(pos) = 0;
	cPatchs = repmat(cPatchs, 1, 4)';
	colormap(axHandle_, flip(gray));
	hd = patch(axHandle_, xPatchs, yPatchs, cPatchs); hold(axHandle_, 'on');
	
	%%Show Streamlines
	for jj=1:length(layer_1_StreamlinePool_)
		hold(axHandle_, 'on');
		plot(axHandle_, layer_1_StreamlinePool_(jj).phyCoordList(:,1), layer_1_StreamlinePool_(jj).phyCoordList(:,2), ...
			'-', 'color', [252 141 98]/255, 'LineWidth', lw);
	end	
	for jj=1:length(layer_2_StreamlinePool_)
		hold(axHandle_, 'on');
		plot(axHandle_, layer_2_StreamlinePool_(jj).phyCoordList(:,1), layer_2_StreamlinePool_(jj).phyCoordList(:,2), ...
			'-', 'color', [102 194 165]/255, 'LineWidth', lw);
	end	
	for jj=1:length(layer_3_StreamlinePool_)
		hold(axHandle_, 'on');
		plot(axHandle_, layer_3_StreamlinePool_(jj).phyCoordList(:,1), layer_3_StreamlinePool_(jj).phyCoordList(:,2), ...
			'-', 'color', [0 0 1]/1, 'LineWidth', lw); %%[141 160 203]/255
	end	
	%%Show silhouette
	if 0
		edgeIndices = meshStruct_.eNodMat(:, [1 2 2 1  2 3 3 2  3 4 4 3  4 1 1 4])';
		edgeIndices = reshape(edgeIndices(:), 4, 4*meshStruct_.numElements);	
		tmp = zeros(meshStruct_.numNodes,1); 
		tmp(meshStruct_.nodesOnBoundary) = 1;
		tmp = tmp(edgeIndices); 
		tmp = sum(tmp,1);
		boundaryEleEdges = edgeIndices(:,find(4==tmp));
		domainSilhouette_.vertices = nodeCoords_;
		domainSilhouette_.faces = boundaryEleEdges';
		
		hold(axHandle_, 'on');
		hdSilho = patch(axHandle_, domainSilhouette_);
		set(hdSilho, 'FaceColor', 'None', 'EdgeColor', [0.5 0.5 0.5], 'LineWidth', 2);
	end
	set(hd, 'FaceColor', 'flat', 'FaceAlpha', 1, 'EdgeColor', 'None');
	axis(axHandle_, 'equal'); axis(axHandle_, 'tight'); axis(axHandle_, 'off');
end

