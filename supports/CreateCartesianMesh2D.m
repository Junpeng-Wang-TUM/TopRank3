function meshInfo = CreateCartesianMesh2D(voxelizedModel)	
	%    __ x
	%   / 
	%  -y         
	%		4--------3
	%	    |		 |		
	%		|		 |
	%		1--------2
	%	rectangular element	

	%% initialize mesh
	meshInfo = MeshStruct();
	[ny, nx] = size(voxelizedModel);
	meshInfo.resX = nx;	meshInfo.resY = ny;
	meshInfo.boundingBox = [0 0; [nx ny]];
	meshInfo.eleSize = [meshInfo.boundingBox(2,:) - meshInfo.boundingBox(1,:)] ./ [nx ny];
	
	%% identify solid&void elements
	meshInfo.eleMapBack = find(1==voxelizedModel);
	meshInfo.numElements = length(meshInfo.eleMapBack);
	meshInfo.eleMapForward = zeros(nx*ny,1);	
	meshInfo.eleMapForward(meshInfo.eleMapBack) = (1:meshInfo.numElements)';
	
	%% discretize
	nodenrs = reshape(1:(nx+1)*(ny+1), 1+ny, 1+nx);
	eNodVec = reshape(nodenrs(1:end-1,1:end-1)+1, nx*ny, 1);
	meshInfo.eNodMat = repmat(eNodVec(meshInfo.eleMapBack),1,4);
	tmp = [0 ny+[1 0] -1]; 
	for ii=1:4
		meshInfo.eNodMat(:,ii) = meshInfo.eNodMat(:,ii) + repmat(tmp(ii), meshInfo.numElements,1);
	end
	meshInfo.nodMapBack = unique(meshInfo.eNodMat);
	meshInfo.numNodes = length(meshInfo.nodMapBack);
	meshInfo.numDOFs = meshInfo.numNodes*2;
	meshInfo.nodMapForward = zeros((nx+1)*(ny+1),1);
	meshInfo.nodMapForward(meshInfo.nodMapBack) = (1:meshInfo.numNodes)';
	for ii=1:4
		meshInfo.eNodMat(:,ii) = meshInfo.nodMapForward(meshInfo.eNodMat(:,ii));
	end
	meshInfo.eDofMat = [2*meshInfo.eNodMat-1 2*meshInfo.eNodMat];
	meshInfo.eDofMat = meshInfo.eDofMat(:, [1 5 2 6 3 7 4 8]);
	meshInfo.nodeCoords = zeros((nx+1)*(ny+1),2);
	xSeed = meshInfo.boundingBox(1,1):(meshInfo.boundingBox(2,1)-meshInfo.boundingBox(1,1))/nx:meshInfo.boundingBox(2,1);
	ySeed = meshInfo.boundingBox(2,2):(meshInfo.boundingBox(1,2)-meshInfo.boundingBox(2,2))/ny:meshInfo.boundingBox(1,2);		
	meshInfo.nodeCoords(:,1) = reshape(repmat(xSeed, ny+1, 1), (nx+1)*(ny+1), 1);
	meshInfo.nodeCoords(:,2) = repmat(ySeed, 1, nx+1)';

	%% identify boundary info.
	meshInfo.numNodsAroundEleVec = zeros(meshInfo.numNodes,1);
	for ii=1:meshInfo.numElements
		iNodes = meshInfo.eNodMat(ii,:);
		meshInfo.numNodsAroundEleVec(iNodes,:) = meshInfo.numNodsAroundEleVec(iNodes) + 1;
	end
	meshInfo.nodesOnBoundary = find(meshInfo.numNodsAroundEleVec<4);
	allNodes = zeros(meshInfo.numNodes,1);
	allNodes(meshInfo.nodesOnBoundary) = 1;	
	tmp = zeros(meshInfo.numElements,1);
	for ii=1:4
		tmp = tmp + allNodes(meshInfo.eNodMat(:,ii));
	end
	meshInfo.elementsOnBoundary = find(tmp>0);	
	
	%% compute element centroids
	eleCentX = meshInfo.nodeCoords(meshInfo.nodMapBack,1); eleCentX = eleCentX(meshInfo.eNodMat);
	eleCentY = meshInfo.nodeCoords(meshInfo.nodMapBack,2); eleCentY = eleCentY(meshInfo.eNodMat);	
	meshInfo.eleCentroidList = [sum(eleCentX,2) sum(eleCentY,2)]/4;
	meshInfo.boundaryNodeCoords = meshInfo.nodeCoords(meshInfo.nodMapBack(meshInfo(1).nodesOnBoundary),:);
end

function val = MeshStruct()
	val = struct(...
		'resX',							0,	...
		'resY',							0,	...
		'boundingBox',					[],	...
		'eleSize',						[],	...
		'numElements',					0,	...		
		'numNodes',						0,	...
		'numDOFs',						0,	...
		'nodeCoords',					0,	...
		'eleCentroidList',				0,	...
		'eNodMat',						0,	...
		'eDofMat',						0,	...
		'freeDOFs',						0,	...
		'fixedDOFs',					0,	...
		'eleMapBack',					[],	...
		'eleMapForward',				[],	...
		'nodMapBack',					[],	...
		'nodMapForward',				[],	...
		'nodesOnBoundary',				[],	...
		'boundaryNodeCoords',			[],	...
		'numNodsAroundEleVec',			[],	...
		'elementsOnBoundary',			[],	...
		'elementUpwardMap',				[]	...
	);
end