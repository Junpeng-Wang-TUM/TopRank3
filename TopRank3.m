%% This code is created for performing the 2D homogenization-based topology optimization subject to multiple loading conditions
%% Author: Junpeng Wang
%% Copyright (c) All Rights Reserved.
%% Version: 2023-11-10
%% E-mail: junpeng.wang@tum.de
function [optiAlphaList, optiThetaList, optiLog, cValuesPerLoadHist] = TopRank3(mesh_, fixingCond_, loadingCond_, matProp_, optiSet_, outPath)
	%%1. SETTINGS
	%%1.1 Material Properties
	%%Only 3==RN is supported for optimization, others can only be used for FEA. To adapt to optimization involving a higner rank,
	%%works in derivative computation needs to be extended, typically on the computation of derivatives of "M" in Rank-N
	RN = 3; RN = max(2, RN); 
	E0 = matProp_.modulus;
	Emin = matProp_.minModulus;
	v0 = matProp_.PossionRatio;
	%%1.2 Geometric Modeling
	nelX = mesh_.resX;
	nelY = mesh_.resY;
	numElements = mesh_.numElements;
	eDofMat = mesh_.eDofMat;
	eleIndexMapForward = mesh_.eleMapForward;
	eleMapBack = mesh_.eleMapBack;
	eleCentroidList = mesh_.eleCentroidList;
	eleSize = mesh_.eleSize(1);
	%%1.3 Apply for Boundary Condition
	[F, freeDOFs] = ApplyBoundaryCondition(mesh_, fixingCond_, loadingCond_);

	%%1.4 Optimization Description
	move = 0.01;
	thetaMove = 1.0*pi/180;
	Vf = optiSet_.Vf;
	wqList = optiSet_.perLoadContri;
	rMin = 2.6;
	nLoop = optiSet_.nloop;
	passiveEles = optiSet_.passiveEles;	
	
	OPT_OriStartingGuess = optiSet_.OPT_OriStartingGuess; 
	OPT_Theta3Only = optiSet_.OPT_Theta3Only; 
	W = optiSet_.W; %%weight factor for compliance item, 1-W for orientation regularization item
	alphaUpperBound = optiSet_.alphaBounds(2);
	alphaLowerBound = optiSet_.alphaBounds(1);
	if 3~=RN, OPT_Theta3Only = 0; end

	%%2. Setup
	%%2.1 FEA
	strainMatrixGP = ObtainStrainMatrix();
	Stensor4thSolid = HookeLaw4th(E0, v0); 
	Stensor4thVoid = HookeLaw4th(Emin, v0); 
	[eK2gKrowMap, eK2gKcolMap, ~] = find(ones(8));
	iK = eDofMat(:,eK2gKrowMap)'; iK = iK(:);
	jK = eDofMat(:,eK2gKcolMap)'; jK = jK(:);

	%%2.2 Rank-N Theory
	Stensor2ndSolid = Tensor4thToTensor2nd(Stensor4thSolid);
	Stensor2ndVoid = Tensor4thToTensor2nd(Stensor4thVoid);
	C0 = inv(Stensor2ndSolid-Stensor2ndVoid);
	S0Plus = E0 / (1-v0^2);	
	
	%%2.3 Filter and Orientation Regularization
	[H, Hs] = BuildDensityFilter(rMin, nelX, nelY, numElements, eleIndexMapForward, eleCentroidList, eleSize);
	eleAdjacencyMap = BuildElementAdjacency(numElements, nelX, nelY, eleIndexMapForward, eleMapBack);
	
	
	%%2.4 Evaluate Compliance of Fully Solid Domain
	[U, ~, ce] = PerformFEA(Stensor4thSolid, F, freeDOFs, strainMatrixGP, iK, jK, numElements, eDofMat);
	cSolid = sum(ce,1)*wqList(:);
	disp(['Compliance of Fully Solid Domain: ' sprintf('%16.6e',cSolid)]);	
	
	%%2.5 Ini. for Starting Guess
	alphaList = repmat(Vf/RN, numElements, RN);
	if 3==RN
		switch OPT_OriStartingGuess
			case 'SGu'
				iniTheta3 = pi/2; %%Normal of Layer-3
				thetaList = OrientationStartingGuess_Uniform(numElements, iniTheta3);
			case 'SGp'
				thetaList = OrientationStartingGuess_PrincipalStress(numElements, eDofMat, U, Stensor4thSolid);
			otherwise
				error('Wrong Option for Starting Guess Strategy!');				
		end
	else
        thetaTemp = linspace(0,pi,RN+1); 
		thetaList = repmat(thetaTemp(1:end-1), numElements, 1); %%Initialization by User
		warning('Not Rank-3 Laminate! Optimization Functionality not Supported!');	
	end
	alphaList(passiveEles,:) = ones(size(passiveEles,1), RN);
	
	%%3. Preparation for Optimization Loop
	%%3.1 Setup 	
	alphaTildeList = alphaList;	
	thetaTildeList = thetaList;
	optiAlphaList = alphaTildeList;
	optiThetaList = thetaTildeList;
	rhoRN = 1-prod(1-alphaTildeList,2);
	%%3.2 Ini. for Orientation Regularization
	%%Initial Values of Different Items in Objective Function 
	Mmat = ComputeMmat(alphaTildeList, thetaTildeList, v0, RN);
	Stensor4thRankN = AeembleStensor4thRankN(Mmat, Stensor2ndSolid, rhoRN, C0, S0Plus);
	[~, ~, ce] = PerformFEA(Stensor4thRankN, F, freeDOFs, strainMatrixGP, iK, jK, numElements, eDofMat);
	c0Total = sum(ce*wqList(:),1); 
	[p0Total, ~] = GetLayerOrientationDeviationMetric(thetaTildeList, eleAdjacencyMap, OPT_Theta3Only);
	p0Total = max(p0Total,1);
	%%3.3 Ini. for Optimizer
	if 1==OPT_Theta3Only
		xold1_theta = thetaList(:,3);
	else
		xold1_theta = thetaList(:);		
	end
	xold2_theta = xold1_theta;
	low = 0;
	upp = 0;
	loop = 0;

	%%3.4 Ini. for Output
	oHist = [];
	cHist = [];
	pHist = [];
	volHist = [];
	cValuesPerLoadHist = [];
	designVariablesRank3_ = [alphaList thetaList];
	designVariablesRank3_ = ProjectTheta0to360(designVariablesRank3_);
	fileName = sprintf(strcat(outPath, 'designVariablesRank3-It-%d.mat'), loop);
	save(fileName, 'designVariablesRank3_');

	%%4. Optimization Loop
	while loop < nLoop
		loop = loop + 1;
		%%4.1 FEA, Objective, and Sensitivity Analysis
		Mmat = ComputeMmat(alphaTildeList, thetaTildeList, v0, RN);
		Stensor4thRankN = AeembleStensor4thRankN(Mmat, Stensor2ndSolid, rhoRN, C0, S0Plus);
		[~, Ue2, ce] = PerformFEA(Stensor4thRankN, F, freeDOFs, strainMatrixGP, iK, jK, numElements, eDofMat);
		cValuesPerLoad = sum(ce, 1);
		ce0 = cValuesPerLoad*wqList(:);
		[pe0, ~] = GetLayerOrientationDeviationMetric(thetaTildeList, eleAdjacencyMap, OPT_Theta3Only);
		f0val = W*ce0/c0Total + (1-W)*pe0/p0Total;
		
		[df0dAlpha, df0dTheta] = Def0DeX(alphaTildeList, thetaTildeList, Mmat, Ue2, wqList, strainMatrixGP, ...
									C0, S0Plus, v0, c0Total, p0Total, OPT_Theta3Only, W, RN, eleAdjacencyMap);

		fval = sum(1-prod(1-alphaTildeList,2))/(numElements*Vf) - 1;
		dfdAlpha = zeros(size(alphaTildeList));	
		allRanks = 1:RN;
		for jj=1:RN
			restRanks = setdiff(allRanks, jj);
			dfdAlpha(:,jj) = prod(1-alphaTildeList(:,restRanks),2);
		end
		dfdAlpha = dfdAlpha/(numElements*Vf);
		dfdTheta = 0.0e-12*ones(size(df0dTheta)); %% Constant Small Values	

		%% 4.2 Filtering
		df0dAlpha = H*(df0dAlpha./Hs);
		dfdAlpha = H*(dfdAlpha./Hs);
		df0dTheta = H*(df0dTheta./Hs);

		%% 4.3 Update of Design Variables by Optimizer
		%% 4.3.1 Alpha 
		l1 = 0; l2 = 1e5;
		while (l2-l1)/(l1+l2) > 1e-3
			lmid = 0.5*(l2+l1);
			alphaListNew = max(alphaLowerBound, max(alphaList-move, min(alphaUpperBound, min(alphaList+move, alphaList.*real(sqrt(-df0dAlpha./dfdAlpha/lmid))))));
			tmp = alphaListNew;
			alphaTildeList = H * tmp ./ Hs;								
			alphaTildeList(passiveEles,:) = ones(size(passiveEles,1), RN);		
			fval = sum(1-prod(1-alphaTildeList,2))/(numElements*Vf) - 1;
			if fval > 0, l1 = lmid; else l2 = lmid; end
		end
		changeAlpha = max([abs(alphaListNew(:)-alphaList(:))]);
		alphaList = alphaListNew;
		rhoRN = 1-prod(1-alphaTildeList,2);
		optiAlphaList = alphaTildeList;
		%% 4.3.2 Theta
		m = 1;
		iter = loop;	
		a0 = 1;
		a = zeros(m,1);     
		c_ = ones(m,1)*1000;
		d = zeros(m,1);		
		if 1==OPT_Theta3Only
			xval_MMA_theta = thetaList(:,3);	
		else		
			xval_MMA_theta = thetaList(:);			
		end				
		xold1_MMA_theta = xold1_theta;
		xold2_MMA_theta = xold2_theta;
		df0dx_MMA_theta = df0dTheta(:);
		df0dx2_MMA_theta = zeros(size(df0dx_MMA_theta));
		dfdx_MMA_theta = dfdTheta(:);
		dfdx2_MMA_theta = zeros(size(dfdx_MMA_theta));
		n = numel(xval_MMA_theta);
		xmin_MMA_theta = xval_MMA_theta;
		xmax_MMA_theta = xval_MMA_theta;
		
		xmax_MMA_theta = min(4*pi, xmax_MMA_theta+thetaMove);
		xmin_MMA_theta = max(-4*pi, xmin_MMA_theta-thetaMove);			
		
		[xmma_MMA_theta,~,~,~,~,~,~,~,~,low,upp] = ...
			mmasub(m,n,iter,xval_MMA_theta,xmin_MMA_theta,xmax_MMA_theta,xold1_MMA_theta,xold2_MMA_theta,...
				f0val,df0dx_MMA_theta,df0dx2_MMA_theta,fval,dfdx_MMA_theta',dfdx2_MMA_theta',low,upp,a0,a,c_,d);
		
		xold1_theta = xval_MMA_theta;
		xold2_theta = xold1_MMA_theta;
		
		if 1==OPT_Theta3Only
			thetaListNew = reshape(xmma_MMA_theta, numElements, 1);
			changeTheta = max([abs(thetaListNew(:)-thetaList(:,3))]);
			thetaList = [thetaListNew+2*pi/3 thetaListNew+pi/3 thetaListNew];
		else
			thetaListNew = reshape(xmma_MMA_theta, numElements, RN);
			changeTheta = max([abs(thetaListNew(:)-thetaList(:))]);
			thetaList = thetaListNew;			
		end
		thetaTildeList = thetaList;
		thetaTildeList = H * thetaTildeList ./ Hs;							
		optiThetaList = thetaTildeList;
		change = max([changeAlpha changeTheta]);
		
		%% 4.4 Print Results
		disp([' It.: ' sprintf('%4i',loop) ' Obj.: ' sprintf('%10.4e',f0val) ' C0: ' sprintf('%10.4e',ce0) ' P0: ' sprintf('%10.4e',pe0) ...
				' Vol.: ' sprintf('%6.3f',sum(rhoRN)/numElements) ' Cons.: ' sprintf('%10.4e ',fval)]);
		
		%% 4.5 Plot and Write history
		oHist(loop,1) = f0val;
		cHist(loop,1) = ce0;
		pHist(loop,1) = pe0;
		volHist(loop,1) = sum(rhoRN)/numElements;
		cValuesPerLoadHist(loop,:) = cValuesPerLoad;
		designVariablesRank3_ = [alphaTildeList thetaTildeList];
		designVariablesRank3_ = ProjectTheta0to360(designVariablesRank3_);
		fileName = sprintf(strcat(outPath, 'designVariablesRank3-It-%d.mat'), loop);
		save(fileName, 'designVariablesRank3_');
		xPhys = zeros(nelX*nelY,1);
		xPhys(eleMapBack) = rhoRN;
        figure(2);clf;
		colormap(gray); 
		imagesc(-reshape(xPhys, nelY, nelX), [-1 0]); axis('equal'); axis('off'); drawnow;		
	end
	optiLog = [oHist volHist cHist pHist];
	
	%% 5. Show Opti. Result
	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%FEA Related
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F_, freeDOFs_] = ApplyBoundaryCondition(mesh_, fixingCond_, loadingCond_)	
	%%Loading
	nodDOF = 2;
	loadingCond_ = loadingCond_(:);
	numLoadSteps = numel(loadingCond_);
	validLoadSteps = [];
	for ii=1:numLoadSteps
		if ~isempty(loadingCond_(ii).arr), validLoadSteps(end+1,1) = ii; end
	end
	loadingCond_ = loadingCond_(validLoadSteps);
	numLoadSteps = numel(loadingCond_);

	if numLoadSteps>0
		F_ = sparse(mesh_(1).numDOFs, numLoadSteps);
		for ii=1:numLoadSteps
			iFVec = sparse(mesh_(1).numNodes, nodDOF);
			iFVec(loadingCond_(ii).arr(:,1),:) = loadingCond_(ii).arr(:,2:end);
			F_(:,ii) = reshape(iFVec',mesh_(1).numDOFs,1);
		end		
	else
		warning('No Loads!'); return;
	end
	
	%%Fixing
	if ~isempty(fixingCond_(1).arr)
		numFixations = numel(fixingCond_);
		flexibleArr = struct('arr', []);
		freeDOFs_ = repmat(flexibleArr, numFixations, 1);
		for ii=1:numFixations
			fixedDOFs = nodDOF*double(fixingCond_(ii).arr(:,1));
			fixedDOFs = fixedDOFs - (nodDOF-1:-1:0);
			fixedDOFs = reshape(fixedDOFs', numel(fixedDOFs), 1);	
			fixingState = fixingCond_(ii).arr(:,2:end)';
			realFixedDOFs = fixedDOFs(find(1==fixingState(:)));
			freeDOFs = setdiff((1:int32(mesh_(1).numDOFs))',realFixedDOFs);		
			freeDOFs_(ii).arr = freeDOFs;			
		end
	else
		warning('No Constraint!'); return;			
	end		
end

function Bmat = ObtainStrainMatrix()
	q3 = sqrt(3)/3;
	Bmat = [  	-1-q3 0 1+q3 0 1-q3 0 -1+q3 0; 0 -1-q3 0 -1+q3 0 1-q3 0 1+q3; -1-q3 -1-q3 -1+q3 1+q3 1-q3 1-q3 1+q3 -1+q3; 
				-1-q3 0 1+q3 0 1-q3 0 -1+q3 0; 0 -1+q3 0 -1-q3 0 1+q3 0 1-q3; -1+q3 -1-q3 -1-q3 1+q3 1+q3 1-q3 1-q3 -1+q3;
				-1+q3 0 1-q3 0 1+q3 0 -1-q3 0; 0 -1+q3 0 -1-q3 0 1+q3 0 1-q3; -1+q3 -1+q3 -1-q3 1-q3 1+q3 1+q3 1-q3 -1-q3;
				-1+q3 0 1-q3 0 1+q3 0 -1-q3 0; 0 -1-q3 0 -1+q3 0 1-q3 0 1+q3; -1-q3 -1+q3 -1+q3 1-q3 1-q3 1+q3 1+q3 -1-q3]/2;
end

function Stensor = HookeLaw4th(E, v)
	%%Hooke's law
	Stensor = E/(1-v^2)*[1 v 0; v 1 0; 0 0 (1-v)/2];
end

function [eleKarr, eleKmat] = ComputeElementStiffMatrix(S, B)
	%%Integral(B' * D * B)
	w = ones(1,12)*0.25; %%0.25 being specific to unit element
	eleD = zeros(12);
	for ii=1:4
		index = (ii-1)*3+1:ii*3;
		eleD(index,index) = S;
	end
	eleKmat = B'*(eleD.*w)*B;
	eleKarr = eleKmat(:);
end

function [U, Ue2, ce] = PerformFEA(Stensor4th, F, freeDOFs, strainMatrixGP, iK, jK, numElements, eDofMat)
	%%Solving
	if 1 == size(Stensor4th,3)
		[Ke, ~] = ComputeElementStiffMatrix(Stensor4th, strainMatrixGP);
		Ks = repmat(Ke,1,numElements);
	elseif numElements == size(Stensor4th,3)
		Ks = zeros(8*8, numElements);
		for ii=1:numElements
			iD = Stensor4th(:,:,ii);
			[Ks(:,ii), ~] = ComputeElementStiffMatrix(iD, strainMatrixGP);       
		end		
	else
		error('Wrong Input of Elasticity Tensor!');
	end
	K = sparse(iK, jK, Ks(:));
	U = zeros(size(F));
	
	if numel(freeDOFs) == size(F,2) %%Concurrently varying loads and fixations
		for ii=1:size(F,2)
			iFreeDofs = freeDOFs(ii).arr;
			U(iFreeDofs,ii) = K(iFreeDofs,iFreeDofs)\F(iFreeDofs,ii);
		end		
	elseif 1==numel(freeDOFs) && size(F,2)>1 %%Constant fixations, only the loads vary
		iFreeDofs = freeDOFs(1).arr;
		U(iFreeDofs,:) = K(iFreeDofs,iFreeDofs)\F(iFreeDofs,:);
	else
		error('Inconsistent Boundary Conditions!')
    end
	
	%%PostP
	ce = zeros(numElements, size(U,2));
	Ue2 = zeros(64, numElements, size(U,2));
	for jj=1:size(U,2)
		ithU = U(:,jj);
		ithUe = ithU(eDofMat');
		ithUe2 = reshape(bsxfun(@times, reshape(ithUe,[1 8 numElements]), reshape(ithUe, [8 1 numElements])), 64, numElements);
		ce(:,jj) = reshape(sum(Ks.*ithUe2,1),numElements,1);
		Ue2(:,:,jj) = ithUe2;
	end	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Rank-N related
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mmat = ComputeMmat(alphaList, thetaList, v0, RN)
	rhoRN = 1-prod(1-alphaList,2);
	numElements = size(alphaList,1);
	PnList = alphaList;
	for ii=RN-1:-1:1
		PnList(:,ii) = prod(1-alphaList(:,ii+1:RN),2) .* alphaList(:,ii);
	end
	PnList = PnList ./ rhoRN;
	
	h = 4*(1-v0);
	l1 = cos(2*thetaList);
	l2 = sin(2*thetaList);
	l3 = cos(4*thetaList);
	l4 = sin(4*thetaList);
	m1 = sum(PnList.*l1, 2); m1 = reshape(m1, 1, 1, numElements);
	m2 = sum(PnList.*l2, 2); m2 = reshape(m2, 1, 1, numElements);
	m3 = sum(PnList.*l3, 2); m3 = reshape(m3, 1, 1, numElements);
	m4 = sum(PnList.*l4, 2); m4 = reshape(m4, 1, 1, numElements);
	
	Mmat = zeros(3,3,numElements);
	Mmat(1,1,:) = (3-v0 - (1+v0)*m3) / h;
	Mmat(1,2,:) = -(1+v0)*m4 / h;
	Mmat(1,3,:) = m1 /2;
	Mmat(2,2,:) = (3-v0 + (1+v0)*m3) / h;
	Mmat(2,3,:) = m2 /2;
	Mmat(3,3,:) = 1/2;
	Mmat(2,1,:) = Mmat(1,2,:);
	Mmat(3,1,:) = Mmat(1,3,:);
	Mmat(3,2,:) = Mmat(2,3,:);	
end

function Stensor4thRankN = AeembleStensor4thRankN(Mmat, Stensor2ndSolid, rhoRN, C0, S0Plus)
	Stensor4thRankN = zeros(size(Mmat));
	numElements = size(Mmat,3);
	for ii=1:numElements
		iRho = rhoRN(ii);
		iStensor2ndRankN = Stensor2ndSolid - (1-iRho)*inv(C0-iRho/S0Plus*Mmat(:,:,ii));
		Stensor4thRankN(:,:,ii) = Tensor2ndToTensor4th(iStensor2ndRankN);
	end
end

function A4th = Tensor2ndToTensor4th(A2nd)
	%% A2nd
    %%	A11		A12		A13
    %%			A22		A23
    %%	syms			A33
	%% outA
    %%	A1111		A1122		A1112
    %%				A2222		A2221
    %%	syms					A1212	
    A4th = zeros(size(A2nd));
	A4th(1,1,:) = 0.5*(A2nd(1,1,:)+A2nd(3,3,:))+A2nd(1,3,:);
	A4th(1,2,:) = -0.5*(A2nd(1,1,:)-A2nd(3,3,:));
	A4th(1,3,:) = 0.5*(A2nd(1,2,:)+A2nd(2,3,:));
	A4th(2,2,:) = 0.5*(A2nd(1,1,:)+A2nd(3,3,:))-A2nd(1,3,:);
	A4th(2,3,:) = -0.5*(A2nd(1,2,:)-A2nd(2,3,:));
	A4th(3,3,:) = 0.5*A2nd(2,2,:);
	A4th(2,1,:) = A4th(1,2,:);
	A4th(3,1,:) = A4th(1,3,:);
	A4th(3,2,:) = A4th(2,3,:);	
end

function A2nd = Tensor4thToTensor2nd(A4th)
	%% A4th
    %%	A1111		A1122		A1112
    %%				A2222		A2221
    %%	syms					A1212
	%% A2nd
    %%	A11		A12		A13
    %%			A22		A23
    %%	syms			A33	
	A2nd = zeros(size(A4th));
	A2nd(1,1,:) = 0.5*(A4th(1,1,:)+A4th(2,2,:))-A4th(1,2,:);
	A2nd(1,2,:) = A4th(1,3,:)-A4th(2,3,:);
	A2nd(1,3,:) = 0.5*(A4th(1,1,:)-A4th(2,2,:));
	A2nd(2,2,:) = 2*A4th(3,3,:);
	A2nd(2,3,:) = A4th(1,3,:)+A4th(2,3,:);
	A2nd(3,3,:) = 0.5*(A4th(1,1,:)+A4th(2,2,:))+A4th(1,2,:);
	A2nd(2,1,:) = A2nd(1,2,:);
	A2nd(3,1,:) = A2nd(1,3,:);
	A2nd(3,2,:) = A2nd(2,3,:);	
end

function alphaAndThetaO = ProjectTheta0to360(alphaAndThetaI)
	alphaAndThetaO = alphaAndThetaI;
	thetaAng = alphaAndThetaI(:,4:6)/pi*180;
	thetaAng = mod(thetaAng,360);
	alphaAndThetaO(:,4:6) = thetaAng/180*pi;	
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Optimization related
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [H, Hs] = BuildDensityFilter(rMin, nelX, nelY, numElements, eleIndexMapForward, eleCentroidList, eleSize)
	p = 2;
	iH = zeros(numElements*(2*(ceil(rMin)-1)+1)^p, 1, 'int32');
	jH = iH;
	iIndex = 0;

	%%	1	4	7
	%%	2	5*	8
	%%	3	6	9	
	for ii = 1:nelX
		for jj = 1:nelY
			e1 = eleIndexMapForward((ii-1)*nelY+jj);
			if e1
				for ii2 = max(ii-(ceil(rMin)-1),1):min(ii+(ceil(rMin)-1),nelX)
					for jj2 = max(jj-(ceil(rMin)-1),1):min(jj+(ceil(rMin)-1),nelY)
						e2 = eleIndexMapForward((ii2-1)*nelY+jj2);
						if e2
							iIndex = iIndex + 1;
							iH(iIndex) = e1;
							jH(iIndex) = e2;									
						end
					end
				end
			end
		end
	end
	iH = iH(1:iIndex,:);
	jH = jH(1:iIndex,:);
	
	% eleSize = 1;
	sH = rMin*eleSize - vecnorm(eleCentroidList(iH,:)-eleCentroidList(jH,:),2,2);
	sH(sH<0) = 0;	
	H = sparse(iH, jH, sH, numElements, numElements);
	Hs = sum(H,2);	
end

function eleAdjacencyMap = BuildElementAdjacency(numElements, nelX, nelY, eleMapForward, eleMapBack)
	eleAdjacencyMap = struct('adjacentEles', []);
	eleAdjacencyMap = repmat(eleAdjacencyMap, numElements, 1);
	eleMapBackTmp = double(eleMapBack);
	for ii=1:numElements
		iEleGlobalIdx = eleMapBackTmp(ii);
		adjEleLeft = iEleGlobalIdx-nelY;
		adjEleRight = iEleGlobalIdx+nelY;
		adjEleUp = iEleGlobalIdx-1;
		adjEleDown = iEleGlobalIdx+1;
		if 1==mod(iEleGlobalIdx, nelY) %%Top Row
			adjEleUp = [];
		end
		if 0==mod(iEleGlobalIdx, nelY) %%Bottom Row
			adjEleDown = [];
		end
		if iEleGlobalIdx<=nelY %%First Column
			adjEleLeft = [];
		end
		if iEleGlobalIdx>(nelX-1)*nelY %%Last Column
			adjEleRight = [];
		end
		adjacentElesGlobal = [adjEleLeft adjEleRight adjEleUp adjEleDown]';
		adjacentElesLocal = eleMapForward(adjacentElesGlobal);
		adjacentElesLocal(0==adjacentElesLocal) = [];
		eleAdjacencyMap(ii).adjacentEles = adjacentElesLocal;
	end
end

function [pTotal, pe] = GetLayerOrientationDeviationMetric(thetaList, eleAdjacencyMap, OPT_Theta3Only)
	numElements = size(thetaList,1);
	pe = zeros(numElements,1);
	coefi = 6; %% 2: k*pi (min), k*pi/2 (max); 6: k*pi/3 (min), k*pi/6 (max)
	for ii=1:numElements
		iTheta = thetaList(ii,3);
		iThetasAdjEles = thetaList(eleAdjacencyMap(ii).adjacentEles,3);
		pe(ii) = sum(0.5 - 0.5 * cos(coefi*(iTheta-iThetasAdjEles)));
	end		
	if ~OPT_Theta3Only
		pe = repmat(pe,1,3);
		for ii=1:numElements
			iTheta2 = thetaList(ii,2);
			iTheta1 = thetaList(ii,1);
			iAdjacentEles = eleAdjacencyMap(ii).adjacentEles;
			iThetas2AdjEles = thetaList(iAdjacentEles,2);
			pe(ii,2) = sum(0.5 - 0.5 * cos(coefi*(iTheta2-iThetas2AdjEles)));
			iThetas1AdjEles = thetaList(iAdjacentEles,1);
			pe(ii,1) = sum(0.5 - 0.5 * cos(coefi*(iTheta1-iThetas1AdjEles)));
		end				
	end	
	pTotal = sum(sum(pe));
end

function thetaList = OrientationStartingGuess_Uniform(numElements, iniTheta3)
	thetaList = repmat(iniTheta3+[2*pi/3 pi/3 0], numElements, 1);
end

function thetaList = OrientationStartingGuess_PrincipalStress(numElements, eDofMat, U, Stensor)
	strainMatrixMid = [-1 0 1 0 1 0 -1 0; 0 -1 0 -1 0 1 0 1; -1 -1 -1 1 1 1 1 -1]/2;		
	numLoads = size(U,2);
	psPool = zeros(numElements,4*numLoads);	
	if 1==size(Stensor,3), Stensor = repmat(Stensor,1,1,numElements); end
	for ii=1:numLoads
		for jj=1:numElements
			relativeDOFsIndex = eDofMat(jj,:);
			iUe = U(relativeDOFsIndex,ii);
			iEleStress = Stensor(:,:,jj) * (strainMatrixMid*iUe);
			devTerm = sqrt( (iEleStress(1)-iEleStress(2))^2/4 + iEleStress(3)^2 );
			ampMajor = (iEleStress(1)+iEleStress(2))/2 + devTerm;
			ampMinor = (iEleStress(1)+iEleStress(2))/2 - devTerm;
			oriMajor = atan((ampMajor-iEleStress(1)) ./ iEleStress(3));
			oriMinor = oriMajor + pi/2;
			if 0==sum(iUe), ampMajor = 1; oriMajor = 0; ampMinor = 1; oriMinor = pi/2; end					
			psPool(jj,4*(ii-1)+1:4*ii) = [ampMinor oriMinor ampMajor oriMajor];
		end
	end
	psAmps = abs(psPool(:,1:2:end));
	psDir = psPool(:,2:2:end);
	[~,principalDirPos] = max(psAmps,[],2);
	principalDir = zeros(numElements,1);
	for jj=1:numElements
		principalDir(jj) = psDir(jj,principalDirPos(jj));
	end
	principalDir = principalDir + pi/2;
	thetaList = repmat(principalDir, 1, 3) + repmat([2*pi/3 pi/3 0], numElements, 1);	
end

function [df0dAlpha, df0dTheta] = Def0DeX(alphaList, thetaList, Mmat, Ue2, ww, B, C0, S0Plus, v0, c0Total, p0Total, OPT_Theta3Only, ...
									W, RN, eleAdjacencyMap)
	numElements = size(alphaList,1);
	eye3 = eye(3);
	zeros12 = zeros(12,12);
	
	%%1. Preparations
	%%Rho, dRho_dAlpha
	rhoRN = 1-prod(1-alphaList,2);
	dRho_dAlpha = zeros(numElements, RN);
	allRanks = 1:RN;
	for ii=1:RN
		restRanks = setdiff(allRanks, ii);
		dRho_dAlpha(:,ii) = prod(1-alphaList(:,restRanks),2);
	end
	
	%%Y, YY
	A0 = 1/S0Plus;
	CAM = repmat(C0,1,1,numElements) - A0*reshape(rhoRN,1,1,numElements) .* Mmat;
	Y = zeros(size(CAM));
	YY = Y;
	for ii=1:numElements
		Y(:,:,ii) = eye3 / CAM(:,:,ii);
		YY(:,:,ii) = Y(:,:,ii) / CAM(:,:,ii); %% dot or times WARNING! WARNING! WARNING!
	end

	%%RhoList (Rho_1, 2, 3), PnList, dRhoList_dAlpha
	RhoList = alphaList;
	for ii=RN-1:-1:1
		RhoList(:,ii) = prod(1-alphaList(:,ii+1:RN),2) .* alphaList(:,ii);
	end
	PnList = RhoList ./ rhoRN;
	
	% dRhoList_dAlpha = zeros(RN,RN,numElements);
	%%						dRho1_dAlpha1		dRho2_dAlpha1		dRho3_dAlpha1
	%%	dRhoList_dAlpha = 	dRho1_dAlpha2		dRho2_dAlpha2		dRho3_dAlpha2
	%%						dRho1_dAlpha3		dRho2_dAlpha3		dRho3_dAlpha3
	dP_dAlpha = zeros(RN,RN,numElements);
	%%					dP1_dAlpha1		dP2_dAlpha1		dP3_dAlpha1
	%%	dP_dAlpha = 	dP1_dAlpha2		dP2_dAlpha2		dP3_dAlpha2
	%%					dP1_dAlpha3		dP2_dAlpha3		dP3_dAlpha3
	tmpCons = 1 ./ rhoRN ./ rhoRN;
	
	dRho1_dAlpha1 = (1-alphaList(:,3)).*(1-alphaList(:,2));
	dP1_dAlpha1 = tmpCons .* (rhoRN .* dRho1_dAlpha1 - RhoList(:,1) .* dRho_dAlpha(:,1));
	dP_dAlpha(1,1,:) = reshape(dP1_dAlpha1,1,1,numElements);
	%%dRho2_dAlpha1 = [0];
	dP2_dAlpha1 = tmpCons .* (-RhoList(:,2) .* dRho_dAlpha(:,1));
	dP_dAlpha(1,2,:) = reshape(dP2_dAlpha1,1,1,numElements);
	%%dRho3_dAlpha1 = [0];
	dP3_dAlpha1 = tmpCons .* (-RhoList(:,3) .* dRho_dAlpha(:,1));
	dP_dAlpha(1,3,:) = reshape(dP3_dAlpha1,1,1,numElements);
	dP_dAlpha1 = [dP1_dAlpha1 dP2_dAlpha1 dP3_dAlpha1];
	
	dRho1_dAlpha2 = (alphaList(:,3)-1).*alphaList(:,1);
	dP1_dAlpha2 = tmpCons .* (rhoRN .* dRho1_dAlpha2 - RhoList(:,1) .* dRho_dAlpha(:,2));
	dP_dAlpha(2,1,:) = reshape(dP1_dAlpha2,1,1,numElements);
	dRho2_dAlpha2 = 1-alphaList(:,3);
	dP2_dAlpha2 = tmpCons .* (rhoRN .* dRho2_dAlpha2 - RhoList(:,2) .* dRho_dAlpha(:,2));
	dP_dAlpha(2,2,:) = reshape(dP2_dAlpha2,1,1,numElements);
	%%dRho3_dAlpha2 = [0];
	dP3_dAlpha2 = tmpCons .* (-RhoList(:,3) .* dRho_dAlpha(:,2));
	dP_dAlpha(2,3,:) = reshape(dP3_dAlpha2,1,1,numElements);
	dP_dAlpha2 = [dP1_dAlpha2 dP2_dAlpha2 dP3_dAlpha2];
	
	dRho1_dAlpha3 = (alphaList(:,2)-1).*alphaList(:,1);
	dP1_dAlpha3 = tmpCons .* (rhoRN .* dRho1_dAlpha3 - RhoList(:,1) .* dRho_dAlpha(:,3));
	dP_dAlpha(3,1,:) = reshape(dP1_dAlpha3,1,1,numElements);
	dRho2_dAlpha3 = -alphaList(:,2);
	dP2_dAlpha3 = tmpCons .* (rhoRN .* dRho2_dAlpha3 - RhoList(:,2) .* dRho_dAlpha(:,3));
	dP_dAlpha(3,2,:) = reshape(dP2_dAlpha3,1,1,numElements);
	%%dRho3_dAlpha3 = [1];
	dP3_dAlpha3 = tmpCons .* (rhoRN - RhoList(:,3) .* dRho_dAlpha(:,3));
	dP_dAlpha(3,3,:) = reshape(dP3_dAlpha3,1,1,numElements);
	dP_dAlpha3 = [dP1_dAlpha3 dP2_dAlpha3 dP3_dAlpha3];
	
	%%dm1_dAlpha1, dm2_dAlpha1, dm3_dAlpha1, dm4_dAlpha1; by page
	dmList_dAlpha = zeros(numElements,4,RN);
	l1 = cos(2*thetaList);
	l2 = sin(2*thetaList);
	l3 = cos(4*thetaList);
	l4 = sin(4*thetaList);
	dmList_dAlpha(:,1,1) = sum(dP_dAlpha1 .* l1, 2); %dm1_dAlpha1
	dmList_dAlpha(:,2,1) = sum(dP_dAlpha1 .* l2, 2); %dm2_dAlpha1
	dmList_dAlpha(:,3,1) = sum(dP_dAlpha1 .* l3, 2); %dm3_dAlpha1
	dmList_dAlpha(:,4,1) = sum(dP_dAlpha1 .* l4, 2); %dm4_dAlpha1

	dmList_dAlpha(:,1,2) = sum(dP_dAlpha2 .* l1, 2); %dm1_dAlpha2
	dmList_dAlpha(:,2,2) = sum(dP_dAlpha2 .* l2, 2); %dm2_dAlpha2
	dmList_dAlpha(:,3,2) = sum(dP_dAlpha2 .* l3, 2); %dm3_dAlpha2
	dmList_dAlpha(:,4,2) = sum(dP_dAlpha2 .* l4, 2); %dm4_dAlpha2

	dmList_dAlpha(:,1,3) = sum(dP_dAlpha3 .* l1, 2); %dm1_dAlpha3
	dmList_dAlpha(:,2,3) = sum(dP_dAlpha3 .* l2, 2); %dm2_dAlpha3
	dmList_dAlpha(:,3,3) = sum(dP_dAlpha3 .* l3, 2); %dm3_dAlpha3
	dmList_dAlpha(:,4,3) = sum(dP_dAlpha3 .* l4, 2); %dm4_dAlpha3
	
	%%2. dM
	%%2.1 dM_dAlpha
	h0 = (1+v0) / 4/(1-v0);
	dM_dAlpha = struct('arr', []); dM_dAlpha.arr = zeros(3,3,numElements); dM_dAlpha = repmat(dM_dAlpha, 1, RN);
	for ii=1:RN
		dM_dAlpha(ii).arr(1,1,:) = reshape(-h0 * dmList_dAlpha(:,3,ii),1,1,numElements);	
		dM_dAlpha(ii).arr(1,2,:) = reshape(-h0 * dmList_dAlpha(:,4,ii),1,1,numElements); 
		dM_dAlpha(ii).arr(1,3,:) = reshape(0.5 * dmList_dAlpha(:,1,ii),1,1,numElements);	
		dM_dAlpha(ii).arr(2,1,:) = dM_dAlpha(ii).arr(1,2,:);	
		dM_dAlpha(ii).arr(2,2,:) = reshape(h0 * dmList_dAlpha(:,3,ii),1,1,numElements);	
		dM_dAlpha(ii).arr(2,3,:) = reshape(0.5 * dmList_dAlpha(:,2,ii),1,1,numElements);
		dM_dAlpha(ii).arr(3,1,:) = dM_dAlpha(ii).arr(1,3,:);	
		dM_dAlpha(ii).arr(3,2,:) = dM_dAlpha(ii).arr(2,3,:);
	end
	%%2.2 dM_dTheta
	dl1 = -2*sin(2*thetaList);
	dl2 = 2*cos(2*thetaList);
	dl3 = -4*sin(4*thetaList);
	dl4 = 4*cos(4*thetaList);	
	dmList_dTheta = zeros(numElements,4,RN);
	dmList_dTheta(:,1,1) = PnList(:,1).*dl1(:,1); %dm1_dTheta1
	dmList_dTheta(:,2,1) = PnList(:,1).*dl2(:,1); %dm2_dTheta1
	dmList_dTheta(:,3,1) = PnList(:,1).*dl3(:,1); %dm3_dTheta1
	dmList_dTheta(:,4,1) = PnList(:,1).*dl4(:,1); %dm4_dTheta1

	dmList_dTheta(:,1,2) = PnList(:,2).*dl1(:,2); %dm1_dTheta2
	dmList_dTheta(:,2,2) = PnList(:,2).*dl2(:,2); %dm2_dTheta2
	dmList_dTheta(:,3,2) = PnList(:,2).*dl3(:,2); %dm3_dTheta2
	dmList_dTheta(:,4,2) = PnList(:,2).*dl4(:,2); %dm4_dTheta2

	dmList_dTheta(:,1,3) = PnList(:,3).*dl1(:,3); %dm1_dTheta3
	dmList_dTheta(:,2,3) = PnList(:,3).*dl2(:,3); %dm2_dTheta3
	dmList_dTheta(:,3,3) = PnList(:,3).*dl3(:,3); %dm3_dTheta3
	dmList_dTheta(:,4,3) = PnList(:,3).*dl4(:,3); %dm4_dTheta3
	
	dM_dTheta = struct('arr', []); dM_dTheta.arr = zeros(3,3,numElements);
	if 1==OPT_Theta3Only
		for ii=1:RN
			dM_dTheta.arr(1,1,:) = dM_dTheta.arr(1,1,:) + reshape(-h0 * dmList_dTheta(:,3,ii),1,1,numElements);	
			dM_dTheta.arr(1,2,:) = dM_dTheta.arr(1,2,:) + reshape(-h0 * dmList_dTheta(:,4,ii),1,1,numElements); 
			dM_dTheta.arr(1,3,:) = dM_dTheta.arr(1,3,:) + reshape(0.5 * dmList_dTheta(:,1,ii),1,1,numElements);	
			dM_dTheta.arr(2,1,:) = dM_dTheta.arr(1,2,:);	
			dM_dTheta.arr(2,2,:) = dM_dTheta.arr(2,2,:) + reshape(h0 * dmList_dTheta(:,3,ii),1,1,numElements);	
			dM_dTheta.arr(2,3,:) = dM_dTheta.arr(2,3,:) + reshape(0.5 * dmList_dTheta(:,2,ii),1,1,numElements);
			dM_dTheta.arr(3,1,:) = dM_dTheta.arr(1,3,:);	
			dM_dTheta.arr(3,2,:) = dM_dTheta.arr(2,3,:);			
		end	
	else	
		dM_dTheta = repmat(dM_dTheta, 1, RN);
		for ii=1:RN
			dM_dTheta(ii).arr(1,1,:) = reshape(-h0 * dmList_dTheta(:,3,ii),1,1,numElements);	
			dM_dTheta(ii).arr(1,2,:) = reshape(-h0 * dmList_dTheta(:,4,ii),1,1,numElements); 
			dM_dTheta(ii).arr(1,3,:) = reshape(0.5 * dmList_dTheta(:,1,ii),1,1,numElements);	
			dM_dTheta(ii).arr(2,1,:) = dM_dTheta(ii).arr(1,2,:);	
			dM_dTheta(ii).arr(2,2,:) = reshape(h0 * dmList_dTheta(:,3,ii),1,1,numElements);	
			dM_dTheta(ii).arr(2,3,:) = reshape(0.5 * dmList_dTheta(:,2,ii),1,1,numElements);
			dM_dTheta(ii).arr(3,1,:) = dM_dTheta(ii).arr(1,3,:);	
			dM_dTheta(ii).arr(3,2,:) = dM_dTheta(ii).arr(2,3,:);
		end		
	end
	
	%%3. dY
	%%3.1 dY_dAlpha
	dY_dAlpha = struct('arr', []); dY_dAlpha.arr = zeros(3,3,numElements); dY_dAlpha = repmat(dY_dAlpha, 1, RN);
	for ii=1:RN
		for jj=1:numElements
			dY_dAlpha(ii).arr(:,:,jj) = YY(:,:,jj) * dM_dAlpha(ii).arr(:,:,jj);
		end
		dY_dAlpha(ii).arr = A0*dY_dAlpha(ii).arr;
	end
	%%3.2 dY_dTheta
	dY_dTheta = struct('arr', []); dY_dTheta.arr = zeros(3,3,numElements);
	if 1==OPT_Theta3Only
		for jj=1:numElements
			dY_dTheta.arr(:,:,jj) = YY(:,:,jj) * dM_dTheta.arr(:,:,jj);
		end
		dY_dTheta.arr = A0*dY_dTheta.arr;
	else		
		dY_dTheta = repmat(dY_dTheta, 1, RN);
		for ii=1:RN
			for jj=1:numElements
				dY_dTheta(ii).arr(:,:,jj) = YY(:,:,jj) * dM_dTheta(ii).arr(:,:,jj);
			end
			dY_dTheta(ii).arr = A0*dY_dTheta(ii).arr;			
		end		
	end
	
	%%4. dS
	%%4.1 dS_dAlpha
	dS_dAlpha = struct('arr', []); dS_dAlpha.arr = zeros(3,3,numElements); dS_dAlpha = repmat(dS_dAlpha, 1, RN);
	for ii=1:RN
		dS_dAlpha(ii).arr = times(Y, reshape(dRho_dAlpha(:,ii),1,1,numElements)) - times(dY_dAlpha(ii).arr, reshape(1-rhoRN,1,1,numElements));
		%% Transfer 2nd-order Tensor to 4th-order Tensor
		dS_dAlpha(ii).arr = Tensor2ndToTensor4th(dS_dAlpha(ii).arr);
	end
	%%4.2 dS_dTheta
	dS_dTheta = struct('arr', []); dS_dTheta.arr = zeros(3,3,numElements); 
	if 1==OPT_Theta3Only
		dS_dTheta.arr = times(dY_dTheta.arr, reshape(rhoRN-1,1,1,numElements));
		%% Transfer 2nd-order Tensor to 4th-order Tensor
		dS_dTheta.arr = Tensor2ndToTensor4th(dS_dTheta.arr);		
	else
		dS_dTheta = repmat(dS_dTheta, 1, RN);
		for ii=1:RN
			dS_dTheta(ii).arr = times(dY_dTheta(ii).arr, reshape(rhoRN-1,1,1,numElements));
			%% Transfer 2nd-order Tensor to 4th-order Tensor
			dS_dTheta(ii).arr = Tensor2ndToTensor4th(dS_dTheta(ii).arr);
		end			
	end	
	
	%%5. dKe
	%%5.1 dKe_dAlpha
	tmp_dK_dAlpha = struct('arr', []); tmp_dK_dAlpha.arr = zeros(8,8,numElements); tmp_dK_dAlpha = repmat(tmp_dK_dAlpha, 1, RN);	
	for ii=1:RN
		for jj=1:numElements
			idS_dAlpha = zeros12;
			for kk=1:4
				idx = (kk-1)*3+1:kk*3;
				idS_dAlpha(idx,idx) = dS_dAlpha(ii).arr(:,:,jj);
			end
			iDe = B'*idS_dAlpha*B;%%Integral
			tmp_dK_dAlpha(ii).arr(:,:,jj) = iDe/4;
		end
	end
	%%5.2 dKe_dTheta
	tmp_dK_dTheta = struct('arr', []); tmp_dK_dTheta.arr = zeros(8,8,numElements);
	if 1==OPT_Theta3Only
		for jj=1:numElements
			idS_dTheta = zeros12;
			for kk=1:4
				idx = (kk-1)*3+1:kk*3;
				idS_dTheta(idx,idx) = dS_dTheta.arr(:,:,jj);
			end
			iDe = B'*idS_dTheta*B;			
			tmp_dK_dTheta.arr(:,:,jj) = iDe/4;
		end		
	else
		tmp_dK_dTheta = repmat(tmp_dK_dTheta, 1, RN);
		for ii=1:RN
			for jj=1:numElements
				idS_dTheta = zeros12;
				for kk=1:4
					idx = (kk-1)*3+1:kk*3;
					idS_dTheta(idx,idx) = dS_dTheta(ii).arr(:,:,jj);
				end
				iDe = B'*idS_dTheta*B;
				tmp_dK_dTheta(ii).arr(:,:,jj) = iDe/4;
			end			
		end		
	end
	
	%%6. Assemble
	%%6.1 Width
	df0dAlpha = zeros(size(alphaList));	
	dK_dAlpha = ones(64, numElements, RN);
	for jj=1:RN
		dK_dAlpha(:,:,jj) = reshape(tmp_dK_dAlpha(jj).arr, 64, numElements);
	end
	for jj=1:RN
		for kk=1:size(Ue2,3)
			df0dAlpha(:,jj) = df0dAlpha(:,jj) + ww(kk)*reshape(sum(dK_dAlpha(:,:,jj).*Ue2(:,:,kk),1),numElements,1);
		end
	end
	df0dAlpha = W*df0dAlpha/c0Total;
	%%6.2 Orientation
	if 1==OPT_Theta3Only
		df0dTheta = zeros(size(thetaList,1),1);
		dK_dTheta = reshape(tmp_dK_dTheta.arr, 64, numElements);
		for kk=1:size(Ue2,3)
			df0dTheta = df0dTheta + ww(kk)*reshape(sum(dK_dTheta.*Ue2(:,:,kk),1),numElements,1);
		end		
	else
		df0dTheta = zeros(size(thetaList));
		dK_dTheta = ones(64, numElements, RN);
		for jj=1:RN
			dK_dTheta(:,:,jj) = reshape(tmp_dK_dTheta(jj).arr, 64, numElements);
		end	
		for jj=1:RN
			for kk=1:size(Ue2,3)
				df0dTheta(:,jj) = df0dTheta(:,jj) + ww(kk)*reshape(sum(dK_dTheta(:,:,jj).*Ue2(:,:,kk),1),numElements,1);
			end
		end			
	end
	df0dTheta = W/c0Total*df0dTheta;	
	
	%%7. Regularization
	coefi = 6; %% 2: k*pi (min), k*pi/2 (max); 6: k*pi/3 (min), k*pi/6 (max)
	dp0dTheta = zeros(size(df0dTheta));
	if W<1
		if OPT_Theta3Only
			for ii=1:numElements
				ijTheta = thetaList(ii,3);
				ijAdjacentEles = eleAdjacencyMap(ii).adjacentEles;
				ijAdjacentThetas = thetaList(ijAdjacentEles,3);
				% dp0dTheta(ii) = 2*sum(sin(coefi*ijTheta - coefi*ijAdjacentThetas));
				dp0dTheta(ii) = sum(coefi*0.5*sin(coefi*ijAdjacentThetas - coefi*ijTheta));
			end		
		else
			for jj=1:RN
				for ii=1:numElements
					ijTheta = thetaList(ii,jj);
					ijAdjacentEles = eleAdjacencyMap(ii).adjacentEles;
					ijAdjacentThetas = thetaList(ijAdjacentEles,jj);
					% dp0dTheta(ii,jj) = 2*sum(sin(2*ijTheta - 2*ijAdjacentThetas));
					dp0dTheta(ii,jj) = sum(coefi*0.5*sin(coefi*ijAdjacentThetas - coefi*ijTheta));
				end			
			end		
		end
		dp0dTheta = (1-W)/p0Total*dp0dTheta;
	end
	df0dTheta = df0dTheta + dp0dTheta;

	df0dAlpha = -df0dAlpha;
	df0dTheta = -df0dTheta;	
end