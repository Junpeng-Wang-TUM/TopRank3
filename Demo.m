%% ----TopRank3----
%% Author: Junpeng Wang
%% Copyright (c) All Rights Reserved.
%% Version: 2023-11-10
%% E-mail: junpeng.wang@tum.de

%%This repository was created for the paper "Design and Optimization of Functionally-graded Lattice Structures for Multiple Loading Conditions" 
%%	by Junpeng Wang, Rüdiger Westermann, Xifeng Gao and Jun Wu. 
%%which was submitted to the journal "xxx", and also available in arXiv (xxxx.xxxxx)

clear; clc;
addpath('./supports');

%%1. Problem Description
if 1 %%User Tailored Model
	%%1.1 Design Domain
	mesh_ = CreateCartesianMesh2D(ones(40, 80));
	
	%%1.2 Boundary Conditions
	ndC = mesh_.nodeCoords(mesh_.nodMapBack(mesh_.nodesOnBoundary),:);
	fixingCond_ = struct('arr', []);
	fixedNodes = mesh_.nodesOnBoundary(find(mesh_.boundingBox(1,1)==ndC(:,1)));
	fixingCond_.arr = [fixedNodes ones(length(fixedNodes), 2)];
	
	loadingCond_ = struct('arr', []); loadingCond_ = repmat(loadingCond_,2,1);
	force = [1 0; 0 -1];
	loadedCoords = [mesh_.boundingBox(2,1) mesh_.boundingBox(2,2)/2; mesh_.boundingBox(2,1)/2 0]; 
	[~,loadedNodes] = min(vecnorm(loadedCoords(1,:)-ndC,2,2));
	loadingCond_(1).arr = [mesh_.nodesOnBoundary(loadedNodes) force(1,:)];
	[~,loadedNodes] = min(vecnorm(loadedCoords(2,:)-ndC,2,2));
	loadingCond_(2).arr = [mesh_.nodesOnBoundary(loadedNodes) force(2,:)];
else %%Import the external model involved in the associated paper
	fileName = ('./data/femur_R150.txt');
	[mesh_, fixingCond_, loadingCond_] = ImportExternalFEAmodel(fileName);
end

%%1.3 Material Properties
matProp_.modulus = 1.0;
matProp_.PossionRatio = 0.3;
matProp_.minModulus = 1.0e-6 * matProp_.modulus;

%%1.4 Optimization Settings
optiSet_.Vf = 0.5; %% Target Volume Fraction
optiSet_.perLoadContri = ones(1, numel(loadingCond_))/numel(loadingCond_);
optiSet_.OPT_Theta3Only = 1; %%1: (theta3, theta2 = theta3+pi/3, theta1 = theta3+2*pi/3), 0: (All Independent Thetas)
optiSet_.OPT_OriStartingGuess = 'SGp'; %% 'SGu' (Uniform), 'SGp' (Principal Stress-guided)
optiSet_.W = 0.5; %%Weight Factor for Compliance Item, 1-W for Orientation Regularization Item
optiSet_.alphaBounds = [1.0e-1 0.5]; %% Lower and Upper Bounds for Layer Widths (Alpha)
optiSet_.nloop = 300; %% Number of Iterations to Conduct
optiSet_.passiveEles = SetBoundaryElementsBePassive(mesh_, 0); %%Set Domain Boundary be Solid, the Second Argument Decides the Width of the Solid Boundary
ShowProblemDescription(mesh_, fixingCond_, loadingCond_);

%%2. Run Optimization
tStart = tic;
outPath = './out/'; if ~exist(outPath, 'dir'), mkdir(outPath); end
[alphaList_, thetaList_, optiLog_, cPerLoad_] = TopRank3(mesh_, fixingCond_, loadingCond_, matProp_, optiSet_, outPath);
disp(['Performing Optimization Costs ' sprintf('%10.3g',toc(tStart)) 's']);

%%3. Visualization
ShowRank3Laminate([alphaList_ thetaList_], mesh_, 0.5);
ShowRank3byStreamlines(alphaList_, thetaList_, mesh_, 30);