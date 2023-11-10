function ShowProblemDescription(meshInfo, fixingCond, loadingCond)
	figure;
	%%Design Domain
	xPatchs = meshInfo.nodeCoords(meshInfo.nodMapBack,1); xPatchs = xPatchs(meshInfo.eNodMat');
	yPatchs = meshInfo.nodeCoords(meshInfo.nodMapBack,2); yPatchs = yPatchs(meshInfo.eNodMat');
	cPatchs = zeros(size(yPatchs));
	hdDesignDomain = patch(xPatchs, yPatchs, cPatchs); hold on;
	xlabel('X'); ylabel('Y');
	set(hdDesignDomain, 'FaceColor', [0.5 0.5 0.5], 'FaceAlpha', 1, 'EdgeColor', 'None');
	
	%%Fixations
	for ii=1:numel(fixingCond)
		tarNodeCoord = meshInfo.nodeCoords(meshInfo.nodMapBack(fixingCond(ii).arr(:,1)),:);
		hold('on'); hdFixation = plot(tarNodeCoord(:,1), tarNodeCoord(:,2), 'x', 'Color', [0.0 0.0 0.0], 'LineWidth', 2, 'MarkerSize', 10);			
	end
	
	%%Loads
	lB = 0.2;
	uB = 1.0;
	numLoadSteps = numel(loadingCond);
	amps = struct('arr', [], 'scalingFac', 1); amps = repmat(amps, numLoadSteps, 1);
	for ii=1:numLoadSteps
		amps(ii).arr = vecnorm(loadingCond(ii).arr(:,2:end),2,2)';
	end
	maxAmp = max([amps.arr]);
	minAmp = min([amps.arr]);
	if abs(minAmp-maxAmp)/(minAmp+maxAmp)>0.1
		if minAmp/maxAmp>lB/uB, lB = minAmp/maxAmp; end
		for ii=1:numLoadSteps
			amps(ii).scalingFac = lB + (uB-lB)*(amps(ii).arr'-minAmp)/(maxAmp-minAmp);
		end
	end	 	
	for ii=1:numLoadSteps		
		loadingDirVec = loadingCond(ii).arr(:,2:end)./amps(ii).arr' .* amps(ii).scalingFac;
		coordLoadedNodes = meshInfo.nodeCoords(meshInfo.nodMapBack(loadingCond(ii).arr(:,1)),:);
		amplitudesF = mean(meshInfo.boundingBox(2,:)-meshInfo.boundingBox(1,:))/5 * loadingDirVec;
		hold('on'); hdLoads = quiver(coordLoadedNodes(:,1), coordLoadedNodes(:,2), amplitudesF(:,1), ...
			amplitudesF(:,2), 0, 'LineWidth', 2, 'MaxHeadSize', 1, 'MaxHeadSize', 1); 			
	end
	
	axis('equal'); axis('tight'); axis('on');
	set(gca, 'FontName', 'Times New Roman', 'FontSize', 20);	
end
