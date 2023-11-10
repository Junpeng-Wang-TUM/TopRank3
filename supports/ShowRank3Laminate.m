function ShowRank3Laminate(rank3LaminateParas, meshInfo, lgthScale)
	figure;
	if 6~=size(rank3LaminateParas,2), warning('Only Works with Rank-3!'); return; end
	alphaList = rank3LaminateParas(:,1:3);
	thetaList = rank3LaminateParas(:,4:6);
	densityField = 1-prod(1-alphaList,2);
	xPatchs = meshInfo.nodeCoords(meshInfo.nodMapBack,1); xPatchs = xPatchs(meshInfo.eNodMat');
	yPatchs = meshInfo.nodeCoords(meshInfo.nodMapBack,2); yPatchs = yPatchs(meshInfo.eNodMat');
	cPatchs = densityField(:)';
	cPatchs = repmat(cPatchs, 4, 1);	
	colormap(flip(gray));
	hCells = patch(xPatchs, yPatchs, cPatchs); hold('on');
	ctrs = meshInfo.eleCentroidList;
	
	ang1 = thetaList(:,1);
	ang2 = thetaList(:,2);
	ang3 = thetaList(:,3);
	if 0==lgthScale
		lgthScale = 0.3;
		amp1 = 1; 
		amp2 = 1; 
		amp3 = 1; 	
    else
        amp1 = alphaList(:,1).*(1-alphaList(:,2)).*(1-alphaList(:,3));
        amp2 = alphaList(:,2).*(1-alphaList(:,3));
        amp3 = alphaList(:,3);
		amp1 = max(amp1,1.0e-6); amp2 = max(amp2,1.0e-6); amp3 = max(amp3,1.0e-6);
	end
	hDir11 = quiver(ctrs(:,1), ctrs(:,2), cos(ang1).*amp1, sin(ang1).*amp1, lgthScale); hold('on');
	%hDir12 = quiver(ctrs(:,1), ctrs(:,2), -cos(ang1).*amp1, -sin(ang1).*amp1, lgthScale); hold('on');
	hDir21 = quiver(ctrs(:,1), ctrs(:,2), cos(ang2).*amp2, sin(ang2).*amp2, lgthScale); hold('on');
	%hDir22 = quiver(ctrs(:,1), ctrs(:,2), -cos(ang2).*amp2, -sin(ang2).*amp2, lgthScale); hold('on');
	hDir31 = quiver(ctrs(:,1), ctrs(:,2), cos(ang3).*amp3, sin(ang3).*amp3, lgthScale); hold('on');
	%hDir32 = quiver(ctrs(:,1), ctrs(:,2), -cos(ang3).*amp3, -sin(ang3).*amp3, lgthScale); hold('on');		
	
	set([hDir11], 'LineWidth', 2, 'Color', [252 141 98]/255);
	set([hDir21], 'LineWidth', 2, 'Color', [102 194 165]/255);
	set([hDir31], 'LineWidth', 2, 'Color', [0 0 255]/255);
	set(hCells, 'FaceColor', 'flat', 'FaceAlpha', 1, 'EdgeColor', 'None');
	axis('equal'); axis('tight'); axis('off');
end
	