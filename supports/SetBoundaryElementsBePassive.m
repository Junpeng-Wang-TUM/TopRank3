function passiveEles = SetBoundaryElementsBePassive(meshInfo, numLayers)
	%%numLayers -> 0: none passive elements
	passiveEles = [];
	if numLayers>0
		index = 1;
		while index<=numLayers
			if 1==index
				passiveEles = double(meshInfo.elementsOnBoundary);
			else
				passiveEles = RelateAdjacentElements(meshInfo, passiveEles);
			end
			index = index + 1;
		end		
	end
end

function oEleList = RelateAdjacentElements(meshInfo, iEleList)
	iEleListMapBack = meshInfo.eleMapBack(iEleList);
	%%	1	4	7
	%%	2	5*	8
	%%	3	6	9
	resX = meshInfo.resX;
	resY = meshInfo.resY;			
	[eleX, eleY] = NodalizeDesignDomain([resX-1 resY-1], [1 1; resX resY]);
	eleX = eleX(iEleListMapBack);
	eleY = eleY(iEleListMapBack);
	tmpX = [eleX-1 eleX-1 eleX-1  eleX eleX eleX  eleX+1 eleX+1 eleX+1]; tmpX = tmpX(:);
	tmpY = [eleY+1 eleY eleY-1  eleY+1 eleY eleY-1  eleY+1 eleY eleY-1]; tmpY = tmpY(:);
	xNegative = find(tmpX<1); xPositive = find(tmpX>resX);
	yNegative = find(tmpY<1); yPositive = find(tmpY>resY);
	allInvalidEles = unique([xNegative; xPositive; yNegative; yPositive]);
	tmpX(allInvalidEles) = []; tmpY(allInvalidEles) = [];
	oEleListMapBack = resY*(tmpX-1) + resY-tmpY + 1;
	oEleList = meshInfo.eleMapForward(oEleListMapBack);
	oEleList(oEleList<1) = []; oEleList = unique(oEleList);
end