function [mesh_, fixingCond_, loadingCond_] = ImportExternalFEAmodel(fileName)
	fid = fopen(fileName, 'r');
	tmp = fscanf(fid, '%s %s', 2);
	domainType = fscanf(fid, '%s', 1);
	tmp = fscanf(fid, '%s', 1);
	tmp = fscanf(fid, '%d %d', 3);
	nelx = tmp(1); nely = tmp(2);
	tmp = fscanf(fid, '%s %s', 2);
	numEles = fscanf(fid, '%d', 1);
	eleList = fscanf(fid, '%d', [1 numEles])';
 	allEles = zeros(nelx*nely,1); allEles(eleList) = 1;
	mesh_ = CreateCartesianMesh2D(reshape(allEles, nely, nelx));           
	tmp = fscanf(fid, '%s %s', 2);
	numFixedNodes = fscanf(fid, '%d', 1);
	fixingCond_ = struct('arr', []);
	if numFixedNodes>0
		fixingCond_(1).arr = fscanf(fid, '%d %d %d', [3 numFixedNodes])';
		fixingCond_(1).arr(:,1) = mesh_.nodMapForward(fixingCond_(1).arr(:,1));
	end
	tmp = fscanf(fid, '%s %s', 2);
	numLoadedNodes = fscanf(fid, '%d', 1);	
	if numLoadedNodes>0
		loadingCond_(1).arr = fscanf(fid, '%d %e %e', [3 numLoadedNodes])';	
		loadingCond_(1).arr(:,1) = double(mesh_.nodMapForward(loadingCond_(1).arr(:,1)));
	end
	tmp1 = fscanf(fid, '%s %s %s', 3);
	if ~isempty(tmp1)
		if strcmp(tmp1, 'additionalboundaryconditions:')
			numAdditionalBCs = fscanf(fid, '%d', 1);
			for ii=1:numAdditionalBCs
				%%Additional Loads
				tmp = fscanf(fid, '%s %s', 2);
				tmp1 = fscanf(fid, '%d %d', 2);
				numLoadedNodes = tmp1(1);
				iLoadIdx = tmp1(2);
				tmp1 = fscanf(fid, '%d %e %e', [3 numLoadedNodes])';
				tmp1(:,1) = double(mesh_.nodMapForward(tmp1(:,1)));
				loadingCond_(iLoadIdx).arr = tmp1;
				%%Additional Fixations
				tmp = fscanf(fid, '%s %s', 2);
				tmp1 = fscanf(fid, '%d %d', 2);
				numFixedNodes = tmp1(1);
				iFixationIdx = tmp1(2);
				tmp1 = fscanf(fid, '%d %d %d', [3 numFixedNodes])';
				tmp1(:,1) = mesh_.nodMapForward(tmp1(:,1));
				fixingCond_(iFixationIdx).arr = tmp1;
			end
		end
	end
	fclose(fid);	
end