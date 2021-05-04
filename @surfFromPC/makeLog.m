function makeLog(obj,filename)
%makeLog prints a log of settings to a file

fileID = fopen(filename,'w');

fprintf(fileID,"Options log for rbfpum surface\r\n\r\n");

if ~isempty(fieldnames(obj.denoiseOpts))
	fprintf(fileID,"denoise nbrs = %d\r\n", obj.denoiseOpts.NumNeighbors);
	fprintf(fileID,"denoise thresh = %g\r\n", obj.denoiseOpts.Threshold);
end

if ~isempty(fieldnames(obj.downsampleOpts))
	fprintf(fileID,"downsample method = %s\r\n", obj.downsampleOpts.method);
	fprintf(fileID,"downsample param = %g\r\n", obj.downsampleOpts.param);
	fprintf(fileID,"\r\n");
end

fprintf(fileID,"preprocessed point count = %d\r\n", obj.pc.Count);
fprintf(fileID,"\r\n");

fprintf(fileID,"normal est nbrs = %d\r\n", obj.normalEstOpts.pcaNbrs);
fprintf(fileID,"normal est coarse grid = %g\r\n", obj.normalEstOpts.coarseGrid);
fprintf(fileID,"normal est graph nbrs = %d\r\n", obj.normalEstOpts.graphNbrs);
fprintf(fileID,"\r\n");

fprintf(fileID,"h = %d\r\n", obj.offSurfPtOpts.h);
fprintf(fileID,"L = %g\r\n", obj.offSurfPtOpts.L);
fprintf(fileID,"\r\n");

fprintf(fileID,"Nmin = %d\r\n", obj.octreeOpts.Nmin);
fprintf(fileID,"Nmax = %d\r\n", obj.octreeOpts.Nmax);
fprintf(fileID,"expand = %g\r\n", obj.octreeOpts.expand);
fprintf(fileID,"\r\n");

fprintf(fileID,"rbf = %s\r\n", func2str(obj.fitOpts.rbf));
fprintf(fileID,"polydegree = %d\r\n", obj.fitOpts.polydegree);
fprintf(fileID,"regFn = %s\r\n", func2str(obj.fitOpts.regFn));

if length(obj.fitOpts.rho) == 1
	fprintf(fileID,"GCV = false\r\n");
	fprintf(fileID,"rho = %g\r\n", obj.fitOpts.rho);
else
	fprintf(fileID,"GCV = true\r\n");
	fprintf(fileID,"mean rho = %g\r\n", mean(obj.fitOpts.rho));
end

fclose(fileID);

end

