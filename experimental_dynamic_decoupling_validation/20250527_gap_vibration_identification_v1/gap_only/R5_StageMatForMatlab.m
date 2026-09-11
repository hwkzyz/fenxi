function stagedFile=R5_StageMatForMatlab(sourceFile,tag)
assert(isfile(sourceFile),'R5:MissingMatFile','MAT file not found: %s',sourceFile);
if nargin<2||isempty(tag),tag='r5';end
[~,~,ext]=fileparts(sourceFile);
% Use a unique ASCII staging name so concurrent MATLAB workers cannot
% overwrite each other's HDF5 file.
stagedFile=[tempname(tempdir) ext];
[ok,msg]=copyfile(sourceFile,stagedFile,'f'); assert(ok,'R5:StageFailed','Cannot stage MAT file %s: %s',sourceFile,msg);
end
