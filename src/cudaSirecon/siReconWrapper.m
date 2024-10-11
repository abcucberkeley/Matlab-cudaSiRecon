function siReconWrapper(inFol,inN,outFol,otfF,configF,chunkSize,overlap, bak, ndir, np)
addpath(genpath('/clusterfs/nvme/matthewmueller/Matlab-cudaSiRecon/src/cudaSirecon/'));
cd '/clusterfs/nvme/matthewmueller/PetaKit5D';
setup;
gpuDevice;

outFullPath = [inFol sprintf('/%s/%s_recon.tif',outFol,inN)];
mkdir([inFol '/' outFol]);
if ~exist(outFullPath,'file')
    % GL: this cudaSiReconChunk.m currently (2023/11/20) messes with z and
    % is only suitable for LLS-SIM
    out = cudaSireconChunk(inFol,inN,otfF,configF,outFullPath,'chunkSize',chunkSize, ...
                           'overlap', overlap, 'bak', bak,'ndirs', ndir, 'nphases', np);
    writetiff(out,outFullPath);
else
    disp('LLS-SIM reconstructed file exists!!!')
end

end