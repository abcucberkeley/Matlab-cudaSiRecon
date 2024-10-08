function makeOTF(psfPath,nphases,angle,background,xyres,zres)
    [pathstr, fsname, ext] = fileparts(psfPath);
    otfFol = [pathstr '/OTFs/'];
    if  ~exist(otfFol, 'dir')
        mkdir(otfFol);
        fileattrib(otfFol, '+w', 'g');
    end
    system(['makeotf ' psfPath ' ' otfFol fsname  '_otf' ext ...,
        ' -nphases ' sprintf('%d',nphases) ' -angle ' sprintf('%.5f',angle) ...,
        ' -background ' sprintf('%d',background) ...,
        ' -xyres ' sprintf('%.5f',xyres)  ' -zres ' sprintf('%.5f',zres) ]);
end