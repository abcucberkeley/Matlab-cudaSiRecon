function [recIm] = cudaSireconChunk(inFol, inN, otfF,configF,outFullPath,varargin)

%export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CONDA_PREFIX/lib
%deskew=-32.45

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('inN', @ischar);
ip.addRequired('otfF', @ischar);
ip.addRequired('configF', @ischar);
ip.addRequired('outFullPath', @ischar);
ip.addParameter('nphases', 5, @isnumeric);
ip.addParameter('ndirs', 1, @isnumeric);
ip.addParameter('chunk', true, @islogical);
ip.addParameter('smooth', false, @islogical);
ip.addParameter('overlap', 128, @isnumeric);
ip.addParameter('chunkSize', [320,320,0], @isnumeric);
ip.addParameter('bak', 99, @isnumeric);
ip.addParameter('padpx', 0, @isnumeric);
ip.addParameter('occupancyRatio', 0.0, @isnumeric);
ip.parse(inN, otfF, configF, outFullPath, varargin{:});

pr = ip.Results;

nphases = pr.nphases;
ndirs = pr.ndirs;
chunk = pr.chunk;
smooth = pr.smooth;
ol = pr.overlap;
s = pr.chunkSize;
bak = pr.bak;
padpx = pr.padpx;
occupancyRatio = pr.occupancyRatio;
zf = 2;

im = readtiff_matlab([inFol '/' inN '.tif']);
[sy, sx,sz] = size(im);

% mask and interpolate
mask = im(:,:,1:nphases:end);
mask(mask<0) = 0;
mask = logical(mask);
imask = imresize3(mask,[ceil(sy*zf),ceil(sx*zf), sz/(nphases*ndirs)], 'nearest');

clear mask

im = single(im);
im = im - bak;
im(im<0) = 0;
otf = readtiff_matlab(otfF);

cd /clusterfs/nvme/matthewmueller/cudasireconMex/src/cudaSirecon/;

if(s(1) > sx)
    s(1) = sx;
end
if(s(2) > sy)
    s(2) = sy;
end
if(s(1) == 0)
    s(1) = sx;
end
if(s(2) == 0)
    s(2) = sy;
end

nyc = floor(sy/(s(2)-ol));
nxc = floor(sx/(s(1)-ol));

ymin = zeros(nyc,1);
ymax = zeros(nyc,1);

ymin_out = zeros(nyc,1);
ymax_out = zeros(nyc,1);


xmin = zeros(nxc,1);
xmax = zeros(nxc,1);

xmin_out = zeros(nxc,1);
xmax_out = zeros(nxc,1);

a = xmin;
edgey = xmin;
edgex = edgey;

nn = 0;
for h = 1:nxc
    for j = 1:nyc
        nn = nn + 1;
        xmin(nn) = (1+((h-1)*(s(1)-ol)));
        ymin(nn) = (1+((j-1)*(s(2)-ol)));
        if j==1
            ymin(nn) = 1;
            ymin_out(nn) = 1;
            edgey(nn) = 1;
        end

        if h==1
            xmin(nn) = 1;
            xmin_out(nn) = 1;
            edgex(nn) = 1;
        end

        if j < nyc
            %             ymax(nn) = (s(2)*j-(ol*(j)));
            %             ymax_out(nn) = (s(2)*j-(ol*(j)))*zf;
            ymax(nn) = (s(2)*j-(ol*(j-1)));
            ymax_out(nn) = (s(2)*j-(ol*(j-1)))*zf;
        else
            ymax(nn) = sy;
            ymax_out(nn) = sy*zf;
            edgey(nn) = 1;
        end

        if h < nxc
            %             xmax(nn) = (s(1)*h-(ol*(h)));
            %             xmax_out(nn) = (s(1)*h-(ol*(h)))*zf;
            xmax(nn) = (s(1)*h-(ol*(h-1)));
            xmax_out(nn) = (s(1)*h-(ol*(h-1)))*zf;
        else
            xmax(nn) = sx;
            xmax_out(nn) = sx*zf;
            edgex(nn) = 1;
        end


    end
end
a = 1:nn;

recIm = zeros(ceil(sy*zf),ceil(sx*zf), sz/(nphases*ndirs),'single');

if nn > 1
    for rr = a(logical(edgey+edgex))
        [csy, csx, csz] = size(im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:));
        ck_im = zeros(csy+padpx*2, csx+padpx*2, csz,'single');
        ck_im(padpx+1:csy+padpx, padpx+1:csx+padpx, :) = im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:);
        OR = nnz(im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:))/numel(im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:))
        if OR>=occupancyRatio
            %         outImg= cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
            %             configF, im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:), otf);
            outImg= cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
                configF, ck_im, otf);
            outImg = outImg(2*padpx+1:csy*2+2*padpx, 2*padpx+1:csx*2+2*padpx, :);

            if ymin(rr) ==1 && xmin(rr) ==1 && ymax(rr) < sy && xmax(rr) < sx % top left edge
                recIm(ymin(rr)*2-1:ymax(rr)*2-ol/2, xmin(rr)*2-1:xmax(rr)*2-ol/2,:) = outImg(1:end-ol/2, 1:end-ol/2,:);
            elseif ymin(rr) ==1 && xmin(rr) >1 && ymax(rr) < sy && xmax(rr) == sx % top right edge
                recIm(ymin(rr)*2-1:ymax(rr)*2-ol/2, xmin(rr)*2-1+ol/2:xmax(rr)*2, :) = outImg(1:end-ol/2, ol/2+1:end, :);
            elseif ymin(rr) >1 && xmin(rr) ==1 && ymax(rr) == sy && xmax(rr) < sx % bottom left edge
                recIm(ymin(rr)*2-1+ol/2:ymax(rr)*2, xmin(rr)*2-1:xmax(rr)*2-ol/2, :) = outImg(ol/2+1:end, 1:end-ol/2, :);
            elseif ymin(rr) >1 && xmin(rr) >1 && ymax(rr) == sy && xmax(rr) == sx % bottom right edge
                recIm(ymin(rr)*2-1+ol/2:ymax(rr)*2, xmin(rr)*2-1+ol/2:xmax(rr)*2, :) = outImg(ol/2+1:end, ol/2+1:end, :);
            elseif ymin(rr) ==1 && xmin(rr) >1 && ymax(rr) < sy && xmax(rr) < sx % top middle
                recIm(ymin(rr)*2-1:ymax(rr)*2-ol/2, xmin(rr)*2-1+ol/2:xmax(rr)*2-ol/2, :) = outImg(1:end-ol/2, ol/2+1:end-ol/2, :);
            elseif ymin(rr) >1 && xmin(rr) ==1 && ymax(rr) < sy && xmax(rr) < sx % left middle
                recIm(ymin(rr)*2-1+ol/2:ymax(rr)*2-ol/2, xmin(rr)*2-1:xmax(rr)*2-ol/2, :) = outImg(ol/2+1:end-ol/2, 1:end-ol/2, :);
            elseif ymin(rr) >1 && xmin(rr) >1 && ymax(rr) < sy && xmax(rr) == sx % right middle
                recIm(ymin(rr)*2-1+ol/2:ymax(rr)*2-ol/2, xmin(rr)*2-1+ol/2:xmax(rr)*2, :) = outImg(ol/2+1:end-ol/2, ol/2+1:end, :);
            elseif ymin(rr) >1 && xmin(rr) >1 && ymax(rr) ==sy && xmax(rr) < sx % bottom middle
                recIm(ymin(rr)*2-1+ol/2:ymax(rr)*2, xmin(rr)*2-1+ol/2:xmax(rr)*2-ol/2, :) = outImg(ol/2+1:end, ol/2+1:end-ol/2, :);
            elseif ymin(rr) ==1 && xmin(rr) ==1 && ymax(rr) < sy && xmax(rr) == sx % top center edge
                recIm(ymin(rr)*2-1:ymax(rr)*2-ol/2, xmin(rr)*2-1:xmax(rr)*2,:) = outImg(1:end-ol/2, 1:end,:);
            elseif ymin(rr) >1 && xmin(rr) ==1 && ymax(rr) < sy && xmax(rr) == sx % middle center edge
                recIm(ymin(rr)*2-1+ol/2:ymax(rr)*2-ol/2, xmin(rr)*2-1:xmax(rr)*2,:) = outImg(ol/2+1:end-ol/2, 1:end,:);
            elseif ymin(rr) >1 && xmin(rr) ==1 && ymax(rr) == sy && xmax(rr) == sx % bottom center edge
                recIm(ymin(rr)*2-1+ol/2:ymax(rr)*2, xmin(rr)*2-1:xmax(rr)*2,:) = outImg(1+ol/2:end, 1:end,:);
            end
        end
    end

    for rr = a(~logical(edgey+edgex))
        [csy, csx, csz] = size(im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:));
        ck_im = zeros(csy+padpx*2, csx+padpx*2, csz,'single');
        ck_im(padpx+1:csy+padpx, padpx+1:csx+padpx, :) = im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:);
        OR = nnz(im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:))/numel(im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:))
        if OR>=occupancyRatio
            %         outImg= cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
            %             configF, im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:), otf);

            outImg= cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
                configF, im(ymin(rr):ymax(rr),xmin(rr):xmax(rr),:), otf);
            outImg = outImg(2*padpx+1:csy*2+2*padpx, 2*padpx+1:csx*2+2*padpx, :);

            recIm(ymin(rr)*2-1+ol/2:ymax(rr)*2-ol/2, xmin(rr)*2-1+ol/2:xmax(rr)*2-ol/2, :) = outImg(ol/2+1:end-ol/2, ol/2+1:end-ol/2, :);
        end
    end
else
    recIm= cudaSireconDriverMex('cudasirecon',inFol,inN,otfF,'-c',...,
        configF, im, otf);
end

% mask to zero
recIm(~imask) = 0;
end