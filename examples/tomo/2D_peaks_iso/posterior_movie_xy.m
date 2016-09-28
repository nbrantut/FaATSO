function posterior_movie_xy(giffile, speed, Vrange, nsamples, mpostfile, CMpostfile, Nx, Ny, Nz, h)
%POSTERIOR_MOVIE_XY(giffile, speed, Vrange, nsamples, mpostfile, CMpostfile, Nx, Ny, Nz, h)
%
%This function is used to display a sort of "movie" of example velocity
%structures from the posterior pdf in the model space, given by a
%multivariate gaussian around the inversion's "best" (i.e., mean) result.
%The movie is then exported into an animated gif file.
%
%input:
%    giffile:      name of .gif file in which movide is saved
%    speed:        speed of movie
%    Vrange:       range of velocities to display using the colormap
%    nsamples:     number of samples to extract from pdf
%    mpostfile:    name of file with posterior mean m
%    CMpostfile:   name of file with posterior covariance
%    Nx:           grid size in x
%    Ny:           grid size in y
%    Nz:           grid size in z
%    h:            grid spacing


%coordinates
x = h*(0:Nx-1);
y = h*(0:Ny-1);

%read mpost vector:
mpost = dlmread(mpostfile);


%read CMpost matrix:
fid = fopen(CMpostfile,'r');
Cmpost = fread(fid,[length(mpost) length(mpost)],'float64');
fclose(fid);

% get samples from multivariate gaussian:
msamples = mvnrnd(mpost,Cmpost,nsamples);

%prepare array strcutures to convert "m" to more usable V0, E0, etc arrays.
ms = repmat(struct('V0',{}, 'E0', {}, 'N', {}, 'E', {}, 'D', {}, 't0', {}),nsamples, 1);

for k=1:nsamples
    ms(k) = m2fields(msamples(k,:), Nx, Ny, Nz);
end

%% display the movie:
figure;

%store colormap:
cmap = colormap;


%plot first model:
k=1;
pcolor(x,y,ms(k).V0');

%set some details of how it should look:
set(gca,'NextPlot','Replacechildren');
axis equal;
xlim([x(1) x(end)]);
ylim([y(1) y(end)]);
colorbar;
caxis(Vrange);

%store frame into indexed bitmap:
F = getframe;
[RGB,~] = frame2im(F);
X = rgb2ind(RGB, cmap);

%prepare arrays to make animated gif:
GIF = zeros([size(X) 1 nsamples]);
%store first frame:
GIF(:,:,1,1) = X;

%pause:
pause(1/speed);

%repeat that for all samples:
for k=2:nsamples
    pcolor(h*(0:Nx-1), h*(0:Ny-1),ms(k).V0');
    
    F = getframe;
    [RGB,~] = frame2im(F);
    X = rgb2ind(RGB, cmap);
    GIF(:,:,1,k) = X;
    
    pause(1/speed)
end

%export into animated gif:
imwrite(GIF,cmap,giffile,'gif', 'DelayTime', 1/speed, 'LoopCount',Inf)