function show_movie_2d_xy(x,y,V,T, ts)
%SHOW_MOVIE_2D_XY(x,y,V,T,ts)
%
%This function dusplays a new fiugure with a colorplot of the V field, and
%a sequence of wavefront propagation based on the arrival times given in T.
%
%input:
%    x:     x coordinates (1d array)
%    y:     y coordinates (1d array)
%    V:     velocity field, size (Nx,Ny,1)
%    T:     arrival times, size (Nx,Ny,1)
%    ts:    time step between successive frames


tt = linspace(0,max(max(T(isfinite(T)))),60);

figure;

pcolor(x,y,squeeze(V(:,:,1))');
caxis([floor(min(min(V))) ceil(max(max(V)))]);
axis equal;

xlim([x(1) x(end)]);
ylim([y(1) y(end)]);
xlabel('{\itx} (mm)');
ylabel('{\ity} (mm)');

hold on;

[~,h2] = contour(x,y,squeeze(T(:,:,1))',tt(1),'k', 'linewidth',0.5);
[~,h1] = contour(x,y,squeeze(T(:,:,1))',[tt(1) tt(1)],'k', 'linewidth',2);

for k=2:length(tt)
    delete(h1);
    delete(h2);
    [~,h2] = contour(x,y,squeeze(T(:,:,1))',tt(1:3:k-1),'k', 'linewidth',0.5);
    [~,h1] = contour(x,y,squeeze(T(:,:,1))',[tt(k) tt(k)],'k', 'linewidth',2);
    
    drawnow;
    
    pause(ts);
    
end
disp('done');