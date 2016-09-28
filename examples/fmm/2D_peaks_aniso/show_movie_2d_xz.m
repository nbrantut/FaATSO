function show_movie_2d_xz(x,z,V,E,T, ts)
%SHOW_MOVIE_2D_XZ(x,z,V,E,T,ts)
%
%This function displays a new fiugure with a colorplot of the V field, and
%a sequence of wavefront propagation based on the arrival times given in T.
%
%input:
%    x:     x coordinates (1d array)
%    z:     z coordinates (1d array)
%    V:     velocity field, size (Nx,1,Nz)
%    E:     anisotropy, size (Nx,1,Nz)
%    T:     arrival times, size (Nx,1,Nz)
%    ts:    time step between successive frames


tt = linspace(0,max(max(T(isfinite(T)))),60);

figure;

pcolor(x,z,squeeze(V(:,1,:))');
caxis([floor(min(min(V))) ceil(max(max(V)))]);
axis equal;

xlim([x(1) x(end)]);
ylim([z(1) z(end)]);
xlabel('{\itx} (mm)');
ylabel('{\itz} (mm)');

hold on;

plotaniso(x,z,E,10)

[~,h2] = contour(x,z,squeeze(T(:,1,:))',tt(1),'k', 'linewidth',0.5);
[~,h1] = contour(x,z,squeeze(T(:,1,:))',[tt(1) tt(1)],'k', 'linewidth',2);

for k=2:length(tt)
    delete(h1);
    delete(h2);
    [~,h2] = contour(x,z,squeeze(T(:,1,:))',tt(1:3:k-1),'k', 'linewidth',0.5);
    [~,h1] = contour(x,z,squeeze(T(:,1,:))',[tt(k) tt(k)],'k', 'linewidth',2);
    
    drawnow;
    
    pause(ts);
    
end
disp('done');