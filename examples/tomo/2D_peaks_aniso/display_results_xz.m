%display_results_xz
%
%This routine is used to extract the inversion results from the output
%files of the tomo code, and plot them in a synoptic way to visually
%determine how the inversion method performed.

%import posterior V structure
[Nx,Ny,Nz,h,vout] = import_data([folder  'v0_out.bin']);
%import posterior E structure
[Nx,Ny,Nz,h,Eout] = import_data([folder  'E0_out.bin']);

%import posterior event locations
events_out = dlmread([folder  'events_out.txt']);

%plot results
figure;
subplot 321
pcolor(h*(0:Nx-1), h*(0:Nz-1),exp(squeeze(vout))');
axis equal;
xlim([0 5.4]);
ylim([0 5.4]);
caxis([2.9 3.2])

hold on
plot(sensors(:,1), sensors(:,3),'v','markersize',8);
plot(events(:,1), events(:,3),'k+');
plot(events_out(:,1), events_out(:,3),'ro');
title('Inversion (mean model)');

subplot 322
pcolor(h*(0:Nx-1), h*(0:Nz-1),squeeze(Vh)');
axis equal;
xlim([0 5.4]);
ylim([0 5.4]);
caxis([2.9 3.2]);
hold on
plot(sensors(:,1), sensors(:,3),'v','markersize',8);
plot(events(:,1), events(:,3),'k+');
plot(events_out(:,1), events_out(:,3),'ro');
title('True model');

subplot 323
pcolor(h*(0:Nx-1), h*(0:Nz-1),squeeze(Eout)');
axis equal;
xlim([0 5.4]);
ylim([0 5.4]);
caxis([-0.15 0.15])

hold on
plot(sensors(:,1), sensors(:,3),'v','markersize',8);
plot(events(:,1), events(:,3),'k+');
plot(events_out(:,1), events_out(:,3),'ro');
title('Inversion (mean model)');

subplot 324
pcolor(h*(0:Nx-1), h*(0:Nz-1),squeeze(E)');
axis equal;
xlim([0 5.4]);
ylim([0 5.4]);
caxis([-0.15 0.15])
hold on
plot(sensors(:,1), sensors(:,3),'v','markersize',8);
plot(events(:,1), events(:,3),'k+');
plot(events_out(:,1), events_out(:,3),'ro');
title('True model');

subplot(3,2,[5 6])
r = dlmread([folder 'residuals.txt']);
semilogy(r,'o-')
xlabel('iteration');
ylabel('residual');