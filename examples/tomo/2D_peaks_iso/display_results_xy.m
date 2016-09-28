%display_results_xy
%
%This routine is used to extract the inversion results from the output
%files of the tomo code, and plot them in a synoptic way to visually
%determine how the inversion method performed.

%import posterior V structure
[Nx,Ny,Nz,h,vout] = import_data([folder  'v0_out.bin']);

%import posterior event locations
events_out = dlmread([folder  'events_out.txt']);

%plot results
figure;
subplot 221
pcolor(h*(0:Nx-1), h*(0:Ny-1),exp(vout)');
axis equal;
xlim([0 5.4]);
ylim([0 5.4]);
caxis([2.9 3.2])

hold on
plot(sensors(:,1), sensors(:,2),'v','markersize',8);
plot(events(:,1), events(:,2),'k+');
plot(events_out(:,1), events_out(:,2),'ro');
title('Inversion (mean model)');

subplot 222
pcolor(h*(0:Nx-1), h*(0:Ny-1),Vh');
axis equal;
xlim([0 5.4]);
ylim([0 5.4]);
caxis([2.9 3.2]);
hold on
plot(sensors(:,1), sensors(:,2),'v','markersize',8);
plot(events(:,1), events(:,2),'k+');
plot(events_out(:,1), events_out(:,2),'ro');
title('True model');

subplot(2,2,[3 4])
r = dlmread([folder 'residuals.txt']);
semilogy(r,'o-')
xlabel('iteration');
ylabel('residual');