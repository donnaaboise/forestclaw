yrbcolormap;
setviews;

parms = read_vars();

fprintf('qmin = %24.8e\n', qmin)
fprintf('qmax = %24.8e\n', qmax)

showpatchborders(1:10);
setpatchborderprops('linewidth',1)
clim([0,1])

hideslices();
showslices('z',1);

axis([-1,1,-1,1,-1,1]*parms.scale(1));
% axis image;

plot_sphere = false;
if plot_sphere
    hold on;
    M = 32;
    dt = 2*pi/M;
    dp = pi/M;
    te = linspace(0,2*pi,M+1);
    pe = linspace(-pi,pi,M+1);
    
    [tem,pem] = meshgrid(te,pe);
    zem = 0*tem;
    R = parms.maxelev*zem + parms.scale(1);
    xem = R.*cos(pem).*cos(tem);
    yem = R.*cos(pem).*sin(tem);
    zem = R.*sin(pem);

    
    tc = te(1:end-1) + dt/2;
    pc = pe(1:end-1) + dp/2;
    [tcm,pcm] = meshgrid(tc,pc);
    zcm = 0*tcm;
    R = parms.maxelev*zcm + parms.scale(1);
    xcm = R.*cos(pcm).*cos(tcm);
    ycm = R.*cos(pcm).*sin(tcm);
    zcm = R.*sin(pcm);
        
    h = patch(surf2patch(xem,yem,zem,0*xcm));
    
    set(h,'edgecolor','k');
    set(h,'facecolor',[1,1,1]);
    hold off

    view(3)
    daspect([1 1 1]);   
    set(gca,'clipping','off');
else
    axis image
end




shg
