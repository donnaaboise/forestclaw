yrbcolormap;
setviews;

parms = read_vars();

showpatchborders(1:10);
setpatchborderprops('linewidth',1)
caxis([0,1])

hideslices();
showslices('z',1);

axis([-1,1,-1,1,-1,1]*parms.scale(1));


daspect([1 1 1]);
view(vtop);
%view(3);

shg
