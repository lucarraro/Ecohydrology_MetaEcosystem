function plot_custom(A,y,path,kol,lty,draw_patch,hh)
if nargin==5; draw_patch=0; hh=0; end
if nargin==6; hh=0; end
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if strcmp(lty,'-.')
    semilogx(A(path), mean(y(path,:),2), 'Marker', '.', 'LineStyle', '-', 'Color',kol);
else
    semilogx(A(path), mean(y(path,:),2), lty, 'Color',kol);
end
if hh; hold on; end
if draw_patch
    patch([A(path); flip(A(path))],[quantile(y(path,:),0.1,2); flip(quantile(y(path,:),0.9,2))], ...
        kol,'FaceAlpha',0.5, 'EdgeColor','none')
end


end