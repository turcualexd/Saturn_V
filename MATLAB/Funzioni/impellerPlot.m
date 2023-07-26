function  impellerPlot(r1, r2, rsh, beta2, rot) 

whitebg('k') 
set(gcf, 'color', 'k') 
set(gcf, 'InvertHardCopy', 'off'); 
ang([0,0], r2, [0, 2*pi], 'w-'); 

hold on ang([0,0], r1, [0, 2*pi], 'w-');
text(-r2, r2+r2/10, datestr(clock)); 
i=0; 
for teta = 0:0.01:2*pi 
    i=i+1; x2(i) = r2*cos(teta);
    y2(i) = r2*sin(teta); 
end 
j=0; 
i=0; 
ii=0; 

for ra=20:1:120 
    i=i+1; 
    delta(i) = inter(ra, x2, y2, beta2); 
    r(i) = ra; 
end 
c = delta(delta~=0); 
pos=find(delta==min(c)); 
rFin=r(pos); 
j=0; 

for teta = 0:0.01:pi j=j+1; xa(j) = -rFin+rFin*cos(teta); 
    ya(j) = rFin*sin(teta); 
end 

[a, b] = intersections(x2, y2, xa, ya); 

for j=1:(pi/0.01) 
    if xa(j)>a 
        xf(j) = xa(j); 
        yf(j) = ya(j); 
    end 
end 

h = plot(xf,yf, 'b'); 
xlabel('[mm]'); 
ylabel('[mm]');

rotate(h, [0,0,1], rot); %chiamare la funzione con ciclo per far ruotare 
% 

for teta = 0:0.01:2*pi i=i+1;
    xsh(i) = rsh*cos(teta);
    ysh(i) = rsh*sin(teta); 
end 
area(xsh,ysh, 'ShowBaseLine', 'off', 'FaceColor', 'w'); 
daspect([1,1,1]); 

end


