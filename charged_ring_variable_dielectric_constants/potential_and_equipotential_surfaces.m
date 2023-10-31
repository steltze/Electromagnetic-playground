function charged_ring_potential_isodynamic_surfaces(L,h,er1,er2)

    %============
    %ask7_delta(0.5, 0.5,1,5);
    %ask7_delta(0.5, 0.5,5,1);
    %============
    purple         = [125  46 143]/255;
    
    Npoints=100;
    cont=[0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.2, 0.3, 0.4, 0.5];
    
    xx=-2:4/Npoints:2;
    zz1=-2:2/Npoints:0;
    [X1,Z1]=meshgrid(xx,zz1);
    
    Phi_array1=(2/(4*pi*(er1+er2))).*log(factors(X1,Z1,L,h));
    
    zz2=0:2/Npoints:2;
    zz2(abs(zz2-0.5)<=0.001)=[];
    [X2,Z2]=meshgrid(xx,zz2);
    
    Phi_array2=(1/(4*pi*er1)).*(log(factors(X2,Z2,L,h))+((er1-er2)/(er1+er2)).*log(factors(X2,Z2,L,-h)));
    
    M1=max(max(Phi_array1));
    M2=max(max(Phi_array2));
    M=max(M1,M2);
    
    % -------------------
    %    Plot Potential
    % -------------------
    figure(1);
    surface(X1,Z1,Phi_array1), shading interp
    hold on
    surface(X2,Z2,Phi_array2), shading interp
    hold on
    
    xlabel('x (m)','Fontsize',12,'FontWeight','bold')
    ylabel('z (m)','Fontsize',12,'FontWeight','bold')
    title(['Potential \Phi (x,z)'], ...
           'Fontsize',10,'FontWeight','bold','Color','b')
    axis equal
    axis equal
    caxis([0 M])
    colorbar
    
    %-------------------
    %Isodynamic surfaces
    %-------------------
    figure(2);
    [CS1,H1] = contour(X1,Z1,Phi_array1,cont,'Linewidth',1,'Color','b');
    clabel(CS1,H1,cont);
    hold on
    
    [CS2,H2] = contour(X2,Z2,Phi_array2,cont,'Linewidth',1,'Color','b');
    clabel(CS2,H2,cont);
    hold on
    
    xlabel('x (m)','Fontsize',12,'FontWeight','bold')
    ylabel('z (m)','Fontsize',12,'FontWeight','bold')
    title(['Isodynamic surfaces of Potential \Phi (x,z)'], ...
           'Fontsize',10,'FontWeight','bold','Color','b')
    axis equal 
    