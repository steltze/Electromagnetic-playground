function charged_ring_electric_field(L,h,er1,er2)
    %-------------------
    %Electric Field
    %-------------------
    
    %z>=0
    Ex_array1=(1/(4*pi*er1)).*(integral(@(x0)Ex_comp(x0,X2,Z2,h),-L,L,'ArrayValued',true,'RelTol',1e-12,'AbsTol',1e-12)+((er1-er2)/(er1+er2)).*integral(@(x0)Ex_comp(x0,X2,Z2,-h),-L,L,'ArrayValued',true,'RelTol',1e-12,'AbsTol',1e-12));
    Ez_array1=(1/(4*pi*er1)).*(integral(@(x0)Ez_comp(x0,X2,Z2,h),-L,L,'ArrayValued',true,'RelTol',1e-12,'AbsTol',1e-12)+((er1-er2)/(er1+er2)).*integral(@(x0)Ez_comp(x0,X2,Z2,-h),-L,L,'ArrayValued',true,'RelTol',1e-12,'AbsTol',1e-12));
    
    %z<=0
    Ex_array2=(2/(4*pi*(er1+er2))).*integral(@(x0)Ex_comp(x0,X1,Z1,h),-L,L,'ArrayValued',true,'RelTol',1e-12,'AbsTol',1e-12);
    Ez_array2=(2/(4*pi*(er1+er2))).*integral(@(x0)Ez_comp(x0,X1,Z1,h),-L,L,'ArrayValued',true,'RelTol',1e-12,'AbsTol',1e-12);
    
    zz3=1/(Npoints):2/Npoints:2;
    [X3,Z3]=meshgrid(xx,zz3);
        
    
    figure(3);
    LL1=sqrt(Ex_array1.^2+Ez_array1.^2);
    LL2=sqrt(Ex_array2.^2+Ez_array2.^2);
    
    p4 = quiver(X3,Z3,Ex_array1./LL1,Ez_array1./LL1,0.5);
    set(p4,'Color', purple);
    hold on
    p5 = quiver(X1,Z1,Ex_array2./LL2,Ez_array2./LL2,0.5);
    set(p5,'Color', purple);
    hold on
    
    contour(X1,Z1,Phi_array1,cont,'Linewidth',1,'Color','b');
    hold on
    contour(X2,Z2,Phi_array2,cont,'Linewidth',1,'Color','b');
    hold on
    
    title(['Normalized Electric Field E (x,z)'], ...
           'Fontsize',10,'FontWeight','bold','Color','b')
    xlabel('x (m)','Fontsize',12,'FontWeight','bold')
    ylabel('z (m)','Fontsize',12,'FontWeight','bold')
    
    axis equal
    
    %============
    function fact = factors(x,z,L,h)
    fact = (L-x+sqrt((L-x).^2+(z-h).^2))./(-L-x+sqrt((L+x).^2+(z-h).^2));
    end
    %============
    % ==================================
    function dEx = Ex_comp (x0,x,z,h)
        
    R = sqrt((x-x0).^2+(z-h).^2);
    dEx  =  (x-x0)./(R.^3);
    end
    % ==================================
    % ==================================
    function dEz = Ez_comp (x0,x,z,h)
        
    R = sqrt((x-x0).^2+(z-h).^2);
    dEz  =  (z-h)./(R.^3);
    end
    % ==================================
    end