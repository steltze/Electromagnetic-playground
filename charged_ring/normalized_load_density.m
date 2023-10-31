function charged_ring_load_density(a,d,h)   
    %----------------------------------
    %calculation of load density
    %----------------------------------
     
    %Circular ring
    theta_r    = 0 : 2*pi/(5*Npoints) : 2*pi;
    xx_r   = 1+ a*cos(theta_r);
    zz_r   = a*sin(theta_r);
     
    zz=-2*h:4*h/Npoints:2*h;
    [X2,Z]=meshgrid(xx,zz);
    Sarray = - integral(@(theta)sig(theta,X2,Z,d,h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12) ...
             - integral(@(theta)sig(theta,X2,Z,d,-h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)...
             + integral(@(theta)sig(theta,X2,Z,-d,-h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)...
             + integral(@(theta)sig(theta,X2,Z,-d,h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12);
     
    Sig_max=max(max(Sarray));
     
    cont2 = [-10,-7.5,-5,-4,-3,-2,-1,-0.25,-0.1];
     
    figure(5);
    [CS2,H2] = contour(X2,Z,Sarray,cont2,'Linewidth',1,'Color','b');
    clabel(CS2,H2,cont2);
    hold on
    cr = plot(xx_r,zz_r,'Linewidth',3);
    set(cr,'Color', red);
     
    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('x','Fontsize',12,'FontWeight','bold');
    ylabel('z','Fontsize',12,'FontWeight','bold');
    title(['Contours of normalized load density \sigma (x,y)/(\lambda\alpha/4\pi\epsilon_0) C/(V m^3) or F/m^3']);
     
    axis equal
    grid on
    % ==================================
    function Phi = potential(theta,x,y,d,h,a)
     
    R = sqrt((x-d).^2 + (y-h).^2 + a.^2 - 2*a*(x-d)*cos(theta));
    Phi=1./R;
    end
    % ==================================
    % ==================================
    function dEx = Ex_comp (theta,x,y,d,h,a)
        
    R = sqrt((x-d).^2 + (y-h).^2 + a.^2 - 2*a*(x-d)*cos(theta));
    dEx  =  (x-d-a*cos(theta))./(R.^3);
    end
    % ==================================
    % ==================================
    function dEy = Ey_comp (theta,x,y,d,h,a)
        
    R = sqrt((x-d).^2 + (y-h).^2 + a.^2 - 2*a*(x-d)*cos(theta));
    dEy  =  (y-h)./(R.^3);
    end
    % ==================================
    % ==================================
    function S = sig(theta,x,z,d,h,a)
     
    R = sqrt((x-d).^2 + h.^2+ z.^2 + a.^2 - 2*a*((x-d)*cos(theta)+z*sin(theta)));
    S=1./(R.^3);
    end
    % ==================================
    end