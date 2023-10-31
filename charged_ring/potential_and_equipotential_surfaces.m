function charged_ring_potential_isodynamic_surfaces(a,d,h)   
 
    % --------------------------------
    % Example:
    % ask7_delta(0.1,1,1)
    % --------------------------------
    white  = [1 1 1];
    black  = [0 0 0];
    red            = [192   0   0]/255;
    green          = [119 172  48]/255; 
    darker_green   = [ 26 108  28]/255;
    purple         = [125  46 143]/255;
    magenta_light  = [255 102 255]/255;
     
    Npoints=250;
    factor=1;
    epsilon0 = 8.85418782e-12;
    % ----------------------------------------------------
    % potential inside the spherical cavity
    % ----------------------------------------------------
     
    xx=0:2*d/Npoints:2*d;
    yy=0:2*h/Npoints:2*h;
     
    [X,Y]=meshgrid(xx,yy);
            
    Parray=integral(@(theta)potential(theta,X,Y,d,h,a),0,2*pi,'ArrayValued',true,'RelTol',1e-12,'AbsTol',1e-12) - integral(@(theta)potential(theta,X,Y,d,-h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)+integral(@(theta)potential(theta,X,Y,-d,-h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)-integral(@(theta)potential(theta,X,Y,-d,h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12);
         
    xr(1) = d-a;
    yr(1) = h;
    xr(2) = d+a;
    yr(2) = h;
    Parray(isinf(Parray))=2.346956877521396e+02; %antikatastash tou apeirou
    Parray=a*Parray;
    Pot_max = max(max(Parray))/factor;
     
    figure(1);
    surface(X,Y,Parray/factor), shading interp
    hold on
    pr = plot(xr,yr,'Linewidth',2);
    set(pr, 'Color', red);
    set(pr,'ZData',Pot_max+1+zeros(size(xr)))   
     
    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('x (m)','Fontsize',12,'FontWeight','bold');
    ylabel('y (m)','Fontsize',12,'FontWeight','bold');
    title(['Normalized Potential \Phi (x,y)/(\lambda\alpha/4\pi\epsilon_0)']);
          
     
    axis equal
    caxis([0 8])
    colorbar
     
    %--------------------------------------
    %Isodynamic surfaces
    %--------------------------------------
    cont = [0.01, 0.1, 0.2, 0.4, 0.8, 1.0, 1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 7.5];
     
    figure(2);
    [CS,H] = contour(X,Y,Parray,cont,'Linewidth',1,'Color','b');
    clabel(CS,H,cont);
    hold on
     
    pr = plot(xr,yr,'Linewidth',2);
    set(pr, 'Color', red);
    set(pr,'ZData',Pot_max+1+zeros(size(xr)))   
     
    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('x (m)','Fontsize',12,'FontWeight','bold');
    ylabel('y (m)','Fontsize',12,'FontWeight','bold');
    title(['Isodynamic surfaces of Potential \Phi (x,y)/(\lambda\alpha/4\pi\epsilon_0)']);
    axis equal