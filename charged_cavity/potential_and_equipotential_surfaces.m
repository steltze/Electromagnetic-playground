function [Y,Z,Farray_1,Farray_3] = charged_cavity_potential_isodynamic_surfaces(a,b,d,L,D,Npoints)
    % --------------------------------
    % Example:
    % [A,B,C,D]=ask6_gamma(10,5,3,3,2.5,100)
    % --------------------------------
     
    Theta=L/d;
     
    white  = [1 1 1];
    black  = [0 0 0];
    red            = [192   0   0]/255;
    green          = [119 172  48]/255; 
    darker_green   = [ 26 108  28]/255;
    purple         = [125  46 143]/255;
    magenta_light  = [255 102 255]/255;
    
    % ----------------------------------------------------
    % potential inside the spherical cavity
    % ----------------------------------------------------
    
    A=0:5/Npoints:5;
    u=0:2*pi/Npoints:2*pi;
    radius1=A(:);
    Z=D+radius1*cos(u);
    Y=radius1*sin(u);
    
    factor=1; 
    
    d_mirror = b^2./d;
    Farray_a = 3*integral(@(theta)potential(theta,0,Y,Z,d,D),-Theta,Theta,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12);
    Farray_b = (-5)*integral(@(theta)potential(theta,0,Y,Z,d_mirror,D),-Theta,Theta,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12);
    Farray_1 = Farray_a + Farray_b + 1;
    
    % ----------------------------------------------------
    % potential inside and outside of the conductive sphere
    % ----------------------------------------------------
    y3min = -15;
    y3max = 15;
    z3min = -15;
    z3max = 15;
    
    y3 = y3min:(y3max-y3min)/Npoints:y3max;
    z3 = z3min:(z3max-z3min)/Npoints:z3max;
    [Y3,Z3]=meshgrid(y3,z3);
    
    Farray_3=zeros(length(y3),length(z3));
    for i=1:length(y3)
        for j=1:length(z3)
        if sqrt(y3(i)^2+(z3(j))^2) < a
            Farray_3(i,j) = 1;
        else
            Farray_3(i,j) = a./(sqrt(y3(i).^2 + z3(j).^2));
        end
        end
    end
    
    % Load charge coordinates
    theta_2  = -Theta : (2*Theta)/(Npoints) : Theta;
    yy_ch = d*sin(theta_2);
    zz_ch = d*cos(theta_2)+D;
    
    % Spherical cavity coordinates
    theta_3    = 0 : 2*pi/(5*Npoints) : 2*pi;
    yy_cav   = b*sin(theta_3);
    zz_cav   = b*cos(theta_3)+D;
    
    % Conductive sphere coordinates
    yy_s   = a*sin(theta_3);
    zz_s   = a*cos(theta_3);
    
    % -----------------------------
    %    Plot Potential
    % -----------------------------
    Phi_max = max(max(Farray_1));

    
    figure(1);
    surface(Y3,Z3,Farray_3), shading interp
    hold on
    surface(Y,Z,Farray_1), shading interp
    hold on
    
    p1 = plot(yy_ch,zz_ch,'Linewidth',3);
    set(p1,'Color', red);
    hold on
    p2 = plot(yy_s,zz_s,'Linewidth',2);
    set(p2,'Color', black);
    hold on
    p3 = plot(yy_cav,zz_cav,'Linewidth',2);
    set(p3,'Color', black);
    
    set(p1,'ZData',1.5*Phi_max+1+zeros(size(zz_ch)))

    set(p2,'ZData',1.5*Phi_max+1+zeros(size(zz_s)))
    
    set(p3,'ZData',1.5*Phi_max+1+zeros(size(zz_cav)))  

    set(gca,'Fontsize',12,'Fontweight','bold')
    xlabel('y (m)','Fontsize',12,'FontWeight','bold')
    ylabel('z (m)','Fontsize',12,'FontWeight','bold')
    title(['Normalized Potential \Phi (x,y)/(\lambda/4\pi\epsilon_0)'], ...
        'Fontsize',10,'FontWeight','bold','Color','b')
    axis equal
    caxis([0 6])
    colorbar
    
    %--------------------------------------
    %Isodynamic surfaces
    %--------------------------------------
    cont = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.25, 2, 3, 5]; 
    
    figure(2);
    [CS1,H1] = contour(Y,Z,Farray_1,cont,'Linewidth',1,'Color','b');
    clabel(CS1,H1,cont);
    hold on
    
    [CS2,H2] = contour(y3,z3,Farray_3,cont,'Linewidth',1,'Color','b');
    clabel(CS2,H2,cont);
    hold on
    
    p1 = plot(yy_ch,zz_ch,'Linewidth',3);
    set(p1,'Color', red);
    hold on
    p2 = plot(yy_s,zz_s,'Linewidth',2);
    set(p2,'Color', black);
    hold on
    p3 = plot(yy_cav,zz_cav,'Linewidth',2);
    set(p3,'Color', black);
    
    set(p1,'ZData',1.5*Phi_max+1+zeros(size(zz_ch))) 
    
    set(p2,'ZData',1.5*Phi_max+1+zeros(size(zz_s)))

    set(p3,'ZData',1.5*Phi_max+1+zeros(size(zz_cav)))  

    set(gca,'Fontsize',12,'Fontweight','bold')
    xlabel('y (m)','Fontsize',12,'FontWeight','bold')
    ylabel('z (m)','Fontsize',12,'FontWeight','bold')
    title(['Isodynamic surfaces of Normalized Potential \Phi (x,y)/(\lambda/4\pi\epsilon_0)'], ...
        'Fontsize',10,'FontWeight','bold','Color','b')
    axis equal  
    
    
    % =====================================
    function Phi = potential(theta,x,y,z,d,D)
    
    R    = (x.^2 + (y-d*sin(theta)).^2 + (z-d*cos(theta)-D).^2 ).^(1/2);
    Phi=1./R;
    end
    % =====================================
 
end