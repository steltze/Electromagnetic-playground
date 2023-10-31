function charged_cavity_electric_field(a,b,d,L,D)
 
    % --------------------------------
    % Example:
    % ask6_delta(10,5,3,3,2.5)
    % --------------------------------
    Theta=L/d;
    d_mirror = b^2./d;
     
    white  = [1 1 1];
    black  = [0 0 0];
    red            = [192   0   0]/255;
    green          = [119 172  48]/255; 
    darker_green   = [ 26 108  28]/255;
    purple         = [125  46 143]/255;
    magenta_light  = [255 102 255]/255;
     
     
    Npoints=100;
     
    ymin = -15;
    ymax = 15;
    zmin = -15;
    zmax = 15;
     
    yy = ymin:(ymax-ymin)/Npoints:ymax;
    zz = zmin:(zmax-zmin)/Npoints:zmax;
    [Cy,Cz,Farray1,Farray2]=ask6_gamma(10,5,3,3,2.5,Npoints);
     
    [Y,Z]=meshgrid(yy,zz); 
     
    for i=1:length(yy)
        for j=1:length(zz)               
            if sqrt(Y(i,j).^2+(Z(i,j)-D).^2) < b
                Eyarray(i,j) = 3*integral(@(theta)Ey_comp(theta,0,Y(i,j),Z(i,j),d,D),-Theta,Theta,'RelTol',1e-6,'AbsTol',1e-12)-5*integral(@(theta)Ey_comp(theta,0,Y(i,j),Z(i,j),d_mirror,D),-Theta,Theta,'RelTol',1e-6,'AbsTol',1e-12);
                Ezarray(i,j) = 3*integral(@(theta)Ez_comp(theta,0,Y(i,j),Z(i,j),d,D),-Theta,Theta,'RelTol',1e-6,'AbsTol',1e-12)-5*integral(@(theta)Ez_comp(theta,0,Y(i,j),Z(i,j),d_mirror,D),-Theta,Theta,'RelTol',1e-6,'AbsTol',1e-12);
            elseif sqrt(Y(i,j).^2+(Z(i,j)).^2) < a
                Eyarray(i,j)=0;
                Ezarray(i,j)=0;
            else Eyarray(i,j)=(a*Y(i,j))./(((Y(i,j)).^2 + (Z(i,j)).^2).^(3/2));
                 Ezarray(i,j)=(a*Z(i,j))./(((Y(i,j)).^2 + (Z(i,j)).^2).^(3/2));
            end
        end
    end
     
    % Load charge coordinates
    theta2  = -Theta : (2*Theta)/(Npoints) : Theta;
    yy_ch = d*sin(theta2);
    zz_ch = d*cos(theta2)+D;
     
    % Spherical cavity coordinates
    theta3    = 0 : 2*pi/(5*Npoints) : 2*pi;
    yy_cav   = b*sin(theta3);
    zz_cav   = b*cos(theta3)+D;
     
    % Conductive sphere coordinates
    yy_s   = a*sin(theta3);
    zz_s   = a*cos(theta3);
     
    cont = [0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1, 1.05, 1.10, 1.25, 2, 3, 5]; 
    figure(3);
    LL = ((Eyarray).^2 + (Ezarray).^2).^(1/2);
    quiver(Y,Z,Eyarray./LL,Ezarray./LL,0.5);
    hold on  
     
    contour(Cy,Cz,Farray1,cont,'Linewidth',1,'Color','b');
    hold on
    contour(yy,zz,Farray2,cont,'Linewidth',1,'Color','b');
    hold on
     
    p1 = plot(yy_ch,zz_ch,'Linewidth',3);
    set(p1,'Color', red);
    hold on
    p2 = plot(yy_s,zz_s,'Linewidth',2);
    set(p2,'Color', black);
    hold on
    p3 = plot(yy_cav,zz_cav,'Linewidth',2);
    set(p3,'Color', black);
     
    set(gca,'Fontsize',12,'Fontweight','bold')
    xlabel('y (m)','Fontsize',12,'FontWeight','bold')
    ylabel('z (m)','Fontsize',12,'FontWeight','bold')
    title(['Normalized Electric Field E (x,y)/(\lambda/4\pi\epsilon_0)'], ...
           'Fontsize',10,'FontWeight','bold','Color','b')
    axis equal
     
    figure(4);
     
    hs = streamslice(Y,Z,Eyarray,Ezarray,2);
    set(hs,'Color','m','Linewidth',1.0);
    hold on
     
    contour(Cy,Cz,Farray1,cont,'Linewidth',1,'Color','b');
    hold on
     
    contour(yy,zz,Farray2,cont,'Linewidth',1,'Color','b');
    hold on
     
    p1 = plot(yy_ch,zz_ch,'Linewidth',3);
    set(p1,'Color', red);
    hold on
    p2 = plot(yy_s,zz_s,'Linewidth',2);
    set(p2,'Color', black);
    hold on
    p3 = plot(yy_cav,zz_cav,'Linewidth',2);
    set(p3,'Color', black);
     
    set(gca,'Fontsize',12,'Fontweight','bold')
    xlabel('y (m)','Fontsize',12,'FontWeight','bold')
    ylabel('z (m)','Fontsize',12,'FontWeight','bold')
    title(['Normalized Electric Field E (x,y)/(\lambda/4\pi\epsilon_0)'], ...
           'Fontsize',10,'FontWeight','bold','Color','b')
    axis equal
     
    % =====================================
    function dEy = Ey_comp (theta,x,y,z,d,D)
    R    = (x.^2 + (y-d*sin(theta)).^2 + (z-d*cos(theta)-D).^2 ).^(1/2);
    dEy  =  (y-d*sin(theta))./(R.^3);
    end
    % =====================================
     
    % =====================================
    function dEz = Ez_comp (theta,x,y,z,d,D)
    R    = (x.^2 + (y-d*sin(theta)).^2 + (z-d*cos(theta)-D).^2 ).^(1/2);
    dEz  =  (z-d*cos(theta)-D)./(R.^3);
    end
    % =====================================
    end