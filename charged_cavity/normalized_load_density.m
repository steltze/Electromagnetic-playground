function charged_cavity_load_density(a,b,d,L,D)
    % --------------------------------
    % Example:
    % ask6_epsilon(10,5,3,3,2.5)
    % --------------------------------
    Npoints=250;
    Theta=L/d;
    d_mirror = b^2./d;
     
    u=0:2*pi/Npoints:2*pi;
     
    Y=b*sin(u);
    Z=b*cos(u)+D;
    Eyarray = 3*integral(@(theta)Ey_comp(theta,0,Y,Z,d,D),-Theta,Theta,...
    'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)-5*integral(@(theta)Ey_comp(theta,0,Y,Z,d_mirror,D),-Theta,Theta,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12);
    Ezarray = 3*integral(@(theta)Ez_comp(theta,0,Y,Z,d,D),-Theta,Theta,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)-5*integral(@(theta)Ez_comp(theta,0,Y,Z,d_mirror,D),-Theta,Theta,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12);
     
    sigma = zeros(length(u));
    for i=1:length(u)
        sigma(i)=-(sin(u(i))*Eyarray(i)+cos(u(i))*Ezarray(i))./4*pi;
    end
     
    plot(u*180/pi,sigma);
    hold on
    set(gca,'Fontsize',12,'Fontweight','bold')
    xlabel('Polar Angle Theta (degrees)','Fontsize',12,'FontWeight','bold')
    ylabel('Normalized Surface Charge sigma(y,z)/lamda (1/m)','Fontsize',12,'FontWeight','bold')
    title(['Induces Surface Charge sigma(y,z)/lamda'], ...
           'Fontsize',10,'FontWeight','bold','Color','b')
    grid on
    % ==================================
    function dEy = Ey_comp (u,x,y,z,d,D)
    R    = (x.^2 + (y-d*sin(u)).^2 + (z-d*cos(u)-D).^2 ).^(1/2);
    dEy  =  (y-d*sin(u))./(R.^3);
    end
    % ==================================
     
    % ==================================
    function dEz = Ez_comp (u,x,y,z,d,D)
    R    = (x.^2 + (y-d*sin(u)).^2 + (z-d*cos(u)-D).^2 ).^(1/2);
    dEz  =  (z-d*cos(u)-D)./(R.^3);
    end
    % ==================================
     
    end