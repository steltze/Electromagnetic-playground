
function current_distributions(d,h,a,I)  
    %===========
    %Current distributions
    %===========

    %-----------
    % x=0
    %-----------
    y = 0:4/Npoints:4;
    [Y,Z]=meshgrid(y,z);

    for i=1:length(y)
        for j=1:length(z)       
            R1=sqrt(d^2+(y(i)-h)^2+(z(j))^2);
            R3=sqrt(d^2+(y(i)+h)^2+(z(j))^2);
            
            K1_z(j,i) = (((3*m0*((y(i)-h)^2)/(R1^2))-m0)/(R1^3)+((-3)*m0*((y(i)+h)^2)/(R3^2)+m0)/(R3^3))/(2*pi);
            K1_y(j,i) = (-6)*m0*(y(i)-h)*z(j)/((R1^5)*4*pi) + 6*m0*(y(i)+h)*z(j)/((R3^5)*4*pi);	
        end
    end
        
    figure(4);
    LL = sqrt((K1_y).^2 + (K1_z).^2);
    quiver(Y,Z,K1_y./LL,K1_z./LL,0.5);
    hold on

    hs4 = streamslice(Y,Z,K1_y,K1_z,1);
    set(hs4,'Color','m','Linewidth',1.0);
    hold on

    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('y (m)','Fontsize',12,'FontWeight','bold');
    ylabel('z (m)','Fontsize',12,'FontWeight','bold');
    title(['Current distribution on plane x=0m']);
    axis equal

            
    %-----------
    % y=0
    %-----------
    for i=1:length(x)
        for j=1:length(z)       
            R1=sqrt((x(i)-d)^2+h^2+(z(j))^2);
            R2=sqrt((x(i)+d)^2+h^2+(z(j))^2);
            
            K2_z(j,i) = h*3*m0*((x(i)-d)/(2*pi*(R1^5)) +(x(i)+d)/(2*pi*(R2^5)));
            K2_x(j,i) = -3*m0*h*z(j)*(1/(4*pi*(R1^5))+1/(4*pi*(R2^5)));
        end
    end

    figure(5);
    LL = sqrt((K2_x).^2 + (K2_z).^2);
    quiver(x,z,K2_x./LL,K2_z./LL,0.5);
    hold on

    hs4 = streamslice(x,z,K2_x,K2_z,1);
    set(hs4,'Color','m','Linewidth',1.0);
    hold on

    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('x (m)','Fontsize',12,'FontWeight','bold');
    ylabel('z (m)','Fontsize',12,'FontWeight','bold');
    title(['Current distribution on plane y=0m']);
    axis equal

    end