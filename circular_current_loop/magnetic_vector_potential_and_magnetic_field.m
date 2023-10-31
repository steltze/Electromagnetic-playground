function magnetic_vector_potential_and_magnetic_field(d,h,a,I)
    % --------------------------------
    % Example:
    % ask8(2,1,0.1,1);
    % --------------------------------
    red            = [192   0   0]/255;
    Npoints=100;
    m0=I*pi*a^2;
    %μ0 = 4*pi*10^(-7);
    
    %===========
    %Magnetic Vector Potential
    %===========
    
    %-----------
    %x-z plane, y=1m
    
    x=0:4/Npoints:4;
    z=-2:4/Npoints:2;
    
    [X,Z] = meshgrid(x,z);
    
    for i=1:length(x)
        for j=1:length(z)
            R1=sqrt((x(i)-d)^2+(1-h)^2+(z(j))^2);
            R2=sqrt((x(i)+d)^2+(1-h)^2+(z(j))^2);
            R3=sqrt((x(i)-d)^2+(1+h)^2+(z(j))^2);
            R4=sqrt((x(i)+d)^2+(1+h)^2+(z(j))^2);
            
            Az(j,i) = m0*(-(x(i)-d)/(R1^3)-(x(i)+d)/(R2^3)+(x(i)-d)/(R3^3)+(x(i)+d)/(R4^3))/(4*pi);
            Ax(j,i) = m0*(z(j)/(R1^3)+z(j)/(R2^3)-z(j)/(R3^3)-z(j)/(R4^3))/(4*pi);
        end
    end
    
    %Line current coordinates
    theta=0:2*pi/Npoints:2*pi;
    l_curr_x=a*sin(theta)+d;
    l_curr_z=a*cos(theta);
    
    figure(1);
    LL = sqrt((Ax).^2 + (Az).^2);
    quiver(X,Z,Ax./LL,Az./LL,0.5);
    hold on
    
    %p1 = plot(l_curr_x,l_curr_z,'Linewidth',2);
    %set(p1,'Color', red);
    %hold on
    
    hs1 = streamslice(X,Z,Ax,Az,1);
    set(hs1,'Color','m','Linewidth',1.0);
    hold on
    
    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('x (m)','Fontsize',12,'FontWeight','bold');
    ylabel('z (m)','Fontsize',12,'FontWeight','bold');
    title(['Magnetic Vector Potential on x-z plane, y=1m']);
    axis equal
    
    %-----------
    %x-y plane, z=2m
    
    y=h/Npoints:2/Npoints:2*h+h/Npoints;
    [X,Y] = meshgrid(x,y);
    
    for i=1:length(x)
        for j=1:length(y)
            R1=sqrt((x(i)-d)^2+(y(j)-h)^2+4);
            R2=sqrt((x(i)+d)^2+(y(j)-h)^2+4);
            R3=sqrt((x(i)-d)^2+(y(j)+h)^2+4);
            R4=sqrt((x(i)+d)^2+(y(j)+h)^2+4);
            
            %Az(j,i) = m0*(-(x(i)-d)/(R1^3)-(x(i)+d)/(R2^3)+(x(i)-d)/(R3^3)+(x(i)+d)/(R4^3))/(4*pi);
            Axx(j,i) = m0*(2/(R1^3)+2/(R2^3)-2/(R3^3)-2/(R4^3))/(4*pi);
            Ayy(j,i) = 0;
        end
    end
    
    figure(2);
    LL = abs(Axx);
    quiver(x,y,Axx./LL, Ayy,0.5);
    hold on
    
    hs2 = streamslice(x,y,Axx,Ayy,1);
    set(hs2,'Color','m','Linewidth',1.0);
    hold on
    
    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('x (m)','Fontsize',12,'FontWeight','bold');
    ylabel('y (m)','Fontsize',12,'FontWeight','bold');
    title(['Magnetic Vector Potential on x-y plane, z=2m']);
    axis equal
    
    %===========
    %Magnetic Field
    %===========
    %x(abs(x-d)<=0.1)=[];
    
    for i=1:length(x)
        for j=1:length(y)       
            R1=sqrt((x(i)-d)^2+(y(j)-h)^2);
            R2=sqrt((x(i)+d)^2+(y(j)-h)^2);
            R3=sqrt((x(i)-d)^2+(y(j)+h)^2);
            R4=sqrt((x(i)+d)^2+(y(j)+h)^2);
            
            H1_x(j,i)=3*m0*(y(j)-h)*(x(i)-d)/(4*pi*R1^5);
            H2_x(j,i)=3*m0*(y(j)-h)*(x(i)+d)/(4*pi*R2^5);
            H3_x(j,i)=-3*m0*(y(j)+h)*(x(i)-d)/(4*pi*R3^5);
            H4_x(j,i)=-3*m0*(y(j)+h)*(x(i)+d)/(4*pi*R4^5);
            
            H1_y(j,i)=((3*m0*(y(j)-h)^2/R1^2)-m0)/(4*pi*R1^3);
            H2_y(j,i)=((3*m0*(y(j)-h)^2/R2^2)-m0)/(4*pi*R2^3);
            H3_y(j,i)=((-3*m0*(y(j)+h)^2/R3^2)+m0)/(4*pi*R3^3);
            H4_y(j,i)=((-3*m0*(y(j)+h)^2/R4^2)+m0)/(4*pi*R4^3);
        end
    end
    
    H_x = H1_x + H2_x + H3_x + H4_x;
    H_y = H1_y + H2_y + H3_y + H4_y;
    
    %Line current coordinates
    l_curr_x=d-a:2*a/Npoints:d+a;
    l_curr_y=h*ones(1,length(l_curr_x));
    
    figure(3);
    LL = sqrt((H_x).^2 + (H_y).^2);
    quiver(X,Y,H_x./LL,H_y./LL,0.5);
    hold on
    
    p2 = plot(l_curr_x,l_curr_y,'Linewidth',2);
    set(p2,'Color', red);
    hold on
    
    hs3 = streamslice(X,Y,H_x,H_y,1);
    set(hs3,'Color','m','Linewidth',1.0);
    hold on
    
    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('x (m)','Fontsize',12,'FontWeight','bold');
    ylabel('y (m)','Fontsize',12,'FontWeight','bold');
    title(['Magnetic Field on x-y plane, z=0m']);
    axis equal