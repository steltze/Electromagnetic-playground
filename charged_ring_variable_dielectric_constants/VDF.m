function [VDF,current,Phi_cond,Rg,Phi_000]=ask8(N, a, L, h, D, sigma0, I)

    %============
    %[VDF,I,Phi_cond,Rg,Phi_000]=ask8(3,0.005,0.99,1.00,1.5,1/160,500);
    %============
    
    Npoints=200;
    Dx = D/N;
    Dz = L/N;
    x = -D/2:Dx:D/2;
    z = -h-a:-Dz:-h-a-L;
    Phivector = ones(3*N,1);
    X = -5: 10/Npoints : 5;
    Farray = zeros([1 length(X)]);
    
    for i=1:N
       xx(i) = 0.5*( x(i+1) + x(i) );
    end
    
    for i=1:N
       zz(i) = 0.5*( z(i+1) + z(i) );
    end
    
    for i=1:N
      for j=1:N
         B1(i,j) = Dz/(4*pi*sigma0)*(1/(sqrt((xx(i)+D/2-a)^2+(h+zz(j))^2))+1/(sqrt((xx(i)+D/2-a)^2+(-h+zz(j))^2)));
         C1(i,j) = Dz/(4*pi*sigma0)*(1/(sqrt((xx(i)-D/2+a)^2+(h+zz(j))^2))+1/(sqrt((xx(i)-D/2+a)^2+(-h+zz(j))^2)));
         A2(i,j) = Dx/(4*pi*sigma0)*(1/(sqrt((xx(j)+D/2-a)^2+(h+zz(i))^2))+1/(sqrt((xx(j)+D/2-a)^2+(-h+zz(i))^2))); 
         C2(i,j) = Dz/(4*pi*sigma0)*(1/(sqrt((D-2*a)^2+(zz(i)-zz(j))^2))+1/(sqrt((D-2*a)^2+(zz(j)+zz(i))^2)));
         A3(i,j) = Dx/(4*pi*sigma0)*(1/(sqrt((xx(j)-D/2+a)^2+(h+zz(i))^2))+1/(sqrt((xx(j)-D/2+a)^2+(-h+zz(i))^2)));
         B3(i,j) = Dz/(4*pi*sigma0)*(1/(sqrt((D-2*a)^2+(zz(i)-zz(j))^2))+1/(sqrt((D-2*a)^2+(zz(j)+zz(i))^2)));
      end
    end
    
    R1 = (Dx/2)+sqrt(a^2+(Dx/2)^2);
    R2 = (-Dx/2)+sqrt(a^2+(Dx/2)^2);
    R3 = (Dz/2)+sqrt(a^2+(Dz/2)^2);
    R4 = (-Dz/2)+sqrt(a^2+(Dz/2)^2);
    
    for i=1:N
      for j=1:N
        if i==j
           A1(i,j) = 1/(4*pi*sigma0)*(log(R1/R2)+Dx/(sqrt((xx(i)-xx(j))^2+4*h*h)));
           B2(i,j) = 1/(4*pi*sigma0)*(log(R3/R4)+Dz/(abs((zz(i)+zz(j)))));
           C3(i,j) = 1/(4*pi*sigma0)*(log(R3/R4)+Dz/(abs((zz(i)+zz(j)))));
        else
           A1(i,j) = Dx/(4*pi*sigma0)*(1/(abs(xx(i)-xx(j)))+1/(sqrt((xx(i)-xx(j))^2+4*h*h)));
           B2(i,j) = Dz/(4*pi*sigma0)*(1/(abs(zz(i)-zz(j)))+1/(abs(zz(i)+zz(j))));
           C3(i,j) = Dz/(4*pi*sigma0)*(1/(abs(zz(i)-zz(j)))+1/(abs(zz(i)+zz(j))));	
        end
      end  
    end
    
    VDF = [A1 B1 C1; A2 B2 C2; A3 B3 C3];
    %Φ0=1
    current = VDF\Phivector; 
    flag=0;
    
    for i=1:N
    flag = flag + current(i)*Dx + (current(i+N)+current(i+2*N))*Dz;
    end
    
    current = I*current/(flag); %current distribution on π-shaped conductor
    
    %===============
    %Beta
    %===============
    %Φ0=1
    Phi_cond = I/flag; %potential on π-shaped conductor
    Rg = 1/flag; %ground resistance 
    
    %-------------
    %Potential on -2m<=x<=2m ,y=z=0
    %-------------
    
    p=0;
    for i=1:length(X)
             for j=1:N
                 R1 = sqrt((X(i)-xx(j))^2+h^2);
                 R2 = sqrt((X(i)+D/2-a)^2+zz(j)^2);
                 R3 = sqrt((X(i)-D/2+a)^2+zz(j)^2);
                 Farray(i) = Farray(i) + current(j)*Dx/R1 + (current(j+N)/R2+current(j+2*N)/R3)*Dz;
             
             end
         if X(i)==0
             gege=i; %point x=y=z=0
         end
    end
    
    Farray = 2*Farray/(4*pi*sigma0);
    figure(1);
    plot(X, Farray, 'linewidth',2);
    
    xlabel('x (m)','Fontsize',14);
    ylabel('Potential on -5m\leq x \leq5m line (V)','Fontsize',14);
    
    grid on
    
    Phi_000 = Farray(gege); %Potential on point x=y=z=0
    
    end