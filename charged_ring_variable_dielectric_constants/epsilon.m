function [charge1] = ask7_epsilon(L, a, N, h, V, er1, er2)

    %στην έξοδο στη μεταβλητή charge εμφανίζεται η κατανομή φορτίου λ/ε0
    
    %============
    %[charge] = ask7_epsilon(0.5, 0.0025, 50, 0.5, 1, 1, 5);
    %[charge] = ask7_epsilon(0.5, 0.0025, 50, 0.5, 1, 5, 1);
    %============
    
    
    Npoints = 200;
    Length = 2*L;
    Dx     = Length/N;
    x      = -L:Dx:L;
    bvector = ones(N,1);
    X = -2: 4/Npoints : 2;
    Farray = zeros([1 length(X)]);
    
    for i=1:N
       xx(i) = 0.5*( x(i+1) + x(i) );
    end
    
    %-------------
    for i=1:N
      for j=1:N
        if i==j        
           R1 = (Dx/2)+sqrt((Dx/2)^2+a^2);
           R2 =(-Dx/2)+sqrt((Dx/2)^2+a^2);
           R3 = 2*h;
           BB(i,j) = log(R1/R2)+(Dx*(er1-er2))/((er1+er2)*R3);
           
        else
           R1 = abs(xx(i)-xx(j));       
           R3 = sqrt((R1)^2 + 4*h^2);
           BB(i,j) = Dx*(1/R1 + (er1-er2)/((er1+er2)*R3));
        end
      end
    end
    
    charge1 = BB\bvector;
    charge1 = 4*pi*er1*V*charge1;
    totalc1 = sum(charge1(:,1))*Dx;
    charge1 = charge1/totalc1;
    
    %-----------
    %line charge density plot
    %-----------
    figure(3);
    plot(xx,charge1,'linewidth',2);
    
    xlabel('x (m)','Fontsize',14);
    ylabel('Line charge (V)','Fontsize',14);
    
    grid on
    
    %----------- 
    %Potential on line  -2m<=x<=2m ,y=z=0
    %-----------
    
    for i=1:length(X)
         for j=1:N
             R1 = sqrt((X(i) - xx(j))^2 + h^2);
             R2 = (er1-er2)/(er1+er2);
             Farray(i) = Farray(i) + charge1(j)*Dx*(1+R2)/R1;
        end
    end
    
    Farray = Farray/(4*pi*er1);
    
    figure(4);
    plot(X, Farray, 'linewidth',2);
    
    xlabel('x (m)','Fontsize',14);
    ylabel('Potential on -2m\leq x \leq2m line (V)','Fontsize',14);
    
    grid on
    
    end