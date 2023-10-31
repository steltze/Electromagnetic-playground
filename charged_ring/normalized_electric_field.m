function charged_ring_electric_field(a,d,h)   
    %----------------------------------
    %calculation of electric field
    %----------------------------------
     
    Ex_array = integral(@(theta)Ex_comp(theta,X,Y,d,h,a),0,2*pi,'ArrayValued',true,'RelTol',1e-12,'AbsTol',1e-12) - integral(@(theta)Ex_comp(theta,X,Y,d,-h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)+integral(@(theta)Ex_comp(theta,X,Y,-d,-h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)-integral(@(theta)Ex_comp(theta,X,Y,-d,h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12);
    Ey_array = integral(@(theta)Ey_comp(theta,X,Y,d,h,a),0,2*pi,'ArrayValued',true,'RelTol',1e-12,'AbsTol',1e-12) - integral(@(theta)Ey_comp(theta,X,Y,d,-h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)+integral(@(theta)Ey_comp(theta,X,Y,-d,-h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12)-integral(@(theta)Ey_comp(theta,X,Y,-d,h,a),0,2*pi,'ArrayValued',true,'RelTol',0,'AbsTol',1e-12);
     
    figure(3);
    LL = sqrt((Ex_array).^2 + (Ey_array).^2);
    quiver(X,Y,Ex_array./LL,Ey_array./LL,0.5);
    hold on
     
    hs = streamslice(X,Y,Ex_array,Ey_array,2);
    set(hs,'Color','m','Linewidth',1.0);
    hold on
     
    pr = plot(xr,yr,'Linewidth',4);
    set(pr, 'Color', red);
    set(pr,'ZData',Pot_max+1+zeros(size(xr)))   
     
    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('x (m)','Fontsize',12,'FontWeight','bold');
    ylabel('y (m)','Fontsize',12,'FontWeight','bold');
    title(['Normalized Electric Field E(x,y)/(\lambda\alpha/4\pi\epsilon_0) (1/m^2)']);
    axis equal
     
    figure(4);
    hs = streamslice(X,Y,Ex_array,Ey_array,2);
    set(hs,'Color','m','Linewidth',1.0);
    hold on
     
    contour(X,Y,Parray,cont,'Linewidth',1,'Color','b');
    %clabel(CS,H,cont);
    pr = plot(xr,yr,'Linewidth',4);
    set(pr, 'Color', red);
    set(pr,'ZData',Pot_max+1+zeros(size(xr)))   
     
    set(gca,'Fontsize',12,'Fontweight','bold');
    xlabel('x (m)','Fontsize',12,'FontWeight','bold');
    ylabel('y (m)','Fontsize',12,'FontWeight','bold');
    title(['Normalized Electric Field E(x,y)/(\lambda\alpha/4\pi\epsilon_0) (1/m^2)']);
    axis equal