load('polarisation_array_V03.mat');
rtt=polarr(:,4);
pi=1;
ska={};
reduced = [];
ipi = 1;
slopes=[];
dslopes=[];
for i=1:142
    %disp(24*60*60*(rtt(i)-rtt(pi)))
if floor(24*6000*(rtt(i)-rtt(pi)))~=0
    figure(2)
    v = rand(3,1);
    x=polarr(pi:i-1,3);
    y=polarr(pi:i-1,1);
    plot(x,y,"Marker",".","MarkerSize",30,"LineStyle",'--',"DisplayName",datestr(rtt(pi)),'Color',v)
    %[p,s] = polyfit(x,y,1);
    %g = fittype('a*x','independent','x')
    %[fitresult, gof1, out1] = fit(x,y,g);
    [fitresult, gof1, out1] = fit(x,y,'poly1');
    coeffs = coeffvalues(fitresult);
    confints = confint(fitresult,0.5);
    dm = confints(2,1)-confints(1,1); % uncertainty of the slope
    %ska{ipi}=s;
    y_est = feval(fitresult,x);
    % Add trend line to plot
    hold on
    plot(x,y_est,'LineWidth',2,"DisplayName",datestr(rtt(pi)),'Color',v)
    %figure(3)
    %x=r;
    %x=polarr(:,15);
    %plot(mean(x(pi:i)),c(1),'Marker','.','MarkerSize',20,"DisplayName",datestr(rtt(pi)))
    slopes(ipi)=coeffs(1);
    dslopes(ipi)=dm;
    for j=1:15
        
        reduced(ipi,j)=mean(polarr(pi:i,j));
    end
    reduced(ipi,16)=coeffs(1);
    ipi=ipi+1;
hold on
pi=i;
end
end
legend()
hold off