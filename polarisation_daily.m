pi=1;
ska={};
reduced = [];
ipi = 1;
slopes=[];
for i=1:142
if floor(10*(rtt(i)-rtt(pi)))~=0
    figure(2)
    v = rand(3,1);
    x=polarr(pi:i-1,3);
    y=polarr(pi:i-1,1);
    plot(x,y,"Marker",".","MarkerSize",30,"LineStyle",'--',"DisplayName",datestr(rtt(pi)),'Color',v)
    [c,s] = polyfit(x,y,1);
    ska{ipi}=s;
    y_est = polyval(c,x);
    % Add trend line to plot
    hold on
    plot(x,y_est,'LineWidth',2,"DisplayName",datestr(rtt(pi)),'Color',v)
    %figure(3)
    %x=r;
    %x=polarr(:,15);
    %plot(mean(x(pi:i)),c(1),'Marker','.','MarkerSize',20,"DisplayName",datestr(rtt(pi)))
    slopes(ipi)=c(1);
    for j=1:15
        
        reduced(ipi,j)=mean(polarr(pi:i,j));
    end
    ipi=ipi+1
hold on
pi=i;
end
end
legend()
hold off