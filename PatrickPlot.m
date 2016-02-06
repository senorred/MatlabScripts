
SS = PatrickData(:,1);
cv2 = PatrickData(:,2);
Radius = PatrickData(:,3)./2;
Radius10 = Radius .* 10;
range = max(Radius10) - min(Radius10);

SS = SS .* (Radius.^3)*4/3*pi();

Radius10 = Radius10 - min(Radius10);
Radius10 = round(Radius10);
c = colormap(jet(range));

hold on

for i = 1:length(SS)
    if Radius10(i) == 0
        Radius10(i) = 1;
    end
    plot(SS(i),cv2(i),'marker','.','markersize',12,'color',c(Radius10(i),:))
    
end

set(gca,'fontsize',20);
%axis([100 1000000 .001 10])
set(gca,'XScale','log');
set(gca,'YScale','log');
xlabel('<p>','FontSize',20)
ylabel('cv2','FontSize',20)
title('Abundance vs. cv2')