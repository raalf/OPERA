function fcnPLOTCOEFF(valTIMESTEP, matCOEFF_HSTRY)

linestyles = {'--';'-.';'-';':'};
markers = {'o';'x';'s';'^';'*';'d';'v';'>';'<';'p';'h'};
colors = {'k';'b';'r';'m'};

hFig21 = figure(21);
clf(21);
hold on
for i = 1:5
    pFig1 = plot(1:valTIMESTEP, abs(reshape(sum((matCOEFF_HSTRY(:,i,2:end) - matCOEFF_HSTRY(:,i,1:end-1))./matCOEFF_HSTRY(:,i,1:end-1),1),[],1,1)).*100,'-ok');
    pFig1.LineStyle = linestyles{1+mod(i,4),:};
    pFig1.Marker = markers{1+mod(i,11),:};
    pFig1.Color = colors{1+mod(i,4),:};
end
hold off
ylim([0, inf]);
grid minor
box on
hold off
xlabel('Timestep','FontSize',15)
ylabel('Percent Change in Value','FontSize',15)
legend('A_1','A_2','B_1','B_2','C_3','Location','NorthEast')
set(gca, 'YScale', 'log')

drawnow;
end

