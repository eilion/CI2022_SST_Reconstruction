function saveFigure(outputFile,CURVE)

close all;

path = ['Outputs/',outputFile,'/results.mat'];
AA = load(path);
data = AA.data;

clear AA;

P = length(data);

TITLE = cell(P,1);
for p = 1:P
    tt = data(p).name;
    tt(tt=='_') = '-';
    TITLE{p} = tt;
end


fig = figure;
for p = 1:P
    subplot(3,3,p);
    hold on;
    title(TITLE{p},'FontSize',12);
    
    N = size(data(p).Y,1);
    
    for n = 1:N
        plot([data(p).T(n),data(p).T(n)],[data(p).LB_INDIV(n),data(p).UB_INDIV(n)],'c','LineWidth',2);
    end
    plot(data(p).T,data(p).MD_INDIV,'*b','LineWidth',2);
    
    xlim([0,800]);
    if p == 1
        ylim([-1,18]);
    elseif p == 2
        ylim([6,23]);
    elseif p == 3
        ylim([-0.5,17.5]);
    elseif p == 4
        ylim([5,22]);
    elseif p == 5
        ylim([4,21]);
    elseif p == 6
        ylim([8,23]);
    elseif p == 7
        ylim([20,32]);
    elseif p == 8
        ylim([19,32]);
    elseif p == 9
        ylim([16.5,31.5]);    
    end
    xlabel('time (ka)','FontSize',12);
    ylabel('SST (°C)','FontSize',12);
    
end

set(fig,'Position',[20 20 900 600]);
movegui(fig,'center');
drawnow;

path = ['Outputs/',outputFile,'/INDIV.fig'];
savefig(fig,path);



fig = figure;
for p = 1:P
    subplot(3,3,p);
    hold on;
    title(TITLE{p},'FontSize',12);
    
    N = size(data(p).Y,1);
    
    for n = 1:N
        plot([data(p).T(n),data(p).T(n)],[data(p).LB_INDIV(n),data(p).UB_INDIV(n)],'c','LineWidth',2);
    end
    
    xx = [data(p).T0;flipud(data(p).T0)];
    yy = [data(p).LB0;flipud(data(p).UB0)];
    patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.5,'EdgeColor','none');
    % plot(data(p).T0,data(p).LB0,':k','LineWidth',2);
    plot(data(p).T0,data(p).MD0,'--k','LineWidth',2);
    % plot(data(p).T0,data(p).UB0,':k','LineWidth',2);
    
    xlim([0,800]);
    if p == 1
        ylim([-1,18]);
    elseif p == 2
        ylim([6,23]);
    elseif p == 3
        ylim([-0.5,17.5]);
    elseif p == 4
        ylim([5,22]);
    elseif p == 5
        ylim([4,21]);
    elseif p == 6
        ylim([8,23]);
    elseif p == 7
        ylim([20,32]);
    elseif p == 8
        ylim([19,32]);
    elseif p == 9
        ylim([16.5,31.5]);    
    end
    xlabel('time (ka)','FontSize',12);
    ylabel('SST (°C)','FontSize',12);
    
end

set(fig,'Position',[20 20 900 600]);
movegui(fig,'center');
drawnow;

path = ['Outputs/',outputFile,'/Summary.fig'];
savefig(fig,path);




for p = 1:P
    fig = figure;
    
    h = zeros(3,1);
    subplot(3,1,1:2);
    hold on;
    title(TITLE{p},'FontSize',16);
    
    N = size(data(p).Y,1);
    
    for n = 1:N
        h(2) = plot([data(p).T(n),data(p).T(n)],[data(p).LB_INDIV(n),data(p).UB_INDIV(n)],'c','LineWidth',2);
    end
    
    xx = [data(p).T0;flipud(data(p).T0)];
    yy = [data(p).LB0;flipud(data(p).UB0)];
    patch(xx,yy,1,'FaceColor','k','FaceAlpha',0.25,'EdgeColor','none');
    plot(data(p).T0,data(p).LB0,':k','LineWidth',2);
    plot(data(p).T0,data(p).UB0,':k','LineWidth',2);
    
    plot(data(p).T,data(p).MD_INDIV,'*b','LineWidth',2);
    h(1) = plot(data(p).T,(data(p).Y-0.044)/0.033,'xr','LineWidth',2);
    
    h(4) = plot(data(p).T0,data(p).MEAN0+data(p).BIAS,'g','LineWidth',2);
    h(3) = plot(data(p).T0,data(p).MD0,'k','LineWidth',2);
    
    xlim([data(p).T0(1),data(p).T0(end)]);
    
    MIN = min(data(p).LB0);
    MIN = min(MIN,min(data(p).LB_INDIV));
    
    MAX = max(data(p).UB0);
    MAX = max(MAX,max(data(p).UB_INDIV));
    
    ylim([MIN-0.05*(MAX-MIN),MAX+0.05*(MAX-MIN)]);
    
    xlabel('time (ka)','FontSize',12);
    ylabel('SST (°C)','FontSize',12);
    
    legend(h,{'MULLER',CURVE,'GPST','SHAKUN'},'Location','EastOutside');
    
    
    h = zeros(3,1);
    subplot(3,1,3);
    hold on;
    title('Deviations of the Point Estimates','FontSize',16);
    
    h(4) = plot([data(p).T0(1),data(p).T0(end)],[0,0],'--g','LineWidth',2);
    
    h(2) = plot(data(p).T,data(p).MD_INDIV-data(p).MEAN-data(p).BIAS,'*b','LineWidth',2);
    h(1) = plot(data(p).T,(data(p).Y-0.044)/0.033-data(p).MEAN-data(p).BIAS,'xr','LineWidth',2);
    h(3) = plot(data(p).T,data(p).MD-data(p).MEAN-data(p).BIAS,'*k','LineWidth',2);
    
    % plot(data(p).T0,mean((data(p).Y-0.044)/0.033-data(p).MEAN-data(p).BIAS).*ones(size(data(p).T0)),'m','LineWidth',2);
    % plot(data(p).T0,mean(data(p).MD_INDIV-data(p).MEAN-data(p).BIAS).*ones(size(data(p).T0)),'c','LineWidth',2);
    
    xlim([data(p).T0(1),data(p).T0(end)]);
    
    xlabel('time (ka)','FontSize',12);
    ylabel('SST (°C)','FontSize',12);
    
    legend(h,{'MULLER',CURVE,'GPST','SHAKUN'},'Location','EastOutside');
    
    set(fig,'Position',[20 20 900 800]);
    movegui(fig,'center');
    drawnow;
    
    path = ['Outputs/',outputFile,'/',data(p).name,'.fig'];
    savefig(fig,path);
end


dd = lines(7);
cc = zeros(P,3);
for k = 1:3
    cc(:,k) = interp1(linspace(0,1,7)',dd(:,k),linspace(0,1,9)');
end
fig = figure;
subplot(4,1,1);
hold on;
title(['Deviations from the SST Prior (',CURVE,')'],'FontSize',16);
plot([data(1).T0(1),data(1).T0(end)],[0,0],'--k','LineWidth',2);
for p = 1:1
    plot(data(p).T,data(p).MD-data(p).MEAN-data(p).BIAS,'*','Color',cc(p,:),'LineWidth',2);
    plot(data(p).T,data(p).MD_INDIV-data(p).MEAN-data(p).BIAS,'x','Color',cc(p,:),'LineWidth',1);
    % IND = (data(p).T0>=data(p).T(1))&(data(p).T0<=data(p).T(end));
    % plot(data(p).T0,data(p).MD0-data(p).MEAN0-data(p).BIAS,':','Color',cc(p,:),'LineWidth',1);
    % h(p) = plot(data(p).T0(IND),data(p).MD0(IND)-data(p).MEAN0(IND)-data(p).BIAS,'Color',cc(p,:),'LineWidth',2);
end

xlim([data(1).T0(1),data(1).T0(end)]);
ylim([-6.5 6.5]);

xlabel('time (ka)','FontSize',12);
ylabel('SST (°C)','FontSize',12);

subplot(4,1,2);
hold on;
plot([data(1).T0(1),data(1).T0(end)],[0,0],'--k','LineWidth',2);
for p = 6:6
    plot(data(p).T,data(p).MD-data(p).MEAN-data(p).BIAS,'*','Color',cc(p,:),'LineWidth',2);
    plot(data(p).T,data(p).MD_INDIV-data(p).MEAN-data(p).BIAS,'x','Color',cc(p,:),'LineWidth',1);
    % IND = (data(p).T0>=data(p).T(1))&(data(p).T0<=data(p).T(end));
    % plot(data(p).T0,data(p).MD0-data(p).MEAN0-data(p).BIAS,':','Color',cc(p,:),'LineWidth',1);
    % h(p) = plot(data(p).T0(IND),data(p).MD0(IND)-data(p).MEAN0(IND)-data(p).BIAS,'Color',cc(p,:),'LineWidth',2);
end

xlim([data(1).T0(1),data(1).T0(end)]);
ylim([-4 4]);

xlabel('time (ka)','FontSize',12);
ylabel('SST (°C)','FontSize',12);

subplot(4,1,3);
hold on;
plot([data(1).T0(1),data(1).T0(end)],[0,0],'--k','LineWidth',2);
for p = P-1:P-1
    plot(data(p).T,data(p).MD-data(p).MEAN-data(p).BIAS,'*','Color',cc(p,:),'LineWidth',2);
    plot(data(p).T,data(p).MD_INDIV-data(p).MEAN-data(p).BIAS,'x','Color',cc(p,:),'LineWidth',1);
    % IND = (data(p).T0>=data(p).T(1))&(data(p).T0<=data(p).T(end));
    % plot(data(p).T0,data(p).MD0-data(p).MEAN0-data(p).BIAS,':','Color',cc(p,:),'LineWidth',1);
    % h(p) = plot(data(p).T0(IND),data(p).MD0(IND)-data(p).MEAN0(IND)-data(p).BIAS,'Color',cc(p,:),'LineWidth',2);
end

xlim([data(1).T0(1),data(1).T0(end)]);
ylim([-4 4]);

xlabel('time (ka)','FontSize',12);
ylabel('SST (°C)','FontSize',12);

subplot(4,1,4);
hold on;
plot([data(1).T0(1),data(1).T0(end)],[0,0],'--k','LineWidth',2);
for p = P:P
    plot(data(p).T,data(p).MD-data(p).MEAN-data(p).BIAS,'*','Color',cc(p,:),'LineWidth',2);
    plot(data(p).T,data(p).MD_INDIV-data(p).MEAN-data(p).BIAS,'x','Color',cc(p,:),'LineWidth',1);
    % IND = (data(p).T0>=data(p).T(1))&(data(p).T0<=data(p).T(end));
    % plot(data(p).T0,data(p).MD0-data(p).MEAN0-data(p).BIAS,':','Color',cc(p,:),'LineWidth',1);
    % h(p) = plot(data(p).T0(IND),data(p).MD0(IND)-data(p).MEAN0(IND)-data(p).BIAS,'Color',cc(p,:),'LineWidth',2);
end

xlim([data(1).T0(1),data(1).T0(end)]);
ylim([-2.5 2.5]);

% legend(h,TITLE,'Location','EastOutside');

xlabel('time (ka)','FontSize',12);
ylabel('SST (°C)','FontSize',12);

set(fig,'Position',[20 20 900 1000]);
movegui(fig,'center');
drawnow;

path = ['Outputs/',outputFile,'/Deviations.fig'];
savefig(fig,path);



end