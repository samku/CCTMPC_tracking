% NOTE: this script must be called from the main one. this separation has
% been done in order to clean up the main code (as the plotting part is
% less relevant compared to the core code. For this reason, a simple check
% is done to avoid running this script directly. Obviously we're not
% responsible for a different use (that it's up to you.)
db_check = dbstack;
if length(db_check) < 3 && ...
        ~any(strcmp(db_check(end).name,{'LTI_2D','evaluateCode'}))
    error("Script was called directly. Execute the example one first.")
end

% % create figure
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 4*[4,2]; % [width,height]

% % plot feasible regions
han_CCTMPC = feasRegionCCTMPC.plot('alpha',0.25,'Linewidth',0.25,'color',[0.466 0.674 0.188],'EdgeColor',[0.466 0.674 0.188]);
hold on
feasRegionETMPC.plot('alpha',1,'Linewidth',0.1,'color',[1 1 1]); % white background
han_ETMPC = feasRegionETMPC.plot('alpha',0.25,'Linewidth',0.25,'color',[0.7 0.4667 0.35],'EdgeColor',[0.7 0.4667 0.35]);
feasRegionRTMPC.plot('alpha',1,'Linewidth',0.1,'color',[1 1 1]); % white background
han_RTMPC = feasRegionRTMPC.plot('alpha',0.25,'Linewidth',0.25,'color',[0.256 0.633 0.8],'EdgeColor',[0.256 0.633 0.8]);

% % plot state constraints
rectangle('Position',[min(X.V(:,1))-2, max(X.V(:,2)),20,5],'FaceColor',[0.7,0.7,0.7],'LineStyle','none')
rectangle('Position',[max(X.V(:,1)), min(X.V(:,2))-2,5,10],'FaceColor',[0.7,0.7,0.7],'LineStyle','none')
rectangle('Position',[min(X.V(:,1))-2, min(X.V(:,2))-5,20,5],'FaceColor',[0.7,0.7,0.7],'LineStyle','none');
rectangle('Position',[min(X.V(:,1))-2, min(X.V(:,2))-5,2,10],'FaceColor',[0.7,0.7,0.7],'LineStyle','none');
X.plot('alpha',0,'LineWidth',1.5,'EdgeColor',[0.5,0.5,0.5],'xdimension',2)

% % plot closed-loop CCPolytopes
for t = 1:r_change-1
    Polyhedron(ccPoly.F,OCP_y{t}(:,1)).plot('alpha',1,'Linewidth',0.1,'color',[1,1,1]); % white background
    han_F0_1stRef = Polyhedron(ccPoly.F,OCP_y{t}(:,1)).plot('alpha',0.7,'Linewidth',0.1,'color',[0.929 0.694 0.125]);
end
for t = r_change:N_mpc
    Polyhedron(ccPoly.F,OCP_y{t}(:,1)).plot('alpha',1,'Linewidth',0.1,'color',[1,1,1]);
    han_F0_2ndRef = Polyhedron(ccPoly.F,OCP_y{t}(:,1)).plot('alpha',0.7,'Linewidth',0.1,'color',[0.85 0.325 0.098]);
end

% % plot state evolution
% first reference
plot(x_sys(1,1:r_change)',x_sys(2,1:r_change)','color','black','LineWidth',1)
for t = 1:r_change
    scatter(x_sys(1,t)',x_sys(2,t)',12,'filled',"MarkerFaceColor",[0,0,0]) % black contour
    scatter(x_sys(1,t)',x_sys(2,t)',10,'filled',"MarkerFaceColor",[0 0.447 0.741])
end
% second reference
plot(x_sys(1,r_change:end)',x_sys(2,r_change:end)','--','color','black','LineWidth',1)
for t = r_change:N_mpc
    scatter(x_sys(1,t)',x_sys(2,t)',12,'filled',"MarkerFaceColor",[0,0,0]) % black contour
    scatter(x_sys(1,t)',x_sys(2,t)',10,'filled',"MarkerFaceColor",[0.635 0.078 0.184])
end
%initial condition
scatter(x_sys(1,1),x_sys(2,1),30, "filled","MarkerFaceColor",[0,0,0])
scatter(x_sys(1,1),x_sys(2,1),20, "filled","MarkerFaceColor",[0 0.447 0.741])
text(x_sys(1,1)-0.4,x_sys(2,1)-0.15,'$x_0$','Interpreter','latex','FontSize',12)
% second starting point
scatter(x_sys(1,r_change)',x_sys(2,r_change)',25,'filled',"MarkerFaceColor",[0 0 0])
scatter(x_sys(1,r_change)',x_sys(2,r_change)',20,'filled',"MarkerFaceColor",[0.635 0.078 0.184])

% % add labels and axes limit
xlabel("$x_1$",'Interpreter','latex'); ylabel("$x_2$",'Interpreter','latex')
xlim([min(X.V(:,1))-0.2, max(X.V(:,1)+0.2)]); ylim([min(X.V(:,2))-0.2, max(X.V(:,2))+0.2])

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

% % add first legend inside the plot for CCPolytopes
han_leg1 = legend(Ax1,[han_F0_1stRef, han_F0_2ndRef],{'$X(y_0^*(x_t,5))$','$X(y_0^*(x_t,\textrm{-}5))$'}, ...
    'Interpreter','latex','Location','northeast');
han_leg1.Box = 'off'; han_leg1.FontSize = 12;
han_leg1.Position = han_leg1.Position + [-0.01, -0.01, 0,0];

% % add second legend outside the plot for feasible regions comparison
Ax2 = copyobj(Ax1,gcf); % copy graphical object
delete( get(Ax2,'Children') )   % delete its children
Ax2.Visible = 'off';
han_leg2 = legend(Ax2,[han_CCTMPC, han_ETMPC, han_RTMPC],{'CCTMPC','ETMPC','RTMPC'}, ...
    'Interpreter','latex','Location','northoutside','Orientation','horizontal');
han_leg2.Box = 'off'; han_leg2.FontSize = 12;
han_leg2.Position = han_leg2.Position + [0, 0.06, 0,0];

% % minimize white borders around plot
set(Ax2,'LooseInset', max(get(Ax2,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/LTI_2D.pdf','pdf')