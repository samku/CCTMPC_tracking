% NOTE: this script must be called from the main one. this separation has
% been done in order to clean up the main code (as the plotting part is
% less relevant compared to the core code. For this reason, a simple check
% is done to avoid running this script directly. Obviously we're not
% responsible for a different use (that it's up to you.)
db_check = dbstack;
if length(db_check) < 3 && ...
        ~any(strcmp(db_check(end).name,{'LateralDynamics','evaluateCode'}))
    error("Script was called directly. Execute the example one first.")
end

% % create figure
fig = figure("Renderer","painters","Units","centimeters");
fig.Position(3:4) = 3*[4,2]; % [width,height]
hold on; grid on

% draw a two-lane road
center = round(max(X.V(:,1))+min(X.V(:,1)),5)/2;
roadwidth = 3;
yline(-roadwidth+center,"LineWidth",2);
hold on
yline(roadwidth+center,"LineWidth",2);
yline(center,"LineWidth",1.25,"LineStyle","--");

% draw reference
stairs(longDyn.pos, r_sys, "Color",[0.85 0.325 0.098],"LineWidth",1)
% draw e_y
plot(longDyn.pos, x_sys(1,1:end-1),"Color",[0 0.447 0.741], "LineWidth",1)

% % draw closed-loop e_y projection tubes
% y_MPC_0 (e_y)
p_yMPC = fill([longDyn.pos,longDyn.pos(end:-1:1)], [min(projVerts.ey_y),max(projVerts.ey_y(:,end:-1:1))],...
    [0 0.447 0.741],'FaceAlpha',0.075,'EdgeColor',[0 0.447 0.741],"LineStyle",':');
% ys_MPC (e_y)
p_ysMPC = fill([longDyn.pos,longDyn.pos(end:-1:1)], [min(projVerts.ey_ys),max(projVerts.ey_ys(:,end:-1:1))],...
    [0.466 0.674 0.188],"FaceAlpha",0.075,"EdgeColor",[0.466 0.674 0.188],"LineStyle",':');
% y_rci (e_y)
p_rci = fill([longDyn.pos,longDyn.pos(end:-1:1)], [min(projVerts.ey_rci),max(projVerts.ey_rci(:,end:-1:1))],...
    [0.85 0.325 0.098],"FaceAlpha",0.125,"EdgeColor",[0.85 0.325 0.098],"LineStyle",':');

% % draw admissible and current facing e_psi
epsi_init = 1; epsi_inter = 15;
% get position
x_quiv = longDyn.pos(epsi_init:epsi_inter:end-1);
y_quiv = x_sys(1,epsi_init:epsi_inter:end-2);

% % e_psi ccPoly projection
% up
u_quiv = cos(max(projVerts.epsi_y(:,epsi_init:epsi_inter:end-1),[],1));
v_quiv = sin(max(projVerts.epsi_y(:,epsi_init:epsi_inter:end-1),[],1));
h_up = quiver(x_quiv,y_quiv,u_quiv,v_quiv, 'Color',[0.466 0.674 0.188]);
h_up.Head.LineStyle = 'solid';
set(h_up,'AutoScale','on', 'AutoScaleFactor',0.25, 'LineWidth',0.75)
% down
u_quiv = cos(min(projVerts.epsi_y(:,epsi_init:epsi_inter:end-1),[],1));
v_quiv = sin(min(projVerts.epsi_y(:,epsi_init:epsi_inter:end-1),[],1));
h_down = quiver(x_quiv,y_quiv,u_quiv,v_quiv, 'Color',[0.466 0.674 0.188]);
h_down.Head.LineStyle = 'solid';
set(h_down,'AutoScale','on', 'AutoScaleFactor',0.25, 'LineWidth',0.75)
% current e_psi
th_quiv = x_sys(3,epsi_init:epsi_inter:end-2);
u_quiv = cos(th_quiv); v_quiv = sin(th_quiv);
h_curr = quiver(x_quiv,y_quiv,u_quiv,v_quiv, 'Color',[0.929 0.694 0.125]);
h_curr.Head.LineStyle = 'solid';
set(h_curr,'AutoScale','on', 'AutoScaleFactor',0.25, 'LineWidth',0.75)

% % add labels and axes limit
xlabel("$x\ [m]$",'Interpreter','latex'); ylabel("$y\ [m]$",'Interpreter','latex')
xlim([0, max(longDyn.pos)]); ylim([-roadwidth+center, roadwidth+center])

% % bring grid in front of everything
Ax1 = gca;
Ax1.YGrid = 'on'; Ax1.Layer = 'top'; Ax1.GridAlpha = 0.05;
set(Ax1, 'TickLabelInterpreter', 'latex');

% % add legend inside the plot
lgd_han = legend([p_yMPC,p_ysMPC,p_rci,h_up,h_curr],...
    '$\mathrm{proj}_{e_y}(X(y_0^*(x_t,r_t))$',...
    '$\mathrm{proj}_{e_y}(X(y_{\mathrm{s}}^*(x_t,r_t))$',...
    '$\mathrm{proj}_{e_y}(X(y_{\mathrm{o}}(r_t)))$',...
    '$\mathrm{proj}_{e_\psi}(X(y_0^*(x_t,r_t))$','$e_\psi(t)$',...
    'Interpreter','latex','Location','southeast','FontSize',10,'LineWidth',0.02);
lgd_han.Position(1:2) = [0.6,0.1925];

% % minimize white borders around plot
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.01)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/LateralDynamics.pdf','pdf')