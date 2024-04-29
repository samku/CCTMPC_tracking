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
fig.Position(3:4) = 3*[2,1.75]; % [width,height]


% impose a solver tolerance
solv_tol = 1e-5;
% minor cosmetic edits
Lyap_cost = max(Lyap_cost, solv_tol/10);
Lyap_cost(Lyap_cost < solv_tol*2) = solv_tol/10;

% plot Lyapunov cost
semilogy(0:N_mpc-1, Lyap_cost,'.-','MarkerSize',7, 'Color',[0 0.447 0.741],'LineWidth',0.5);
hold on

% mask points under solver tolerance
rectangle('Position',[0,1e-15,N_mpc,solv_tol-1e-15],'FaceColor',[0.7,0.7,0.7],'LineStyle','none')
grid on
xl = xline(r_change-1,'--');
yl = yline(solv_tol);

% % add labels and axes limit
xlabel("$t$",'Interpreter','latex');    set(gca, 'TickLabelInterpreter', 'latex')
xlim([0,N_mpc-1]);  ylim([solv_tol/1.5, max(Lyap_cost)*1.5])

% % minimize white borders around plot
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/Lyapunov_2D_LTI.pdf','pdf')