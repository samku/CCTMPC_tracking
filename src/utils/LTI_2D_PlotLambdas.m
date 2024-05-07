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
fig.Position(3:4) = 3*[4,1.75]; % [width,height]

% plot data until reference change
l_MPC_mat = cell2mat(l_MPC);
plot(0:r_change-2, l_MPC_mat(:,1:r_change-1))
hold on;    grid on

% add labels and axes limit
xlabel("$t$",'Interpreter','latex');    ylabel("$\lambda$",'Interpreter','latex');
set(gca, 'TickLabelInterpreter', 'latex')
xlim([0,r_change-2])

% % minimize white borders around plot
set(gca,'LooseInset', max(get(gca,'TightInset'), 0.02)) % remove border from axis
set(fig,'PaperPositionMode','Auto','PaperUnits',...
    'centimeters','PaperSize',fig.Position(3:4)) % resize pdf page

fig.Renderer = 'painters'; % a way to force saving in vector graphics
% % save the plot as PDF file
saveas(fig, '../figures/LTI_2D_Lambdas.pdf','pdf')