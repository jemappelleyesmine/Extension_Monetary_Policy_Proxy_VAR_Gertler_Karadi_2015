function no_fig=plot_figure_sep_save(VAR,VARChol,VARbs,VARCholbs,nCol,nRow,M,switch_extern)

fig=figure(M);
set(gcf,'DefaultAxesFontSize',VAR.fontsize);
set(gcf,'DefaultTextFontSize',VAR.fontsize);
set(fig,'Position',[100 100 1000 1200]);  % Larger figure for 4 rows
M_new=M;
xmin    =   1;
xmax    =   VAR.irhor;
ymin    =   NaN(1,VAR.n);
ymax    =   NaN(1,VAR.n);
for nvar=1:VAR.n
    %Create axes
    ymin_temp =  min(min([VARCholbs.irsL(:,nvar); VARbs.irsL(:,VARChol.chol_order(nvar))]));
    ymax_temp =  max(max([VARCholbs.irsH(:,nvar); VARbs.irsH(:,VARChol.chol_order(nvar))]));
    ydiff = abs(ymax_temp-ymin_temp);
    ymin(1,VARChol.chol_order(nvar)) =  floor(10*(ymin_temp-ydiff*0.025))/10;
    ymax(1,VARChol.chol_order(nvar)) =  ceil(10*(ymax_temp+ydiff*0.025))/10;
end;

for nvar=1:VAR.n
    
    subplot(nRow,nCol,2*nvar-1); plot(1:VAR.irhor,VAR.irs(:,nvar),'LineWidth',2,'Color','k'); title(VAR.select_vars_label(1,nvar));  axis manual; axis([xmin xmax ymin(1,nvar) ymax(1,nvar)]); ylabel('%');
    hold on
    subplot(nRow,nCol,2*nvar-1); plot(1:VAR.irhor,VARbs.irsH(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    hold on
    subplot(nRow,nCol,2*nvar-1); plot(VARbs.irsL(:,nvar),'LineWidth',1,'Color','k','LineStyle','--');
    set(gca,'YGrid','off','XGrid','on');
    grid
    if nvar==1
        legend('External Instruments');
    end;
    
end;
for nvar=1:VAR.n
    subplot(nRow,nCol,2*VARChol.chol_order(nvar)); plot(1:VARChol.irhor,VARChol.irs(:,nvar),'LineWidth',2,'Color','r');  title(VAR.select_vars_label(1,VARChol.chol_order(nvar))); axis manual; axis([xmin xmax ymin(1,VARChol.chol_order(nvar)) ymax(1,VARChol.chol_order(nvar))]);  set(gcf, 'Color', 'w'); ylabel('%');
    hold on
    subplot(nRow,nCol,2*VARChol.chol_order(nvar)); plot(1:VARChol.irhor,VARCholbs.irsH(:,nvar),'LineWidth',1,'Color','r','LineStyle','--');
    hold on
    subplot(nRow,nCol,2*VARChol.chol_order(nvar)); plot(VARCholbs.irsL(:,nvar),'LineWidth',1,'Color','r','LineStyle','--');
    hold on

    set(gca,'YGrid','off','XGrid','on');
    grid
    if VARChol.chol_order(nvar)==1
        legend('Cholesky');
    end;
end;

% Add statistics at bottom
h=axes('Position',[0,0,1,1],'Xlim',[0 1],'Ylim',[0 1]);
set(h,'Visible','off');
string_F_m=['First stage regression: F: ' num2str(VAR.F_m(1),'%2.2f')];
text('Position',[0.1 0.02],'string',string_F_m, 'FontSize', VAR.fontsize);
string_F_m_rob=['robust F: ' num2str(VAR.F_m_rob(1),'%2.2f')];
text('Position',[0.35 0.02],'string',string_F_m_rob, 'FontSize', VAR.fontsize);
string_R2_m = ['R2: ' num2str(VAR.R2_m(1)*100,'%1.2f') '%'];
text('Position',[0.55 0.02],'string',string_R2_m, 'FontSize', VAR.fontsize);
string_R2adj_m = ['Adjusted R2: ' num2str(VAR.R2adj_m(1)*100,'%1.2f') '%'];
text('Position',[0.7 0.02],'string',string_R2adj_m, 'FontSize', VAR.fontsize);

% Save the figure BEFORE closing
if isfield(VAR, 'figure_name')
    % Save as PDF
    pdf_file = [VAR.figure_name '.pdf'];
    print(pdf_file, '-dpdf', '-fillpage');
    
    % Also save as PNG for easy viewing
    png_file = [VAR.figure_name '.png'];
    print(png_file, '-dpng', '-r300');
    
    % Check if files were created
    if exist(pdf_file, 'file')
        fprintf('    ✓ Saved PDF: %s\n', pdf_file);
    else
        fprintf('    ✗ Failed to save PDF: %s\n', pdf_file);
    end
    
    if exist(png_file, 'file')
        fprintf('    ✓ Saved PNG: %s\n', png_file);
    else
        fprintf('    ✗ Failed to save PNG: %s\n', png_file);
    end
end

% Now close
close(fig);
M_new=M_new-1;

no_fig=M_new-M+1;
