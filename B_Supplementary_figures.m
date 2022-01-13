%% Supplementary Figures
% Script created by: Tamika Bassman, 01/12/2021

close all

%% Fig. S1

clear

load('./param/fig_colors.mat') 
load('./param/figS1_param.mat')

% Plot T_LB, T_UB
figure(in.fig_num); hold on;
F(1) = subplot(2,2,1); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);

% Read and plot target spectra
for gm = 1:length(in.fold_names)
    foldname = in.fold_names{gm};
    load(['./Ground_Motion/',foldname,'/GMSel_Param.mat']);
    
    plot(target_x,target_y,'LineStyle',in.line_types{gm},'Color',...
        in.color_rgbs{gm},'Marker',in.markers{gm},'linewidth',...
        in.line_widths(gm),'DisplayName',in.disp_names{gm},...
        'HandleVisibility',in.handles{gm});
    
    clearvars -except in rgb F
end

% One plot of all the target spectra
xlabel('Period (s)');
ylabel('Spectral Acceleration (g)');
lgd1 = legend('Location','northeast');
lgd1.FontSize = in.lgd_font_size;
xlim([0.1 10]);
ylim([0.05 4]);
xticks([0.1:0.1:1,2:1:10]);
xticklabels({'0.1','0.2','','','0.5','','','','','1','2','','','5','',...
    '','','','10'});
yticks([0.01:0.01:0.1,0.2:0.1:1.0,2.0:1.0:4.0]);
yticklabels({'0.01','0.02','','','0.05','','','','','0.1','0.2','','',...
    '0.5','','','','','1.0','2.0','','4.0'});
set(gca,'YScale','log','XScale','log');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Read UHS-based EDPs
fpath = ['./Analysis_Results/',in.UHSfile,'/'];
UHS_sdr    = mean(csvread([fpath,'sdr.csv']),2);
UHS_strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
UHS_pfa    = mean(csvread([fpath,'pfa.csv']),2);
max_UHS_sdr = max(UHS_sdr);
max_UHS_strshr = max(UHS_strshr);
max_UHS_pfa = max(UHS_pfa);

fpath = ['./Analysis_Results/',in.UHS75file,'/'];
UHS75_sdr    = mean(csvread([fpath,'sdr.csv']),2);
UHS75_strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
UHS75_pfa    = mean(csvread([fpath,'pfa.csv']),2);
max_UHS75_sdr = max(UHS75_sdr);
max_UHS75_strshr = max(UHS75_strshr);
max_UHS75_pfa = max(UHS75_pfa);

% Read CMS-based EDPs
all_sdrs = zeros(43,length(in.fold_names));
all_strshrs = zeros(43,length(in.fold_names));
all_pfas = zeros(43,length(in.fold_names));

max_sdrs = zeros(length(in.fold_names),1);
max_strshrs = zeros(length(in.fold_names),1);
max_pfas = zeros(length(in.fold_names),1);

for gm = 1:length(in.fold_names)
    foldname = in.fold_names{gm};
    fpath = ['./Analysis_Results/',foldname,'/'];
    
    sdr    = mean(csvread([fpath,'sdr.csv']),2);
    strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
    pfa    = mean(csvread([fpath,'pfa.csv']),2);
    
    all_sdrs(:,gm) = sdr;
    all_strshrs(:,gm) = strshr;
    all_pfas(:,gm) = pfa;
    
    max_sdrs(gm) = max(sdr);
    max_strshrs(gm) = max(strshr);
    max_pfas(gm) = max(pfa);
    
    clear sdr strshr pfa
end

% Read response spectra data and compute Sa-EDP correlations
all_Sa_pairs  = cell(length(in.Ts_corr),1);
all_IDR_pairs = cell(length(in.Ts_corr),1);
all_SHR_pairs = cell(length(in.Ts_corr),1);
all_PFA_pairs = cell(length(in.Ts_corr),1);
x_Sa  = 1; % counters for filling the above cells
x_EDP = 1;

% For each CMS-based target spectrum
for tgt = 1:length(in.fold_names)-2 % Omit the UHS-based spectra at the end
    foldname = in.fold_names{tgt};
    load(['./Ground_Motion/',foldname,'/GMSel_Param.mat']);
    
    % Get T,Sa values of the 40 scaled selected response spectra
    Ts = knownPer;
    Sas = SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2));
    
    % For each record associated with this target
    for gm = 1:size(Sas,1)
        % Get and store interpolated Sa values at Ts_corr
        Sas_corr = interp1(Ts,Sas(gm,:),in.Ts_corr);
        for t = 1:length(in.Ts_corr)
            all_Sa_pairs{t}(x_Sa) = Sas_corr(t);
        end
        x_Sa = x_Sa + 1;
    end
    
    fpath = ['./Analysis_Results/',foldname,'/'];
    
    % Get and store EDPs for this target
    sdr = csvread([fpath,'sdr.csv']);
    strshr = csvread([fpath,'strshr.csv'])*4.4482216;
    pfa = csvread([fpath,'pfa.csv']);
    sdrs    = max(sdr(:,1:40),[],1);
    strshrs = max(strshr(:,1:40),[],1);
    pfas    = max(pfa(:,1:40),[],1);
    
    for t = 1:length(in.Ts_corr)
        all_IDR_pairs{t}(x_EDP:x_EDP+39) = sdrs;
        all_SHR_pairs{t}(x_EDP:x_EDP+39) = strshrs;
        all_PFA_pairs{t}(x_EDP:x_EDP+39) = pfas;
    end
    
    x_EDP = x_EDP + 40;
end

% Compute and store correlations between Sa and EDP values
IDR_Sa_corrs = zeros(length(in.Ts_corr),1);
SHR_Sa_corrs = zeros(length(in.Ts_corr),1);
PFA_Sa_corrs = zeros(length(in.Ts_corr),1);
% At each period, there is a unique set of Sas, but identical set of EDPs
for t = 1:length(in.Ts_corr)
    IDR_Sa_corrs(t) = corr(log(all_Sa_pairs{t}(:)),...
        log(all_IDR_pairs{t}(:)));
    SHR_Sa_corrs(t) = corr(log(all_Sa_pairs{t}(:)),...
        log(all_SHR_pairs{t}(:)));
    PFA_Sa_corrs(t) = corr(log(all_Sa_pairs{t}(:)),...
        log(all_PFA_pairs{t}(:)));
end

% Plot T_LB, T_UB; UHS-based EDPs
figure(in.fig_num); hold on;
F(2) = subplot(2,2,2); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
plot([0 10],max_UHS_sdr*ones(1,2),':k','LineWidth',1.25,'DisplayName',...
    'UHS');
plot([0 10],max_UHS75_sdr*ones(1,2),':','LineWidth',1.25,'Color',...
    rgb.gray,'DisplayName','UHS75');

F(3) = subplot(2,2,3); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
plot([0 10],max_UHS_strshr*ones(1,2),':k','LineWidth',1.25,...
    'HandleVisibility','off');
plot([0 10],max_UHS75_strshr*ones(1,2),':','LineWidth',1.25,'Color',...
    rgb.gray,'HandleVisibility','off');

F(4) = subplot(2,2,4); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
plot([0 10],max_UHS_pfa*ones(1,2),':k','LineWidth',1.25,...
    'HandleVisibility','off');
plot([0 10],max_UHS75_pfa*ones(1,2),':','LineWidth',1.25,'Color',...
    rgb.gray,'HandleVisibility','off');

% Plot governing EDP value for each CMS-based target (i.e., for each T*)
for gm = 1:length(in.fold_names)-2
    subplot(2,2,2); hold on;
    plot(in.periods(gm),max_sdrs(gm),'Marker',in.markers2{gm},'LineWidth',...
        2,'MarkerSize',7.5,'MarkerFaceColor',in.color_rgbs{gm},...
        'MarkerEdgeColor',in.color_rgbs{gm},'HandleVisibility','off');
    subplot(2,2,3); hold on;
    plot(in.periods(gm),max_strshrs(gm),'Marker',in.markers2{gm},...
        'LineWidth',2,'MarkerSize',7.5,'MarkerFaceColor',...
        in.color_rgbs{gm},'MarkerEdgeColor',in.color_rgbs{gm},...
        'HandleVisibility','off');
    subplot(2,2,4); hold on;
    plot(in.periods(gm),max_pfas(gm),'Marker',in.markers2{gm},'LineWidth',...
        2,'MarkerSize',7.5,'MarkerFaceColor',in.color_rgbs{gm},...
        'MarkerEdgeColor',in.color_rgbs{gm},'HandleVisibility','off');
end

% Plot correlations and format plots
subplot(2,2,2); hold on;
plot(20,10,'s','MarkerFaceColor',rgb.gray,'MarkerEdgeColor',rgb.gray,...
    'DisplayName','CMS');
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it IDR}_{max}');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0.005 0.035]);
yticklabels({'0.010','0.020','0.030'});
lgd1 = legend('Location','Southeast');
lgd1.FontSize = in.lgd_font_size;
yyaxis right;
plot(in.Ts_corr,IDR_Sa_corrs,in.corr_line_type,'Color',in.corr_color,...
    'LineWidth',in.line_width,'DisplayName','\rho');
ylabel('\rho_{{\itIDR, Sa}({\itT})}');
ylim([0 1.5]);
yticklabels({'0','0.5','1.0',''});
ax = gca; ax.YColor = rgb.black;

subplot(2,2,3); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it V}_{base,max} (MN)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([10000 50000]);
yticklabels({'10','20','30','40','50'});
yyaxis right;
plot(in.Ts_corr,SHR_Sa_corrs,in.corr_line_type,'Color',in.corr_color,...
    'LineWidth',in.line_width);
ylabel('\rho_{{\itV}_{base}, {\itSa}({\itT})}');
ylim([0 1.5]);
yticklabels({'0','0.5','1.0',''});
ax = gca; ax.YColor = rgb.black;

subplot(2,2,4); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it RA}_{max} (g)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0.2 1.6]);
yticklabels({'0.5','1.0','1.5'});
yyaxis right;
plot(in.Ts_corr,PFA_Sa_corrs,in.corr_line_type,'Color',in.corr_color,...
    'LineWidth',in.line_width);
ylabel('\rho_{{\itRA, Sa}({\itT})}');
ylim([0 1.5]);
yticklabels({'0','0.5','1.0',''});
ax = gca; ax.YColor = rgb.black;

% Format sizing
set(gcf,'Position',[0,0,850,750])
set(F(2), 'position', [0.115, 0.63, 0.3, 0.34]);
set(F(3), 'position', [0.6, 0.63, 0.3, 0.34]);
set(F(4), 'position', [0.115, 0.13, 0.3, 0.34]);
set(F(1), 'position', [0.6, 0.13, 0.3, 0.34]);

annotation('textbox',[0.25,0.05,0,0],'string','(c)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.75,0.05,0,0],'string','(d)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.25,0.55,0,0],'string','(a)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.75,0.55,0,0],'string','(b)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')

% a = figure(in.fig_num);
% set(a,'Units','Inches');
% pos = get(a,'Position');
% set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(a,'figS1.pdf','-dpdf','-r0')
% savefig('figS1.fig');

%% Fig. S2

clear;

load('./param/fig_colors.mat') 
load('./param/figS2_param.mat')

% Read CMS-based EDPs
all_sdrs = zeros(43,length(in.fold_names));
all_strshrs = zeros(43,length(in.fold_names));
all_pfas = zeros(43,length(in.fold_names));

max_sdrs = zeros(length(in.fold_names),1);
max_strshrs = zeros(length(in.fold_names),1);
max_pfas = zeros(length(in.fold_names),1);

for gm = 1:length(in.fold_names)
    foldname = in.fold_names{gm};
    fpath = ['./Analysis_Results/',foldname,'/'];
    
    sdr    = mean(csvread([fpath,'sdr.csv']),2);
    strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
    pfa    = mean(csvread([fpath,'pfa.csv']),2);
    
    all_sdrs(:,gm) = sdr;
    all_strshrs(:,gm) = strshr;
    all_pfas(:,gm) = pfa;
    
    max_sdrs(gm) = max(sdr);
    max_strshrs(gm) = max(strshr);
    max_pfas(gm) = max(pfa);
    
    clear sdr strshr pfa
end

% Read response spectra data and compute Sa-EDP correlations
all_Sa_pairs  = cell(length(in.Ts_corr),1);
all_IDR_pairs = cell(length(in.Ts_corr),1);
all_SHR_pairs = cell(length(in.Ts_corr),1);
all_PFA_pairs = cell(length(in.Ts_corr),1);
x_Sa  = 1; % counters for filling the above cells
x_EDP = 1;

% For each CMS-based target spectrum
for tgt = 1:length(in.fold_names)-2 % Omit the UHS-based spectra at the end
    foldname = in.fold_names{tgt};
    load(['./Ground_Motion/',foldname,'/GMSel_Param.mat']);
    
    % Get T,Sa values of the 40 scaled selected response spectra
    Ts = knownPer;
    Sas = SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2));
    
    % For each record associated with this target
    for gm = 1:size(Sas,1) % N
        % Get and store interpolated Sa values at Ts_corr
        Sas_corr = interp1(Ts,Sas(gm,:),in.Ts_corr);
        for t = 1:length(in.Ts_corr)
            all_Sa_pairs{t}(x_Sa) = Sas_corr(t);
        end
        x_Sa = x_Sa + 1;
    end
    
    fpath = ['./Analysis_Results/',foldname,'/'];
    
    % Get and store EDPs for this target
    sdr = csvread([fpath,'sdr.csv']);
    strshr = csvread([fpath,'strshr.csv'])*4.4482216;
    pfa = csvread([fpath,'pfa.csv']);
    sdrs    = max(sdr(:,1:40),[],1);
    strshrs = max(strshr(:,1:40),[],1);
    pfas    = max(pfa(:,1:40),[],1);
    
    for t = 1:length(in.Ts_corr)
        all_IDR_pairs{t}(x_EDP:x_EDP+39) = sdrs;
        all_SHR_pairs{t}(x_EDP:x_EDP+39) = strshrs;
        all_PFA_pairs{t}(x_EDP:x_EDP+39) = pfas;
    end
    
    x_EDP = x_EDP + 40;
end

% Compute and store correlations between Sa and EDP values
IDR_Sa_corrs = zeros(length(in.Ts_corr),1);
SHR_Sa_corrs = zeros(length(in.Ts_corr),1);
PFA_Sa_corrs = zeros(length(in.Ts_corr),1);
% At each period, there is a unique set of Sas, but identical set of EDPs
for t = 1:length(in.Ts_corr)
    IDR_Sa_corrs(t) = corr(log(all_Sa_pairs{t}(:)),...
        log(all_IDR_pairs{t}(:)));
    SHR_Sa_corrs(t) = corr(log(all_Sa_pairs{t}(:)),...
        log(all_SHR_pairs{t}(:)));
    PFA_Sa_corrs(t) = corr(log(all_Sa_pairs{t}(:)),...
        log(all_PFA_pairs{t}(:)));
end

% Plot T_LB, T_UB
figure(in.fig_num); hold on;
F(2) = subplot(2,2,2); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);

F(3) = subplot(2,2,3); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);

F(4) = subplot(2,2,4); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);

% Plot governing EDP value for each CMS-based target (i.e., for each T*)
for gm = 1:length(in.fold_names)-2
    subplot(2,2,2); hold on;
    plot(in.periods(gm),max_sdrs(gm),'Marker',in.markers2{gm},'LineWidth',...
        2,'MarkerSize',7.5,'MarkerFaceColor',in.color_rgbs{gm},...
        'MarkerEdgeColor',in.color_rgbs{gm},'HandleVisibility','off');
    subplot(2,2,3); hold on;
    plot(in.periods(gm),max_strshrs(gm),'Marker',in.markers2{gm},...
        'LineWidth',2,'MarkerSize',7.5,'MarkerFaceColor',...
        in.color_rgbs{gm},'MarkerEdgeColor',in.color_rgbs{gm},...
        'HandleVisibility','off');
    subplot(2,2,4); hold on;
    plot(in.periods(gm),max_pfas(gm),'Marker',in.markers2{gm},'LineWidth',...
        2,'MarkerSize',7.5,'MarkerFaceColor',in.color_rgbs{gm},...
        'MarkerEdgeColor',in.color_rgbs{gm},'HandleVisibility','off');
end

% Format plots
subplot(2,2,2); hold on;
plot(20,10,'s','MarkerFaceColor',rgb.gray,'MarkerEdgeColor',rgb.gray,...
    'DisplayName','CMS');
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it IDR}_{max}');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0 0.002]);
yticklabels({'0','0.001','0.002'});
lgd1 = legend('Location','Northeast');
lgd1.FontSize = in.lgd_font_size;
yyaxis right;
plot(in.Ts_corr,IDR_Sa_corrs,in.corr_line_type,'Color',in.corr_color,...
    'LineWidth',in.line_width,'DisplayName','\rho');
ylabel('\rho_{{\itIDR, Sa}({\itT})}');
ylim([0 2]);
yticklabels({'0','1.0',''});
ax = gca; ax.YColor = rgb.black;

subplot(2,2,3); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it V}_{base,max} (MN)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0 8000]);
yticklabels({'0','','4','','8'});
yyaxis right;
plot(in.Ts_corr,SHR_Sa_corrs,in.corr_line_type,'Color',in.corr_color,...
    'LineWidth',in.line_width);
ylabel('\rho_{{\itV}_{base}, {\itSa}({\itT})}');
ylim([0 2]);
yticklabels({'0','1.0',''});
ax = gca; ax.YColor = rgb.black;

subplot(2,2,4); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it RA}_{max} (g)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0 0.2]);
yticklabels({'0','0.1','0.2'});
yyaxis right;
plot(in.Ts_corr,PFA_Sa_corrs,in.corr_line_type,'Color',in.corr_color,...
    'LineWidth',in.line_width);
ylabel('\rho_{{\itRA, Sa}({\itT})}');
ylim([0 2]);
yticklabels({'0','1.0',''});
ax = gca; ax.YColor = rgb.black;

% Format sizing
set(gcf,'Position',[0,0,850,750])
set(F(2), 'position', [0.12, 0.63, 0.3, 0.34]);
set(F(3), 'position', [0.59, 0.63, 0.3, 0.34]);
set(F(4), 'position', [0.12, 0.13, 0.3, 0.34]);

annotation('textbox',[0.25,0.05,0,0],'string','(c)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.25,0.55,0,0],'string','(a)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.75,0.55,0,0],'string','(b)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')

% a = figure(in.fig_num);
% set(a,'Units','Inches');
% pos = get(a,'Position');
% set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
% print(a,'figS2.pdf','-dpdf','-r0')
% savefig('figS2.fig');

%% Fig. S3

clear;

load('./param/fig_colors.mat') 
load('./param/figS3_param.mat')

% Read UHS-based EDPs
fpath = ['./Analysis_Results/',in.UHSfile,'/'];
UHS_sdr    = mean(csvread([fpath,'sdr.csv']),2);
UHS_strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
UHS_pfa    = mean(csvread([fpath,'pfa.csv']),2);
max_UHS_sdr = max(UHS_sdr);
max_UHS_strshr = max(UHS_strshr);
max_UHS_pfa = max(UHS_pfa);

fpath = ['./Analysis_Results/',in.UHS75file,'/'];
UHS75_sdr    = mean(csvread([fpath,'sdr.csv']),2);
UHS75_strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
UHS75_pfa    = mean(csvread([fpath,'pfa.csv']),2);
max_UHS75_sdr = max(UHS75_sdr);
max_UHS75_strshr = max(UHS75_strshr);
max_UHS75_pfa = max(UHS75_pfa);

fpath = ['./Analysis_Results/',in.UHS90file,'/'];
UHS90_sdr    = mean(csvread([fpath,'sdr.csv']),2);
UHS90_strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
UHS90_pfa    = mean(csvread([fpath,'pfa.csv']),2);
max_UHS90_sdr = max(UHS90_sdr);
max_UHS90_strshr = max(UHS90_strshr);
max_UHS90_pfa = max(UHS90_pfa);

% Read CMS-based EDPs
all_sdrs = zeros(43,length(in.fold_names));
all_strshrs = zeros(43,length(in.fold_names));
all_pfas = zeros(43,length(in.fold_names));

max_sdrs = zeros(length(in.fold_names),1);
max_strshrs = zeros(length(in.fold_names),1);
max_pfas = zeros(length(in.fold_names),1);

for gm = 1:length(in.fold_names)
    N = in.Ns(gm);
    foldname = in.fold_names{gm};
    fpath = ['./Analysis_Results/',foldname,'/'];
    
    sdr    = mean(csvread([fpath,'sdr.csv']),2);
    strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
    pfa    = mean(csvread([fpath,'pfa.csv']),2);
    
    all_sdrs(:,gm) = sdr;
    all_strshrs(:,gm) = strshr;
    all_pfas(:,gm) = pfa;
    
    max_sdrs(gm) = max(sdr);
    max_strshrs(gm) = max(strshr);
    max_pfas(gm) = max(pfa);
    
    clear sdr strshr pfa
end

% Read response spectra data and compute Sa-EDP correlations
all_Sa_pairs  = cell(length(in.Ts_corr),length(in.corr_fold_names));
all_IDR_pairs = cell(length(in.Ts_corr),length(in.corr_fold_names));
all_SHR_pairs = cell(length(in.Ts_corr),length(in.corr_fold_names));
all_PFA_pairs = cell(length(in.Ts_corr),length(in.corr_fold_names));

% For each set of target spectra
N = 40;
for gmset = 1:length(in.corr_fold_names)
    fold_names = in.corr_fold_names{gmset};
    x_Sa  = 1; % counters for filling the above cells
    x_EDP = 1;
    
    % For each CMS-based target spectrum
    for tgt = 1:length(fold_names)
        
        foldname = fold_names{tgt};
        load(['./Ground_Motion/',foldname,'/GMSel_Param.mat']);
        
        % get scaled response spectra associated with this target
        Ts = knownPer;
        Sas = SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2));
                
        % For each record associated with this target
        for gm = 1:size(Sas,1)
            % Get and store interpolated Sa values at Ts_corr
            Sas_corr = interp1(Ts,Sas(gm,:),in.Ts_corr);
            for t = 1:length(in.Ts_corr)
                all_Sa_pairs{t,gmset}(x_Sa) = Sas_corr(t);
            end
            x_Sa = x_Sa + 1;
        end
        
        fpath = ['./Analysis_Results/',foldname,'/'];

        % Get and store EDPs for this target
        sdr    = csvread([fpath,'sdr.csv']);
        strshr = csvread([fpath,'strshr.csv'])*4.4482216;
        pfa    = csvread([fpath,'pfa.csv']);
        sdrs    = max(sdr(:,1:N),[],1);
        strshrs = max(strshr(:,1:N),[],1);
        pfas    = max(pfa(:,1:N),[],1);
        
        for t = 1:length(in.Ts_corr)
            all_IDR_pairs{t,gmset}(x_EDP:x_EDP+(N-1)) = sdrs;
            all_SHR_pairs{t,gmset}(x_EDP:x_EDP+(N-1)) = strshrs;
            all_PFA_pairs{t,gmset}(x_EDP:x_EDP+(N-1)) = pfas;
        end
        
        x_EDP = x_EDP + N;
    end
end

% Compute and store correlations between Sa and EDP values
IDR_Sa_corrs = zeros(length(in.Ts_corr),length(in.corr_fold_names));
SHR_Sa_corrs = zeros(length(in.Ts_corr),length(in.corr_fold_names));
PFA_Sa_corrs = zeros(length(in.Ts_corr),length(in.corr_fold_names));
% At each period, there is a unique set of Sas, but identical set of EDPs
for s = 1:length(in.corr_fold_names)
    for t = 1:length(in.Ts_corr)
        IDR_Sa_corrs(t,s) = corr(log(all_Sa_pairs{t,s}(:)),...
            log(all_IDR_pairs{t,s}(:)));
        SHR_Sa_corrs(t,s) = corr(log(all_Sa_pairs{t,s}(:)),...
            log(all_SHR_pairs{t,s}(:)));
        PFA_Sa_corrs(t,s) = corr(log(all_Sa_pairs{t,s}(:)),...
            log(all_PFA_pairs{t,s}(:)));
    end
end

% Plot T_LB, T_UB; UHS-based EDPs
figure(in.fig_num); hold on;
F(1) = subplot(1,3,1); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
plot([0 10],max_UHS_sdr*ones(1,2),':k','LineWidth',1.25,...
    'HandleVisibility','off');
plot([0 10],max_UHS75_sdr*ones(1,2),':','LineWidth',1.25,'Color',...
    rgb.gray,'HandleVisibility','off');

F(2) = subplot(1,3,2); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
plot([0 10],max_UHS_strshr*ones(1,2),':k','LineWidth',1.25,...
    'DisplayName','UHS');
plot([0 10],max_UHS75_strshr*ones(1,2),':','LineWidth',1.25,'Color',...
    rgb.gray,'DisplayName','UHS75');

F(3) = subplot(1,3,3); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
plot([0 10],max_UHS_pfa*ones(1,2),':k','LineWidth',1.25,...
    'HandleVisibility','off');
plot([0 10],max_UHS75_pfa*ones(1,2),':','LineWidth',1.25,'Color',...
    rgb.gray,'HandleVisibility','off');

% Plot governing EDP value for each CMS-based target (i.e., for each T*)
for gm = 1:length(in.fold_names)-2
    subplot(1,3,1); hold on;
    plot(in.periods(gm),max_sdrs(gm),'Marker',in.markertypes2{gm},...
        'LineWidth',in.line_widths(gm),'MarkerSize',...
        in.markersizes(gm),'MarkerFaceColor',in.color_rgbs{gm},...
        'MarkerEdgeColor',in.color_rgbs2{gm},'HandleVisibility','off');
    subplot(1,3,2); hold on;
    plot(in.periods(gm),max_strshrs(gm),'Marker',in.markertypes2{gm},...
        'LineWidth',in.line_widths(gm),'MarkerSize',...
        in.markersizes(gm),'MarkerFaceColor',in.color_rgbs{gm},...
        'MarkerEdgeColor',in.color_rgbs2{gm},'HandleVisibility','off');
    subplot(1,3,3); hold on;
    plot(in.periods(gm),max_pfas(gm),'Marker',in.markertypes2{gm},...
        'LineWidth',in.line_widths(gm),'MarkerSize',...
        in.markersizes(gm),'MarkerFaceColor',in.color_rgbs{gm},...
        'MarkerEdgeColor',in.color_rgbs2{gm},'HandleVisibility','off');
end

% Plot correlations and format plots
subplot(1,3,1); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it IDR}_{max}');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0.005 0.03]);
yticklabels({'','0.01','','0.02','','0.03'});
yyaxis right; hold on;
for s = 1:length(in.corr_fold_names)
    plot(in.Ts_corr,IDR_Sa_corrs(:,s),in.corrlinetype{s},'Color',...
        in.corrcolor{s},'LineWidth',in.linewidth,'HandleVisibility','off');
end
ylabel('\rho_{{\itIDR, Sa}({\itT})}');
ylim([0 2]);
yticklabels({'0','0.5','1.0','',''});
ax = gca; ax.YColor = rgb.black;

subplot(1,3,2); hold on;
plot(20,10,'s','MarkerFaceColor',rgb.gray,'MarkerEdgeColor',rgb.gray,...
    'DisplayName','CMS, theoretical');
plot(20,10,'s','MarkerSize',8.5,'MarkerFaceColor','None',...
    'MarkerEdgeColor',rgb.brick,'LineWidth',1.5,'DisplayName',...
    'CMS, ASCE 7-22');
plot(20,10,'d','MarkerSize',7.5,'MarkerFaceColor','None',...
    'MarkerEdgeColor',rgb.brick,'LineWidth',1.5,'DisplayName',...
    'CMS, C&A 2014');
plot(20,10,'s','MarkerFaceColor','None','MarkerEdgeColor',rgb.black,...
    'LineWidth',0.5,'MarkerSize',5,'DisplayName',...
    'CMS, ASCE 7-22, sets of 11');
plot(20,10,'-','Color',rgb.gray,'LineWidth',1,'DisplayName',...
    '\rho, CMS, theoretical');
plot(20,10,'-','Color',rgb.brick,'LineWidth',1,'DisplayName',...
    '\rho, CMS, ASCE 7-22');
plot(20,10,'--','Color',rgb.brick,'LineWidth',1,'DisplayName',...
    '\rho, CMS, C&A 2014');
plot(20,10,'.','Color',rgb.medblue,'MarkerSize',15,'DisplayName',...
    '{\itT}* = \{1s,5s\}');
plot(20,10,'.','Color',rgb.brick,'MarkerSize',15,'DisplayName',...
    '{\itT}* = \{0.15s,9s\}');
plot(20,10,'.','Color',rgb.dandelion,'MarkerSize',15,'DisplayName',...
    '{\itT}* = \{3s,8s\}');
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it V}_{base,max} (MN)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([10000 45000]);
yticklabels({'','15','','25','','35','','45'});
lgd1 = legend('Location','SouthOutside','NumColumns',1);
lgd1.Position = [0.695,0.253,0.1,0.1];
lgd1.FontSize = in.lgd_font_size;
yyaxis right;
for s = 1:length(in.corr_fold_names)
    plot(in.Ts_corr,SHR_Sa_corrs(:,s),in.corrlinetype{s},'Color',...
        in.corrcolor{s},'LineWidth',in.linewidth,'HandleVisibility','off');
end
ylabel('\rho_{{\itV}_{base}, {\itSa}({\itT})}');
ylim([0 2]);
yticklabels({'0','0.5','1.0','',''});
ax = gca; ax.YColor = rgb.black;

subplot(1,3,3); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it RA}_{max} (g)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0.2 1.4]);
yticklabels({'0.2','','0.6','','1.0','','1.4'});
yyaxis right;
for s = 1:length(in.corr_fold_names)
    plot(in.Ts_corr,PFA_Sa_corrs(:,s),in.corrlinetype{s},'Color',...
        in.corrcolor{s},'LineWidth',in.linewidth,'HandleVisibility','off');
end
ylabel('\rho_{{\itRA, Sa}({\itT})}');
ylim([0 2]);
yticklabels({'0','0.5','1.0','',''});
ax = gca; ax.YColor = rgb.black;

% Format sizing
set(gcf,'Position',[0,0,850,750])
set(F(1), 'position', [0.105, 0.63, 0.3, 0.34]);
set(F(2), 'position', [0.59, 0.63, 0.3, 0.34]);
set(F(3), 'position', [0.105, 0.13, 0.3, 0.34]);

annotation('textbox',[0.25,0.05,0,0],'string','(c)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.25,0.55,0,0],'string','(a)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.734,0.55,0,0],'string','(b)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')

% a = figure(in.fig_num); 
% set(a,'Units','Inches');
% pos = get(a,'Position');
% set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
%     [pos(3), pos(4)])
% print(a,'figS3.pdf','-dpdf','-r0')
% savefig('figS3.fig');

%% Fig. S4-6

clear
figs = {'figS4','figS5','figS6'};

for I = 1:3
    
    load('./param/fig_colors.mat') 
    load(['./param/',figs{I},'_param.mat'])
    
    % Read UHS-based EDPs
    fpath = ['./Analysis_Results/',in.UHSfile,'/'];
    UHS_sdr    = mean(csvread([fpath,'sdr.csv']),2);
    UHS_strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
    UHS_pfa    = mean(csvread([fpath,'pfa.csv']),2);
    max_UHS_sdr = max(UHS_sdr);
    max_UHS_strshr = max(UHS_strshr);
    max_UHS_pfa = max(UHS_pfa);
    
    fpath = ['./Analysis_Results/',in.UHS90file,'/'];
    UHS90_sdr    = mean(csvread([fpath,'sdr.csv']),2);
    UHS90_strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
    UHS90_pfa    = mean(csvread([fpath,'pfa.csv']),2);
    max_UHS90_sdr = max(UHS90_sdr);
    max_UHS90_strshr = max(UHS90_strshr);
    max_UHS90_pfa = max(UHS90_pfa);
    
    fpath = ['./Analysis_Results/',in.UHS75file,'/'];
    UHS75_sdr    = mean(csvread([fpath,'sdr.csv']),2);
    UHS75_strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
    UHS75_pfa    = mean(csvread([fpath,'pfa.csv']),2);
    max_UHS75_sdr = max(UHS75_sdr);
    max_UHS75_strshr = max(UHS75_strshr);
    max_UHS75_pfa = max(UHS75_pfa);
    
    % Read CMS-based EDPs
    all_sdrs = zeros(43,length(in.fold_names));
    all_strshrs = zeros(43,length(in.fold_names));
    all_pfas = zeros(43,length(in.fold_names));
    
    max_sdrs = zeros(length(in.fold_names),1);
    max_strshrs = zeros(length(in.fold_names),1);
    max_pfas = zeros(length(in.fold_names),1);
    
    for gm = 1:length(in.fold_names)
        foldname = in.fold_names{gm};
        N = in.Ns(gm);
        fpath = ['./Analysis_Results/',foldname,'/'];
        
        sdr    = mean(csvread([fpath,'sdr.csv']),2);
        strshr = mean(csvread([fpath,'strshr.csv']),2)*4.4482216;
        pfa    = mean(csvread([fpath,'pfa.csv']),2);
        
        all_sdrs(:,gm) = sdr;
        all_strshrs(:,gm) = strshr;
        all_pfas(:,gm) = pfa;
        
        max_sdrs(gm) = max(sdr);
        max_strshrs(gm) = max(strshr);
        max_pfas(gm) = max(pfa);
        
        clear sdr strshr pfa
    end
    
    % Read response spectra data and compute Sa-EDP correlations
    all_Sa_pairs  = cell(length(in.Ts_corr),length(in.corr_fold_names));
    all_IDR_pairs = cell(length(in.Ts_corr),length(in.corr_fold_names));
    all_SHR_pairs = cell(length(in.Ts_corr),length(in.corr_fold_names));
    all_PFA_pairs = cell(length(in.Ts_corr),length(in.corr_fold_names));
    
    % For each set of target spectra
    for gmset = 1:length(in.corr_fold_names)
        N = 11;
        if gmset == 1
            N = 40;
        end
        c_fold_names = in.corr_fold_names{gmset};
        x_Sa  = 1; % counters for filling the above cells
        x_EDP = 1;
        
        % For each CMS-based target spectrum
        for tgt = 1:length(c_fold_names)
            
            foldname = c_fold_names{tgt};
            load(['./Ground_Motion/',foldname,'/GMSel_Param.mat']);
            
            % Get T,Sa values of the 40 scaled selected response spectra
            Ts = knownPer;
            Sas = SaKnown(IMs.recID,:).*repmat(IMs.scaleFac,1,size(SaKnown,2));
            
            % For each record associated with this target
            for gm = 1:size(Sas,1) % N
                % Get and store interpolated Sa values at Ts_corr
                Sas_corr = interp1(Ts,Sas(gm,:),in.Ts_corr);
                for t = 1:length(in.Ts_corr)
                    all_Sa_pairs{t,gmset}(x_Sa) = Sas_corr(t);
                end
                x_Sa = x_Sa + 1;
            end
            
            fpath = ['./Analysis_Results/',foldname,'/'];
            
            % Get and store EDPs for this target
            sdr = csvread([fpath,'sdr.csv']);
            strshr = csvread([fpath,'strshr.csv'])*4.4482216;
            pfa = csvread([fpath,'pfa.csv']);
            sdrs    = max(sdr(:,1:N),[],1);
            strshrs = max(strshr(:,1:N),[],1);
            pfas    = max(pfa(:,1:N),[],1);
            
            for t = 1:length(in.Ts_corr)
                all_IDR_pairs{t,gmset}(x_EDP:x_EDP+(N-1)) = sdrs;
                all_SHR_pairs{t,gmset}(x_EDP:x_EDP+(N-1)) = strshrs;
                all_PFA_pairs{t,gmset}(x_EDP:x_EDP+(N-1)) = pfas;
            end
            
            x_EDP = x_EDP + N;
        end
    end
    
    % Compute and store correlations between Sa and EDP values
    IDR_Sa_corrs = zeros(length(in.Ts_corr),length(in.corr_fold_names));
    SHR_Sa_corrs = zeros(length(in.Ts_corr),length(in.corr_fold_names));
    PFA_Sa_corrs = zeros(length(in.Ts_corr),length(in.corr_fold_names));
    % At each period, there is a unique set of Sas, but identical set of EDPs
    for s = 1:length(in.corr_fold_names)
        for t = 1:length(in.Ts_corr)
            IDR_Sa_corrs(t,s) = corr(log(all_Sa_pairs{t,s}(:)),...
                log(all_IDR_pairs{t,s}(:)));
            SHR_Sa_corrs(t,s) = corr(log(all_Sa_pairs{t,s}(:)),...
                log(all_SHR_pairs{t,s}(:)));
            PFA_Sa_corrs(t,s) = corr(log(all_Sa_pairs{t,s}(:)),...
                log(all_PFA_pairs{t,s}(:)));
        end
    end
    
    % Plot T_LB, T_UB; UHS-based EDPs
    figure(in.fig_num); hold on;
    F(2) = subplot(2,2,2); hold on;
    xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
        0.75);
    xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
        'LineWidth',0.75);
    plot([0 10],max_UHS_sdr*ones(1,2),':k','LineWidth',1.25,...
        'HandleVisibility','off');
    plot([0 10],max_UHS75_sdr*ones(1,2),':','LineWidth',1.25,'Color',...
        rgb.gray,'HandleVisibility','off');
    
    F(3) = subplot(2,2,3); hold on;
    xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
        0.75);
    xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
        'LineWidth',0.75);
    plot([0 10],max_UHS_strshr*ones(1,2),':k','LineWidth',1.25,...
        'DisplayName','UHS');
    plot([0 10],max_UHS75_strshr*ones(1,2),':','LineWidth',1.25,'Color',...
        rgb.gray,'DisplayName','UHS75');
    
    F(4) = subplot(2,2,4); hold on;
    xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
        0.75);
    xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
        'LineWidth',0.75);
    plot([0 10],max_UHS_pfa*ones(1,2),':k','LineWidth',1.25,...
        'HandleVisibility','off');
    plot([0 10],max_UHS75_pfa*ones(1,2),':','LineWidth',1.25,'Color',...
        rgb.gray,'HandleVisibility','off');
    
    % Plot governing EDP value for each CMS-based target (i.e., for each T*)
    for gm = 1:length(in.fold_names)-2
        subplot(2,2,2); hold on;
        plot(in.periods(gm),max_sdrs(gm),'Marker',in.markers2{gm},...
            'LineWidth',in.line_widths(gm),'MarkerSize',in.marker_sizes(gm),...
            'MarkerFaceColor',in.color_rgbs{gm},'MarkerEdgeColor',...
            in.color_rgbs2{gm},'HandleVisibility','off');
        subplot(2,2,3); hold on;
        plot(in.periods(gm),max_strshrs(gm),'Marker',in.markers2{gm},...
            'LineWidth',in.line_widths(gm),'MarkerSize',in.marker_sizes(gm),...
            'MarkerFaceColor',in.color_rgbs{gm},'MarkerEdgeColor',...
            in.color_rgbs2{gm},'HandleVisibility','off');
        subplot(2,2,4); hold on;
        plot(in.periods(gm),max_pfas(gm),'Marker',in.markers2{gm},...
            'LineWidth',in.line_widths(gm),'MarkerSize',in.marker_sizes(gm),...
            'MarkerFaceColor',in.color_rgbs{gm},'MarkerEdgeColor',...
            in.color_rgbs2{gm},'HandleVisibility','off');
    end
    
    % Format plots
    subplot(2,2,2); hold on;
    xticks(0:9); xlim([0 9.5]);
    xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it IDR}_{max}');
    set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
    ylim([0.005 0.03]);
    yticklabels({'0.01','0.02','0.03'});
    yyaxis right; hold on;
    for s = 1:length(in.corr_fold_names)
        plot(in.Ts_corr,IDR_Sa_corrs(:,s),in.corr_line_type{s},'Color',...
            in.corr_color{s},'LineWidth',in.line_width(s),...
            'HandleVisibility','off');
    end
    ylabel('\rho_{{\itIDR, Sa}({\itT})}');
    ylim([0 2]);
    yticklabels({'0','1.0',''});
    ax = gca; ax.YColor = rgb.black;
    
    subplot(2,2,3); hold on;
    xticks(0:9); xlim([0 9.5]);
    plot(20,10,'s','MarkerFaceColor',rgb.gray,'MarkerEdgeColor',...
        rgb.gray,'DisplayName','CMS, theoretical');
    plot(20,10,'s','MarkerSize',8.5,'MarkerFaceColor','None',...
        'MarkerEdgeColor',in.corr_color{end},'LineWidth',1.5,'DisplayName',...
        'CMS, ASCE 7-22');
    plot(20,10,'d','MarkerSize',7.5,'MarkerFaceColor','None',...
        'MarkerEdgeColor',in.corr_color{end},'LineWidth',1.5,'DisplayName',...
        'CMS, C&A 2014');
    plot(20,10,'s','MarkerFaceColor','None','MarkerEdgeColor',...
        rgb.black,'LineWidth',0.5,'MarkerSize',5,'DisplayName',...
        'CMS, ASCE 7-22, sets of 11');
    plot(20,10,'-','Color',rgb.gray,'LineWidth',1,'DisplayName',...
        '\rho, CMS, theoretical');
    plot(20,10,'-','Color',in.corr_color{end},'LineWidth',...
        in.line_width(end),'DisplayName','\rho, CMS, ASCE 7-22, sets of 11');
    xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it V}_{base,max} (MN)');
    set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
    ylim([10000 45000]);
    yticklabels({'10','20','30','40'});
    lgd1 = legend('Location','SouthOutside','NumColumns',1);
    lgd1.Position = [0.706,0.308,0.1,0.1];
    lgd1.FontSize = in.lgd_font_size;
    yyaxis right;
    for s = 1:length(in.corr_fold_names)
        plot(in.Ts_corr,SHR_Sa_corrs(:,s),in.corr_line_type{s},'Color',...
            in.corr_color{s},'LineWidth',in.line_width(s),...
            'HandleVisibility','off');
    end
    ylabel('\rho_{{\itV}_{base}, {\itSa}({\itT})}');
    ylim([0 2]);
    yticklabels({'0','1.0',''});
    ax = gca; ax.YColor = rgb.black;
    
    subplot(2,2,4); hold on;
    xticks(0:9); xlim([0 9.5]);
    xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it RA}_{max} (g)');
    set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
    ylim([0.2 1.4001]);
    yticklabels({'0.2','','0.6','','1.0','','1.4'});
    yyaxis right;
    for s = 1:length(in.corr_fold_names)
        plot(in.Ts_corr,PFA_Sa_corrs(:,s),in.corr_line_type{s},'Color',...
            in.corr_color{s},'LineWidth',in.line_width(s),...
            'HandleVisibility','off');
    end
    ylabel('\rho_{{\itRA, Sa}({\itT})}');
    ylim([0 2]);
    yticklabels({'0','1.0',''});
    ax = gca; ax.YColor = rgb.black;
    
    
    set(gcf,'Position',[0,0,850,750])
    
    set(F(2), 'position', [0.105, 0.63, 0.3, 0.34]);
    set(F(3), 'position', [0.59, 0.63, 0.3, 0.34]);
    set(F(4), 'position', [0.105, 0.13, 0.3, 0.34]);
    
    annotation('textbox',[0.25,0.05,0,0],'string','(c)','FontName',...
        in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
    annotation('textbox',[0.25,0.55,0,0],'string','(a)','FontName',...
        in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
    annotation('textbox',[0.734,0.55,0,0],'string','(b)','FontName',...
        in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
    
%     a = figure(in.fig_num);
%     set(a,'Units','Inches');
%     pos = get(a,'Position');
%     set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
%     print(a,[figs{I},'.pdf'],'-dpdf','-r0')
%     savefig([figs{I},'.fig']);
    
    clearvars -except figs
    
end