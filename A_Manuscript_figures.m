%% Manuscript Figures
% Script created by: Tamika Bassman, 01/12/2021

close all

%% Fig. 1

clear

load('./param/fig_colors.mat') 
load('./param/fig1_param.mat')
load('./param/fig1_median_spect.mat')

% Plot T_LB, T_UB
figure(in.fig_num); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);

% Read and plot target spectra
for gm = 1:length(in.fold_names)
    foldname = in.fold_names{gm};
    load(['./Ground_Motion/',foldname,'/gm_sel_param.mat']);
    
    figure(in.fig_num); hold on;
    plot(target_x,target_y,'LineStyle',in.line_types{gm},'Color',...
        in.color_rgbs{gm},'Marker',in.markers{gm},'linewidth',...
        in.line_widths(gm),'DisplayName',in.disp_names{gm},...
        'HandleVisibility',in.handles{gm});
    
    clearvars -except in rgb med
end

% Plot median spectrum from ground-motion model from T*=5s disagg data
plot(med.gmpe_T,med.gmpe_sa,'LineStyle','-','Color',rgb.black,...
    'linewidth',0.5,'DisplayName',['Predicted median spectrum' newline...
    'for CMS with {\itT*}=5s'],'HandleVisibility','on');

% Format plot
h = figure(in.fig_num); 
hold on;
xlabel('Period (s)');
ylabel('Spectral Acceleration (g)');
lgd1 = legend('Location','southwest');
lgd1.FontSize = 10;
xlim([0.1 10]);
set(gca,'YScale','log','XScale','log');
xticks([0.1:0.1:1,2:1:10]);
ylim([0.05 4]);
xticklabels({'0.1','0.2','','','0.5','','','','','1','2','','','5','',...
    '','','','10'});
yticks([0.01:0.01:0.1,0.2:0.1:1.0,2.0:1.0:4.0]);
yticklabels({'0.01','0.02','','','0.05','','','','','0.1','0.2','','',...
    '0.5','','','','','1.0','2.0','','4.0'});
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);

% Annotate additional concepts
% T*
annotation('textarrow',[0.52,0.52],[0.75,0.625],'string','{\itT*}=1s',...
    'FontName',in.font_name,'FontSize',in.font_size,...
    'HorizontalAlignment','center','HeadStyle','vback2','HeadLength',...
    in.ahead_size,'HeadWidth',in.ahead_size)
annotation('textarrow',[0.79,0.79],[0.435,0.3075],'string','{\itT*}=5s',...
    'FontName',in.font_name,'FontSize',in.font_size,...
    'HorizontalAlignment','center','HeadStyle','vback2','HeadLength',...
    in.ahead_size,'HeadWidth',in.ahead_size)

% T_LB, T_UB
annotation('textarrow',[0.43,0.39],[0.88,0.88],'string',...
    '{\itT}_{LB}=0.46s','FontName',in.font_name,'FontSize',...
    in.font_size,'HorizontalAlignment','center','HeadStyle','vback2',...
    'HeadLength',in.ahead_size,'HeadWidth',in.ahead_size)
annotation('textarrow',[0.835,0.875],[0.88,0.88],'string',...
    '{\itT}_{UB}=8.4s','FontName',in.font_name,'FontSize',in.font_size,...
    'HorizontalAlignment','center','HeadStyle','vback2','HeadLength',...
    in.ahead_size,'HeadWidth',in.ahead_size)

% Epsilon double arrows
annotation('doublearrow',[0.79,0.79],[0.145,0.295],'HeadStyle','vback2',...
    'Head1Length',in.ahead_size,'Head1Width',in.ahead_size,...
    'Head2Length',in.ahead_size,'Head2Width',in.ahead_size)

% Epsilon labels
annotation('textarrow',[0.7,0.79],[0.18,0.23],'string','({\itT*}=5s)',...
    'FontName',in.font_name,'FontSize',in.font_size,...
    'HorizontalAlignment','center','HeadStyle','vback3','HeadLength',...
    in.ahead_size,'HeadWidth',in.ahead_size)
annotation('textbox',[0.59 0.21 0 0],'string','\boldmath$\varepsilon$',...
    'FontName',in.font_name,'FontSize',in.font_size+1,...
    'HorizontalAlignment','center','interpreter','latex')

% % Save figure
% set(h,'Units','Inches');
% pos = get(h,'Position');
% set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
%     [pos(3), pos(4)]);
% print(h,'fig1.pdf','-dpdf','-r0');
% savefig('fig1.fig');

%% Fig. 3

clear

load('./param/fig_colors.mat') 
load('./param/fig3_param.mat')
load('./Ground_Motion/MCE_CMS_T1_v1/gm_sel_param.mat');

figure(in.fig_num); hold on
F(1) = subplot(1,4,1); hold on

% Plot T_LB, T_UB
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);

% Plot selected scaled spectra
plot(knownPer,sel_spectra,'Color',rgb.gray,'HandleVisibility','off');
plot(100,100,'Color',rgb.gray,'DisplayName','Scaled Selected Spectra');

% Plot target spectrum
plot(target_x,target_y,'LineStyle','-','Color',rgb.black,...
    'linewidth',2.5,'DisplayName','Target');

% Plot mean of selected scaled spectra
plot(meansel_x,meansel_y,'LineStyle','--','Color',rgb.dandelion,...
    'linewidth',1.5,'DisplayName','Mean of Selected Spectra');

% Format plot
xlabel('Period (s)');
ylabel('Spectral Acceleration (g)');
lgd1 = legend('Location','northeast');
lgd1.FontSize = in.lgd_font_size;
xlim([0.1 10]);
ylim([0.05 4]);
set(gca,'YScale','log','XScale','log');
xticks([0.1:0.1:1,2:1:10]);
xticklabels({'0.1','0.2','','','0.5','','','','','1','2','','','5','',...
    '','','','10'});
yticks([0.01:0.01:0.1,0.2:0.1:1.0,2.0:1.0:4.0]);
yticklabels({'0.01','0.02','','','0.05','','','','','0.1','0.2','','',...
    '0.5','','','','','1.0','2.0','','4.0'});
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);

% Read NLTHA results (EDPs) - drift (sdr), story shear (strshr), floor 
% acceleration (pfa)
foldname = in.fold_names{1};
fpath = ['./Analysis_Results/',foldname,'/'];

all_sdrs = csvread([fpath,'sdr.csv']);
all_strshrs = csvread([fpath,'strshr.csv'])*4.4482216;
all_pfas = csvread([fpath,'pfa.csv']);

sdr    = mean(all_sdrs,2);
strshr = mean(all_strshrs,2);
pfa    = mean(all_pfas,2);

max_sdr = max(sdr);
id_max_sdr = find(sdr==max_sdr,1);
max_strshr = max(strshr);
id_max_shr = find(strshr==max_strshr,1);
max_pfa = max(pfa);
id_max_pfa = find(pfa==max_pfa,1);

% Plot EDPs
% Drift
figure(in.fig_num);
hold on;
F(2) = subplot(1,4,2); hold on;
plot(100,100,'Color',rgb.gray,'DisplayName',['Individual' newline...
    'Analyses']);
plot(all_sdrs,1:43,'LineStyle',in.line_types{1},...
    'LineWidth',0.75,'HandleVisibility','off','Color',rgb.gray);
plot(sdr,1:43,'LineStyle',in.line_types{1},...
    'LineWidth',1.5,'DisplayName','Mean','Color',rgb.black);
plot(max_sdr,id_max_sdr,'Marker','s','MarkerFaceColor',rgb.dandelion,...
    'MarkerEdgeColor',rgb.dandelion,...
    'DisplayName','Max(Mean)','Color',[1 1 1]);
xticks(0:0.005:0.04); xlim([0 0.04]);
xticklabels({'0','','0.01','','0.02','','0.03','','0.04'});
yticks([0:5:40,43]); ylim([1 43]);
yticklabels({'0','5','10','15','20','25','30','35','40','43'});
xlabel('{\it IDR}'); ylabel('Story');
lgd2 = legend('Location','Southeast');
lgd2.FontSize = in.lgd_font_size;
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);

% Story shear
hold on;
F(3) = subplot(1,4,3); hold on;
plot(all_strshrs,1:43,'LineStyle',in.line_types{1},...
    'LineWidth',0.75,'HandleVisibility','off','Color',rgb.gray);
plot(strshr,1:43,'LineStyle',in.line_types{1},...
    'LineWidth',1.5,'DisplayName','Mean','Color',rgb.black);
plot(max_strshr,id_max_shr,'LineStyle',in.line_types{1},'Marker','s',...
    'MarkerFaceColor',rgb.dandelion,'DisplayName','Mean','Color',...
    rgb.dandelion);
xticks(5*(0:1000:10000));
xticklabels({'0','','10','','20','','30','','40','','50',''});
xlim([0 50000]);
yticks([0:5:40,43]); ylim([1 43]);
yticklabels({'0','5','10','15','20','25','30','35','40','43'});
xlabel('{\it V}_{story} (MN)'); % ylabel('Story');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);

% Floor acceleration
hold on;
F(4) = subplot(1,4,4); hold on;
plot(all_pfas,1:43,'LineStyle',in.line_types{1},...
    'LineWidth',0.75,'HandleVisibility','off','Color',rgb.gray);
plot(pfa,1:43,'LineStyle',in.line_types{1},...
    'LineWidth',1.5,'DisplayName','Mean','Color',rgb.black);
plot(max_pfa,id_max_pfa,'Marker','s','MarkerFaceColor',rgb.dandelion,...
    'DisplayName','Max(Mean)','Color',rgb.dandelion);
xticks(0:0.2:1.6); xlim([0 1.6]);
xticklabels({'0','','0.4','','0.8','','1.2','','1.6'});
yticks([0:5:40,43]); ylim([1 43]);
yticklabels({'0','5','10','15','20','25','30','35','40','43'});
xlabel('{\it FA} (g)'); % ylabel('Story');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);

% Format plot sizing
set(F(1), 'position', [0.065, 0.21, 0.3, 0.75]);
set(F(2), 'position', [0.425, 0.21, 0.15, 0.75]);
set(F(3), 'position', [0.625, 0.21, 0.15, 0.75]);
set(F(4), 'position', [0.825, 0.21, 0.15, 0.75]);

a = figure(in.fig_num); 
set(a,'Position',[0,0,1150,420])
annotation('textbox',[0.215,0.07,0,0],'string','(a)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.5,0.07,0,0],'string','(b)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.7,0.07,0,0],'string','(c)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.9,0.07,0,0],'string','(d)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')

% set(a,'Units','Inches');
% pos = get(a,'Position');
% set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
%     [pos(3), pos(4)]);
% print(a,'fig3.pdf','-dpdf','-r0');
% savefig('fig3.fig');

%% Fig. 4

clear;

load('./param/fig_colors.mat') 

% (a) Load theoretical procedure target spectra plot parameters
load('./param/fig4_1_param.mat')

figure(in.fig_num); hold on; 
F(1) = subplot(1,3,1); hold on;

% Plot T_LB, T_UB
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);

% Read and plot target spectra
for gm = 1:length(in.fold_names)
    foldname = in.fold_names{gm};
    load(['./Ground_Motion/',foldname,'/gm_sel_param.mat']);

    plot(target_x,target_y,'LineStyle',in.line_types{gm},'Color',...
        in.color_rgbs{gm},'Marker',in.markers{gm},'linewidth',...
        in.line_widths(gm),'DisplayName',in.disp_names{gm},...
        'HandleVisibility',in.handles{gm});
    
    clearvars -except in rgb F
end

% Format plot
xlabel('Period (s)');
ylabel('Spectral Acceleration (g)');
lgd1 = legend('Location','Northeast');
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

% (b) Load ASCE 7-22 procedure target spectra
load('./param/fig4_2_param.mat')

% Plot T_LB, T_UB
figure(in.fig_num); hold on; 
F(2) = subplot(1,3,2); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);

% Read and plot target spectra
for gm = 1:length(in.fold_names)
    foldname = in.fold_names{gm};
    load(['./Ground_Motion/',foldname,'/gm_sel_param.mat']);

    plot(target_x,target_y,'LineStyle',in.line_types{gm},'Color',...
        in.color_rgbs{gm},'Marker',in.markers{gm},'linewidth',...
        in.line_widths(gm),'DisplayName',in.disp_names{gm},...
        'HandleVisibility',in.handles{gm});
    
    clearvars -except in rgb F
end

% Format plot
xlabel('Period (s)');
ylabel('Spectral Acceleration (g)');
plot(20,10,'-','Color',rgb.medblue,'LineWidth',1,'DisplayName',...
    '{\itT}*=\{1s,5s\}');
plot(20,10,'-','Color',rgb.brick,'LineWidth',1,'DisplayName',...
    '{\itT}*=\{0.15s,9s\}');
plot(20,10,'-','Color',rgb.dandelion,'LineWidth',1,'DisplayName',...
    '{\itT}*=\{3s,8s\}');
plot(20,10,'--','Color',rgb.black,'LineWidth',0.75,'DisplayName','UHS75');
plot(20,10,'-','Color',rgb.black,'LineWidth',0.75,'DisplayName','UHS');
lgd2 = legend('Location','northeast');
lgd2.FontSize = in.lgd_font_size;
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


% (b) Load C&A-14 procedure target spectra
load('./param/fig4_3_param.mat');

% Plot T_LB, T_UB
figure(in.fig_num); hold on; 
F(3) = subplot(1,3,3); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off','LineWidth',...
    0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);

% Read and plot target spectra
for gm = 1:length(in.fold_names)
    foldname = in.fold_names{gm};
    load(['./Ground_Motion/',foldname,'/gm_sel_param.mat']);

    plot(target_x,target_y,'LineStyle',in.line_types{gm},'Color',...
        in.color_rgbs{gm},'Marker',in.markers{gm},'linewidth',...
        in.line_widths(gm),'HandleVisibility','off');
    
    clearvars -except in rgb F
end

% Format plot
xlabel('Period (s)');
ylabel('Spectral Acceleration (g)');
plot(20,10,'-','Color',rgb.medblue,'LineWidth',1,'DisplayName',...
    '{\itT}*=\{1s,5s\}');
plot(20,10,'-','Color',rgb.brick,'LineWidth',1,'DisplayName',...
    '{\itT}*=\{0.15s,9s\}');
plot(20,10,'-','Color',rgb.dandelion,'LineWidth',1,'DisplayName',...
    '{\itT}*=\{3s,8s\}');
plot(20,10,'-.','Color',rgb.black,'LineWidth',0.75,'DisplayName','UHS90');
plot(20,10,'-','Color',rgb.black,'LineWidth',0.75,'DisplayName','UHS');
lgd3 = legend('Location','northeast');
lgd3.FontSize = in.lgd_font_size;
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

% Format plot sizing
set(gcf,'Position',[0,0,1100,350])
set(F(1), 'Position', [0.08, 0.23, 0.23, 0.73]);
set(F(2), 'Position', [0.41, 0.23, 0.23, 0.73]);
set(F(3), 'Position', [0.74, 0.23, 0.23, 0.73]);

annotation('textbox',[0.196,0.085,0,0],'string','(a)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.526,0.085,0,0],'string','(b)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')
annotation('textbox',[0.857,0.085,0,0],'string','(c)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center')

% a = figure(in.fig_num); 
% set(a,'Units','Inches');
% pos = get(a,'Position');
% set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
%     [pos(3), pos(4)]);
% print(a,'fig4.pdf','-dpdf','-r0');
% savefig('fig4.fig');

%% Fig. 5

clear;

load('./param/fig_colors.mat') 
load('./param/fig5_param.mat')

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
    load(['./Ground_Motion/',foldname,'/gm_sel_param.mat']);
    
    % Get T,Sa values of the 40 scaled selected response spectra
    Ts = knownPer;
    Sas = sel_spectra;
        
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
F(1) = subplot(1,3,1); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
plot([0 10],max_UHS_sdr*ones(1,2),':k','LineWidth',1.25,...
    'DisplayName','UHS');
plot([0 10],max_UHS75_sdr*ones(1,2),':','LineWidth',1.25,...
    'Color',rgb.gray,'DisplayName','UHS75');

F(2) = subplot(1,3,2); hold on;
xline(0.463,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
xline(2.0*4.22,'--','Color',rgb.gray,'HandleVisibility','off',...
    'LineWidth',0.75);
plot([0 10],max_UHS_strshr*ones(1,2),':k','LineWidth',1.25,...
    'HandleVisibility','off');
plot([0 10],max_UHS75_strshr*ones(1,2),':','LineWidth',1.25,...
    'Color',rgb.gray,'HandleVisibility','off');

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
    plot(in.periods(gm),max_sdrs(gm),'Marker',in.markertypes{gm},...
        'LineWidth',2,'MarkerSize',7.5,'MarkerFaceColor',...
        in.color_rgbs{gm},'MarkerEdgeColor',in.color_rgbs{gm},...
        'HandleVisibility','off');
    subplot(1,3,2); hold on;
    plot(in.periods(gm),max_strshrs(gm),'Marker',in.markertypes{gm},...
        'LineWidth',2,'MarkerSize',7.5,'MarkerFaceColor',...
        in.color_rgbs{gm},'MarkerEdgeColor',in.color_rgbs{gm},...
        'HandleVisibility','off');
    subplot(1,3,3); hold on;
    plot(in.periods(gm),max_pfas(gm),'Marker',in.markertypes{gm},...
        'LineWidth',2,'MarkerSize',7.5,'MarkerFaceColor',...
        in.color_rgbs{gm},'MarkerEdgeColor',in.color_rgbs{gm},...
        'HandleVisibility','off');
end

% Plot correlations and format plots
subplot(1,3,1); hold on;
plot(20,10,'s','MarkerFaceColor',rgb.gray,'MarkerEdgeColor',rgb.gray,...
    'DisplayName','CMS');
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it IDR}_{max}');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0.005 0.03]);
yticklabels({'','0.010','','0.020','','0.030'});
lgd1 = legend('Location','northeast');
lgd1.FontSize = in.lgd_font_size;
yyaxis right;
plot(in.Ts_corr,IDR_Sa_corrs,'Color',rgb.gray,'LineWidth',...
    1,'DisplayName','\rho');
ylabel('\rho_{{\itIDR, Sa}({\itT})}');
ylim([0 1.5]);
yticklabels({'0','0.5','1.0',''});
ax = gca; ax.YColor = rgb.black;

subplot(1,3,2); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it V}_{base,max} (MN)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([10000 45000]);
yticklabels({'','15','','25','','35','','45'});
yyaxis right;
plot(in.Ts_corr,SHR_Sa_corrs,'Color',rgb.gray,'LineWidth',1);
ylabel('\rho_{{\itV}_{base}, {\itSa}({\itT})}');
ylim([0 1.5]);
yticklabels({'0','0.5','1.0',''});
ax = gca; ax.YColor = rgb.black;

subplot(1,3,3); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT} or {\itT}* (s)'); ylabel('{\it RA}_{max} (g)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0.2 1.4]);
yticklabels({'0.2','','0.6','','1.0','','1.4'});
yyaxis right;
plot(in.Ts_corr,PFA_Sa_corrs,'Color',rgb.gray,'LineWidth',1);
ylabel('\rho_{{\itRA, Sa}({\itT})}');
ylim([0 1.5]);
yticklabels({'0','0.5','1.0',''});
ax = gca; ax.YColor = rgb.black;

% Format sizing
set(gcf,'Position',[0,0,1300,350])
set(F(1), 'Position', [0.09, 0.25, 0.19, 0.7]);
set(F(2), 'Position', [0.405, 0.25, 0.19, 0.7]);
set(F(3), 'Position', [0.735, 0.25, 0.19, 0.7]);

annotation('textbox',[0.1875,0.085,0,0.01],'string','(a)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment',...
    'center','EdgeColor','none')
annotation('textbox',[0.503,0.085,0,0.01],'string','(b)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment',...
    'center','EdgeColor','none')
annotation('textbox',[0.833,0.085,0,0.01],'string','(c)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment',...
    'center','EdgeColor','none')

% a = figure(in.fig_num); 
% set(a,'Units','Inches');
% pos = get(a,'Position');
% set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
%     [pos(3), pos(4)]);
% print(a,'fig5.pdf','-dpdf','-r0');
% savefig('fig5.fig');

%% Fig. 6

clear;

load('./param/fig_colors.mat') 
load('./param/fig6_param.mat')

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
for gm = 1:length(in.fold_names)
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

% Format plots
subplot(1,3,1); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT}* (s)'); ylabel('{\it IDR}_{max}');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0.005 0.03]);
yticklabels({'','0.010','','0.020','','0.030'});

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
plot(20,10,'.','Color',rgb.medblue,'MarkerSize',15,'DisplayName',...
    '{\itT}*=\{1s,5s\}');
plot(20,10,'.','Color',rgb.brick,'MarkerSize',15,'DisplayName',...
    '{\itT}*=\{0.15s,9s\}');
plot(20,10,'.','Color',rgb.dandelion,'MarkerSize',15,'DisplayName',...
    '{\itT}*=\{3s,8s\}');
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT}* (s)'); ylabel('{\it V}_{base,max} (MN)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([10000 45000]);
yticklabels({'','15','','25','','35','','45'});
lgd1 = legend('Location','SouthOutside','NumColumns',4);
lgd1.Position = [0.399,0.06,0.28,0.15];
lgd1.FontSize = in.lgd_font_size;

subplot(1,3,3); hold on;
xticks(0:9); xlim([0 9.5]);
xlabel('{\itT}* (s)'); ylabel('{\it RA}_{max} (g)');
set(gca, 'FontName', in.font_name,'FontSize',in.font_size);
ylim([0.2 1.4]);
yticklabels({'0.2','','0.6','','1.0','','1.4'});

set(gcf,'Position',[0,0,1100,500])
set(F(1), 'position', [0.1, 0.41, 0.24, 0.54]);
set(F(2), 'position', [0.42, 0.41, 0.24, 0.54]);
set(F(3), 'position', [0.745, 0.41, 0.24, 0.54]);

annotation('textbox',[0.225,0.28,0,0.01],'string','(a)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center',...
    'EdgeColor','none')
annotation('textbox',[0.545,0.28,0,0.01],'string','(b)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center',...
    'EdgeColor','none')
annotation('textbox',[0.87,0.28,0,0.01],'string','(c)','FontName',...
    in.font_name,'FontSize',in.font_size,'HorizontalAlignment','center',...
    'EdgeColor','none')

% a = figure(in.fig_num); 
% set(a,'Units','Inches');
% pos = get(a,'Position');
% set(a,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',...
%     [pos(3), pos(4)]);
% print(a,'fig6.pdf','-dpdf','-r0');
% savefig('fig6.fig');
