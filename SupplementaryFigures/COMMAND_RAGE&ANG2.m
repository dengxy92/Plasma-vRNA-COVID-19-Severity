clear all
format long
%pathname=fileparts('/Users/xiaoyandeng/Desktop/ViralLoadsFitting/monolix_modified/Plot_319/');
%par1 = readtable('dv_427.xlsx','Sheet','all');
par1 = readtable('viralelimination.xlsx');

data1 = readtable('all_mean.csv');
data2 = readtable('women_mean.csv');
data3 = readtable('men_mean.csv');
data4 = readtable('moderate_mean.csv');
data5= readtable('severe_mean.csv');
data6 = readtable('critical_mean.csv');
data7 = readtable('alive_mean.csv');
data8 = readtable('deceased_mean.csv');

rage1 = table2array(data1(:,6));
ang1 = table2array(data1(:,9));
time1 = table2array(data1(:,2));
virus1 = table2array(data1(:,3));
std1 = table2array(data1(:,4));
std1_r =  table2array(data1(:,7));
std1_a =  table2array(data1(:,10));


rage2 = table2array(data2(:,6));
ang2 = table2array(data2(:,9));
time2 = table2array(data2(:,2));
virus2 = table2array(data2(:,3));
std2 = table2array(data2(:,4));
std2_r = table2array(data2(:,7));
std2_a = table2array(data2(:,10));

rage3 = table2array(data3(:,6));
ang3 = table2array(data3(:,9));
time3 = table2array(data3(:,2));
virus3 = table2array(data3(:,3));
std3 = table2array(data3(:,4));
std3_r = table2array(data3(:,7));
std3_a = table2array(data3(:,10));

rage4 = table2array(data4(:,6));
ang4 = table2array(data4(:,9));
time4 = table2array(data4(:,2));
virus4 = table2array(data4(:,3));
std4 = table2array(data4(:,4));
std4_r = table2array(data4(:,7));
std4_a = table2array(data4(:,10));

rage5 = table2array(data5(:,6));
ang5 = table2array(data5(:,9));
time5 = table2array(data5(:,2));
virus5 = table2array(data5(:,3));
std5 = table2array(data5(:,4));
std5_r = table2array(data5(:,7));
std5_a = table2array(data5(:,10));

rage6 = table2array(data6(:,6));
ang6 = table2array(data6(:,9));
time6 = table2array(data6(:,2));
virus6 = table2array(data6(:,3));
std6 = table2array(data6(:,4));
std6_r = table2array(data6(:,7));
std6_a = table2array(data6(:,10));

rage7 = table2array(data7(:,6));
ang7 = table2array(data7(:,9));
time7 = table2array(data7(:,2));
virus7 = table2array(data7(:,3));
std7 = table2array(data7(:,4));
std7_r = table2array(data7(:,7));
std7_a = table2array(data7(:,10));

rage8 = table2array(data8(:,6));
ang8 = table2array(data8(:,9));
time8 = table2array(data8(:,2));
virus8 = table2array(data8(:,3));
std8 = table2array(data8(:,4));
std8_r = table2array(data8(:,7));
std8_r(4) = std8_r(4)-0.25;
std8_a = table2array(data8(:,10));



% rage
[R1 P1] = corr(rage1,virus1,'Type','Pearson');
[R2 P2] = corr(rage2,virus2,'Type','Pearson');
[R3 P3] = corr(rage3,virus3,'Type','Pearson');
[R4 P4]= corr(rage4,virus4,'Type','Pearson');
[R5 P5] = corr(rage5,virus5,'Type','Pearson');
[R6 P6] = corr(rage6,virus6,'Type','Pearson');
[R7 P7] = corr(rage7,virus7,'Type','Pearson');
[R8 P8] = corr(rage8,virus8,'Type','Pearson');

R = round([R1 R2 R3 R4 R5 R6 R7 R8],5);
P = round([P1 P2 P3 P4 P5 P6 P7 P8],5);

RP =[R' P'];


%%  Simulate viral loads and calculate mean and SE of dynamics

% Initial condition and parameter guesses --------------------------------

    p.T0 = 1.27;
    p.p = 420; %420;              % production rate of new virions (virions/cell/day)
    p.I0 = 0; % Initial amount of infectious virus
    p.d_I = 0.1; 
    p.t_inf = 0;
    p.bet = 0.18;
%all
% all
dv_patient = table2array(par1(1:8,2));
V0_patient = table2array(par1(1:8,1));
%p.V0 = 0.00603484006202506;
solutions_patient = zeros(8,1000);
for i = 1:8
        p.d_V = dv_patient(i);
        p.V0 = V0_patient(i);
        p.IC = [p.T0,p.I0,p.V0];

        [sol,p] = simulation_virus_model_with_delay_no_tinf(p,[0,31]);
        
        t = linspace(0,31,1000);
        curves = deval(sol,t,3);

        solutions_patient(i,:) = curves;

end

  

%
pathname=fileparts('/Users/dengxiaoyan/Desktop/ViralLoadsFittingCode/Paper Code_reversion/Matlab code/Plot/');
t = linspace(0,31,1000);
%%
labelR0 = string(R);
labelR = strcat ('r ='," ", labelR0);
labelP0 = string(P);
labelP = strcat ('p ='," ", labelP0);
labelN = string({'All','Female','Male','Moderate','Severe','Critical','Alive','Deceased'});
labelRAGE = strcat(labelN," ",'('," ",labelR,','," ",labelP," ",')');
%% rage severity
figure(1)
yyaxis left
plot(t, log10(solutions_patient(1,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time1, rage1, std1_r,'o','MarkerEdgeColor','#787878','MarkerFaceColor','#787878','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#787878');
hold off
ylim([1.5 4.5])
xlim([0 31])
xlabel('Time from symptom onset (days)')
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
%legend([h2,h, h1],{'Data','Standard Error','Prediction'},'Location','southeast')
ax = gca;
ax.YAxis(1).Color = [37,37,37]/255;
ax.YAxis(2).Color = '#787878';
box off
%title('Full cohort')
figfile2 = fullfile(pathname,'all_patients_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7.5 4.5]);
saveas(gcf,figfile2,'pdf')
%% female
figure(2)
yyaxis left
plot(t,log10(solutions_patient(2,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time2, rage2, std2_r,'o','MarkerEdgeColor','#B22222','MarkerFaceColor','#B22222','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#B22222');
hold off
ylim([1.5 4.5])
xlim([0 31])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#B22222';
box off
figfile2 = fullfile(pathname,'women_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7.5 4.5]);
saveas(gcf,figfile2,'pdf')

%% male
figure(3)
yyaxis left
plot(t,log10(solutions_patient(3,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time3, rage3, std3_r,'o','MarkerEdgeColor','#1874CD','MarkerFaceColor','#1874CD','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#1874CD');
hold off
ylim([1.5 4.5])
xlim([0 31])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#1874CD';
box off
figfile2 = fullfile(pathname,'men_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7.5 4.5]);
saveas(gcf,figfile2,'pdf')

%% moderate
figure(4)
yyaxis left
plot(t,log10(solutions_patient(4,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255)
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time4, rage4, std4_r,'o','MarkerEdgeColor','#1874CD','MarkerFaceColor','#1874CD','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#1874CD');
hold off
ylim([1.5 4.5])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color = [37,37,37]/255;
ax.YAxis(2).Color = '#1874CD';
xlim([0 31])
box off
figfile2 = fullfile(pathname,'moderate_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7.5 4.5]);
saveas(gcf,figfile2,'pdf')

%%
figure(5)
yyaxis left
coplot(t,log10(solutions_patient(5,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time5, rage5, std5_r,'o','MarkerEdgeColor','#CD9B1D','MarkerFaceColor','#CD9B1D','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#CD9B1D');
hold off
ylim([1.5 4.5])
xlim([0 31])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#CD9B1D';
box off
figfile2 = fullfile(pathname,'severe_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7.5 4.5]);
saveas(gcf,figfile2,'pdf')


%% critical
figure(6)
yyaxis left
plot(t,log10(solutions_patient(6,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time6, rage6, std6_r,'o','MarkerEdgeColor','#B22222','MarkerFaceColor','#B22222','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#B22222');
hold off
ylim([1 4.5])
xlim([0 31])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#B22222';
box off
figfile2 = fullfile(pathname,'critical_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7.5 4.5]);
saveas(gcf,figfile2,'pdf')


%% alive
figure(7)
yyaxis left
plot(t,log10(solutions_patient(7,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time7, rage7, std7_r,'d','MarkerEdgeColor','#1874CD','MarkerFaceColor','#1874CD','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#1874CD');
hold off
ylim([1.5 4.5])
xlim([0 31])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#1874CD';
box off
figfile2 = fullfile(pathname,'alive_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7.5 4.5]);
saveas(gcf,figfile2,'pdf')
%% deceased
figure(8)
yyaxis left
plot(t,log10(solutions_patient(8,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time8, rage8, std8_r,'d','MarkerEdgeColor','#B22222','MarkerFaceColor','#B22222','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#B22222');
hold off
ylim([0 4.5])
xlim([0 31])
ylabel('RAGE (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#B22222';
box off
figfile2 = fullfile(pathname,'deceased_rage');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca, 'LooseInset', get(gca,'TightInset'))
set(gcf, 'PaperPosition', [0 0 7.5 4.5]); 
set(gca,'FontSize',6)
set(gcf, 'PaperSize', [7.5 4.5]);
saveas(gcf,figfile2,'pdf')

%% ang02
% ang
[R1 P1] = corr(ang1,virus1,'Type','Pearson');
[R2 P2] = corr(ang2,virus2,'Type','Pearson');
[R3 P3] = corr(ang3,virus3,'Type','Pearson');
[R4 P4]= corr(ang4,virus4,'Type','Pearson');
[R5 P5] = corr(ang5,virus5,'Type','Pearson');
[R6 P6] = corr(ang6,virus6,'Type','Pearson');
[R7 P7] = corr(ang7,virus7,'Type','Pearson');
[R8 P8] = corr(ang8,virus8,'Type','Pearson');

R_a = round([R1 R2 R3 R4 R5 R6 R7 R8],5);
P_a = round([P1 P2 P3 P4 P5 P6 P7 P8],5);

RP_a =[R_a' P_a'];


%%
labelR0 = string(R);
labelR = strcat ('r ='," ", labelR0);
labelP0 = string(P);
labelP = strcat ('p ='," ", labelP0);
labelN = string({'All','Female','Male','Moderate','Severe','Critical','Alive','Deceased'});
labelang = strcat(labelN," ",'('," ",labelR,','," ",labelP," ",')');

figure(9)
set(gca, 'LooseInset', get(gca,'TightInset'))
%tiledlayout(4,2,'TileSpacing','Compact','Padding','Compact');
% ang severity
%nexttile
subplot(4,2,1)
yyaxis left
plot(t, log10(solutions_patient(1,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time1, rage1, std1_a,'o','MarkerEdgeColor','#787878','MarkerFaceColor','#787878','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#787878');
hold off
ylim([1.5 4.5])
xlim([0 31])
xlabel('Time from symptom onset (days)')
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
%legend([h2,h, h1],{'Data','Standard Error','Prediction'},'Location','southeast')
ax = gca;
ax.YAxis(1).Color = [37,37,37]/255;
ax.YAxis(2).Color = '#787878';
box off
set(gca,'FontSize',6)

% female
%nexttile
subplot(4,2,2)
yyaxis left
plot(t,log10(solutions_patient(2,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time2, rage2, std2_a,'o','MarkerEdgeColor','#B22222','MarkerFaceColor','#B22222','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#B22222');
hold off
ylim([1.5 4.5])
xlim([0 31])
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#B22222';
box off
set(gca,'FontSize',6)

subplot(4,2,3)
% male
%nexttile
yyaxis left
plot(t,log10(solutions_patient(3,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time3, rage3, std3_a,'o','MarkerEdgeColor','#1874CD','MarkerFaceColor','#1874CD','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#1874CD');
hold off
ylim([1.5 4.5])
xlim([0 31])
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#1874CD';
box off
set(gca,'FontSize',6)

subplot(4,2,4)
% moderate
%nexttile
yyaxis left
plot(t,log10(solutions_patient(4,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255)
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time4, rage4, std4_a,'^','MarkerEdgeColor','#1874CD','MarkerFaceColor','#1874CD','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#1874CD');
hold off
ylim([1.5 4.5])
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color = [37,37,37]/255;
ax.YAxis(2).Color = '#1874CD';
xlim([0 31])
box off
set(gca,'FontSize',6)

subplot(4,2,5)
%nexttile
yyaxis left
plot(t,log10(solutions_patient(5,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time5, rage5, std5_a,'^','MarkerEdgeColor','#CD9B1D','MarkerFaceColor','#CD9B1D','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#CD9B1D');
hold off
ylim([1.5 4.5])
xlim([0 31])
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#CD9B1D';
box off
set(gca,'FontSize',6)

subplot(4,2,6)
% critical
%nexttile
yyaxis left
plot(t,log10(solutions_patient(6,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time6, rage6, std6_a,'^','MarkerEdgeColor','#B22222','MarkerFaceColor','#B22222','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#B22222');
hold off
ylim([1 4.5])
xlim([0 31])
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#B22222';
box off
set(gca,'FontSize',6)


subplot(4,2,7)
% alive
%nexttile
yyaxis left
plot(t,log10(solutions_patient(7,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time7, rage7, std7_a,'d','MarkerEdgeColor','#1874CD','MarkerFaceColor','#1874CD','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#1874CD');
hold off
ylim([1.5 4.5])
xlim([0 31])
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#1874CD';
box off
set(gca,'FontSize',6)

subplot(4,2,8)
% deceased
%nexttile
yyaxis left
plot(t,log10(solutions_patient(8,:)),'-','LineWidth',0.5,'Color',[37,37,37]/255);
%h2.Annotation.%legendInformation.IconDisplayStyle = 'off';
hold on
ylim([0 3])
ylabel('Plasma vRNA load (log_{10}(copies/mL))')
yyaxis right
errorbar(time8, rage8, std8_a,'d','MarkerEdgeColor','#B22222','MarkerFaceColor','#B22222','MarkerSize',2,'Capsize',2,'LineWidth',0.25,'Color','#B22222');
hold off
ylim([0 4.5])
xlim([0 31])
ylabel('Ang-2 (log_{10}(pg/mL))','Rotation',270,'VerticalAlignment','bottom')
xlabel('Time from symptom onset (days)')
ax = gca;
ax.YAxis(1).Color =[37,37,37]/255;
ax.YAxis(2).Color ='#B22222';
box off
set(gca,'FontSize',6)


figfile2 = fullfile(pathname,'ang');
set(gcf, 'PaperPositionMode', 'manual'); 
set(gcf, 'PaperUnits', 'centimeters'); 
set(gca,'FontSize',6);
set(gcf, 'PaperPosition', [0 0 18 20]); 
set(gcf, 'PaperSize', [18 20]);
set(gca, 'LooseInset', get(gca,'TightInset'))
saveas(gcf,figfile2,'pdf')
