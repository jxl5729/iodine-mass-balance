%%   Box model for Iodine 
% ---- The unit of time is **Ma** ------
% ---- The unit of amount is **mol** -----
clear all

% -------  SET PARAMETERS -----------
M_ocean = 1.386E21;% kg or liter
%---------INPUT-----------------------
F_water = 3.15E9*1.0E6; %10 ug/L; input_river(freshwater)mol/Ma;lowest value from Milliman et al.1995 and Moran et al. 2002
F_riversed = 6.3E6*1.0E6;%0.05 ppm; input_river(sediment);mol/Ma;lowest value from Lu et al. 2010 and Han et al. 2023. 
F_volc = 7.87E6*1.0E6; %1E3 ton; input_volcanic_emission;mol/Ma;lowest value from Witt et al. 2008 and Calting and Kasting 2017
%----------Input_total = 3.16417E9*1E6 mol/Ma-----------

tmax     = 2000;      % final time in Ma  
maxstep  = 0.1;        % maximum step for integration (Ma)

% ----------- Initial conditions ------------
t0=0;                 % initial time  

I10=5.0E-7; % mol/L;initial amounts of iodine
I20=5.0E-7; % mol/L;initial amounts of iodine
I30=5.0E-7; % mol/L;initial amounts of iodine
I40=5.0E-7; % mol/L;initial amounts of iodine
I50=5.0E-7; % mol/L;initial amounts of iodine
I60=5.0E-7; % mol/L;initial amounts of iodine
I70=5.0E-7; % mol/L;initial amounts of iodine
I80=5.0E-7; % mol/L;initial amounts of iodine
I90=5.0E-7; % mol/L;initial amounts of iodine
I100=5.0E-7; % mol/L;initial amounts of iodine
I110=5.0E-7; % mol/L;initial amounts of iodine
I120=5.0E-7; % mol/L;initial amounts of iodine
I130=5.0E-7; % mol/L;initial amounts of iodine
I140=5.0E-7; % mol/L;initial amounts of iodine
I150=5.0E-7; % mol/L;initial amounts of iodine
I160=5.0E-7; % mol/L;initial amounts of iodine
I170=5.0E-7; % mol/L;initial amounts of iodine
I180=5.0E-7; % mol/L;initial amounts of iodine

T10 = 0; % time

W0 = [I10 I20 I30 I40 I50 I60 I70 I80 I90 I100 I110 I120 I130 I140 I150...
      I160 I170 I180 T10]'; %' means inverse 

I1conc = @(W) W(1);   % Iodine concentration 
I2conc = @(W) W(2);   % Iodine concentration 
I3conc = @(W) W(3);   % Iodine concentration 
I4conc = @(W) W(4);   % Iodine concentration 
I5conc = @(W) W(5);   % Iodine concentration 
I6conc = @(W) W(6);   % Iodine concentration 
I7conc = @(W) W(7);   % Iodine concentration 
I8conc = @(W) W(8);   % Iodine concentration 
I9conc = @(W) W(9);   % Iodine concentration 
I10conc = @(W) W(10); % Iodine concentration
I11conc = @(W) W(11);   % Iodine concentration 
I12conc = @(W) W(12);   % Iodine concentration 
I13conc = @(W) W(13);   % Iodine concentration 
I14conc = @(W) W(14);   % Iodine concentration 
I15conc = @(W) W(15);   % Iodine concentration 
I16conc = @(W) W(16);   % Iodine concentration 
I17conc = @(W) W(17);   % Iodine concentration 
I18conc = @(W) W(18);   % Iodine concentration 


T1 = @(W) W(19);   %Time

% ---- Defining time series -------------
 
% standard Phanerozoic PO2,digitalize from figure 6 of the Mills et al. 2023
%--------time is only from 0 to 500 Ma: 
%--------8 term Fourier transfer. R-square: 0.9621: 
       a0 =       1.042;
       a1 =     -0.3147;
       b1 =      0.3836;
       a2 =        0.14;
       b2 =      0.1284;
       a3 =    -0.04214;
       b3 =      0.1274;
       a4 =     0.03857;
       b4 =    -0.05201;
       a5 =    0.008465;
       b5 =     0.01688;
       a6 =    0.002226;
       b6 =    -0.01917;
       a7 =     0.04496;
       b7 =     0.05274;
       a8 =     0.07475;
       b8 =    -0.02854;
       w =     0.01076;


 PO2 = @(W,t) (a0 + a1*cos(T1(W)*w) + b1*sin(T1(W)*w) + ...
       a2*cos(2*T1(W)*w) + b2*sin(2*T1(W)*w) + a3*cos(3*T1(W)*w) + b3*sin(3*T1(W)*w) + ...
       a4*cos(4*T1(W)*w) + b4*sin(4*T1(W)*w) + a5*cos(5*T1(W)*w) + b5*sin(5*T1(W)*w) + ...
       a6*cos(6*T1(W)*w) + b6*sin(6*T1(W)*w) + a7*cos(7*T1(W)*w) + b7*sin(7*T1(W)*w) + ...
       a8*cos(8*T1(W)*w) + b8*sin(8*T1(W)*w)).*heaviside(540-T1(W)) + 0.01.*heaviside(T1(W)-540);
%   
%-----now take wider range of PO2 from the Phanerozoic 
%-----higher end of PO2, same figure 6 from Mills et al.2023
%------X in Ma, Y in PAL; R-square: 0.9861
%7-term Fourier; Coefficients (with 95% confidence bounds):
       ca0 =      -14.64;
       ca1 =      -12.01;
       cb1 =       27.65;
       ca2 =       16.46;
       cb2 =       18.12;
       ca3 =       16.86;
       cb3 =      -5.525;
       ca4 =       1.007;
       cb4 =      -10.73;
       ca5 =       -4.93;
       cb5 =      -2.673;
       ca6 =      -1.655;
       cb6 =       1.396;
       ca7 =      0.2024;
       cb7 =      0.4816;
       cw =     0.007333;
       
PO2_High = @(W,t)(ca0 + ca1*cos(T1(W)*cw) + cb1*sin(T1(W)*cw) + ...
               ca2*cos(2*T1(W)*cw) + cb2*sin(2*T1(W)*cw) + ca3*cos(3*T1(W)*cw) + cb3*sin(3*T1(W)*cw) + ...
               ca4*cos(4*T1(W)*cw) + cb4*sin(4*T1(W)*cw) + ca5*cos(5*T1(W)*cw) + cb5*sin(5*T1(W)*cw) + ...
               ca6*cos(6*T1(W)*cw) + cb6*sin(6*T1(W)*cw) + ca7*cos(7*T1(W)*cw) + cb7*sin(7*T1(W)*cw))...
               *heaviside(540-T1(W)) + 0.1.*heaviside(T1(W)-540);

%-----lower end of PO2, same figure 6 from Mills et al.2023
%------X in Ma, Y in PAL; R-square: 0.9646
%7-term Fourier; Coefficients (with 95% confidence bounds):
       da0 =       2.714;
       da1 =       2.211;
       db1 =      -3.011;
       da2 =      -1.381;
       db2 =      -3.094;
       da3 =      -2.452;
       db3 =     -0.1268;
       da4 =     -0.9774;
       db4 =       1.327;
       da5 =       0.392;
       db5 =       0.694;
       da6 =      0.4158;
       db6 =     0.01046;
       da7 =    -0.05347;
       db7 =     -0.1213;
       dw =    0.007676;

PO2_Low = @(W,t)(da0 + da1*cos(T1(W)*dw) + db1*sin(T1(W)*dw) + ...
               da2*cos(2*T1(W)*dw) + db2*sin(2*T1(W)*dw) + da3*cos(3*T1(W)*dw) + db3*sin(3*T1(W)*dw) + ...
               da4*cos(4*T1(W)*dw) + db4*sin(4*T1(W)*dw) + da5*cos(5*T1(W)*dw) + db5*sin(5*T1(W)*dw) + ...
               da6*cos(6*T1(W)*dw) + db6*sin(6*T1(W)*dw) + da7*cos(7*T1(W)*dw) + db7*sin(7*T1(W)*dw))...
               *heaviside(540-T1(W)) + 1.0E-3.*heaviside(T1(W)-540);
%
%
%digitalized from figure 3a of Daines et al. 2017
%case 1: default scenario
%-------4 term Fourier; R-square: 0.9996: 
       aa0 =       5.002;
       aa1 =       2.107;
       ab1 =       4.688;
       aa2 =      -1.948;
       ab2 =      0.7423;
       aa3 =     -0.1685;
       ab3 =     -0.8751;
       aa4 =      0.1799;
       ab4 =    -0.02787;
       aw =      0.9655;

F_org= @(W,t) (aa0 + aa1*cos(log10(PO2(W,t))*aw) + ab1*sin(log10(PO2(W,t))*aw) + ...
              aa2*cos(2*log10(PO2(W,t))*aw) + ab2*sin(2*log10(PO2(W,t))*aw) + ...
              aa3*cos(3*log10(PO2(W,t))*aw) + ab3*sin(3*log10(PO2(W,t))*aw) + ...
              aa4*cos(4*log10(PO2(W,t))*aw) + ab4*sin(4*log10(PO2(W,t))*aw))*1.0E18;% mol/Ma
% at t=0 for PO2_standard, F_org = 5.16393E18; 
F_org_High= @(W,t) (aa0 + aa1*cos(log10(PO2_High(W,t))*aw) + ab1*sin(log10(PO2_High(W,t))*aw) + ...
              aa2*cos(2*log10(PO2_High(W,t))*aw) + ab2*sin(2*log10(PO2_High(W,t))*aw) + ...
              aa3*cos(3*log10(PO2_High(W,t))*aw) + ab3*sin(3*log10(PO2_High(W,t))*aw) + ...
              aa4*cos(4*log10(PO2_High(W,t))*aw) + ab4*sin(4*log10(PO2_High(W,t))*aw))*1.0E18;% mol/Ma
%  at t=0 for PO2_high, F_org = 5.57258E18;         
F_org_Low= @(W,t) (aa0 + aa1*cos(log10(PO2_Low(W,t))*aw) + ab1*sin(log10(PO2_Low(W,t))*aw) + ...
              aa2*cos(2*log10(PO2_Low(W,t))*aw) + ab2*sin(2*log10(PO2_Low(W,t))*aw) + ...
              aa3*cos(3*log10(PO2_Low(W,t))*aw) + ab3*sin(3*log10(PO2_Low(W,t))*aw) + ...
              aa4*cos(4*log10(PO2_Low(W,t))*aw) + ab4*sin(4*log10(PO2_Low(W,t))*aw))*1.0E18;% mol/Ma
               
%
%digitalized from figure 3a of Daines et al. 2017
%case 2: low uplift rate and high shale porosity: 
% Four-terms Fourier; R-square: 0.9998:
       ea0 =       2.899;
       ea1 =       2.004;
       eb1 =     -0.1498;
       ea2 =       0.325;
       eb2 =      0.2067;
       ea3 =     -0.1172;
       eb3 =      0.1778;
       ea4 =    -0.06195;
       eb4 =      0.0276;
       ew =        1.21;
       
F_org_uplift_highO2 = @(W,t) (ea0 + ea1*cos(log10(PO2_High(W,t))*ew) + eb1*sin(log10(PO2_High(W,t))*ew) + ...
           ea2*cos(2*log10(PO2_High(W,t))*ew) + eb2*sin(2*log10(PO2_High(W,t))*ew) + ... 
           ea3*cos(3*log10(PO2_High(W,t))*ew) + eb3*sin(3*log10(PO2_High(W,t))*ew) + ...
           ea4*cos(4*log10(PO2_High(W,t))*ew) + eb4*sin(4*log10(PO2_High(W,t))*ew))*1E18; %mol/Ma

F_org_uplift_lowO2 = @(W,t) (ea0 + ea1*cos(log10(PO2_Low(W,t))*ew) + eb1*sin(log10(PO2_Low(W,t))*ew) + ...
           ea2*cos(2*log10(PO2_Low(W,t))*ew) + eb2*sin(2*log10(PO2_Low(W,t))*ew) + ... 
           ea3*cos(3*log10(PO2_Low(W,t))*ew) + eb3*sin(3*log10(PO2_Low(W,t))*ew) + ...
           ea4*cos(4*log10(PO2_Low(W,t))*ew) + eb4*sin(4*log10(PO2_Low(W,t))*ew))*1E18; %mol/Ma

%digitalized from figure 3b of Daines et al. 2017: 
%case 3: larger reduced atmospheric flux from methane:default(black line): 
% Four terms Fourier; R-square: 0.9998
       fa0 =       5.977;
       fa1 =       3.831;
       fb1 =       3.425;
       fa2 =     -0.8922;
       fb2 =       1.558;
       fa3 =     -0.6602;
       fb3 =     -0.3886;
       fa4 =     0.07779;
       fb4 =     -0.1603;
       fw =       1.158;

F_org_reduced_highO2 = @(W,t) (fa0 + fa1*cos(log10(PO2_High(W,t))*fw) + fb1*sin(log10(PO2_High(W,t))*fw) + ...
           fa2*cos(2*log10(PO2_High(W,t))*fw) + fb2*sin(2*log10(PO2_High(W,t))*fw) + ... 
           fa3*cos(3*log10(PO2_High(W,t))*fw) + fb3*sin(3*log10(PO2_High(W,t))*fw) + ...
           fa4*cos(4*log10(PO2_High(W,t))*fw) + fb4*sin(4*log10(PO2_High(W,t))*fw))*1E18; %mol/Ma
       
 
F_org_reduced_lowO2 = @(W,t) (fa0 + fa1*cos(log10(PO2_Low(W,t))*fw) + fb1*sin(log10(PO2_Low(W,t))*fw) + ...
           fa2*cos(2*log10(PO2_Low(W,t))*fw) + fb2*sin(2*log10(PO2_Low(W,t))*fw) + ... 
           fa3*cos(3*log10(PO2_Low(W,t))*fw) + fb3*sin(3*log10(PO2_Low(W,t))*fw) + ...
           fa4*cos(4*log10(PO2_Low(W,t))*fw) + fb4*sin(4*log10(PO2_Low(W,t))*fw))*1E18; %mol/Ma

%digitalized from figure 3b of Daines et al. 2017;
%case 4: larger reduced atmospheric flux from methane:Proterozoic(green) + 'modern' case(blue): 
% 'green' and 'blue' represents their color in the quoted paper.
% apply Daines's 'modern' case to the Phanaerozoic.
% Four terms Fourier; R-square: 0.9997
       ga0 =       6.626;
       ga1 =       2.423;
       gb1 =       1.694;
       ga2 =     -0.3166;
       gb2 =       1.011;
       ga3 =     -0.4866;
       gb3 =     -0.1087;
       ga4 =     0.00697;
       gb4 =     -0.1397;
       gw =        1.22;

F_org_reduced_modern_L = @(W,t) (ga0 + ga1*cos(log10(PO2_Low(W,t))*gw) + gb1*sin(log10(PO2_Low(W,t))*gw) + ...
               ga2*cos(2*log10(PO2_Low(W,t))*gw) + gb2*sin(2*log10(PO2_Low(W,t))*gw) + ...
               ga3*cos(3*log10(PO2_Low(W,t))*gw) + gb3*sin(3*log10(PO2_Low(W,t))*gw) + ...
               ga4*cos(4*log10(PO2_Low(W,t))*gw) + gb4*sin(4*log10(PO2_Low(W,t))*gw))*1.0E18; %mol/Ma
           
F_org_reduced_modern_H = @(W,t) (ga0 + ga1*cos(log10(PO2_High(W,t))*gw) + gb1*sin(log10(PO2_High(W,t))*gw) + ...
               ga2*cos(2*log10(PO2_High(W,t))*gw) + gb2*sin(2*log10(PO2_High(W,t))*gw) + ...
               ga3*cos(3*log10(PO2_High(W,t))*gw) + gb3*sin(3*log10(PO2_High(W,t))*gw) + ...
               ga4*cos(4*log10(PO2_High(W,t))*gw) + gb4*sin(4*log10(PO2_High(W,t))*gw))*1.0E18; %mol/Ma 
%           
% Four terms Fourier; R-square: 0.9998       
       ha0 =       5.121;
       ha1 =       2.354;
       hb1 =       2.017;
       ha2 =     -0.4191;
       hb2 =      0.9675;
       ha3 =     -0.4111;
       hb3 =      -0.166;
       ha4 =     0.01831;
       hb4 =    -0.09697;
       hw =        1.109;
       
F_org_reduced_Pro_L = @(W,t) (ha0 + ha1*cos(log10(PO2_Low(W,t))*hw) + hb1*sin(log10(PO2_Low(W,t))*hw) + ...
               ha2*cos(2*log10(PO2_Low(W,t))*hw) + hb2*sin(2*log10(PO2_Low(W,t))*hw) + ...
               ha3*cos(3*log10(PO2_Low(W,t))*hw) + hb3*sin(3*log10(PO2_Low(W,t))*hw) + ...
               ha4*cos(4*log10(PO2_Low(W,t))*hw) + hb4*sin(4*log10(PO2_Low(W,t))*hw))*1.0E18; % Ma        
       
F_org_reduced_Pro_H = @(W,t) (ha0 + ha1*cos(log10(PO2_High(W,t))*hw) + hb1*sin(log10(PO2_High(W,t))*hw) + ...
               ha2*cos(2*log10(PO2_High(W,t))*hw) + hb2*sin(2*log10(PO2_High(W,t))*hw) + ...
               ha3*cos(3*log10(PO2_High(W,t))*hw) + hb3*sin(3*log10(PO2_High(W,t))*hw) + ...
               ha4*cos(4*log10(PO2_High(W,t))*hw) + hb4*sin(4*log10(PO2_High(W,t))*hw))*1.0E18; % Ma 
           
% use combo expression to make sure the individual case express correctly
F_org_combo_lowO2 = @(W,t)F_org_reduced_modern_L(W,t) *heaviside(540-T1(W))+ F_org_reduced_Pro_L(W,t) *heaviside(T1(W)-540);
F_org_combo_highO2 = @(W,t)F_org_reduced_modern_H(W,t) *heaviside(540-T1(W))+ F_org_reduced_Pro_H(W,t) *heaviside(T1(W)-540);
       
% Use curve fitting tool to find an expression for the I/Ca evolution
% overtime. N = 3397, T= 0-2000 Ma 
%------8-point Gaussian; R-square: 0.6729
%       ba1 =       6.794;
%       bb1 =       75.17;
%       bc1 =        6.43;
%       ba2 =   1.098e+06;
%       bb2 =      -24.71;
%       bc2 =       7.126;
%       ba3 =       1.735;
%       bb3 =       440.5;
%       bc3 =       5.743;
%       ba4 =       5.535;
%       bb4 =        51.3;
%       bc4 =        42.3;
%       ba5 =       11.35;
%       bb5 =       154.1;
%       bc5 =       22.01;
%       ba6 =       3.695;
%       bb6 =       390.6;
%       bc6 =       24.59;
%       ba7 =      0.5757;
%       bb7 =       171.8;
%       bc7 =       51.87;
%       ba8 =      0.5192;
%       bb8 =       583.3;
%       bc8 =       76.48;
%            
%I/Ca ratio at t= 0 is 7.8571E-6(old);      
%F_carb = @(W,t)((ba1*exp(-((T1(W)-bb1)/bc1)^2) + ba2*exp(-((T1(W)-bb2)/bc2)^2) + ...
%                ba3*exp(-((T1(W)-bb3)/bc3)^2) + ba4*exp(-((T1(W)-bb4)/bc4)^2) + ...
%                ba5*exp(-((T1(W)-bb5)/bc5)^2) + ba6*exp(-((T1(W)-bb6)/bc6)^2) + ...
%                ba7*exp(-((T1(W)-bb7)/bc7)^2) + ba8*exp(-((T1(W)-bb8)/bc8)^2))*3.198E13).*heaviside(540-T1(W))+ ...
%                0.5*3.198E13*heaviside(T1(W)-540);
            %assumes I/Ca ratio of 5E-7 in the Pre-cambrian to maximize
            %iodine removal. 
%
%
%--------------------------new parameterization--------------------------------------------
% break the I/Ca record into different session to better captures the curve: 
% at t = 0, I/Ca is 7.64447E-6
%carbonate sedimentary accumulation 3.198E13*1E6 mol/Ma from Milliman 1993 
%F_carb at time zero = (7.64447E-6*1.0E-6*3.198E13)*1.0E6 = 2.4447E8 mol/Ma. 
% 0-100Ma, 6-term Gaussian, R-square=0.5389: 
       ia1 =       6.659;
       ib1 =        74.3;
       ic1 =       5.049;
       ia2 =      0.3742;
       ib2 =        52.5;
       ic2 =       1.484;
       ia3 =   1.808e+13;
       ib3 =      -143.6;
       ic3 =       26.81;
       ia4 =      -1.108;
       ib4 =       41.83;
       ic4 =       3.126;
       ia5 =     -0.8475;
       ib5 =       26.74;
       ic5 =       10.57;
       ia6 =       5.803;
       ib6 =       50.22;
       ic6 =       41.79;
       
F_carb_100 = @(W,t) (ia1*exp(-((T1(W)-ib1)/ic1)^2) + ia2*exp(-((T1(W)-ib2)/ic2)^2) + ...
                   ia3*exp(-((T1(W)-ib3)/ic3)^2) + ia4*exp(-((T1(W)-ib4)/ic4)^2) + ...
                   ia5*exp(-((T1(W)-ib5)/ic5)^2) + ia6*exp(-((T1(W)-ib6)/ic6)^2))*3.198E13;
            
    
% 100-400Ma, 5-term Fourier, R-square=0.5005: 
       ja0 =       1.737;
       ja1 =      -1.832;
       jb1 =      0.7779;
       ja2 =      0.4988;
       jb2 =     -0.2788;
       ja3 =     -0.1669;
       jb3 =     -0.5104;
       ja4 =      0.5299;
       jb4 =      0.5939;
       ja5 =     -0.2794;
       jb5 =     -0.3037;
       jw =     0.02069;

F_carb_400 = @(W,t)(ja0 + ja1*cos(T1(W)*jw) + jb1*sin(T1(W)*jw) + ...
       ja2*cos(2*T1(W)*jw) + jb2*sin(2*T1(W)*jw) + ja3*cos(3*T1(W)*jw) + jb3*sin(3*T1(W)*jw) + ...
       ja4*cos(4*T1(W)*jw) + jb4*sin(4*T1(W)*jw) + ja5*cos(5*T1(W)*jw) + jb5*sin(5*T1(W)*jw))*3.198E13;
 
% 400-600 Ma; four term Guassian, R-square = 0.0777
       ka1 =       1.127;
       kb1 =       443.1;
       kc1 =        29.3;
       ka2 =       0.491;
       kb2 =         580;
       kc2 =       5.132;
       ka3 =      -0.143;
       kb3 =         445;
       kc3 =        2.61;
       ka4 =      0.3436;
       kb4 =         540;
       kc4 =       7.135;
       
F_carb_600 = @(W,t)(ka1*exp(-((T1(W)-kb1)/kc1)^2) + ka2*exp(-((T1(W)-kb2)/kc2)^2) + ...
                    ka3*exp(-((T1(W)-kb3)/kc3)^2) + ka4*exp(-((T1(W)-kb4)/kc4)^2)+1.0)*3.198E13;% 1.0 is not from the original expression
              
% 600-2000 Ma; 4 terms Guassian, R-square = 0.2476
       la1 =      0.8706;
       lb1 =        1565;
       lc1 =        17.9;
       la2 =      0.7757;
       lb2 =        1439;
       lc2 =       8.258;
       la3 =      0.4565;
       lb3 =       601.8;
       lc3 =       61.48;
       la4 =      0.2267;
       lb4 =        1037;
       lc4 =       263.8;
       
F_carb_rest = @(W,t) ((la1*exp(-((T1(W)-lb1)/lc1)^2) + la2*exp(-((T1(W)-lb2)/lc2)^2) + ...
              la3*exp(-((T1(W)-lb3)/lc3)^2) + la4*exp(-((T1(W)-lb4)/lc4)^2))+0.25)*3.198E13;%0.25 are added on top of the original expression
              
F_carb = @(W,t) 6.0*3.198E13*heaviside(50-T1(W))+ F_carb_100(W,t).*((50 <T1(W))&(T1(W)<100)) + F_carb_400(W,t).*((100 <T1(W))&(T1(W)<400)) + ...
    F_carb_600(W,t).*((400<T1(W))&(T1(W)<600)) + F_carb_rest(W,t).*heaviside(T1(W)-600);%fix 0-50 Ma to smooth out spike up increase of iodine
%
%
%
%----------------above F_carb expression as either the Fourier or Gausian expression
%contains too much zero in the Proterozoic. 
% 4 terms polynomial; R-square: 0.4381
%
%       ap1 =   5.345e-12;
%       ap2 =  -2.408e-08;
%       ap3 =   3.701e-05;
%       ap4 =    -0.02249;
%       ap5 =       5.204;
%       
%F_carb = @(W,t) (ap1*T1(W)^4 + ap2*T1(W)^3 + ap3*T1(W)^2 + ap4*T1(W) + ap5)*3.198E13;% mol/Ma
% when t= 0, I/Ca ratio is 5.204E-6
%
%evolution of iodine in the upper continental crust; digitalized from figure 7g of Han
%et al. 2023 GCA   
%-----R-square: 0.2786:       
       p1 =   2.095e-18;  
       p2 =  -6.152e-15;  
       p3 =  -1.322e-11;  
       p4 =   6.534e-08; 
       
F_crust = @(W,t) p1*(T1(W))^3 + p2*(T1(W))^2 + p3*(T1(W)) + p4;

F_water_evo = @(W,t) F_water*(F_crust(W,t)/6.534E-8);
F_riversed_evo = @(W,t) F_riversed*(F_crust(W,t)/6.534E-8);
F_volc_evo = @(W,t) F_volc*(F_crust(W,t)/6.534E-8);

%--------------------set up pO2 vs. I/TOC -----------------------------
% establish the relationship between atmospheric O2 in PAL and seafloor O2
% in uM. Reading results from figure 5,Cole et al. 2022: 
%2 terms Fourier; R-square: 0.9982
%       ma0 =       87.67;
%       ma1 =      -78.72;
%       mb1 =      -38.08;
%       ma2 =      -9.398;
%       mb2 =       16.63;
%       mw =       3.694;
%       
%O2_seafloor_default = @(W,t) (ma0 + ma1*cos(PO2(W,t)*mw) + mb1*sin(PO2(W,t)*mw) + ...
%                              ma2*cos(2*PO2(W,t)*mw) + mb2*sin(2*PO2(W,t)*mw)).*heaviside(450-T1(W))+ ...
%                             0.*heaviside(T1(W)-450);
%                               
%O2_seafloor_LowO2 = @(W,t) (ma0 + ma1*cos(PO2_Low(W,t)*mw) + mb1*sin(PO2_Low(W,t)*mw) + ...
%                            ma2*cos(2*PO2_Low(W,t)*mw) + mb2*sin(2*PO2_Low(W,t)*mw)).*heaviside(450-T1(W))+...
%                             0.*heaviside(T1(W)-450); 
%                             
%O2_seafloor_HighO2 = @(W,t) ma0 + ma1*cos(PO2_High(W,t)*mw) + mb1*sin(PO2_High(W,t)*mw) + ...
%                            ma2*cos(2*PO2_High(W,t)*mw) + mb2*sin(2*PO2_High(W,t)*mw);
%
%--------abandon the above expression as the 2 term fourier handles poorly
%above 1 PAL of pO2: 
% three term polynomials,R-square: 0.997: 
       ap1 =      -76.45;
       ap2 =       196.2;
       ap3 =       70.69;
       ap4 =      -1.947; 
       
O2_seafloor_default = @(W,t) (ap1*PO2(W,t)^3 + ap2*PO2(W,t)^2 + ap3*PO2(W,t) + ap4).*heaviside(540-T1(W))+ ...
                             0.*heaviside(T1(W)-540);

O2_seafloor_LowO2 = @(W,t) (ap1*PO2_Low(W,t)^3 + ap2*PO2_Low(W,t)^2 + ap3*PO2_Low(W,t) + ap4).*heaviside(540-T1(W))+ ...
                             0.*heaviside(T1(W)-540);

O2_seafloor_HighO2 = @(W,t) (ap1*PO2_High(W,t)^3 + ap2*PO2_High(W,t)^2 + ap3*PO2_High(W,t) + ap4).*heaviside(540-T1(W))+ ...
                             0.*heaviside(T1(W)-540);

% ----------now set up the relationship between ambient seafloor O2 and changes in I/TOC:                        
% use compiled modern surface sediment data from Zhou et al. 2017. 
% revised on Nov 6 2023; previous 2 term power expression can't handle full
% anoxia; 
% revised on Oct 10 2023; previous scaling on the y-axis is wrong
% one term gassuain,R-square: 0.5022
       a_50 =        1841;  %(1312, 2369)
       b_50 =       256.9;  %(111.9, 401.8)
       c_50 =       171.2;  %(39.56, 302.8)
       
PO2IC_default_wrongscale =  @(W,t) a_50*exp(-((O2_seafloor_default(W,t)-b_50)/c_50)^2); 
PO2IC_LowO2_wrongscale =  @(W,t) a_50*exp(-((O2_seafloor_LowO2(W,t)-b_50)/c_50)^2);
PO2IC_HighO2_wrongscale =  @(W,t) a_50*exp(-((O2_seafloor_HighO2(W,t)-b_50)/c_50)^2);
%
%For the above expression,we pick the modern baseline to be 1551(umol/mol)of I/OC corresponds
% with 186 uM/kg of abmient seafloor oxygen. The 186 uM/kg derives from Cole et al.
% 2022 corresponds with 1 PAL of O2. 
%
% scale the above wrong scale expression into unitless one, and track the
% relative changes of I/OC ratio against 1 PAL of O2: 

PO2IC_default = @(W,t) PO2IC_default_wrongscale(W,t)/1551.0;%*1E-6; 
PO2IC_LowO2 = @(W,t)   PO2IC_LowO2_wrongscale(W,t)/1551.0; %*1E-6; 
PO2IC_HighO2 = @(W,t)  PO2IC_HighO2_wrongscale(W,t)/1551.0;%*1E-6; 

% now consider biological speciation:
% Y = [PO2IC * (Dt/D0 * K)]^C
% Y shall be the relative changes in I/OC ratio to the base line considered
% all factors; 
% PO2IC(defined above already, and expressed as RO2 in the manuscript) is the relative change of I/OC ratio against 1 PAL of O2(just
% considers the influence of seafloor O2). 
% Dt is the iodine concentration at time t;
% D0 is the present iodine concentration of 5E-7( parameter I10
% previously). 
% K is the mean biological uptake coefficient assembled from van Bergeijk
% et al. 2016. K = 0.57 for 0-450 Ma, K = 0.4 for 450- 1000 Ma, and K = 0.1
% for 1000-2000 Ma. 
% C is the scaling factor, values in theory can be between 0.5 and 1.

K_factor = @(W,t) 0.57*heaviside(540-T1(W)) + 0.4*((540 <T1(W))&(T1(W)<1000)) + ...
    0.1*((1000 <T1(W))&(T1(W)<2000));

%-------C equals to 0(default,uplift,reduced etc.have not been introduced yet,for numbering purpose only) 
PO2IC_default_0_1 = @(W,t) (PO2IC_default(W,t)*(I1conc(W)/5.0E-7 * K_factor(W,t)))^0.5;% Iodine_default_standard PO2
PO2IC_HighO2_0_2 = @(W,t) (PO2IC_HighO2(W,t)*(I2conc(W)/5.0E-7 * K_factor(W,t)))^0.5; % Iodine_default_High PO2
PO2IC_LowO2_0_3 = @(W,t) (PO2IC_LowO2(W,t)*(I3conc(W)/5.0E-7 * K_factor(W,t)))^0.5; % Iodine_default_Low PO2
PO2IC_HighO2_0_4 = @(W,t) (PO2IC_HighO2(W,t)*(I4conc(W)/5.0E-7 * K_factor(W,t)))^0.5;% Iodine_uplift_High PO2
PO2IC_LowO2_0_5 = @(W,t) (PO2IC_LowO2(W,t)*(I5conc(W)/5.0E-7 * K_factor(W,t)))^0.5;% Iodine_uplift_Low PO2
PO2IC_HighO2_0_6 = @(W,t) (PO2IC_HighO2(W,t)*(I6conc(W)/5.0E-7 * K_factor(W,t)))^0.5;% Iodine_reduced_High PO2
PO2IC_LowO2_0_7 = @(W,t) (PO2IC_LowO2(W,t)*(I7conc(W)/5.0E-7 * K_factor(W,t)))^0.5;% Iodine_reduced_Low PO2
PO2IC_HighO2_0_8 = @(W,t) (PO2IC_HighO2(W,t)*(I8conc(W)/5.0E-7 * K_factor(W,t)))^0.5;% Iodine_reduced_Combo_High PO2
PO2IC_LowO2_0_9 = @(W,t) (PO2IC_LowO2(W,t)*(I9conc(W)/5.0E-7 * K_factor(W,t)))^0.5;% Iodine_reduced_Combo_Low PO2
%
%-------C equals to 1
PO2IC_default_1_10 = @(W,t) (PO2IC_default(W,t)*(I10conc(W)/5.0E-7 * K_factor(W,t)))^0.9;% Iodine_default_standard PO2
PO2IC_HighO2_1_11 = @(W,t) (PO2IC_HighO2(W,t)*(I11conc(W)/5.0E-7 * K_factor(W,t)))^0.9; % Iodine_default_High PO2
PO2IC_LowO2_1_12 = @(W,t) (PO2IC_LowO2(W,t)*(I12conc(W)/5.0E-7 * K_factor(W,t)))^0.9; % Iodine_default_Low PO2
PO2IC_HighO2_1_13 = @(W,t) (PO2IC_HighO2(W,t)*(I13conc(W)/5.0E-7 * K_factor(W,t)))^0.9;% Iodine_uplift_High PO2
PO2IC_LowO2_1_14 = @(W,t) (PO2IC_LowO2(W,t)*(I14conc(W)/5.0E-7 * K_factor(W,t)))^0.9;% Iodine_uplift_Low PO2
PO2IC_HighO2_1_15 = @(W,t) (PO2IC_HighO2(W,t)*(I15conc(W)/5.0E-7 * K_factor(W,t)))^0.9;% Iodine_reduced_High PO2
PO2IC_LowO2_1_16 = @(W,t) (PO2IC_LowO2(W,t)*(I16conc(W)/5.0E-7 * K_factor(W,t)))^0.9;% Iodine_reduced_Low PO2
PO2IC_HighO2_1_17 = @(W,t) (PO2IC_HighO2(W,t)*(I17conc(W)/5.0E-7 * K_factor(W,t)))^0.9;% Iodine_reduced_Combo_High PO2
PO2IC_LowO2_1_18 = @(W,t) (PO2IC_LowO2(W,t)*(I18conc(W)/5.0E-7 * K_factor(W,t)))^0.9;% Iodine_reduced_Combo_Low PO2

%----F_input_total at default = 3.16417E9 mol/yr.
%----time dependent I/TOC ratio in sedimentary rocks-----
%----base output for iodine in F_org at time zero should be 2.97229E9 mol/yr(need revision if base ...
%input values or base I/Ca changes) 
%----referenced F_org,read from figure 3: 
% default-standardO2 reads 5.16E12 mol/yr;
% default-HighO2 use 5.58E12 mol/yr;
% uplift-LowO2 and default-LowO2 has 4.98E12 mol/yr;
% uplift-HighO2 has 5.16E12 mol/yr;
% reduced_LowO2(both combo and default) are 8.05E12 mol/yr;
% reduced_combo_HighO2 is 8.67E12 mol/yr;
% reduced_HighO2 is 9.0E12 mol/yr;
%----the I/Ca depletion in OAE2 is about 53%. 

R_C_base_0 = @(W,t)        PO2IC_default_0_1(W,t)*5.746E-4;% apply to defult-standard O2;C=0           
R_C_base_highO2_0 = @(W,t) PO2IC_HighO2_0_2(W,t)*5.3267E-4;% apply to defult-high O2;C=0                                           
R_C_base_lowO2_0 = @(W,t)  PO2IC_LowO2_0_3(W,t)*5.968E-4;% apply to default-low O2; C=0
R_C_uplift_highO2_0 = @(W,t) PO2IC_HighO2_0_4(W,t)*5.76E-4;% apply to uplift-High O2;C=0
R_C_uplift_lowO2_0 = @(W,t)  PO2IC_LowO2_0_5(W,t)*5.968E-4;% apply to uplift-lowO2; C=0 
R_C_reduced_highO2_0 = @(W,t) PO2IC_HighO2_0_6(W,t)*3.30E-4;  % apply to reduced-(default)-High O2;C=0           
R_C_reduced_lowO2_0 = @(W,t) PO2IC_LowO2_0_7(W,t)*3.69E-4; % apply to reduced-(default)-Low O2;C=0                             
R_C_combo_highO2_0 =  @(W,t) PO2IC_HighO2_0_8(W,t)*3.3597E-4; % applt to reduced-combo-High O2;C=0
R_C_combo_lowO2_0 = @(W,t) PO2IC_LowO2_0_9(W,t)*3.69E-4; % apply to reduced-combo-Low O2;C=0 

R_C_base_1 = @(W,t)        PO2IC_default_1_10(W,t)*5.746E-4;% apply to defult-standard O2;C=1           
R_C_base_highO2_1 = @(W,t) PO2IC_HighO2_1_11(W,t)*5.3267E-4;% apply to defult-high O2;C=1                                           
R_C_base_lowO2_1 = @(W,t)  PO2IC_LowO2_1_12(W,t)*5.968E-4;% apply to default-low O2; C=1
R_C_uplift_highO2_1 = @(W,t) PO2IC_HighO2_1_13(W,t)*5.76E-4;% apply to uplift-High O2;C=1
R_C_uplift_lowO2_1 = @(W,t)  PO2IC_LowO2_1_14(W,t)*5.968E-4;% apply to uplift-lowO2; C=1 
R_C_reduced_highO2_1 = @(W,t) PO2IC_HighO2_1_15(W,t)*3.30E-4;  % apply to reduced-(default)-High O2;C=1           
R_C_reduced_lowO2_1 = @(W,t) PO2IC_LowO2_1_16(W,t)*3.69E-4; % apply to reduced-(default)-Low O2;C=1                             
R_C_combo_highO2_1 =  @(W,t) PO2IC_HighO2_1_17(W,t)*3.3597E-4; % applt to reduced-combo-High O2;C=1
R_C_combo_lowO2_1 = @(W,t) PO2IC_LowO2_1_18(W,t)*3.69E-4; % apply to reduced-combo-Low O2;C=1 

% ------------- SETTING UP THE ODEs ------------------
% Iodine_default_standard PO2,C=0
dIdt1  = @(t,W) (F_water_evo(W,t)+ F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_base_0(W,t).*F_org(W,t)- ...
                F_carb(W,t))/M_ocean;
% Iodine_default_High PO2,C=0  
dIdt2  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_base_highO2_0(W,t).*F_org_High(W,t) - ...
                F_carb(W,t))/M_ocean;
% Iodine_default_Low PO2,C=0  
dIdt3  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_base_lowO2_0(W,t).*F_org_Low(W,t) - ...
                F_carb(W,t))/M_ocean;   
% Iodine_uplift_High PO2,C=0
dIdt4  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_uplift_highO2_0(W,t).*F_org_uplift_highO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_uplift_Low PO2,C=0
dIdt5  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_uplift_lowO2_0(W,t).*F_org_uplift_lowO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_reduced_High PO2,C=0
dIdt6  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_reduced_highO2_0(W,t).*F_org_reduced_highO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_reduced_Low PO2,C=0
dIdt7  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_reduced_lowO2_0(W,t).*F_org_reduced_lowO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_reduced_Combo_High PO2,C=0
dIdt8  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_combo_highO2_0(W,t).*F_org_combo_highO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_reduced_Combo_Low PO2,C=0
dIdt9  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_combo_lowO2_0(W,t).*F_org_combo_lowO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
            
            
            
% Iodine_default_standard PO2,C=1
dIdt10  = @(t,W) (F_water_evo(W,t)+ F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_base_1(W,t).*F_org(W,t)- ...
                F_carb(W,t))/M_ocean;
% Iodine_default_High PO2,C=1  
dIdt11  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_base_highO2_1(W,t).*F_org_High(W,t) - ...
                F_carb(W,t))/M_ocean;
% Iodine_default_Low PO2,C=1  
dIdt12  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_base_lowO2_1(W,t).*F_org_Low(W,t) - ...
                F_carb(W,t))/M_ocean;   
% Iodine_uplift_High PO2,C=1
dIdt13  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_uplift_highO2_1(W,t).*F_org_uplift_highO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_uplift_Low PO2,C=1
dIdt14  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_uplift_lowO2_1(W,t).*F_org_uplift_lowO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_reduced_High PO2,C=1
dIdt15  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_reduced_highO2_1(W,t).*F_org_reduced_highO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_reduced_Low PO2,C=1
dIdt16  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_reduced_lowO2_1(W,t).*F_org_reduced_lowO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_reduced_Combo_High PO2,C=1
dIdt17  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_combo_highO2_1(W,t).*F_org_combo_highO2(W,t) - ...
                F_carb(W,t))/M_ocean; 
% Iodine_reduced_Combo_Low PO2,C=1
dIdt18  = @(t,W) (F_water_evo(W,t)+F_riversed_evo(W,t)+F_volc_evo(W,t) - R_C_combo_lowO2_1(W,t).*F_org_combo_lowO2(W,t) - ...
                F_carb(W,t))/M_ocean;             

% Time series
dTdt = @(t,W) 1;


% ----- Call the numerical solver  -----------
dWdt1 = @(t,W) [dIdt1(t,W) dIdt2(t,W) dIdt3(t,W) dIdt4(t,W) dIdt5(t,W) dIdt6(t,W) dIdt7(t,W) dIdt8(t,W) dIdt9(t,W) ...
    dIdt10(t,W) dIdt11(t,W) dIdt12(t,W) dIdt13(t,W) dIdt14(t,W) dIdt15(t,W) dIdt16(t,W) dIdt17(t,W) dIdt18(t,W)...
    dTdt(t,W)]' ; % define the derivative

options = odeset('MaxStep',maxstep);        % set solver options
[T,Y] = ode15s(dWdt1,[t0 tmax],W0,options);% call the solver 1



% -------------- PLOT RESULTS  -------------
Iodinec0 = Y(:,1);  % plot concentrations in Mol/L under standard PO2,default F_org,C=0
Iodinech0 =Y(:,2);  %iodine concentrations in Mol/L under high PO2,default F_org,C=0
Iodinecl0 =Y(:,3);  %iodine concentrations in Mol/L under low PO2,default F_org,C=0
Iodinecuh0 =Y(:,4);  %iodine concentrations in Mol/L under high PO2,low uplift F_org,C=0
Iodinecul0 =Y(:,5);  %iodine concentrations in Mol/L under low PO2,low uplift F_org,C=0
Iodinecrh0 =Y(:,6);  %iodine concentrations in Mol/L under high PO2,F_org under high reduced flux,default,C=0 
Iodinecrl0 =Y(:,7);  %iodine concentrations in Mol/L under low PO2,F_org under high reduced flux,default,C=0
Iodinecch0 =Y(:,8);  %iodine concentrations in Mol/L under high PO2,F_org under high reduced flux,combo,C=0
Iodineccl0 =Y(:,9);  %iodine concentrations in Mol/L under low PO2,F_org under high reduced flux,combo,C=0

Iodinec1 = Y(:,10);  % plot concentrations in Mol/L under standard PO2,default F_org,C=1
Iodinech1 =Y(:,11);  %iodine concentrations in Mol/L under high PO2,default F_org,C=1
Iodinecl1 =Y(:,12);  %iodine concentrations in Mol/L under low PO2,default F_org,C=1
Iodinecuh1 =Y(:,13);  %iodine concentrations in Mol/L under high PO2,low uplift F_org,C=1
Iodinecul1 =Y(:,14);  %iodine concentrations in Mol/L under low PO2,low uplift F_org,C=1
Iodinecrh1 =Y(:,15);  %iodine concentrations in Mol/L under high PO2,F_org under high reduced flux,default,C=1 
Iodinecrl1 =Y(:,16);  %iodine concentrations in Mol/L under low PO2,F_org under high reduced flux,default,C=1
Iodinecch1 =Y(:,17);  %iodine concentrations in Mol/L under high PO2,F_org under high reduced flux,combo,C=1
Iodineccl1 =Y(:,18);  %iodine concentrations in Mol/L under low PO2,F_org under high reduced flux,combo,C=1
%
% Iodine
figure(1)
%for cases of C=0.5
plot(T, Iodinec0,'r','LineWidth',2.0);
hold on;
plot(T, Iodinech0,'g','LineWidth',2.0);
plot(T,Iodinecl0,'b','LineWidth',2.0);
plot(T,Iodinecuh0,'c','LineWidth',2.0);
plot(T,Iodinecul0,'m','LineWidth',2.0);
plot(T,Iodinecrh0,'y','LineWidth',2.0);
plot(T,Iodinecrl0,'k','LineWidth',2.0);
plot(T,Iodinecch0,'Color', [0.5 0.5 0.5],'LineWidth',2.0);
plot(T,Iodineccl0,'Color', [0.7 0.7 0.7],'LineWidth',2.0);

%for cases of C=1
plot(T, Iodinec1,'r--','LineWidth',2.0);
plot(T, Iodinech1,'g--','LineWidth',2.0);
plot(T,Iodinecl1,'b--','LineWidth',2.0);
plot(T,Iodinecuh1,'c--','LineWidth',2.0);
plot(T,Iodinecul1,'m--','LineWidth',2.0);
plot(T,Iodinecrh1,'y--','LineWidth',2.0);
plot(T,Iodinecrl1,'k--','LineWidth',2.0);
plot(T,Iodinecch1,'--','Color', [0.5 0.5 0.5],'LineWidth',2.0);
plot(T,Iodineccl1,'--','Color', [0.7 0.7 0.7],'LineWidth',2.0);
hold off; 
grid off;
set(gca,'XDir','reverse')
ylabel ('Iodine,mol/L' );
xlabel ('Geological time(Ma)');
% Increase the font size of the axes tick labels
ax = gca;
set(ax, 'FontSize', 17, 'LineWidth', 2.0,'FontWeight', 'bold'); % Set the 'LineWidth' property for the axes
legend('default-standardO2-0.5','default-highO2-0.5','default-lowO2-0.5','uplift-highO2-0.5','uplift-lowO2-0.5',...
    'reduced-highO2-0.5','reduced-lowO2-0.5','combo-highO2-0.5','combo-lowO2-0.5',...
    'default-standardO2-1','default-highO2-1','default-lowO2-1','uplift-highO2-1','uplift-lowO2-1',...
    'reduced-highO2-1','reduced-lowO2-1','combo-highO2-1','combo-lowO2-1')

figure(2)
%for cases of C=0.5: 
semilogy(T,Iodinec0,'r','LineWidth',2.0);
hold on; 
semilogy(T,Iodinech0,'g','LineWidth',2.0);
semilogy(T,Iodinecl0,'b','LineWidth',2.0);
semilogy(T,Iodinecuh0,'c','LineWidth',2.0);
semilogy(T,Iodinecul0,'m','LineWidth',2.0);
semilogy(T,Iodinecrh0,'y','LineWidth',2.0);
semilogy(T,Iodinecrl0,'k','LineWidth',2.0);
semilogy(T,Iodinecch0,'Color', [0.5 0.5 0.5],'LineWidth',2.0);
semilogy(T,Iodineccl0,'Color', [0.7 0.7 0.7],'LineWidth',2.0);

%for cases of C=1: 
semilogy(T,Iodinec1,'r--','LineWidth',2.0);
semilogy(T,Iodinech1,'g--','LineWidth',2.0);
semilogy(T,Iodinecl1,'b--','LineWidth',2.0);
semilogy(T,Iodinecuh1,'c--','LineWidth',2.0);
semilogy(T,Iodinecul1,'m--','LineWidth',2.0);
semilogy(T,Iodinecrh1,'y--','LineWidth',2.0);
semilogy(T,Iodinecrl1,'k--','LineWidth',2.0);
semilogy(T,Iodinecch1,'--', 'Color', [0.5 0.5 0.5],'LineWidth',2.0);
semilogy(T,Iodineccl1,'--', 'Color', [0.7 0.7 0.7],'LineWidth',2.0);
hold off;
grid off;
set(gca,'XDir','reverse')
ylabel ('Iodine,mol/L' );
xlabel ('Geological time(Ma) ');
% Increase the font size of the axes tick labels
ax = gca;
set(ax, 'FontSize', 17, 'LineWidth', 2.0,'FontWeight', 'bold'); % Set the 'LineWidth' property for the axes
%legend('default-standardO2-0.5','default-highO2-0.5','default-lowO2-0.5','uplift-highO2-0.5','uplift-lowO2-0.5',...
%    'reduced-highO2-0.5','reduced-lowO2-0.5','combo-highO2-0.5','combo-lowO2-0.5',...
%    'default-standardO2-1','default-highO2-1','default-lowO2-1','uplift-highO2-1','uplift-lowO2-1',...
%    'reduced-highO2-1','reduced-lowO2-1','combo-highO2-1','combo-lowO2-1')

figure(3)
for i=1:size(T,1)
   
    wForg(i) = F_org(Y(i,:),T(i));
    hForg(i) = F_org_High(Y(i,:),T(i));
    lForg(i) = F_org_Low(Y(i,:),T(i));
    uhForg(i) = F_org_uplift_highO2(Y(i,:),T(i));
    ulForg(i) = F_org_uplift_lowO2(Y(i,:),T(i));
    rhForg(i) = F_org_reduced_highO2(Y(i,:),T(i));
    rlForg(i) = F_org_reduced_lowO2(Y(i,:),T(i));
    chForg(i) = F_org_combo_highO2(Y(i,:),T(i));
    clForg(i) = F_org_combo_lowO2(Y(i,:),T(i));
    
end
plot(T,wForg,'k','LineWidth',3.0);
hold on; 
plot(T,hForg,'LineWidth',3.0);
plot(T,lForg,'--','LineWidth',3.0);
plot(T,uhForg,'LineWidth',3.0);
plot(T,ulForg,'--','LineWidth',3.0);
plot(T,rhForg,'LineWidth',3.0);
plot(T,rlForg,'--','LineWidth',3.0);
plot(T,chForg,'LineWidth',3.0);
plot(T,clForg,'--','LineWidth',3.0);
hold off; 
grid off;
set(gca,'XDir','reverse')
ylabel ('Organic Carbon Burial rate, Mol/Ma','FontSize', 17);
xlabel ('Geological time(Ma)','FontSize', 17);
% Increase the font size of the axes tick labels
ax = gca;
set(ax, 'FontSize', 17, 'LineWidth', 2.0,'FontWeight', 'bold'); % Set the 'LineWidth' property for the axes
legend('default-standardO2','default-highO2','default-lowO2','uplift-highO2','uplift-lowO2',...
    'reduced-highO2','reduced-lowO2','combo-highO2','combo-lowO2')

figure(4)
for i=1:size(T,1)
   
    wPO2(i) = PO2(Y(i,:),T(i));
    hwPO2(i) = PO2_High(Y(i,:),T(i));
    lwPO2(i) = PO2_Low(Y(i,:),T(i));
    
end
semilogy(T,wPO2,'LineWidth',3.0);
hold on;
semilogy(T,hwPO2,'LineWidth',3.0);
semilogy(T,lwPO2,'LineWidth',3.0);
hold off; 
grid off;
set(gca,'XDir','reverse')
ylabel ('Atmospheric O2, PAL' );
xlabel ('Geological time(Ma) ','FontSize', 17);
% Increase the font size of the axes tick labels
ax = gca;
set(ax, 'FontSize', 17, 'LineWidth', 2.0,'FontWeight', 'bold'); % Set the 'LineWidth' property for the axes
legend('standard O2','high O2','low O2')

%figure(5)
%for i=1:size(T,1)
%   
%    wR_C_base_0(i) = R_C_base_0(Y(i,:),T(i));
%    wR_C_base_highO2_0(i) = R_C_base_highO2_0(Y(i,:),T(i));
%    wR_C_base_lowO2_0(i) = R_C_base_lowO2_0(Y(i,:),T(i));
%    wR_C_uplift_highO2_0(i) =R_C_uplift_highO2_0(Y(i,:),T(i));
%    wR_C_reduced_highO2_0(i) = R_C_reduced_highO2_0(Y(i,:),T(i));
%    wR_C_reduced_lowO2_0(i) = R_C_reduced_lowO2_0(Y(i,:),T(i));
%    wR_C_combo_highO2_0(i) = R_C_combo_highO2_0(Y(i,:),T(i));
%    wR_C_reduced_combo_lowO2_0(i) = R_C_combo_lowO2_0(Y(i,:),T(i));
%    
%end
%plot(T, wR_C_base_0,'k','LineWidth',2.0);
%hold on;
%plot(T, wR_C_base_highO2_0,'LineWidth',2.0); 
%plot(T,wR_C_base_lowO2_0,'--','LineWidth',2.0);
%plot(T,wR_C_uplift_highO2_0,'LineWidth',2.0);
%plot(T,wR_C_reduced_highO2_0,'LineWidth',2.0);
%plot(T,wR_C_reduced_lowO2_0,'--','LineWidth',2.0);
%plot(T,wR_C_combo_highO2_0,'LineWidth',2.0);
%plot(T,wR_C_reduced_combo_lowO2_0,'--','LineWidth',2.0)
%old off; 
%grid off;
%set(gca,'XDir','reverse')
%ylabel ('I/TOC' );
%xlabel ('Ma ');
%legend('RC-base-0.5','RC-base-highO2-0.5','RC-base-lowO2-0.5','RC-uplift-highO2-0.5','RC-reduced-highO2-0.5',...
%   'RC-reduced-lowO2-0.5','RC-combo-highO2-0.5','RC-combo-lowO2-0.5')
%
%figure(6)
%for i=1:size(T,1)
%   
%    wR_C_base_1(i) = R_C_base_1(Y(i,:),T(i));
%    wR_C_base_highO2_1(i) = R_C_base_highO2_1(Y(i,:),T(i));
%    wR_C_base_lowO2_1(i) = R_C_base_lowO2_1(Y(i,:),T(i));
%    wR_C_uplift_highO2_1(i) =R_C_uplift_highO2_1(Y(i,:),T(i));
%    wR_C_reduced_highO2_1(i) = R_C_reduced_highO2_1(Y(i,:),T(i));
%    wR_C_reduced_lowO2_1(i) = R_C_reduced_lowO2_1(Y(i,:),T(i));
%    wR_C_combo_highO2_1(i) = R_C_combo_highO2_1(Y(i,:),T(i));
%    wR_C_reduced_combo_lowO2_1(i) = R_C_combo_lowO2_1(Y(i,:),T(i));
%    
%end
%plot(T, wR_C_base_1,'k','LineWidth',2.0);
%hold on;
%plot(T, wR_C_base_highO2_1,'LineWidth',2.0);
%plot(T,wR_C_base_lowO2_1,'--','LineWidth',2.0); 
%plot(T,wR_C_uplift_highO2_1,'LineWidth',2.0);
%plot(T,wR_C_reduced_highO2_1,'LineWidth',2.0);
%plot(T,wR_C_reduced_lowO2_1,'--','LineWidth',2.0);
%plot(T,wR_C_combo_highO2_1,'LineWidth',2.0);
%plot(T,wR_C_reduced_combo_lowO2_1,'--','LineWidth',2.0)
%hold off; 
%grid on;
%set(gca,'XDir','reverse')
%ylabel ('I/TOC' );
%xlabel ('Ma ');
%legend('RC-base-0.9','RC-base-highO2-0.9','RC-base-lowO2-0.9','RC-uplift-highO2-0.9','RC-reduced-highO2-0.9',...
%    'RC-reduced-lowO2-0.9','RC-combo-highO2-0.9','RC-combo-lowO2-0.9');
%
%figure(7) 
%for i=1:size(T,1)
%   
%    wPO2IC(i) = PO2IC_default(Y(i,:),T(i));
%    hwPO2IC(i) = PO2IC_HighO2(Y(i,:),T(i));
%    lwPO2IC(i) = PO2IC_LowO2(Y(i,:),T(i));
%    
%end
%plot(T,wPO2IC,'LineWidth',2.0);
%hold on;
%plot(T,hwPO2IC,'LineWidth',2.0);
%plot(T,lwPO2IC,'LineWidth',2.0);
%hold off; 
%grid on;
%set(gca,'XDir','reverse')
%ylabel ('I/TOC ratio adjustment,pO2 dependent only' );
%xlabel ('Ma ');
%legend('standard PO2','high PO2','low PO2')
%
%figure(8)
%for i=1:size(T,1)
%   
%    mwPO2(i) = O2_seafloor_default(Y(i,:),T(i));
%    mhwPO2(i) = O2_seafloor_HighO2(Y(i,:),T(i));
%    mlwPO2(i) = O2_seafloor_LowO2(Y(i,:),T(i));
%    
%end
%plot(T,mwPO2,'LineWidth',2.0);
%hold on;
%plot(T,mhwPO2,'LineWidth',2.0);
%plot(T,mlwPO2,'LineWidth',2.0);
%hold off; 
%grid on;
%set(gca,'XDir','reverse')
%ylabel ('seafloor O2,uM' );
%xlabel ('Ma ');
%legend('standard PO2','high PO2','low PO2')
%
%figure(9)
%for i=1:size(T,1)
%    FCarb(i) = F_carb(Y(i,:),T(i))/3.198E13;
%end
%plot(T,FCarb)
%grid on;
%set(gca,'XDir','reverse')
%ylabel ('I to Ca ratio/1E6' );
%xlabel ('Ma ');
%
%
%
%figure(9)
%for i=1:size(T,1)
%    PO2IC_0_1_plot(i) = PO2IC_default_0_1(Y(i,:),T(i));
%    PO2IC_0_2_plot(i) = PO2IC_HighO2_0_2(Y(i,:),T(i));
%    PO2IC_0_3_plot(i) = PO2IC_LowO2_0_3(Y(i,:),T(i));
%    PO2IC_0_4_plot(i) = PO2IC_HighO2_0_4(Y(i,:),T(i));
%    PO2IC_0_5_plot(i) = PO2IC_LowO2_0_5(Y(i,:),T(i));
%    PO2IC_0_6_plot(i) = PO2IC_HighO2_0_6(Y(i,:),T(i));
%    PO2IC_0_7_plot(i) = PO2IC_LowO2_0_7(Y(i,:),T(i));
%    PO2IC_0_8_plot(i) = PO2IC_HighO2_0_8(Y(i,:),T(i));
%    PO2IC_0_9_plot(i) = PO2IC_LowO2_0_9(Y(i,:),T(i));
%end
%plot(T,PO2IC_0_1_plot,'LineWidth',2.0);
%hold on;
%plot(T,PO2IC_0_2_plot,'LineWidth',2.0);
%plot(T,PO2IC_0_3_plot,'LineWidth',2.0);
%plot(T,PO2IC_0_4_plot,'LineWidth',2.0);
%plot(T,PO2IC_0_5_plot,'LineWidth',2.0);
%plot(T,PO2IC_0_6_plot,'--','LineWidth',2.0);
%plot(T,PO2IC_0_7_plot,'--','LineWidth',2.0);
%plot(T,PO2IC_0_8_plot,'--','LineWidth',2.0);
%plot(T,PO2IC_0_9_plot,'--','LineWidth',2.0);
%hold off;
%set(gca,'XDir','reverse')
%ylabel ('relative Y,all effect considered includes biological speciation,C=0.5' );
%xlabel ('Ma ');
%legend('R1','R2','R3','R4','R5','R6','R7','R8','R9')
%
%figure(10)
%for i = 1: size(T,1)
%  PO2IC_1_10_plot(i) = PO2IC_default_1_10(Y(i,:),T(i));
%    PO2IC_1_11_plot(i) = PO2IC_HighO2_1_11(Y(i,:),T(i));
%    PO2IC_1_12_plot(i) = PO2IC_LowO2_1_12(Y(i,:),T(i));
%    PO2IC_1_13_plot(i) = PO2IC_HighO2_1_13(Y(i,:),T(i));
%    PO2IC_1_14_plot(i) = PO2IC_LowO2_1_14(Y(i,:),T(i));
%    PO2IC_1_15_plot(i) = PO2IC_HighO2_1_15(Y(i,:),T(i));
%    PO2IC_1_16_plot(i) = PO2IC_LowO2_1_16(Y(i,:),T(i));
%    PO2IC_1_17_plot(i) = PO2IC_HighO2_1_17(Y(i,:),T(i));
%    PO2IC_1_18_plot(i) = PO2IC_LowO2_1_18(Y(i,:),T(i));
%end 
%plot(T,PO2IC_1_10_plot,'LineWidth',2.0);
%hold on;
%plot(T,PO2IC_1_11_plot,'LineWidth',2.0);
%plot(T,PO2IC_1_12_plot,'LineWidth',2.0);
%plot(T,PO2IC_1_13_plot,'LineWidth',2.0);
%plot(T,PO2IC_1_14_plot,'LineWidth',2.0);
%plot(T,PO2IC_1_15_plot,'LineWidth',2.0);
%plot(T,PO2IC_1_16_plot,'LineWidth',2.0);
%plot(T,PO2IC_1_17_plot,'LineWidth',2.0);
%plot(T,PO2IC_1_18_plot,'LineWidth',2.0);
%hold off;
%set(gca,'XDir','reverse')
%ylabel ('relative Y,all effect considered includes biological speciation,C=0.9' );
%xlabel ('Ma ');
%legend('R10','R11','R12','R13','R14','R15','R16','R17','R18')