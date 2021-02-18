%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code quantifies the flux of dissolved and gaseous N2O from the
% TEMMS floodplain.
% Written by Corey Wallace on December 21, 2020.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all; close all; clc;

%% important parameters
N2O_weight = 44.013; % g/mol
Henrys_Constant = 2.471e-10; % mol/Pa.L
Gas_Constant = 8.31446261815324e3; % Pa.L/mol.K
Air_Water_Partition = Henrys_Constant*Gas_Constant*298.15;
K0 = 2.5483e-2; % mol/L equilibrium constant at 1 atm and 25C
mole_frac = 3.2622e-7; % based on 326.22 ppb NOAA value
Temp = 298.15; % absolute pressure 25C = 298.15 K
Tot_Press = 1; % total pressure = 1 atm
% Air_Water_Partition = Aqueous Conc. / Gas Conc.

%% Calculate equilibrium N2O concentration in water from Weiss and Price 1980

pH2O = exp(24.4543 - 67.4509*(100/Temp) - 4.8489*log(Temp/100));
F = K0*(Tot_Press - pH2O)*exp(-9.4563/Temp + 0.04739 - 6.427e-05*Temp);
C_eq = mole_frac*F*1000; % mol/m3 in water at 1 atm and 25C

%% Calculate Schmidt Number from Raymond et al., 2012
Kin_Vis = 1.735e-02 + (-5.023e-04*(Temp-273.15)) + (8.598e-06*(Temp-273.15)^2)...
    + (6.805e-08*(Temp-273.15)^3);
Sc_N2O = 2056 + (-137.11*(Temp-273.15)) + (4.317*(Temp-273.15)^2)...
    + (-.0543*(Temp-273.15)^3);
Sc_O2 = 1801 + (-120.1*(Temp-273.15)^2) + (3.782*(Temp-273.15)^2)...
    + (-0.0476*(Temp-273.15)^3);
kO2 = 5; % m/day

% scale kN2O by kO2 measured in other rivers
kN2O = (Sc_N2O/Sc_O2)^-0.5 * kO2; % m/day

%% import mass balance file from PFLOTRAN
readfile = sprintf('TEMMS_chem_R3-mas.dat');
db = importdata(readfile);
header = strsplit(db{1,1},",");
for i = 2:size(db,1)
    db_data(i,:) = strsplit(db{i,1}," ");
end
data = cell2table(db_data(2:end,2:end));
db = str2double(data{:,:});

% differentiate necessary columns
time = db(:,1);
dt(1,1) = time(1,1);
j = 2;
for i = 1:size(time,1)-1
    dt(j,1) = time(i+1,1) - time(i,1);
    j = j + 1;
end
% dt = db(:,2);
% bank_volume = db(:,53); % to compute conc. of N2O -- not needed
bank_NO3_aq = db(:,75); % to check sign convention
bank_N2O_aq = db(:,78);
bank_N2O_g = db(:,79);
top_N2O_g = db(:,111);
gw_boundary_flow = db(:,21);

%% read in stage file used for visualization
readfile = sprintf('stilling_well_stage.txt');
stage_file = importdata(readfile);
stage_time = stage_file.data(:,1)/86400;
stage = (stage_file.data(:,2)-101325)/9806;

%% Compute gas flux across top and bank
top_flux = top_N2O_g.*dt; % gives moles
bank_flux = bank_N2O_g.*dt; % gives moles

% convert from mol to g
% moles x (g/mol) = grams
top_flux = top_flux*N2O_weight; 
bank_flux = bank_flux*N2O_weight;

%% Compute dissolved flux across bank
bank_aq_flux = bank_N2O_aq.*dt; % gives moles
bank_NO3 = bank_NO3_aq.*dt*62.001; % g of nitrate

% compute outgassing based on air-water partition
bank_gas_partition = bank_aq_flux/Air_Water_Partition;

% comvert from mol to g
bank_gas_partition = bank_gas_partition*N2O_weight;

%% Compute flow across the GW boundary
gw_flow = gw_boundary_flow.*dt;

%% Compute the integral of the flux for each day
for j = 1:1
    end_index = find(time==j);
    if j == 1
        bank_gas = bank_flux(1:end_index);
        top_gas = top_flux(1:end_index);
        diss_flux = bank_gas_partition(1:end_index);
        total_bank_gas_flux(j) = trapz(time(1:end_index),bank_gas);
        total_top_gas_flux(j) = trapz(time(1:end_index),top_gas);
        total_diss_flux(j) = trapz(time(1:end_index),diss_flux);
    else
        start_index = find(time==j-1);
        bank_gas = bank_flux(start_index:end_index);
        top_gas = top_flux(start_index:end_index);
        diss_flux = bank_gas_partition(start_index:end_index);
        total_bank_gas_flux(j) = trapz(time(start_index:end_index),bank_gas);
        total_top_gas_flux(j) = trapz(time(start_index:end_index),top_gas);
        total_diss_flux(j) = trapz(time(start_index:end_index),diss_flux);
    end
    total_gas_flux(j) = total_bank_gas_flux(j) + total_top_gas_flux(j);
    total_flux(j) = total_gas_flux(j) + total_diss_flux(j);
end

% this output is in grams. multiply by 1000 to get mg, or divide by 1000 to
% get kg
gas_flux = sum(total_gas_flux*1000) % milligrams
top_gas_flux = sum(total_top_gas_flux) % grams
N2O_flux = sum(total_flux/1000) % kilograms


%% Create plots of flux for visualization
figure
subplot(4,1,1)
plot(stage_time,stage,'k-');
% ylim([0 2]);
% xlim([0 1]);
ylabel('Stage');
subplot(4,1,2)
plot(time,bank_NO3,'k-');
ylabel('Nitrate (g)');
subplot(4,1,3)
plot(time,bank_gas_partition,'k-');
ylabel('Dissolved N_2O (g)');
subplot(4,1,4)
plot(time,bank_flux,'k-');
hold on
plot(time,top_flux,'b-');
ylabel('Gaseous N_2O (g)');
legend('Bank','Top');
% xticklabels({'10/18','10/19','10/20','10/21','10/22','10/23','10/24',...
%     '10/25','10/26','10/27'});
% xticks(0,1,2,3,4);
% xticklabels({'10/18','10/19','10/20','10/21','10/22'});
xlabel('Date');

x = 1:1:10;
figure
subplot(4,1,1)
plot(stage_time,stage,'k-');
ylabel('Stage');
subplot(4,1,2)
plot(x,total_gas_flux,'k-');
ylabel('Gas Flux (g)');
subplot(4,1,3)
plot(x,total_diss_flux/1000,'k-');
ylabel('Diss. Flux (kg)');
subplot(4,1,4)
plot(x,total_flux/1000,'k-');
ylabel('Total Flux (kg)');
% figure
% plot(time,gw_flow);

