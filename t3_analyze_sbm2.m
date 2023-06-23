function t3_analyze_sbm2(filename, opts)
% t3_analyze_sbm2 This function performs analysis on the wave polarisation
% during a type III event recorded in SBM2 mode by SolO
%   Use:
%    Important opts settings:
%       
%    Other opts settings:
%       hrs - how many hours to include in the plot/analysis
%
%   Requirements: MAG data


% List of figures: 
% - TNR + EPD + wave polarisation + angle between B and the Y-Z plane
% - wave polarisation vs. beam energy (color ~ quality ~ angle B vs. Y-Z)


% Loading data
tsbm = tdscdf_load_l2_surv_tswf(filename);
ep0 = tsbm.epoch(1);
rswf = tdscdf_load_l2_surv_rswf(ep0, 1);


% Setting parameters
if exist("opts") %#ok<EXIST>
    if isfield(opts, 'hrs')
        tlen = opts.hrs*3600;
    else
        tlen = 2.5*3600;
    end

    if isfield(opts, 'tnr_comp')
        tnr_comp = opts.tnr_comp;
    else
        tnr_comp = 2;
    end
end


% Detect Langmuir waves
[~,~,recs] = size(tsbm.data);
for rec = 1:recs

% Determine the Langmuir wave polarization



% Subplot setup
tsubplot = @(m,n,p) subtightplot (m,n,p,[0.06 0.06], [0.05 0.07], [0.07 0.07]);
figure(1); clf;
k = 1; nfig = 1; 






end