function t3m_analyze_aug_18_2010()
% T3M_ANALYZE_AUG_18_2010 makes several plots from stereo data from august
% 18 2010.

ster = load_stereo_tds(datenum(2010,8,18),'a');
figure(1)
t3m_stereo_spectrogram(ster);

%figure(2)

%subplot(3,2,1)
%nsamps = ster.n_samps(26);
%plot((1e3/ster.samp_rate(26)).*(1:nsamps),ster(1:nsamps,1,26))

end

