function t3m_stereo_spectrogram(ster)
%T3M_STEREO_SPECTROGRAM Plots all snapshots after each other like triggered
% snapshot spectrograms

for i = 1:size(ster.channel,2)
    [sp,fq,nav] = make_spectrum(ster.data(:,1,i),length(ster.data(:,1,i)),1/ster.samp_rate(i),ster.samp_rate(i)/2);
    sps(:,i) = sp;
end

imagesc(1:79,fq*1e-3,log10(sps))
set(gca, 'YDir', 'normal')
title(sprintf('STEREO %c swaves TDS triggered snapshots spectrogram from %s',upper(ster.sc),datestr(ster.ep(1),'yyyy-mm-dd')))
ylabel('Frequency [kHz]')
xlabel('Number of snapshot')

end

