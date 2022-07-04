function ist3 = load_events(fname)
%LOAD_EVENTS Loads data from the event list containing informarion about in situ type III emissions recorded by SolO.
%   Detailed explanation below

if ~exist('fname','var') || isempty(fname)
    fname = ['events' filesep 'event_list.txt'];
end

T = readtable(fname);           
ist3.year = table2array(T(:,1));             % year of event
ist3.month = table2array(T(:,2));            % month of event
ist3.day = table2array(T(:,3));              % day of beginning of event
ist3.radio_hour = table2array(T(:,4));       % hour when the radio emmision gets detected by TDS
ist3.radio_minute = table2array(T(:,5));     % minute when the radio emmision gets detected by TDS
ist3.epd_nxtday = table2array(T(:,6));       % does the energetic beam arrive on the next day
ist3.epd_hour = table2array(T(:,7));         % hour at which EPD STEP detects an increase in flux
ist3.epd_minute = table2array(T(:,8));       % minute at which EPD STEP detects an increase in flux
ist3.tswf_nxtday = table2array(T(:,9));      % are the triggered snapshots with langmuir waves on the next day
ist3.tswf_indexes = table2array(T(:,10));    % indexes from the tswf cdf file which contain langmuir waves possibly caused by the energetic particles
ist3.tswf_frequency = table2array(T(:,11));  % frequency/frequencies of langmuir waves from the triggered snapshots
ist3.langmuir_nxtday = table2array(T(:,12)); % do the langmuir waves show up on regular snapshots first on the next day
ist3.langmuir_hour = table2array(T(:,13));   % hour when langmuir waves show up on regular snapshots
ist3.langmuir_minute = table2array(T(:,14)); % minute when langmuir waves show up on regular snapshots

end

