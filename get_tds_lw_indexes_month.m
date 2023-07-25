function get_tds_lw_indexes_month(year,month,type)
% get_tds_lw_indexes_month(year, month)
% this function calls the get_tds_lw_indexes function for each day of month
% and runs it in the mode which saves the days into the database.

type = "tswf"; %#ok<NASGU> % TODO: add other types 
for day = 1:eomday(year, month)
    % interferences
    if year == 2021
        if month == 1
            if day>23
                continue
            end
        elseif month == 2
            continue
        elseif month == 3
            if day<13
                continue
            end
        end

    end
    cdf = tdscdf_load_l2_surv_tswf(datenum(year,month,day));
    if isempty(cdf)
        continue
    end
    get_tds_lw_indexes(cdf, type, 1)
end