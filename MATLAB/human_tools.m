EUWS_DateRange = daterange(EUWS);
pythonArray = arrayGen(EUWS_DateRange);

EUWS_DateRange.duration = EUWS_DateRange.endTime - EUWS_DateRange.startTime;
totaltime = sum(years(EUWS_DateRange.duration));

function EUWS_DateRange = daterange(EUWS)

EUWS_DateRange = table(EUWS.time, EUWS.time - hours(36), EUWS.time + hours(36), ...
    'VariableNames', {'time', 'startTime', 'endTime'});

EUWS_DateRange = sortrows(EUWS_DateRange, 'startTime');

% Initialize output table
EUWS_filtered = EUWS_DateRange(1, :);  % Copy the first row

% Loop through the rest of the rows
for i = 2:height(EUWS_DateRange)
    lastRow = EUWS_filtered(end, :);

    % Check if the new row is entirely within the last row in the output table
    if EUWS_DateRange.startTime(i) >= lastRow.startTime && EUWS_DateRange.endTime(i) <= lastRow.endTime
        continue;  % Skip this row, it's entirely within the previous one
    end

    % Check if the new row overlaps with the last row in the output table
    if EUWS_DateRange.startTime(i) <= lastRow.endTime
        % Combine the two timeframes and update the last row in the output table
        EUWS_filtered.endTime(end) = max(EUWS_DateRange.endTime(i), lastRow.endTime);
    else
        % If the new row does not overlap with the last row in the output table, add it to the table
        EUWS_filtered = [EUWS_filtered; EUWS_DateRange(i, :)];
    end
end

EUWS_DateRange = EUWS_filtered;

end


function pythonArray = arrayGen(EUWS_DateRange)

EUWS_DateRange = EUWS_DateRange(1:end-4, :);

% Split date string into parts
[startYear, startMonth, startDay] = datevec(EUWS_DateRange.startTime);
[endYear, endMonth, endDay] = datevec(EUWS_DateRange.endTime);

% Generate input strings for Python
for i = 1:height(EUWS_DateRange)
    year        = startYear(i):endYear(i);
    yearStr     = sprintf('''%04d'',', year);
    yearStr     = yearStr(1:end-1);

    if startMonth(i) <= endMonth(i)
        month = startMonth(i):endMonth(i);
    else
        month = [startMonth(i):12, 1:endMonth(i)];
    end
    monthStr    = sprintf('''%02d'',', month);
    monthStr    = monthStr(1:end-1);

    if startDay(i) <= endDay(i)
        days = startDay(i):endDay(i);
    else
        days1 = startDay(i):eomday(startYear(i), startMonth(i));
        days2 = 1:endDay(i);
        days = [days1, days2];
    end
    dayStr      = sprintf('''%02d'',', days);
    dayStr      = dayStr(1:end-1);

    pythonArray{i} = sprintf('{''year'': [%s], ''month'': [%s], ''day'': [%s]},\n', yearStr, monthStr, dayStr);
end

end
