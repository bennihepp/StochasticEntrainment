function NFkB_TNF_ArnoldTongue_BinarySearch_JobArrayWrapper(n, filename_prefix, ...
    volume, input_periods, input_amplitude_tolerance,Ntrials, population_average)

    threads = 1;
    if ~isempty(getenv('NUM_OF_THREADS'))
        threads = str2double(getenv('NUM_OF_THREADS'));
    end
    if threads > 1
        matlabpool('local', threads);
    end

    try
		NFkB_TNF_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, volume, input_periods, input_amplitude_tolerance,Ntrials, population_average);
	catch e
		display('pause');
		pause(5);
		try
			NFkB_TNF_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, volume, input_periods, input_amplitude_tolerance,Ntrials, population_average);
		catch e
			display('pause');
			pause(5);
			try
				NFkB_TNF_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, volume, input_periods, input_amplitude_tolerance,Ntrials, population_average);
			catch e
				display('pause');
				pause(5);
				NFkB_TNF_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, volume, input_periods, input_amplitude_tolerance,Ntrials, population_average);
			end
        end
    end

    if threads > 1
        matlabpool close;
    end

end
