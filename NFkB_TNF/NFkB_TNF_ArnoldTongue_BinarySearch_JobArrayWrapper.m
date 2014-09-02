function NFkB_TNF_ArnoldTongue_BinarySearch_JobArrayWrapper(n, filename_prefix, ...
    volume, input_periods, input_amplitude_tolerance,Ntrials, population_average)

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

end
