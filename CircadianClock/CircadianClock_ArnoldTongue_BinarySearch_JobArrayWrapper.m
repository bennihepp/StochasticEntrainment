function CircadianClock_ArnoldTongue_BinarySearch_JobArrayWrapper(n, filename_prefix, Ntrials, population_average)

	try
		CircadianClock_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, Ntrials, population_average);
	catch e
		display('pause');
		pause(5);
		try
			CircadianClock_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, Ntrials, population_average);
		catch e
			display('pause');
			pause(5);
			try
				CircadianClock_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, Ntrials, population_average);
			catch e
				display('pause');
				pause(5);
				CircadianClock_ArnoldTongue_BinarySearch_JobArray(n, filename_prefix, Ntrials, population_average);
			end
		end
	end


end

