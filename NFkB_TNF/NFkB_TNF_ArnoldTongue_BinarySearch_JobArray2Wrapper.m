function NFkB_TNF_ArnoldTongue_BinarySearch_JobArray2Wrapper(n, filename_prefix, ...
    volume, input_periods, input_amplitude_tolerance, Ntrials, population_average, entrainment_threshold)

    if nargin < 8
        entrainment_threshold = 0.9;
    end

    NFkB_TNF_ArnoldTongue_BinarySearch_JobArray2(n, filename_prefix, volume, input_periods, input_amplitude_tolerance,Ntrials, population_average, entrainment_threshold);

%     try
% 		NFkB_TNF_ArnoldTongue_BinarySearch_JobArray2(n, filename_prefix, volume, input_periods, input_amplitude_tolerance,Ntrials, population_average);
% 	catch e
% 		display('pause');
% 		pause(5);
% 		try
% 			NFkB_TNF_ArnoldTongue_BinarySearch_JobArray2(n, filename_prefix, volume, input_periods, input_amplitude_tolerance,Ntrials, population_average);
% 		catch e
% 			display('pause');
% 			pause(5);
% 			try
% 				NFkB_TNF_ArnoldTongue_BinarySearch_JobArray2(n, filename_prefix, volume, input_periods, input_amplitude_tolerance,Ntrials, population_average);
% 			catch e
% 				display('pause');
% 				pause(5);
% 				NFkB_TNF_ArnoldTongue_BinarySearch_JobArray2(n, filename_prefix, volume, input_periods, input_amplitude_tolerance,Ntrials, population_average);
% 			end
%         end
%     end

end
