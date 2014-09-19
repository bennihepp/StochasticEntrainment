function VanDerPol_ArnoldTongue_BinarySearch_JobArray2Wrapper(n, filename_prefix, ...
    input_periods, input_amplitude_tolerance, Ntrials, population_average, entrainment_threshold)

    if nargin < 7
        entrainment_threshold = 0.9;
    end

    VanDerPol_ArnoldTongue_BinarySearch_JobArray2(n, filename_prefix, input_periods, input_amplitude_tolerance, Ntrials, population_average, entrainment_threshold);

end

