FOLDERS="output_500_thres_0.75_B_10_d_2 output_500_thres_0.75_B_9_d_2 output_population_500_thres_0.75_B_10_d_2 output_population_500_thres_0.75_B_10_d_2.2 output_population_500_thres_0.75_B_9_d_2"

for FOLDER in $FOLDERS; do
    echo "folder: ${FOLDER}"
    mkdir ${FOLDER}
    scp brutus:BSSE/StochasticEntrainment/VanDerPol/sensitivity/${FOLDER}/VanDerPol*.mat ${FOLDER}
done

