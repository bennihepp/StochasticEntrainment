FOLDERS="output_inf_thres_0.75_B_10_d_1.8 output_inf_thres_0.75_B_10_d_2 output_inf_thres_0.75_B_10_d_2.2 output_inf_thres_0.75_B_11_d_2 output_inf_thres_0.75_B_9_d_2"

for FOLDER in $FOLDERS; do
    echo "folder: ${FOLDER}"
    mkdir ${FOLDER}
    scp euler:BSSE/StochasticEntrainment/VanDerPol/sensitivity/${FOLDER}/VanDerPol*.mat ${FOLDER}
done

