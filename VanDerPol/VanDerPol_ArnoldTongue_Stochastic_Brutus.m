cluster = parcluster('BrutusLSF8h')
matlabpool(cluster, 128)
%cluster = parcluster('BrutusLSF36h')
%matlabpool(cluster, 16)

VanDerPol_ArnoldTongue_Stochastic;

matlabpool close
