#
# Figure 6 growth tests; hexagonal pattern on an ellipsoid with direct viscosity profiles.
#

mkdir -p mu100_st0 mu1000_st0 mu1000_rate4x_st0
cd mu100_st0
../../growth --param ../mu100_st0.txt --niter 72 --step 1 --id 0 >log.txt 2>&1
cd ..
cd mu1000_st0
../../growth --param ../mu1000_st0.txt --niter 72 --step 1 --id 0 >log.txt 2>&1
cd ..
cd mu1000_rate4x_st0
../../growth --param ../mu1000_rate4x_st0.txt --niter 72 --step 1 --id 0 >log.txt 2>&1
cd ..
