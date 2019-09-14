#
# Figure 5 growth tests with hexagonal pattern on an ellipsoid.
#

mkdir -p mu1_st0 mu100_inv_st0
cd mu1_st0
../../growth --param ../mu1_st0.txt --niter 72 --step 1 --id 0 >log.txt 2>&1
cd ..
cd mu100_inv_st0
../../growth --param ../mu100_inv_st0.txt --niter 72 --step 1 --id 0 >log.txt 2>&1
cd ..

# With P2 for reference (very slow):
# mkdir -p mu1_st0_p2
# cd mu1_st0_p2
# ../../growth --param ../mu1_st0_p2.txt --niter 72 --step 1 --id 0 >log.txt 2>&1
# cd ..