#
# Figure 4 viscosity, surface tension tests on a sphere.
#

mkdir -p fig4a_mu1 fig4a_muinv fig4b_mu1 fig4b_muinv
cd fig4a_mu1
../../growth --param ../mu1_st0.txt --niter 68 --step 1 --id 0 >log.txt 2>&1
cd ..
cd fig4a_muinv
../../growth --param ../mu100_inv_st0.txt --niter 68 --step 1 --id 0 >log.txt 2>&1
cd ..
cd fig4b_mu1
../../growth --param ../mu1_st1.txt --niter 68 --step 1 --id 0 >log.txt 2>&1
cd ..
cd fig4b_muinv
../../growth --param ../mu100_inv_st1.txt --niter 68 --step 1 --id 0 >log.txt 2>&1
cd ..

# With P2 for reference (very slow):
# mkdir -p fig4a_muinv_p2 fig4b_muinv_p2
# cd fig4a_muinv_p2
# ../../growth --param ../mu100_inv_st0_p2.txt --niter 68 --step 1 --id 0 >log.txt 2>&1
# cd ..
# cd fig4b_muinv_p2
# ../../growth --param ../mu100_inv_st1_p2.txt --niter 68 --step 1 --id 0 >log.txt 2>&1
# cd ..
