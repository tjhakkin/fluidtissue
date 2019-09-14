#
# Figure 3A, Table 1 surface tension simulations on a cube.
#

mkdir -p 32 64 128 32_P2 64_P2 128_P2
cd 32
../../growth --param ../par_32.txt --niter 28 --step 1 --id 0 >log.txt 2>&1
cd ..
cd 64
../../growth --param ../par_64.txt --niter 56 --step 1 --id 0 >log.txt 2>&1
cd ..
cd 128
../../growth --param ../par_128.txt --niter 112 --step 1 --id 0 >log.txt 2>&1
cd ..
cd 32_P2
../../growth --param ../par_32_p2.txt --niter 28 --step 1 --id 0 >log.txt 2>&1
cd ..
cd 64_P2
../../growth --param ../par_64_p2.txt --niter 56 --step 1 --id 0 >log.txt 2>&1
cd ..
cd 128_P2
# Requires tons of memory!
# ../../growth --param ../par_128_p2.txt --niter 112 --step 1 --id 0 >log.txt 2>&1
cd ..

