#
# Figure 3B, Table 2 Enright's tests.
#

mkdir -p 64 128 256
cd 64
../../leveque --param ../LeVeque_64.txt --niter 48 --step 1 --id 0 >log.txt 2>&1
cd ..
cd 128
../../leveque --param ../LeVeque_128.txt --niter 96 --step 1 --id 0 >log.txt 2>&1
cd ..
cd 256
../../leveque --param ../LeVeque_256.txt --niter 192 --step 1 --id 0 >log.txt 2>&1
cd ..
