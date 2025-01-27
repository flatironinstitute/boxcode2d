rm -f print_testres.txt
./int2-test_l2dpwrouts
./int2-test_lrt2d
./int2-test_pbox2drouts
mv print_testres.txt ../print_testres.txt
rm -f fort.13
