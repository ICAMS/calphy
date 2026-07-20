
cd example_03
echo "testing TS solid"
rm -rf ts_*
calphy_kernel -i input.yaml -k 0 -s True
echo "testing TS liquid"
calphy_kernel -i input.yaml -k 1 -s True
cd ..

cd example_04
rm -rf ts-*
echo "testing Automated Tm"
calphy_kernel -i input.1.yaml -k 0 -s True
cd ..

cd example_07
rm -rf alchemy-*
echo "testing alchemy"
calphy_kernel -i input.1.yaml -k 0 -s True
cd ..

cd example_08
rm -rf tscale-*
echo "testing tscale"
calphy_kernel -i input.1.yaml -k 1 -s True
cd ..

cd example_09
rm -rf pscale-*
echo "testing pscale"
calphy_kernel -i input2.yaml -k 0 -s True
cd ..

cd example_10
rm -rf composition-*
echo "testing composition scaling"
calphy_kernel -i input-composition.yaml -k 0 -s True
cd ..

