cd meshes/
rm *.params
rm *.mesh

cd ..

./meshbuilder > log.txt

cd meshes/
cp *.params ../../Task/meshes/
cp *.mesh ../../Task/meshes/
cd ..