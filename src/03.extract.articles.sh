# deprecated: use zgrep output instead
cd /lustre/scratch117/cellgen/cellgeni/pasham/data/2204.pmc.github/jsons

for i in `ls -1`
do
  echo $i
  rm -f ${i}/*xml
  cd ${i}
  tar -xzvf ../../tars/${i}_json_unicode.tar.gz --files-from files.txt
  cd ..
done
