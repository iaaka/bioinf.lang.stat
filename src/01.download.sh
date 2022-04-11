# download archives
cd /lustre/scratch117/cellgen/cellgeni/pasham/data/2204.pmc.github
mkdir tars
cd tarr
wget https://ftp.ncbi.nlm.nih.gov/pub/wilbur/BioC-PMC/sum
for i in `grep json_unicode sum | cut -d ' ' -f3`; 
do 
  wget https://ftp.ncbi.nlm.nih.gov/pub/wilbur/BioC-PMC/$i; 
done

# check sums
for i in `ls -1 *gz`
do
  md5sum $i
done > check.sums

# find articles that mention github
mkdir ../ghmatches
for i in `ls -1 *gz | sed s/.tar.gz//`
do
  echo $i
  zgrep -a github ${i}.tar.gz > ../ghmatches/${i}.github.txt
done

cd ../ghmatches/
gzip ./*



