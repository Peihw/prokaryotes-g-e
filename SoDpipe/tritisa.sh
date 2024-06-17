#bin/bash

mkdir -p result/prepare/Tritisa/
cd bin/TriTISA/codes/
while read line
	do
	./TriTISA ../../../result/prepare/fna/$line.fna ../../../result/prepare/med/$line.med ../../../result/prepare/Tritisa/$line
	done < ../../../$1
cd ../../../
