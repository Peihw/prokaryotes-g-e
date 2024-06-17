configfile: "config_s.yaml"


rule all:
	input:
		"result/duplication/{sample}.info"


##Assembly pipeline

rule data_clean:
	input:
		fq1=expand("{fq1}",fq1=config["Fastq1"]),
		fq2=expand("{fq2}",fq2=config["Fastq2"]),
	output:
		fq1=temp("result/cleandata/{sample}_1.fq"),
		fq2=temp("result/cleandata/{sample}_2.fq"),
		m1="result/cleandata/{sample}_1.fastq",
		m2="result/cleandata/{sample}_2.fastq"
	shell:
		"""
		./bin/prinseq/prinseq-lite -verbose -fastq {input.fq1} -ns_max_p 10 -min_qual_mean 25 -out_good result/cleandata/{wildcards.sample}_1_good  -out_bad  result/cleandata/{wildcards.sample}_1_bad
		./bin/prinseq/prinseq-lite -verbose -fastq {input.fq2} -ns_max_p 10 -min_qual_mean 25 -out_good result/cleandata/{wildcards.sample}_2_good -out_bad result/cleandata/{wildcards.sample}_2_bad
		awk '/^@/ {{ print $1; n = "" }} !/^@/ {{ print $0}} END {{ printf "%s", n }}' result/cleandata/{wildcards.sample}_1_good.fastq | sed  '/@/ s/$/\/1/'  > {output.fq1}
		awk '/^@/ {{ print $1; n = "" }} !/^@/ {{ print $0}} END {{ printf "%s", n }}' result/cleandata/{wildcards.sample}_2_good.fastq | sed  '/@/ s/$/\/2/'  > {output.fq2}
		./bin/manageread extractpair pairfastq={output.fq1},{output.fq2} outprefix=result/cleandata/{wildcards.sample}
		"""


rule Spades:
	input:
		m1="result/cleandata/{sample}_1.fastq",
		m2="result/cleandata/{sample}_2.fastq"
	output:
		contigs="result/spades_out/{sample}/contigs.fasta",
		scaffolds="result/spades_out/{sample}/scaffolds.fasta"
	params:
		"result/spades_out/{sample}"
	shell:
		"spades.py -1 {input.m1} -2 {input.m2} --phred-offset 33 -o {params} -t 30 -m 1024"


rule Prokka:
	input:
		expand("{fa}",fa=config["Fasta"])
	output:
		"result/prokka_out/{sample}/{sample}.fna",
		"result/prokka_out/{sample}/{sample}.gff",
		"result/prokka_out/{sample}/{sample}.faa",
	params:
		dir="result/prokka_out/{sample}",
		gns=expand("{gns}",gns=config["Genus"]),
		spe=expand("{spe}",spe=config["Species"])
	shell:
		"""
		prokka --compliant --genus {params.gns} --species {params.spe} --cpus 20 --outdir {params.dir} --force --prefix {wildcards.sample}  {input}
		sed -i 's/|/_/g' result/prokka_out/{wildcards.sample}/*
		"""

	

##Signal annotation pipeline

rule data_reformat:
	input:
		gff="result/prokka_out/{sample}/{sample}.gff",
	output:
		ft="result/prokka_out/{sample}/{sample}.ft",
	shell:
		"""
		awk -F'[\\t=;]' -v name={wildcards.sample} '!/^#/{{print$3"\\t""\\t"name.1"\\t""\\t""\\t""\\t"$1"\\t"$4"\\t"$5"\\t"$7"\\t"$16"\\t""\\t""\\t"$24"\\t""\\t"$14"\\t"$10"\\t""\\t""\\t"0}}' {input.gff} > {output.ft} 
		sed -i '1i#' {output.ft}
		sed -i '/^	/d' {output.ft}
		"""

##Tritisa with scaffolds

rule tritisa:
	input:
		ft="result/prokka_out/{sample}/{sample}.ft",
		genome="result/prokka_out/{sample}/{sample}.fna"
	output:
		list=temp("result/prokka_out/{sample}/{sample}.list"),
		tritisa="result/Tritisa/{sample}.tritisa.rec.dat"
	params:
		dir1="result/prepare/fna/",
		dir2="result/prepare/med/",
		dir3="result/Tritisa/{sample}"
	shell:
		"""
		./bin/protisa/codes/Prokaryotes -o {params.dir1} -f {input.genome} -N
		awk -F'\\t' '{{print$7}}' {input.ft} |sort|uniq > {output.list}
		while read line; do ./bin/protisa/codes/Prokaryotes -o {params.dir2} -nc $line -a {input.ft} -M; done < {output.list}
		./tritisa.sh {output.list}
		while read line; do ./bin/protisa/codes/Prokaryotes -o {params.dir3} -nc $line -m prepare/Tritisa/$line.tritisa.dat -a {input.ft} -T; done < {output.list} 
		rm -rf result/prepare/fna/
		rm -rf result/prepare/med/
		"""	
		

rule utr_sequence:
	input:
		rec="result/Tritisa/{sample}.tritisa.rec.dat",
		genome="result/prokka_out/{sample}/{sample}.fna"
	output:
		"result/TISseq/{sample}.tis.fa"
	params:
		dir="result/TISseq/",
		len=expand("{length}",length=config["Length"])
	shell:
		"./bin/extract/extract -o {params.dir} -g {wildcards.sample} -p tis -r {input.rec} -n {input.genome} -l {params.len} -R"
	

rule signal_identification:
	input:
		TISseq="result/TISseq/{sample}.tis.fa",
		genome="result/prokka_out/{sample}/{sample}.fna"
	output:
		"result/signal/model/{sample}.sig.dat"
	params:
		dir="result/signal/model/"
	shell:
		"./bin/TIS/codes/TIS -o {params.dir} -g {wildcards.sample} -f {input.TISseq} -pb {input.genome} -mj 0.52 -mn 0.16 -ep 1.0e-4 -tr 3 -k 9 -r 30 -m 3 -w 8 -M"


rule signal_classification:
	input:
		TISseq="result/TISseq/{sample}.tis.fa",
		model="result/signal/model/{sample}.sig.dat"
	output:
		sig="result/signal/signal/{sample}.sig",
	params:
		ref=expand("{ref}",ref=config["Reference"]),
		dir="result/signal/model/"
	shell:
		"./bin/FP/codes/TIS -o {output.sig} -g {wildcards.sample} -r {params.ref} -q {params.dir} -sd 6 -ta 0.75 -l 5 -C"


rule signal_scanning:
	input:
		tritisa="result/Tritisa/{sample}.tritisa.rec.dat",
		model="result/signal/model/{sample}.sig.dat",
		sig="result/signal/signal/{sample}.sig",
		TISseq="result/TISseq/{sample}.tis.fa"
	output:
		"result/records/{sample}.tis.rec.dat"
	params:
		dir="result/records/"
	shell:
		"./bin/FP/codes/TIS -o {params.dir} -g {wildcards.sample} -rd result/Tritisa/ -q result/signal/model/ -f result/TISseq/ -d {input.sig} -R"


## Promoters detection

rule longer_utr_sequence:
	input:
		rec="result/Tritisa/{sample}.tritisa.rec.dat",
		genome="result/prokka_out/{sample}/{sample}.fna",
		dup="result/duplication/{sample}.dup"
	output:
		seq1="result/TSseq/{sample}.tis.fa",
		seq2="result/TSseq/{sample}.dup.tis.fa",
		list="result/TSseq/{sample}.list"
	params:
		dir="result/TSseq/"
	shell:
		"""
		./bin/extract/extract -o {params.dir} -g {wildcards.sample} -p tis -r {input.rec} -n {input.genome} -l 400 -R
		awk '/^>/&&NR>1 {{print "";}}{{ printf "%s",/^>/ ? $0"%":$0,n="\\n" }}END {{ printf "%s", n }}' {output.seq1} > result/TSseq/{wildcards.sample}.c
		awk -F'[|\\t]' 'NR==FNR{{a[$3]=$0;}}NR!=FNR && a[$2] {{print a[$2]}}' result/TSseq/{wildcards.sample}.c {input.dup} |tr -s "%" "\n" > {output.seq2}
		mkdir -p result/TSseq/{wildcards.sample}.seq
		split -l 2 {output.seq2} result/TSseq/{wildcards.sample}.seq/{wildcards.sample}
		ls -1 result/TSseq/{wildcards.sample}.seq/ > {output.list}
		rm -f result/TSseq/{wildcards.sample}.c
		"""


rule promoters:
	input:
		"result/TSseq/{sample}.list"
	output:
		"result/promoter/{sample}.mp.proseq"
	shell:
		"""
		cd bin/promotech/
		while read line; do python promotech.py -pg -m RF-HOT -f ../../result/TSseq/{wildcards.sample}.seq/$line -o ../../result/promoter/{wildcards.sample}/$line; python promotech.py -g -t 0.6 -i ../../result/promoter/{wildcards.sample}/$line -o ../../result/promoter/{wildcards.sample}/$line; awk -F'\\t' '$5=="+"{{print$0}}' ../../result/promoter/{wildcards.sample}/$line/genome_predictions.csv >> ../../result/promoter/{wildcards.sample}.all.proseq; awk -F'\\t' '$5=="+"{{print$0}}' ../../result/promoter/{wildcards.sample}/$line/genome_predictions.csv | awk -F'\\t'  'BEGIN {{max = 0}} {{if ($4+0 > max+0) {{max=$4 ;content=$0}} }} END {{print content}}' >> ../../result/promoter/{wildcards.sample}.mp.proseq; done < ../../{input}
		cd ../../
		rm -rf result/promoter/{wildcards.sample}
		rm -rf result/TSseq/{wildcards.sample}.seq
		"""


#Gene duplication pipeline

rule duplication_1:
	input:
		cds="result/prokka_out/{sample}/{sample}.faa",
	output:
		clust1="result/duplication/{sample}.clust",
		clust2=temp("result/duplication/{sample}.c2"),
		vfca="result/duplication/{sample}.vfca"
	params:
		idttfd=expand("{idttfd}",idttfd=config["identity"]),
		coverfd=expand("{coverfd}",coverfd=config["coverage"]),
		qid=expand("{qid}",qid=config["query_id"]),
		qc=expand("{qc}",qc=config["query_cover"]),
	shell:
		"""
		./bin/blast/bin/blastclust -i {input.cds} -o {output.clust1} -e F -S {params.idttfd} -L {params.coverfd} -b T
		awk '{{if (NF >1) print $0}}' {output.clust1} | sed 's/ /\\n/g' > {output.clust2}
		diamond blastp --db bin/db/VFDB --query {input.cds} --out result/duplication/temp.vf1 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-5 --id {params.qid} --query-cover {params.qc} --index-chunks 1
		diamond blastp --db bin/db/CARD --query {input.cds} --out result/duplication/temp.ca1 --outfmt 6 --sensitive --max-target-seqs 1 --evalue 1e-5 --id {params.qid} --query-cover {params.qc} --index-chunks 1
		cat result/duplication/temp.vf1 result/duplication/temp.ca1 |awk -F'[\\t@]' '{{print$1"\\t"$3}}' > result/duplication/temp.vfca
		awk -F'\\t' 'NR==FNR{{a[$1]=$2;}}NR!=FNR && a[$1] {{print $0"\\t"a[$1]}}' result/duplication/temp.vfca {output.clust2} > {output.vfca}
		rm -f  result/duplication/temp.*
		"""


rule duplication_2:
	input:
		"result/duplication/{sample}.c2"
	output:
		"result/duplication/{sample}.dup"
	run:
		with open(output[0], "w") as fileout:
			with open(input[0], "r") as filein:
				newmol = True
				ele = 1
				for line in filein:
					newmol = True
					if line == '\n':
						ele = ele + 1
						newmol = False
					if newmol == True:
						fileout.write(str(ele) + '\t' + str(line))	


rule integrate_tis:
	input:
		dup="result/duplication/{sample}.dup",
		tis="result/records/{sample}.tis.rec.dat",
		vfca="result/duplication/{sample}.vfca"
	output:
		fin="result/duplication/{sample}.dup.tis"
	shell:
		"""
		cut -f 1-10,17- {input.tis} >  result/duplication/{wildcards.sample}.tis
		awk -F'\\t' 'NR==FNR{{a[$1]=$2;}}NR!=FNR {{print $0"\\t"a[$2]}}'  {input.vfca} {input.dup} >  result/duplication/{wildcards.sample}.dup.vfca
		awk -F'\\t' 'NR==FNR{{a[$3]=$0;}}NR!=FNR && a[$2] {{print $1"\\t"$3"\\t"a[$2]}}' result/duplication/{wildcards.sample}.tis result/duplication/{wildcards.sample}.dup.vfca > {output.fin}
		header=`head -n 1 result/duplication/{wildcards.sample}.tis`
		sed -i '1i\\Cluster Number\tVF/AMR\t'"$header"'' {output.fin}
		rm -f result/duplication/{wildcards.sample}.tis result/duplication/{wildcards.sample}.dup.vfca
		"""


rule integrate_pro:
	input:
		pro="result/promoter/{sample}.mp.proseq",
		dup="result/duplication/{sample}.dup",
	output:
		fin="result/duplication/{sample}.dup.pro"
	shell:
		"""
		awk -F'[|\\t]' 'NR==FNR{{a[$3]=$3"\\t"$0;}}NR!=FNR && a[$2] {{print $1"\\t"a[$2]}}' {input.pro} {input.dup} > {output.fin}
		sed -i '1i\\Cluster Number\tLocus\tChrom\tStart\tEnd\tScore\tStrand\tPromoter sequence' {output.fin}
		"""


rule integrate_all:
	input:
		pro="result/promoter/{sample}.mp.proseq",
		tis="result/duplication/{sample}.dup.tis"
	output:
		"result/duplication/{sample}.info"
	shell:
		"""
		awk -F'[|\\t]' 'NR==FNR{{a[$3]=$NF;}}NR!=FNR{{print $0"\\t"a[$5]}}' {input.pro} {input.tis} > {output}
		header=`head -n 1 {input.tis}`
		sed -i '1i\\'"$header"'\tPromoter sequence' {output}
		"""
