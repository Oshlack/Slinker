/*==================================================================================================

	     ___                                  ___          ___          ___          ___     
	    /  /\                    ___         /__/\        /__/|        /  /\        /  /\    
	   /  /:/_                  /  /\        \  \:\      |  |:|       /  /:/_      /  /::\   
	  /  /:/ /\   ___     ___  /  /:/         \  \:\     |  |:|      /  /:/ /\    /  /:/\:\  
	 /  /:/ /::\ /__/\   /  /\/__/::\     _____\__\:\  __|  |:|     /  /:/ /:/_  /  /:/~/:/  
	/__/:/ /:/\:\\  \:\ /  /:/\__\/\:\__ /__/::::::::\/__/\_|:|____/__/:/ /:/ /\/__/:/ /:/___
	\  \:\/:/~/:/ \  \:\  /:/    \  \:\/\\  \:\~~\~~\/\  \:\/:::::/\  \:\/:/ /:/\  \:\/:::::/
	 \  \::/ /:/   \  \:\/:/      \__\::/ \  \:\  ~~~  \  \::/~~~~  \  \::/ /:/  \  \::/~~~~ 
	  \__\/ /:/     \  \::/       /__/:/   \  \:\       \  \:\       \  \:\/:/    \  \:\     
	    /__/:/       \__\/        \__\/     \  \:\       \  \:\       \  \::/      \  \:\    
	    \__\/                                \__\/        \__\/        \__\/        \__\/    

	
	Full Slinker Pipeline - superTranscript visualisation of novel splice variants in RNA-Seq data.

	Author: Breon Schmidt
	Version: 0.1
	Date: 6 MAY 2020

==================================================================================================*/

/*-----------------------------------------------------------

	L O A D 

-----------------------------------------------------------*/

// location of this program, i.e. /path/to/program/, default is level up from this script
core = System.getenv("SLINKERDIR")

// Internal tools
load core+"/workflows/tools.groovy"
FLATTEN = core+"/scripts/generate_flattened_gtf.R"
VIS = core+"/scripts/slinker_vis.py"

/*-----------------------------------------------------------

	P I P E L I N E

-----------------------------------------------------------*/

/*===========================================================
	Get user defined variables or define defaults
===========================================================*/

if(!binding.variables.containsKey("results")){
	print("Please enter a destination for the ouput. i.e. -p results='destination/to/results/'")
	System.exit(0)
}

if(!binding.variables.containsKey("gene")){
	print("Please enter a gene name for inspection! Exiting.'")
	System.exit(0)
}

if(!binding.variables.containsKey("case_sample")){
	print("Please enter a sample to use as the case! Exiting.'")
	System.exit(0)
}

if(binding.variables.containsKey("genome")){
	if(genome == 38){ 
		GTF_REF = core+"/references/"+"Homo_sapiens.GRCh38.99.chr.tsl1.gtf"
		GENOME = core+"/references/"+"Homo_sapiens.GRCh38.dna.primary_assembly.fa"
	} else if(genome == 19) {
		GTF_REF = core+"/references/"+"Homo_sapiens.GRCh37.87.chr.tsl1.gtf"
		GENOME = core+"/references/"+"Homo_sapiens.GRCh37.dna.primary_assembly.fa"
	} else {
		print("Using user supplied GTF and FASTA")
	}
} else {

	if(binding.variables.containsKey("fasta")){
		GENOME = fasta
	} else {
		print("Please either set the genome version OR supply a valid reference FASTA path.'")
		System.exit(0)
	}

	if(binding.variables.containsKey("gtf")){
		GTF_REF = gtf
	} else {
		print("Please either set the genome version OR supply a valid reference gtf path.'")
		System.exit(0)
	}
}


if(!binding.variables.containsKey("threads")){
	threads=1
}

if(!binding.variables.containsKey("conservative")){
	conservative=false
}

if(!binding.variables.containsKey("g_mem")){
	threads=4000000000 // 4GB
}

if(!binding.variables.containsKey("a_mem")){
	threads=4000000000 // 4GB
}

if(!binding.variables.containsKey("width")){
	width=25 // 4GB
}

if(!binding.variables.containsKey("height")){
	height=18 // 4GB
}

if(!binding.variables.containsKey("min_support")){
	min_support=10 // 4GB
}

if(!binding.variables.containsKey("format")){
	format="*/%_Aligned.sortedByCoord.out.bam" // 4GB
}

/*===========================================================
	Segment the cases and controls
===========================================================*/

// Expect the case to be identified with -p CASE=<file>
requires case_sample : 'The sample to be treated as the case'
// Controls are everything except the case
controls = args.grep { it != case_sample }
// Cases are a list of samples, just including the single case 
cases = [ case_sample ]


/*===========================================================
	Setup paths
===========================================================*/

results_folder = results + "/" + gene
resources_folder = results_folder + "/resources"
plots_folder = results_folder + "/plots"
temp_folder = results_folder + "/temp/"
temp_seq = temp_folder + "/sequence"
temp_assembly = temp_folder + "/assembly"
temp_gene_info = temp_folder + "/gene_info"


/*===========================================================
	Stages
===========================================================*/

// ### Setup

make_dir = {
	produce("$temp_folder/info.log"){
		def gene_folder = temp_source+"genes/"+gene+"/"

		exec """mkdir -p $resources_folder;
				mkdir -p $plots_folder;
				mkdir -p $temp_seq;
				mkdir -p $temp_assembly;
				mkdir -p $temp_gene_info;				
				touch $temp_folder/info.log"""
	}
	forward inputs.bam
}


// ### Extract reads from gene

get_gene_region = {
	
	output.dir = resources_folder
	def bed_file = temp_gene_info+"/"+gene+".region.bed"
	def gtf_file = gene+".gtf"

	from('.bam') produce(gtf_file) {
		exec """grep '"$gene"' $GTF_REF | awk -F"\\t" '{print \$1"\\t"\$4"\\t"\$5"\\tslinker\\t"\$6"\\t"\$7}' > $bed_file"""
		exec """grep '"$gene"' $GTF_REF > $output.dir/$gtf_file"""
	}

	forward inputs.bam

}

extract_reads = {

	output.dir= temp_gene_info
	def bed_file = output.dir+"/"+gene+".region.bed"

	transform(".bam") to (".gene.bam"){
		exec """$samtools view -@ $threads -Sb -h -L $bed_file $input.bam > $output.gene.bam""", "samtools"
	}

}

// ### Prepare files

qsort_bam = {
	transform('.gene.bam') to('.qsort.bam') {
		output.dir=temp_gene_info
		exec """$samtools sort -n $input.bam > $output.qsort.bam""", "samtools"
	}
}

// ### Genome guided assembly of transcripts

assemble_transcripts = {
	transform('.gene.bam') to('.assembly.gtf')  {
		output.dir=temp_assembly
		if(conservative == "true"){
			exec """$STRINGTIE $input.gene.bam -G $GTF_REF -p $threads -o $output.assembly.gtf -t -c 100 -f 0.3 -v""", "stringtie"
		} else {
			exec """$STRINGTIE $input.gene.bam -G $GTF_REF -p $threads -o $output.assembly.gtf -t -c 100 -f 0.3 -v""", "stringtie"
		}
		
	}
}


merge_original = {

	produce(final_resources + '/assembly.combined.gtf'){
		from("**.assembly.gtf") {
			output.dir=resources_folder
			def gtf = resources_folder+"/"+gene+".gtf"
			exec """$STRINGTIE --merge -G $gtf -p $threads -c 100 -f 0.2 $inputs > $output.dir/assembly.combined.gtf"""
		}
	}

}


// ### Create new superTranscript 

flatten_gtf = {

	from("assembly.combined.gtf") produce (resources_folder + "/flattened.gtf"){
		output.dir=resources_folder
		def out_file = resources_folder + "/flattened.gtf"
		exec """R CMD BATCH --no-restore --no-save 
				"--args input_gtf='$input' output_gtf='$out_file'"
				$FLATTEN $output.dir/generate_flattened_genome_annotation.log;""", "supertranscript"
	}

}

create_st = {

	from("flattened.gtf")  produce(resources_folder + "/st.fasta"){
		output.dir=resources_folder
		def out_file = resources_folder + "/st.fasta"
		exec """gffread $input -g $GENOME -w $out_file -W""", "supertranscript"
	}

}

// ### Align to new reference

get_fastq = {

	transform(".qsort.bam") to (".o.fastq", ".0.fastq", ".all.fastq"){
		output.dir = temp_seq
		exec """$samtools fastq -o $output1 -0 $output2 -s /dev/null -n $input.qsort.bam"""
		exec """cat $output1 $output2 > $output3"""
		//exec """bamToFastq -i $input.qsort.bam -fq $output.R1.fastq -fq2 $output.R2.fastq""", fastq
	}
	
}


star_genome = {

	output.dir = temp_folder+"/genome/"
	from("st.fasta") produce("Genome"){
		exec """$STAR --runMode genomeGenerate 
						 --genomeDir $output.dir
						 --limitGenomeGenerateRAM $g_mem 
						 --runThreadN $threads 
						 --genomeFastaFiles $input
						 --genomeSAindexNbases 4""", "genome_gen"
	}
}


star_align = {

	output.dir = resources_folder

	from(".all.fastq") produce(branch.name+"_Aligned.sortedByCoord.out.bam"){

		def sample_name = branch.name+"_"
		def genome_dir = temp_folder+"/genome/"
		exec """$STAR --genomeDir $genome_dir
					 --runMode alignReads
					 --twopassMode Basic 
					 --alignIntronMin 4 
					 --outSJfilterOverhangMin 12 12 12 12
					 --outSJfilterCountUniqueMin 1 1 1 1
					 --outSJfilterCountTotalMin 1 1 1 1
					 --outSJfilterDistToOtherSJmin 0 0 0 0
					 --scoreGapNoncan 0
					 --scoreGapGCAG 0
					 --scoreGapATAC 0
					 --readFilesIn $input
					 --runThreadN $threads
					 --limitBAMsortRAM $a_mem 
					 --outFileNamePrefix $output.dir/$sample_name
					 --outSAMtype BAM SortedByCoordinate""", "align"
	}

}

star_index = {

	output.dir = resources_folder
	from(".bam") produce(branch.name+"_Aligned.sortedByCoord.out.bam.bai"){
		def sample_name = branch.name
		exec """$samtools index $input""", "samtools"
	}

}


visualise = {

	output.dir = plots_folder
	exec """python $VIS $gene $case_sample $resources_folder $genome $output.dir"""

}

/*===========================================================
	Run
===========================================================*/

run { 
	[ cases: cases ]*[make_dir + get_gene_region] +
	format*[extract_reads + qsort_bam + get_fastq + assemble_transcripts] +
	merge_original + flatten_gtf + create_st + star_genome + 
	format*[extract_reads + qsort_bam + get_fastq + star_align + star_index] +
	[visualise]
}
