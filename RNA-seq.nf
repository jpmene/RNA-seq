#!/usr/bin/env nextflow

/*
*params input 
*/
params.reads = "$baseDir/data/*{1,2}.fastq"
params.genome = "$baseDir/data/GRCh37_region1.fa"
params.annotation = "$baseDir/data/GRCh37_region1.gtf"
params.index = null
//params.featurecount = "/home/jp/featurecount/feauturecount.nf"



//print usage
if (params.help) {
    log.info ''
    log.info 'RNA-seq'
    log.info '-----------------------'
    log.info '.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow RNA-seq.nf --genome PATH/genome.fa --reads "PATH/*{1,2}.fq" '
    log.info ' Or '
    log.info ' nextflow RNA-seq.nf --genome PATH/genome.fa  --index "PATH/genome*" --reads "PATH/*{1,2}.fq" '
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --reads                             File reads paired in fastq [ex : "data/*_{1,2}.{fastq,fq}"'
    log.info '    --genome GENOME_FILE                Reference genome file in fomrat .fa  .'
    log.info '    --index GENOME_INDEX_FILE           Index file with the same prefix of genome reference '
    log.info '                                        (ex : --genome PATH/hg19.fa --index "PATH/hg19*").[Optional but most faster].'   
    exit 1
}


/*
*create a read_pairs by params.reads and genome ref 
*/
genome_file = file(params.genome)
annotation_file = file(params.annotation)
//featurecount= params.featurecount

/*
*Path to the tool trimmomatic (need the adapters file)
*/
Channel
    .fromFilePairs( params.reads )                                             
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }  
    .set { read_pairs }

read_pairs.into{
  read_pairs_qc
  read_pairs_map
} 



/*
* Trimming of read_pairs
*/
/*
process trimming{
    tag{pair_id}

    input:
    set pair_id, file(reads) from read_pairs2 

    output: 
    set pair_id, '*trim.fq'  into read_pairs_trim
    """
    java -jar $trimmomatic/trimmomatic-0.36.jar PE -phred33 ${reads} $pair_id'_1_trim.fq' $pair_id'_1_unp.fq' $pair_id'_2_trim.fq' $pair_id'_2_unp.fq' ILLUMINACLIP:$trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
    """
}
*/

/*
* Control Quality of original read_pairs
*/

process fastQC {
    
    tag{pair_id}
    publishDir "result/RNA-seq/$pair_id/QC/", mode: "move"
    
    input:
    set pair_id, file(reads) from read_pairs_qc 

    output:
    file '*fastqc*' into qc

    """
    fastqc -t  ${task.cpus} ${reads}
    """
}

/*
 * Step 1.0 Builds the genome index required by the mapping process
 */

if(params.index == null){


    process buildIndex {
        
        tag "$genome_file.baseName"
    
        input:
        file genome from genome_file
     
        output:
        file '*' into genome_index
       
        """
        bowtie2-build ${genome} $genome.baseName
        """
    }
}
 else{
    genome_index =  Channel
                .fromPath(params.index).toList()
}



/*
 * Step 1.1 Builds the strascipt index required by the mapping process
 */

process build_transcriptome_index{

    tag "$genome"
    publishDir "result/RNA-seq/$pair_id/bam/", mode: "copy"
    

    input:
    file genome from genome_file 
    file index from genome_index
    file annotation from annotation_file
 
    output:
    file  '*_tr' into transcriptome_index

    """
    tophat2 -p ${task.cpus} --GTF ${annotation} \
    --transcriptome-index=${genome.baseName}_tr $genome.baseName
    """
}

/*
 * Step 2. Maps each read-pair by using Tophat2.1.1 mapper tool
 */
process mapping {

    tag "$pair_id"
    publishDir "result/RNA-seq/$pair_id/bam/", mode: "copy"
    

    input:
    file genome from genome_file 
    file index from genome_index
    file index_tr from transcriptome_index
    set pair_id, file(reads) from read_pairs_map
    file annotation from annotation_file
 
    output:
    set pair_id, "tophat_out" into log_map
    set pair_id, "${pair_id}.bam" into bam

    """
    tophat2 -p ${task.cpus} -r 100 --GTF ${annotation} \
    --transcriptome-index=$index_tr \
    $genome.baseName $reads 

    mv tophat_out/accepted_hits.bam ./${pair_id}.bam
    """
}


bam.into{
    bam_count_gene
    bam_count_exon
    bam_count_transcript  
} 


methods = ['gene', 'exon', 'transcript']

process count {
    tag "$pair_id,$mode"
    publishDir "result/RNA-seq/$pair_id/count", mode: "copy"
    
       
    input:
    set pair_id, file (bam_file) from bam_count_gene
    file annotation from annotation_file
    each mode from methods
     
    output: 

    file "${pair_id}_${mode}" into count_gene 
    file "*.summary" into summary_gene 


    """
    featureCounts -T ${task.cpus} \
    -p -M -O --largestOverlap -s 2 -f \
    -t ${mode} -g ${mode}_id \
    -a ${annotation} \
    -o ${pair_id}_${mode} ${bam_file} \
    """    
}

/*
 * Step 3. Assembles the transcript by using the "cufflinks" tool
 */
 
/* 
process cufflinks {
    cpus 2
    tag "$pair_id"
    publishDir "result/RNA-seq/$pair_id/transcripts", mode: "copy"
       
    input:
    set pair_id, file(bam_file) from bam_transcipt
    file annotation from annotation_file
     
    output:
    file('transcript_*.gtf') into transcripts
 
    """
    cufflinks -q -p $task.cpus -G $annotation $bam_file
    mv transcripts.gtf transcript_${pair_id}.gtf 
    """
}



process list_GTF {
    input: 
    val(x) from transcripts.collect()

    output: 
    file list into list_GTF

    script:
    values = x.collect{"$it"}.join('\\n')
    """
    printf \$'${values}' > list
    """
}


process cuffmerge {
    //cpus 4
    publishDir "result/RNA-seq/$genome_file/merge/", mode: "copy"
       
    input:
    file transcipts from list_GTF
    file annotation from annotation_file
    file genome from genome_file 
     
    output:
    file "merged_asm" into log_merge
    file "merge.gtf" into merge
 
    """
    cuffmerge  -p $task.cpus -s $genome -g $annotation $transcipts 
    mv merged_asm/merged.gtf .
    """
}


process cuffdiff {
    //cpus 4
    publishDir "result/RNA-seq/$pair_id/merge", mode: "copy"
       
    input:
    file gtf_merge from merge
    set pair_id, file(bam_file) from bam_diff
    
    output:
    file('cuffdiff') into diff
 
    """
    cuffdiff  -p $task.cpus $gtf_merge $bam_file
    mkdir cuffdiff
    mv * cuffdiff/
    """
}
*/

workflow.onComplete { 
    println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
