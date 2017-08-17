#!/usr/bin/env nextflow

/*
*params input 
*/
params.reads = "$baseDir/data/*{1,2}.fastq"
params.genome = "$baseDir/data/GRCh37_region1.fa"
params.annotation = "$baseDir/data/GRCh37_region1.gtf"
params.index = null




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
    log.info '    --genome GENOME_FILE                Reference genome file in fomrat fa .'
    log.info '    --annotation ANNOTATION_FILE        Annotation file in format GTF .'
    log.info '    --index GENOME_INDEX_FILE           Index file with the same prefix of genome reference '
    log.info '                                        (ex : --genome PATH/hg19.fa --index "PATH/hg19*").[Optional but most faster].'   
    exit 1
}


/*
*Path of genome file and annotation file   
*/
genome_file = file(params.genome)
annotation_file = file(params.annotation)


/*
*Channel for reads pairs   
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
* Control Quality of original read_pairs by fastQC
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
 * Step 1.0 : Builds the genome index required by the mapping process.
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
 * Step 1.1 : Builds the transcript index, required by the mapping process.
 */

process buildTranscriptomeIndex{

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
* Step 2. : Maps each read-pair by using Tophat2.1.1 mapper tool
*OPTION :
*
*-p     Use this many threads to align reads.
*-r     This is the expected (mean) inner distance between mate pairs.
*--gtf  Annotation file in GTF (see doc of tophat2.1.1)
*--transcriptome-index=     Name of annotation index (see doc of tophat2.1.1)
*--coverage-search      Enables the coverage based search for junctions. Use when coverage 
*    search is disabled by default (such as for reads 75bp or longer), for maximum sensitivity.
*--microexon-search       With this option, the pipeline will attempt to find 
*    alignments incident to micro-exons. Works only for reads 50bp or longer.
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
    --coverage-search --microexon-search \
    $genome.baseName $reads 

    mv tophat_out/accepted_hits.bam ./${pair_id}.bam
    """
}


bam.into{
    bam_count_gene
    bam_count_exon
    bam_count_transcript  
} 

/*
* Step 3. : Count the transcript by using the "featureCounts" tool
*OPTION :
* 
*-T     Number of the threads. 1 by default.
*-p     Count fragments (read pairs) instead of individual reads
*-M     Multi-mapping reads will also be counted. For a multi-mapping read,
*       all its reported alignments will be counted
*-largestOverlap    Assign reads to a meta-feature/feature that has the largest number of overlapping bases.
*-t     Specify feature type in GTF annotation. `exon' by default.
*-g     Specify attribute type in GTF annotation. `gene_id' by default.
*-a     Name of an annotation file. GTF format by default.
*-o     Name of the output file including read counts.
*-O     Assign reads to all their overlapping meta-features
*-s     Perform strand-specific read counting. Possible values:  
*           0 (unstranded), 1 (stranded) and 2 (reversely stranded).
*-f     Perform read counting at feature level (eg. counting reads for exons rather than genes). 
*/

 
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
