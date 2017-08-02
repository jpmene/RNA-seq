#!/usr/bin/env nextflow

/*
*params input 
*/
params.reads = "$baseDir/data/*{1,2}.fastq"
params.genome = "$baseDir/data/GRCh37_region1.fasta"
params.annotation = "/home/jp/featurecount/data/Homo_sapiens.GRCh38.89.gtf"
params.compare = null
params.index = null
params.featurecount = "/home/jp/featurecount/feauturecount.nf"
params.tophat2 = "tophat2"


//print usage
if (params.help) {
    log.info ''
    log.info 'RNA-seq'
    log.info '-----------------------'
    log.info '.'
    log.info ''
    log.info 'Usage: '
    log.info '    nextflow RNA-seq.nf --file_bam_compare ../file.bam --type_data Illumina'
    log.info '            --name_dir toto --annotation file.gtf --path_featurecounts /bin/featureCqounts '
    log.info ''
    log.info 'Options:'
    log.info '    --help                              Show this message and exit.'
    log.info '    --reads                             File reads paired in fastq [ex : "data/*_{1,2}.{fastq,fq}"'
    log.info '    --genome GENOME_FILE                Reference genome file(s).'
    log.info '    --annotation ANNOTATION_FILE        Annotation genome use in GTF .'
    log.info '    --index GENOME_INDEX_FILE           Index file.[Optional but most faster].'   
    exit 1
}


/*
*create a read_pairs by params.reads and genome ref 
*/
genome_file = file(params.genome)
annotation_file = file(params.annotation)
featurecount= params.featurecount
tophat2= params.tophat2
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
    cpus 4
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
 * Step 1. Builds the genome index required by the mapping process
 */

if(params.index == null){


    process buildIndex {
        cpus 4
        tag "$genome_file.baseName"
    
        input:
        file genome from genome_file
     
        output:
        file 'genome.index*' into genome_index
       
        """
        bowtie2-build ${genome} genome.index
        """
    }
}
 else{
    genome_index =  Channel
                .fromPath(params.index).toList()
}

//genome_index.subscribe { println "D: $it" }

/*
 * Step 2. Maps each read-pair by using Tophat2 mapper tool
 */
process mapping {
    tag "$pair_id"
    publishDir "result/RNA-seq/$pair_id/bam/", mode: "copy"
    cpus 4
    input:
    file genome from genome_file 
    file index from genome_index
    set pair_id, file(reads) from read_pairs_map
 
    output:
    set pair_id, "tophat_out" into log_map
    set pair_id, "${pair_id}.bam" into bam

    """
    ${tophat2} -p ${task.cpus} genome.index $reads
    mv tophat_out/accepted_hits.bam ./${pair_id}.bam
    """
}

/*
bam.into{
    bam_count
    bam_transcipt
    bam_diff  
} 



if(params.compare == null){

    process count {
        tag "$pair_id"
        publishDir "result/RNA-seq/$pair_id/", mode: "copy"
           
        input:
        set pair_id, file(bam_file) from bam_count
        file annotation from annotation_file
     
        output:
        set pair_id, "count" into count 

        """
        nextflow ${featurecount} \
        --file_bam ${bam_file} \
        --annotation ${annotation} \
        --name_dir ${pair_id} \
        --path_featurecounts '/bin/featureCqounts' \
        --type_data Illumina \

        mkdir count
        mv result/featureCounts/$pair_id/* count/.
        """
    }
}
else{


    compare_bam =  Channel
                .fromPath(params.compare)
                .map { file -> tuple(file.baseName, file) }

    list_bam = bam_count.combine(compare_bam,by:0)



    process count_and_compare {
        tag "$pair_id"
        publishDir "result/RNA-seq/$pair_id/", mode: "copy"
       
        input:
        set pair_id, file(bam_file) , file (bam_comapre) from list_bam
        file annotation from annotation_file
     
        output:
        set pair_id, "count" into count 

        """
        nextflow ${featurecount} \
        --file_bam ${bam_file} \
        --annotation ${annotation} \
        --name_dir ${pair_id} \
        --path_featurecounts '/bin/featureCqounts' \
        --type_data Illumina \
        --file_bam_compare ${bam_comapre}

        mkdir count
        mv result/featureCounts/$pair_id/* count/.
        """
    }

}


*/

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
    cpus 2
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
    cpus 2
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
