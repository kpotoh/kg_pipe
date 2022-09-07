$HOSTNAME = ""
params.outdir = 'results'  

//* autofill
if ($HOSTNAME == "default"){
	$SINGULARITY_IMAGE = "/export/src/image_pipeline-2.5.sif"
	$SINGULARITY_OPTIONS = "--bind /export"
}
//* autofill

if (!params.species_name){params.species_name = ""} 
if (!params.sequence){params.sequence = ""} 
if (!params.Mt_DB){params.Mt_DB = ""} 
if (!params.number_of_threads){params.number_of_threads = ""} 
if (!params.Substitution_model){params.Substitution_model = ""} 
if (!params.inputparam){params.inputparam = ""} 

Channel.value(params.species_name).into{g_1_species_name_g_49;g_1_species_name_g_98}
Channel.value(params.sequence).set{g_2_sequence_g_98}
Channel.value(params.Mt_DB).into{g_15_commondb_path_g_56;g_15_commondb_path_g_98}
Channel.value(params.number_of_threads).into{g_75_nthreads_g_85;g_75_nthreads_g_101}
Channel.value(params.Substitution_model).set{g_102_model_g_101}


process s1_tblastn {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /query.fasta$/) "query_file/$filename"
	else if (filename =~ /query_out.gb$/) "query_result/$filename"
	else if (filename =~ /hits_(yes|no).txt$/) "hits_report/$filename"
}

input:
 val SPNAME from g_1_species_name_g_98
 val SEQUENCE from g_2_sequence_g_98
 val DB from g_15_commondb_path_g_98

output:
 file "query.fasta"  into g_98_genome
 set val("query_out"), file("query_out.gb")  into g_98_genbank_file_g_97
 file "hits_{yes,no}.txt"  into g_98_outputFileTxt_g_97

"""
printf ">$SPNAME\n$SEQUENCE\n" 1>query.fasta
tblastn -db $DB -db_gencode 2 -num_alignments 500 -query query.fasta -out query_out.gb

if [ `grep -c "No hits found" query_out.gb` -eq 0 ]
then 
	echo "Found hits in the database for given query" > hits_yes.txt
else
	echo "There are no hits in the database for given query" > hits_no.txt
fi

"""
}


process s2_query_result2fasta {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.fasta$/) "s2_out/$filename"
}

input:
 set val(name), file(qout) from g_98_genbank_file_g_97
 file hits_file from g_98_outputFileTxt_g_97

output:
 set val(name), file("${name}.fasta")  into g_97_genomes_g_49

when:
hits_file.toString() == "hits_yes.txt"

script:
"""
mview -in blast -out fasta $qout 1>${name}.fasta
echo "query result converted to fasta"

echo $hits_file
cat $hits_file
"""
}


process s3_extract_outgroup {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_sel.fasta$/) "sel_genomes/$filename"
	else if (filename =~ /${name}_sel.hash$/) "sel_hash/$filename"
}

input:
 val SPNAME from g_1_species_name_g_49
 set val(name), file(query_out_fasta) from g_97_genomes_g_49

output:
 set val("${name}_sel"), file("${name}_sel.fasta")  into g_49_genomes
 set val("${name}_sel"), file("${name}_sel.hash")  into g_49_hash_file_g_56

"""
/export/src/dolphin/scripts/header_sel_mod3.pl $query_out_fasta $SPNAME 1>${name}_sel.fasta 2>${name}_sel.hash

"""
}


process s4_extract_sequences {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.nuc$/) "nuc/$filename"
}

input:
 set val(name), file(hash) from g_49_hash_file_g_56
 val DB from g_15_commondb_path_g_56

output:
 set val(name), file("${name}.nuc")  into g_56_nucleotides_g_99

"""
/export/src/dolphin/scripts/nuc_coding_mod.pl $hash $DB 1>${name}.nuc

"""
}


process s5_drop_dublicates {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_unique.fasta$/) "fasta/$filename"
	else if (filename =~ /report_(yes|no).txt$/) "report/$filename"
}

input:
 set val(name), file(nuc) from g_56_nucleotides_g_99

output:
 set val(name), file("${name}_unique.fasta")  into g_99_genomes_g_100
 file "report_{yes,no}.txt" optional true  into g_99_outputFileTxt_g_100

"""
/export/src/dolphin/scripts/codon_alig_unique.pl $nuc 1>${name}_unique.fasta

"""
}

g_99_outputFileTxt_g_100= g_99_outputFileTxt_g_100.ifEmpty([""]) 


process s6_macse_mulal {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_mulal.fna$/) "mulal/$filename"
	else if (filename =~ /${name}_mulal.faa$/) "mulal/$filename"
}

input:
 set val(name), file(seqs) from g_99_genomes_g_100
 file report from g_99_outputFileTxt_g_100

output:
 set val(name), file("${name}_mulal.fna")  into g_100_nucl_mulal_g_88
 set val(name), file("${name}_mulal.faa")  into g_100_aa_mulal

when:
report.toString() == "report_yes.txt"

script:
"""
java -jar /opt/macse_v2.05.jar -prog alignSequences -gc_def 2 -out_AA ${name}_mulal.faa -out_NT ${name}_mulal.fna -seq $seqs
"""
}


process s7_mulal_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_checked.fna$/) "mulal_nucl_checked/$filename"
}

input:
 set val(name), file(mulal) from g_100_nucl_mulal_g_88

output:
 set val("${name}_checked"), file("${name}_checked.fna")  into g_88_nucl_mulal_g_70, g_88_nucl_mulal_g_72

"""
/export/src/dolphin/scripts/macse2.pl $mulal ${name}_checked.fna

"""
}


process s9_prepare_nucl_for_phylo {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.phy$/) "phylip_mulal/$filename"
}

input:
 set val(name), file(seqs) from g_88_nucl_mulal_g_72

output:
 set val(name), file("${name}.phy")  into g_72_phylip_g_85, g_72_phylip_g_101

"""
java -jar /opt/readseq.jar -a -f Phylip -o ${name}.phy $seqs

"""
}


process s10_phylogeny_recontruction_ml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /RAxML_bestTree.Rax_tree$/) "RAX_tree/$filename"
}

input:
 set val(name), file(mulal) from g_72_phylip_g_85
 val THREADS from g_75_nthreads_g_85

output:
 set val("RAxML_bestTree"), file("RAxML_bestTree.Rax_tree")  into g_85_rax_tree_g_86, g_85_rax_tree_g_87

"""
/opt/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -x 987654 -p 987654 -T $THREADS -N 50 -f a -s $mulal -n Rax_tree -m GTRGAMMAIX

"""
}


process s13_tree_rooting {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_rooted.Rax_tree$/) "tree_rooted/$filename"
}

input:
 set val(name), file(tree) from g_85_rax_tree_g_87

output:
 set val(name), file("${name}_rooted.Rax_tree")  into g_87_rax_tree

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.Rax_tree

"""
}


process s11_12_extract_terminal_branches_lengths {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "branches/$filename"
}

input:
 set val(name), file(tree) from g_85_rax_tree_g_86

output:
 set val(name), file("${name}.branches")  into g_86_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1>${name}.branches

"""
}


process s10_phylogeny_IQTREE2 {

input:
 set val(name), file(mulal) from g_72_phylip_g_101
 val THREADS from g_75_nthreads_g_101
 val MODEL from g_102_model_g_101


"""
iqtree2 -s $mulal -m $MODEL -asr -nt $THREADS --prefix anc
"""
}


process s8_calc_prior_mut_num {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /mutnumbers.tsv$/) "mutnumbers/$filename"
}

input:
 set val(name), file(seqs) from g_88_nucl_mulal_g_70

output:
 file "mutnumbers.tsv"  into g_70_outputFileTSV

"""
/export/src/dolphin/scripts/mutnumbers.pl $seqs 1>mutnumbers.tsv

"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
