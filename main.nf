$HOSTNAME = ""
params.outdir = 'results'  

THREADS = 4

//* autofill
if ($HOSTNAME == "default"){
	$SINGULARITY_IMAGE = "/export/src/image_pipeline-2.7.sif"
	$SINGULARITY_OPTIONS = "--bind /export"
}
//* autofill

if (!params.species_name){params.species_name = ""} 
if (!params.sequence){params.sequence = ""} 
if (!params.Mt_DB){params.Mt_DB = ""} 

Channel.value(params.species_name).into{g_1_species_name_g_49;g_1_species_name_g_125}
Channel.value(params.sequence).set{g_2_sequence_g_125}
Channel.value(params.Mt_DB).into{g_15_commondb_path_g_56;g_15_commondb_path_g_125}


process tblastn {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /report.blast$/) "report/$filename"
}

input:
 val SPNAME from g_1_species_name_g_125
 val SEQUENCE from g_2_sequence_g_125
 val DB from g_15_commondb_path_g_125

output:
 file "query.fasta"  into g_125_genome
 file "report.blast"  into g_125_blast_output_g_126
 file "hits_{yes,no}.txt"  into g_125_outputFileTxt_g_126

script:

nseqs = params.tblastn.nseqs
gencode = params.tblastn.gencode

"""
printf ">$SPNAME\n$SEQUENCE\n" 1>query.fasta
tblastn -db $DB -db_gencode $gencode -num_alignments $nseqs -query query.fasta -out report.blast

if [ `grep -c "No hits found" report.blast` -eq 0 ]
then 
	echo "Found hits in the database for given query" > hits_yes.txt
else
	echo "There are no hits in the database for given query" > hits_no.txt
fi
"""


}


process blast_result2fasta {

input:
 file blast_report from g_125_blast_output_g_126
 file hits_file from g_125_outputFileTxt_g_126

output:
 set val("sequences"), file("sequences.fasta")  into g_126_genomes_g_49

when:
hits_file.toString() == "hits_yes.txt"

script:
"""
mview -in blast -out fasta $blast_report 1>sequences.fasta
"""
}


process extract_outgroup {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_sel.hash$/) "report/$filename"
}

input:
 val SPNAME from g_1_species_name_g_49
 set val(name), file(query_out_fasta) from g_126_genomes_g_49

output:
 set val("${name}_sel"), file("${name}_sel.fasta")  into g_49_genomes
 set val("${name}_sel"), file("${name}_sel.hash")  into g_49_hash_file_g_56

"""
/export/src/dolphin/scripts/header_sel_mod3.pl $query_out_fasta $SPNAME 1>${name}_sel.fasta 2>${name}_sel.hash

"""
}


process extract_sequences {

input:
 set val(name), file(hash) from g_49_hash_file_g_56
 val DB from g_15_commondb_path_g_56

output:
 set val(name), file("${name}.nuc")  into g_56_nucleotides_g_99

"""
/export/src/dolphin/scripts/nuc_coding_mod.pl $hash $DB 1>${name}.nuc

"""
}


process drop_dublicates {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_unique.fasta$/) "sequences/$filename"
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


process macse_mulal {

input:
 set val(name), file(seqs) from g_99_genomes_g_100
 file report from g_99_outputFileTxt_g_100

output:
 set val(name), file("${name}_mulal.fna")  into g_100_nucl_mulal_g_128
 set val(name), file("${name}_mulal.faa")  into g_100_aa_mulal

when:
report.toString() == "report_yes.txt"

script:
"""
java -jar /opt/macse_v2.05.jar -prog alignSequences -gc_def 2 -out_AA ${name}_mulal.faa -out_NT ${name}_mulal.fna -seq $seqs
"""
}


process mulal_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /alignment_checked.fasta$/) "sequences/$filename"
}

input:
 set val(name), file(mulal) from g_100_nucl_mulal_g_128

output:
 set val("alignment_checked"),file("alignment_checked.fasta")  into g_128_nucl_mulal_g_70, g_128_nucl_mulal_g_129

"""
/export/src/dolphin/scripts/macse2.pl $mulal alignment_checked.fasta

"""
}


process convert_alignment_to_phylip {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.phy$/) "sequences/$filename"
}

input:
 set val(name), file(mulal) from g_128_nucl_mulal_g_129

output:
 set val(name), file("${name}.phy")  into g_129_phylip_g_116, g_129_phylip_g_118, g_129_phylip_g_130, g_129_phylip_g_135

"""
java -jar /opt/readseq.jar -a -f Phylip -o ${name}.phy $mulal

"""
}

run_RAXML = params.RAxML_build_tree.run_RAXML
raxml_model = params.RAxML_build_tree.raxml_model
//* @style @condition:{run_RAXML="true", raxml_model}


process RAxML_build_tree {

input:
 set val(name), file(mulal) from g_129_phylip_g_130

output:
 set val("raxml"), file("raxml.nwk")  into g_130_tree_g_109, g_130_tree_g_111

when:
run_RAXML == "true"

script:
"""
/opt/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -x 987654 -p 987654 -T $THREADS -N 50 -f a -s $mulal -n Rax_tree -m $raxml_model
mv RAxML_bestTree.Rax_tree raxml.nwk
"""

}
params.RAxML_model = raxml_model

process extract_terminal_branch_lengths_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "RAxML/$filename"
}

input:
 set val(name), file(tree) from g_130_tree_g_109

output:
 set val(name), file("${name}.branches")  into g_109_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1>${name}.branches

"""
}

run_IQTREE = params.IQTREE_build_tree.run_IQTREE
iqtree_model = params.IQTREE_build_tree.iqtree_model
//* @style @condition:{run_IQTREE="true", iqtree_model}


process IQTREE_build_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /iqtree.log$/) "IQTREE/$filename"
}

input:
 set val(name), file(mulal) from g_129_phylip_g_118

output:
 set val("iqtree"), file("iqtree.nwk")  into g_118_tree_g_131, g_118_tree_g_132
 file "iqtree.log"  into g_118_logFile

when:
run_IQTREE == "true"

script:

"""
iqtree2 -s $mulal -m $iqtree_model -nt $THREADS --prefix phylo
mv phylo.treefile iqtree.nwk
mv phylo.iqtree iqtree.log
"""

}
params.IQTREE_model = iqtree_model

process extract_terminal_branch_lengths_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "IQTREE/$filename"
}

input:
 set val(name), file(tree) from g_118_tree_g_132

output:
 set val(name), file("${name}.branches")  into g_132_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1>${name}.branches

"""
}


process rooting_iqtree_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "IQTREE/$filename"
}

input:
 set val(name), file(tree) from g_118_tree_g_131

output:
 set val("${name}_rooted"), file("*.nwk")  into g_131_tree_g_135

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk

"""
}


process IQTREE_anc_rec {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /iqtree_anc.state$/) "IQTREE/$filename"
	else if (filename =~ /iqtree_anc.log$/) "IQTREE/$filename"
}

input:
 set val(name), file(mulal) from g_129_phylip_g_135
 set val(namet), file(tree) from g_131_tree_g_135

output:
 set val(name), file("iqtree_anc.state")  into g_135_state
 file "iqtree_anc.log"  into g_135_logFile

"""
iqtree2 -s $mulal -m $params.IQTREE_model -asr -nt $THREADS --prefix anc
mv anc.iqtree iqtree_anc.log
mv anc.state iqtree_anc.state
"""
}


process count_prior_mutations {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /mutnumbers.tsv$/) "report/$filename"
}

input:
 set val(name), file(seqs) from g_128_nucl_mulal_g_70

output:
 file "mutnumbers.tsv"  into g_70_outputFileTSV

"""
/export/src/dolphin/scripts/mutnumbers.pl $seqs 1>mutnumbers.tsv

"""
}


process rooting_raxml_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "RAxML/$filename"
}

input:
 set val(name), file(tree) from g_130_tree_g_111

output:
 set val("${name}_rooted"), file("*.nwk")  into g_111_tree_g_116

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk

"""
}


process RAxML_anc_rec {

input:
 set val(namet), file(tree) from g_111_tree_g_116
 set val(name), file(mulal) from g_129_phylip_g_116


"""
/opt/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -T $THREADS -f A -m $params.RAxML_model -s $mulal -t $tree -n ANCESTORS

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
