$HOSTNAME = ""
params.outdir = 'results'  

THREADS = 4

//* autofill
if ($HOSTNAME == "default"){
	$SINGULARITY_IMAGE = "/export/src/image_pipeline-2.8.sif"
	$SINGULARITY_OPTIONS = "--bind /export"
}
//* autofill

if (!params.species_name){params.species_name = ""} 
if (!params.sequence){params.sequence = ""} 
if (!params.Mt_DB){params.Mt_DB = ""} 
if (!params.mode){params.mode = ""} 

Channel.value(params.species_name).set{g_1_species_name_g_49}
g_2_multipleFasta_g_206 = file(params.sequence, type: 'any') 
g_2_multipleFasta_g_212 = file(params.sequence, type: 'any') 
Channel.value(params.Mt_DB).into{g_15_commondb_path_g_205;g_15_commondb_path_g_208}
Channel.value(params.mode).into{g_193_mode_g_206;g_193_mode_g_212}


if (!(params.mode == "single")){
g_2_multipleFasta_g_206.set{g_206_multipleFasta_g_205}
} else {

process query_qc {

input:
 file query from g_2_multipleFasta_g_206
 val mode from g_193_mode_g_206

output:
 file "query_single.fasta"  into g_206_multipleFasta_g_205

when:
params.mode == "single"

script:
"""
if [ `grep -c ">" $query` -eq 1 ]
then 
	mv $query query_single.fasta
else
	echo "Query file must contain only one record when mode == 'single'"
	exit 1
fi
"""
}
}



process tblastn {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /report.blast$/) "report/$filename"
}

input:
 file query from g_206_multipleFasta_g_205
 val DB from g_15_commondb_path_g_205

output:
 file "report.blast"  into g_205_blast_output_g_126
 file "hits_{yes,no}.txt"  into g_205_outputFileTxt_g_126

script:

nseqs = params.tblastn.nseqs
gencode = params.tblastn.gencode
params.gencode = gencode

"""
tblastn -db $DB -db_gencode $gencode -num_alignments $nseqs -query $query -out report.blast

if [ `grep -c "No hits found" report.blast` -eq 0 ]
then 
	echo "Found hits in the database for given query" > hits_yes.txt
else
	echo "There are no hits in the database for given query" > hits_no.txt
fi
"""


}


process mview {

input:
 file blast_report from g_205_blast_output_g_126
 file hits_file from g_205_outputFileTxt_g_126

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
 set val("${name}_sel"), file("${name}_sel.hash")  into g_49_hash_file_g_208

"""
/export/src/dolphin/scripts/header_sel_mod3.pl $query_out_fasta $SPNAME 1>${name}_sel.fasta 2>${name}_sel.hash

"""
}


process extract_sequences {

input:
 set val(name), file(hash) from g_49_hash_file_g_208
 val DB from g_15_commondb_path_g_208

output:
 file "${name}.fasta"  into g_208_multipleFasta_g_207

"""
/export/src/dolphin/scripts/nuc_coding_mod.pl $hash $DB 1>${name}.fasta

"""
}


process aln_qc {

input:
 val mode from g_193_mode_g_212
 file query from g_2_multipleFasta_g_212

output:
 file "query_multiple.fasta"  into g_212_multipleFasta_g_207

script:
"""
if [ `grep -c ">" $query` -gt 1 ]
then 
	mv $query query_multiple.fasta
else
	echo "Query file must contain more than one record when mode == 'multiple'"
	exit 1
fi
"""

}


process drop_dublicates {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /seqs_unique.fasta$/) "sequences/$filename"
	else if (filename =~ /report_(yes|no).txt$/) "report/$filename"
}

input:
 file nucleotides from g_208_multipleFasta_g_207
 file nucleotides from g_212_multipleFasta_g_207

output:
 file "seqs_unique.fasta"  into g_207_multipleFasta_g_209
 file "report_{yes,no}.txt" optional true  into g_207_outputFileTxt_g_209

"""
/export/src/dolphin/scripts/codon_alig_unique.pl $nucleotides 1>seqs_unique.fasta

"""
}

g_207_outputFileTxt_g_209= g_207_outputFileTxt_g_209.ifEmpty([""]) 


process macse {

input:
 file seqs from g_207_multipleFasta_g_209
 file report from g_207_outputFileTxt_g_209

output:
 set val(name), file("${name}_mulal.fna")  into g_209_nucl_mulal_g_128
 set val(name), file("${name}_mulal.faa")  into g_209_aa_mulal

when:
report.toString() == "report_yes.txt"

script:
name = seqs.toString() - '.fasta'

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
 set val(name), file(mulal) from g_209_nucl_mulal_g_128

output:
 set val("alignment_checked"),file("alignment_checked.fasta")  into g_128_nucl_mulal_g_70, g_128_nucl_mulal_g_129, g_128_nucl_mulal_g_151

"""
/export/src/dolphin/scripts/macse2.pl $mulal alignment_checked.fasta

"""
}


process terminal_genomes_states {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /leaves_states.state$/) "IQTREE/$filename"
}

input:
 set val(name), file(mulal) from g_128_nucl_mulal_g_151

output:
 set val("leaves_states"), file("leaves_states.state")  into g_151_state_g_177, g_151_state_g_178

"""
python3 /export/src/mutspec-utils/scripts/alignment2iqtree_states.py $mulal leaves_states.state
"""
}


process convert_alignment_to_phylip {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.phy$/) "tmp/$filename"
}

input:
 set val(name), file(mulal) from g_128_nucl_mulal_g_129

output:
 set val(name), file("${name}.phy")  into g_129_phylip_g_130, g_129_phylip_g_145, g_129_phylip_g_189, g_129_phylip_g_191

"""
java -jar /opt/readseq.jar -a -f Phylip -o ${name}.phy $mulal

"""
}

run_IQTREE = params.IQTREE_build_tree.run_IQTREE
iqtree_model = params.IQTREE_build_tree.iqtree_model
//* @style @condition:{run_IQTREE="true", iqtree_model}


process IQTREE_build_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.log$/) "IQTREE/$filename"
}

input:
 set val(name), file(mulal) from g_129_phylip_g_145

output:
 set val("iqtree"), file("iqtree.nwk")  into g_145_tree_g_131, g_145_tree_g_132
 file "*.log"  into g_145_logFile

when:
run_IQTREE == "true"

errorStrategy 'retry'
maxRetries 3

script:

"""
iqtree2 -s $mulal -m $iqtree_model -nt $THREADS --prefix phylo
mv phylo.treefile iqtree.nwk
mv phylo.iqtree iqtree_report.log
mv phylo.log iqtree.log
"""

}
params.IQTREE_model = iqtree_model

process extract_terminal_branch_lengths_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "IQTREE/$filename"
}

input:
 set val(name), file(tree) from g_145_tree_g_132

output:
 set val(name), file("${name}.branches")  into g_132_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1>${name}.branches

"""
}


process rooting_iqtree_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(tree) from g_145_tree_g_131

output:
 set val("${name}_rooted"), file("*.nwk")  into g_131_tree_g_191

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk

"""
}


process IQTREE_anc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /iqtree_anc_tree.nwk$/) "IQTREE/$filename"
	else if (filename =~ /iqtree_anc.state$/) "IQTREE/$filename"
	else if (filename =~ /.*.log$/) "IQTREE/$filename"
}

input:
 set val(name), file(mulal) from g_129_phylip_g_191
 set val(namet), file(tree) from g_131_tree_g_191

output:
 set val("iqtree_anc_tree"), file("iqtree_anc_tree.nwk")  into g_191_tree_g_177
 set val("iqtree_anc"), file("iqtree_anc.state")  into g_191_state_g_177
 file "*.log"  into g_191_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
iqtree2 -s $mulal -m $params.IQTREE_model -asr -nt $THREADS --prefix anc
mv anc.iqtree iqtree_anc_report.log
# mv anc.state iqtree_anc.state
mv anc.log iqtree_anc.log
mv anc.treefile iqtree_anc_tree.nwk

python3 /export/src/mutspec-utils/scripts/iqtree_states_add_part.py anc.state iqtree_anc.state
"""

}

ms_options = params.mutspec_iqtree.ms_options


process mutspec_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "mutspec_iqtree/$filename"
	else if (filename =~ /.*.log$/) "mutspec_iqtree/$filename"
}

input:
 set val(namet), file(tree) from g_191_tree_g_177
 set val(names1), file(states1) from g_191_state_g_177
 set val(names2), file(states2) from g_151_state_g_177

output:
 file "*.tsv"  into g_177_outputFileTSV
 file "*.log"  into g_177_logFile

"""
python3 /export/src/mutspec-utils/scripts/3.calculate_mutspec.py --tree $tree --states $states1 --states $states2 --gencode $params.gencode $ms_options --outdir mout
mv mout/* .

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

process rooting_raxml_tree {

input:
 set val(name), file(tree) from g_130_tree_g_111

output:
 set val("${name}_rooted"), file("*.nwk")  into g_111_tree_g_189

"""
nw_reroot -l $tree OUTGRP 1>${name}_rooted.nwk

"""
}


process RAxML_anc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /RAxML_nodeLabelledRootedTree.nwk$/) "RAxML/$filename"
	else if (filename =~ /RAxML_marginalAncestralProbabilities.state$/) "RAxML/$filename"
	else if (filename =~ /RAxML_anc_rec.log$/) "RAxML/$filename"
}

input:
 set val(namet), file(tree) from g_111_tree_g_189
 set val(name), file(mulal) from g_129_phylip_g_189

output:
 set val("RAxML_nodeLabelledRootedTree"),file("RAxML_nodeLabelledRootedTree.nwk")  into g_189_tree_g_178
 set val("RAxML_marginalAncestralProbabilities"), file("RAxML_marginalAncestralProbabilities.state")  into g_189_state_g_178
 file "RAxML_marginalAncestralStates.fasta"  into g_189_multipleFasta
 file "RAxML_anc_rec.log"  into g_189_logFile

"""
/opt/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -T $THREADS -f A -m $params.RAxML_model -s $mulal -t $tree -n ANCESTORS
mv RAxML_info.ANCESTORS RAxML_anc_rec.log
mv RAxML_marginalAncestralStates.ANCESTORS RAxML_marginalAncestralStates.fasta
# mv RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.txt
# mv RAxML_nodeLabelledRootedTree.ANCESTORS RAxML_nodeLabelledRootedTree.nwk

python3 /export/src/mutspec-utils/scripts/raxml_states2iqtree_states.py RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.state
python3 /export/src/mutspec-utils/scripts/rename_internal_nodes.py $tree RAxML_nodeLabelledRootedTree.ANCESTORS RAxML_nodeLabelledRootedTree.nwk

"""
}

ms_options = params.mutspec_raxml.ms_options


process mutspec_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "mutspec_raxml/$filename"
	else if (filename =~ /.*.log$/) "mutspec_raxml/$filename"
}

input:
 set val(namet), file(tree) from g_189_tree_g_178
 set val(names1), file(states1) from g_189_state_g_178
 set val(names2), file(states2) from g_151_state_g_178

output:
 file "*.tsv"  into g_178_outputFileTSV
 file "*.log"  into g_178_logFile

"""
python3 /export/src/mutspec-utils/scripts/3.calculate_mutspec.py --tree $tree --states $states1 --states $states2 --gencode $params.gencode $ms_options --outdir mout
mv mout/* .

"""
}


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


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
