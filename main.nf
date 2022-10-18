$HOSTNAME = ""
params.outdir = 'results'  

THREADS = 2
//* params.mode =  "single"  //* @dropdown @options:"single, multiple" @description:"How to run pipeline, either on one record or on multiple records. WARNING: if you run mutational spectrum reconstruction on multiple records, you must be TODO"
//* params.OUTGRP =  "OUTGRP"   //* @input

//* autofill
if ($HOSTNAME == "default"){
	$SINGULARITY_IMAGE = "/export/src/image_pipeline-2.8.sif"
	$SINGULARITY_OPTIONS = "--bind /export"
}
if (params.mode == "multiple"){
	params.OUTGRP = "PASS THE OUTGROUP NAME"
}
else if (params.mode == "single"){
	params.OUTGRP = "OUTGRP"
}

//* autofill




if (!params.species_name){params.species_name = ""} 
if (!params.sequence){params.sequence = ""} 
if (!params.Mt_DB){params.Mt_DB = ""} 
if (!params.gencode){params.gencode = ""} 

Channel.value(params.species_name).set{g_1_species_name_g_49}
g_2_multipleFasta_g_303 = file(params.sequence, type: 'any') 
Channel.value(params.Mt_DB).into{g_15_commondb_path_g_213;g_15_commondb_path_g_308}
Channel.value(params.gencode).into{g_220_gencode_g_222;g_220_gencode_g_308}


process query_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /query_single.fasta$/) "tmp/$filename"
	else if (filename =~ /query_multiple.fasta$/) "tmp/$filename"
}

input:
 file query from g_2_multipleFasta_g_303

output:
 file "query_single.fasta" optional true  into g_303_multipleFasta_g_308
 file "query_multiple.fasta" optional true  into g_303_multipleFasta_g_306

"""
if [ $params.mode == "single" ]; then
	if [ `grep -c ">" $query` -eq 1 ]; then
		if [ $params.OUTGRP != "OUTGRP" ]; then
			echo "If you run 'single' mode you mustn't change the OUTGRP argument and use default one"
			exit 1
		else
			mv $query query_single.fasta
		fi
	else
		echo "Query file must contain only one record when mode == 'single'"
		exit 1
	fi
elif [ $params.mode == "multiple" ]; then
	if [ `grep -c ">" $query` -gt 1 ]; then 
		if [ $params.OUTGRP == "" ]; then
			echo "If you run 'multiple' mode you must pass the OUTGRP argument"
			exit 1
		else
			mv $query query_multiple.fasta
		fi
	else
		echo "Query file must contain more than one record when mode == 'multiple'"
		exit 1
	fi
fi


"""
}

g_303_multipleFasta_g_308= g_303_multipleFasta_g_308.ifEmpty([""]) 


process tblastn {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /report.blast$/) "report/$filename"
	else if (filename =~ /char_numbers.log$/) "report/$filename"
}

input:
 file query from g_303_multipleFasta_g_308
 val DB from g_15_commondb_path_g_308
 val gencode from g_220_gencode_g_308

output:
 file "report.blast"  into g_308_blast_output_g_126
 file "hits_{yes,no}.txt"  into g_308_outputFileTxt_g_126
 file "char_numbers.log"  into g_308_logFile

when:
query.toString() == "query_single.fasta"

script:

nseqs = params.tblastn.nseqs

"""
grep -v  ">" $query | grep -o . | sort | uniq -c | sort -nr > char_numbers.log
if [ `head -n 4 char_numbers.log | grep -Ec "[ACGTacgt]"` -lt 4 ] || [ `grep -Ec "[EFILPQU]" char_numbers.log` -ne 0 ]; then
	tblastn -db $DB -db_gencode $gencode -num_alignments $nseqs -query $query -out report.blast
else
	echo "Query fasta must contain single amino acid sequence"
	exit 1
fi


if [ `grep -c "No hits found" report.blast` -eq 0 ]; then 
	echo "Found hits in the database for given query" > hits_yes.txt
else
	echo "There are no hits in the database for given query" > hits_no.txt
fi
"""


}


process mview {

input:
 file blast_report from g_308_blast_output_g_126
 file hits_file from g_308_outputFileTxt_g_126

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
 set val("${name}_sel"), file("${name}_sel.hash")  into g_49_hash_file_g_213

"""
/export/src/dolphin/scripts/header_sel_mod3.pl $query_out_fasta $SPNAME 1>${name}_sel.fasta 2>${name}_sel.hash

"""
}


process extract_sequences {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /sequences.fasta$/) "tmp/$filename"
}

input:
 set val(name), file(hash) from g_49_hash_file_g_213
 val DB from g_15_commondb_path_g_213

output:
 file "sequences.fasta"  into g_213_multipleFasta_g_277

"""
/export/src/dolphin/scripts/nuc_coding_mod.pl $hash $DB 1>sequences.fasta

"""
}

g_303_multipleFasta_g_306= g_303_multipleFasta_g_306.ifEmpty([""]) 


if (!(params.mode == "multiple")){
g_303_multipleFasta_g_306.set{g_306_multipleFasta_g_277}
g_306_logFile = Channel.empty()
} else {

process qc_aln {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /sequences.fasta$/) "tmp/$filename"
	else if (filename =~ /char_numbers.log$/) "report/$filename"
}

input:
 file query from g_303_multipleFasta_g_306

output:
 file "sequences.fasta"  into g_306_multipleFasta_g_277
 file "char_numbers.log" optional true  into g_306_logFile

when:
params.mode == "multiple"

script:
"""
grep -v  ">" $query | grep -o . | sort | uniq -c | sort -nr > char_numbers.log
if [ `head -n 5 char_numbers.log | grep -Ec "[ACGTacgt]"` -ge 3 ] && [ `grep -Ec "[EFILPQU]" char_numbers.log` -eq 0 ]; then
	mv $query sequences.fasta
else
	echo "Query fasta must contain nucleotides when mode == 'multiple'"
	exit 1
fi
"""

}
}


g_213_multipleFasta_g_277= g_213_multipleFasta_g_277.ifEmpty([""]) 
g_306_multipleFasta_g_277= g_306_multipleFasta_g_277.ifEmpty([""]) 


process drop_dublicates {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /seqs_unique.fasta$/) "sequences/$filename"
	else if (filename =~ /report_(yes|no).txt$/) "report/$filename"
	else if (filename =~ /condensation.log$/) "tmp/$filename"
}

input:
 file nucleotides_single from g_213_multipleFasta_g_277
 file nucleotides_multiple from g_306_multipleFasta_g_277

output:
 file "seqs_unique.fasta"  into g_277_multipleFasta_g_222
 file "report_{yes,no}.txt"  into g_277_outputFileTxt_g_222
 file "condensation.log"  into g_277_logFile

script:
seqs = ""
if (nucleotides_single.toString() == "sequences.fasta") {
	seqs = nucleotides_single
} else if (nucleotides_multiple.toString() == "sequences.fasta") {
	seqs = nucleotides_multiple
}

"""
/export/src/dolphin/scripts/codon_alig_unique.pl $seqs 1>seqs_unique.fasta
echo "$seqs\n$nucleotides_single\n$nucleotides_multiple" > condensation.log

"""

}


process macse {

input:
 file seqs from g_277_multipleFasta_g_222
 file report from g_277_outputFileTxt_g_222
 val gencode from g_220_gencode_g_222

output:
 set val(name), file("${name}_mulal.fna")  into g_222_nucl_mulal_g_128
 set val(name), file("${name}_mulal.faa")  into g_222_aa_mulal

when:
report.toString() == "report_yes.txt"

script:
name = seqs.toString() - '.fasta'

"""
java -jar /opt/macse_v2.05.jar -prog alignSequences -gc_def $gencode -out_AA ${name}_mulal.faa -out_NT ${name}_mulal.fna -seq $seqs
"""
}


process mulal_qc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /alignment_checked.fasta$/) "sequences/$filename"
}

input:
 set val(name), file(mulal) from g_222_nucl_mulal_g_128

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
 set val("leaves_states"), file("leaves_states.state")  into g_151_state_g_304, g_151_state_g_305

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
 set val(name), file("${name}.phy")  into g_129_phylip_g_130, g_129_phylip_g_145, g_129_phylip_g_189, g_129_phylip_g_279

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
 set val("iqtree"), file("iqtree.nwk")  into g_145_tree_g_132, g_145_tree_g_302
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

process rooting_iqtree_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(tree) from g_145_tree_g_302

output:
 set val("${name}_rooted"), file("*.nwk")  into g_302_tree_g_279

"""
nw_reroot -l $tree $params.OUTGRP 1>${name}_rooted.nwk

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
 set val(name), file(mulal) from g_129_phylip_g_279
 set val(namet), file(tree) from g_302_tree_g_279

output:
 set val("iqtree_anc_tree"), file("iqtree_anc_tree.nwk")  into g_279_tree_g_304
 set val("iqtree"), file("iqtree_anc.state")  into g_279_state_g_304
 file "*.log"  into g_279_logFile

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

ms_options = params.mutations_iqtree.ms_options
mnum192 = params.mutations_iqtree.mnum192
proba_min = params.mutations_iqtree.proba_min


process mutations_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "mutspec_tables/$filename"
	else if (filename =~ /.*.log$/) "mutspec_logs/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
}

input:
 set val(namet), file(tree) from g_279_tree_g_304
 set val(label), file(states1) from g_279_state_g_304
 set val(names2), file(states2) from g_151_state_g_304

output:
 file "*.tsv"  into g_304_outputFileTSV
 file "expected_mutations.txt"  into g_304_outputFileTxt
 file "*.log"  into g_304_logFile
 file "*.pdf"  into g_304_outputFilePdf

"""
python3 /export/src/mutspec-utils/scripts/3.collect_mutations.py --tree $tree --states $states1 --states $states2 --gencode $params.gencode $ms_options --outdir mout
mv mout/* .
mv mutations.tsv observed_mutations_${label}.tsv
mv expected_mutations.tsv expected_mutations.txt
mv run.log ${label}_run.log

python3 /export/src/mutspec-utils/scripts/calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.txt -o . --outgrp $params.OUTGRP -m $mnum192 -l $label -p $proba_min -x pdf

"""
}


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

run_RAXML = params.RAxML_build_tree.run_RAXML
raxml_model = params.RAxML_build_tree.raxml_model
//* @style @condition:{run_RAXML="true", raxml_model}


process RAxML_build_tree {

input:
 set val(name), file(mulal) from g_129_phylip_g_130

output:
 set val("raxml"), file("raxml.nwk")  into g_130_tree_g_109, g_130_tree_g_301

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
 set val(name), file(tree) from g_130_tree_g_301

output:
 set val("${name}_rooted"), file("*.nwk")  into g_301_tree_g_189

"""
nw_reroot -l $tree $params.OUTGRP 1>${name}_rooted.nwk

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
 set val(namet), file(tree) from g_301_tree_g_189
 set val(name), file(mulal) from g_129_phylip_g_189

output:
 set val("RAxML_nodeLabelledRootedTree"),file("RAxML_nodeLabelledRootedTree.nwk")  into g_189_tree_g_305
 set val("RAxML"), file("RAxML_marginalAncestralProbabilities.state")  into g_189_state_g_305
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

ms_options = params.mutations_raxml.ms_options
mnum192 = params.mutations_raxml.mnum192
proba_min = params.mutations_raxml.proba_min


process mutations_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "mutspec_tables/$filename"
	else if (filename =~ /expected_mutations.txt$/) "tmp/$filename"
	else if (filename =~ /.*.log$/) "mutspec_logs/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
}

input:
 set val(namet), file(tree) from g_189_tree_g_305
 set val(label), file(states1) from g_189_state_g_305
 set val(names2), file(states2) from g_151_state_g_305

output:
 file "*.tsv"  into g_305_outputFileTSV
 file "expected_mutations.txt"  into g_305_outputFileTxt
 file "*.log"  into g_305_logFile
 file "*.pdf"  into g_305_outputFilePdf

"""
python3 /export/src/mutspec-utils/scripts/3.collect_mutations.py --tree $tree --states $states1 --states $states2 --gencode $params.gencode $ms_options --outdir mout
mv mout/* .
mv mutations.tsv observed_mutations_${label}.tsv
mv expected_mutations.tsv expected_mutations.txt
mv run.log ${label}_run.log

python3 /export/src/mutspec-utils/scripts/calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.txt -o . --outgrp $params.OUTGRP -m $mnum192 -l $label -p $proba_min -x pdf

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
