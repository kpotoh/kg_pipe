$HOSTNAME = ""
params.outdir = 'results'  

THREADS = 1
//* params.mode =  "single"  //* @dropdown @options:"single, multiple" @description:"How to run pipeline, either on one record or on multiple records. WARNING: if you run mutational spectrum reconstruction on multiple records, you must be TODO"
//* params.OUTGRP =  "OUTGRP"   //* @input
//* params.run_pyvolve =  "true"  //* @dropdown @options:"true, false"


//* autofill
if ($HOSTNAME == "default"){
	$SINGULARITY_IMAGE = "/export/src/image_pipeline-2.9.sif"
	$SINGULARITY_OPTIONS = "--bind /export"
}

if ($HOSTNAME == "85.175.149.198"){
	$SINGULARITY_IMAGE = "/home/dolphin/image_pipeline-2.9.sif"
	$SINGULARITY_OPTIONS = "--bind /home/dolphin,/scratch/genkvg/mmseqs2_db_tax-specific:/db"
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
		echo "Query fasta must contain only one record when mode == 'single'"
		exit 1
	fi
elif [ $params.mode == "multiple" ]; then
	if [ `grep -c ">" $query` -gt 1 ]; then 
		if [ $params.OUTGRP == "PASS THE OUTGROUP NAME" ] || [ `grep -c $params.OUTGRP $query` -eq 0 ]; then
			echo "If you run 'multiple' mode you must pass the OUTGRP argument and it must be presented in the fasta file"
			exit 1
		else
			mv $query query_multiple.fasta
		fi
	else
		echo "Query fasta must contain more than one record when mode == 'multiple'"
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
/opt/dolphin/scripts/header_sel_mod3.pl $query_out_fasta $SPNAME 1>${name}_sel.fasta 2>${name}_sel.hash

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
/opt/dolphin/scripts/nuc_coding_mod.pl $hash $DB 1>sequences.fasta

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
	if (filename =~ /char_numbers.log$/) "report/$filename"
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
/opt/dolphin/scripts/codon_alig_unique.pl $seqs 1>seqs_unique.fasta
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

script:
if (report.toString() == "report_no.txt") {println "ERROR: Number of sequences too small to build tree"}

name = seqs.toString() - '.fasta'

"""
if [ -f report_no.txt ]; then 
	cat report_no.txt
	exit 1
fi

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
 set val("alignment_checked"),file("alignment_checked.fasta")  into g_128_nucl_mulal_g_70, g_128_nucl_mulal_g_129, g_128_nucl_mulal_g_151, g_128_nucl_mulal_g_374, g_128_nucl_mulal_g_375

"""
/opt/dolphin/scripts/macse2.pl $mulal alignment_checked.fasta

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
 set val("leaves_states"), file("leaves_states.state")  into g_151_state_g_380, g_151_state_g_381

"""
alignment2iqtree_states.py $mulal leaves_states.state
"""
}


process convert_alignment_to_phylip {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.phy$/) "IQTREE/$filename"
	else if (filename =~ /${name}.phy$/) "RAxML/$filename"
}

input:
 set val(name), file(mulal) from g_128_nucl_mulal_g_129

output:
 set val(name), file("${name}.phy")  into g_129_phylip_g_130, g_129_phylip_g_326, g_129_phylip_g_327, g_129_phylip_g_383

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
	if (filename =~ /iqtree.nwk$/) "tmp/$filename"
	else if (filename =~ /.*.log$/) "report/$filename"
}

input:
 set val(name), file(mulal) from g_129_phylip_g_383

output:
 set val("iqtree"), file("iqtree.nwk")  into g_383_tree_g_315
 file "*.log"  into g_383_logFile

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
run_RAXML = params.RAxML_build_tree.run_RAXML
raxml_model = params.RAxML_build_tree.raxml_model
//* @style @condition:{run_RAXML="true", raxml_model}


process RAxML_build_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /raxml.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(mulal) from g_129_phylip_g_130

output:
 set val("raxml"), file("raxml.nwk")  into g_130_tree_g_317

when:
run_RAXML == "true"

script:
"""
/opt/standard-RAxML-8.2.12/raxmlHPC-PTHREADS-SSE3 -x 987654 -p 987654 -T $THREADS -N 50 -f a -s $mulal -n Rax_tree -m $raxml_model
mv RAxML_bestTree.Rax_tree raxml.nwk
"""

}
params.RAxML_model = raxml_model
quantile = params.shrink_raxml.quantile


process shrink_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_shrinked.nwk$/) "tmp/$filename"
	else if (filename =~ /.*.log$/) "RAxML/$filename"
}

input:
 set val(name), file(tree) from g_130_tree_g_317

output:
 set val("${name}_shrinked"), file("${name}_shrinked.nwk")  into g_317_tree_g_301, g_317_tree_g_109
 file "*.log"  into g_317_logFile

"""
if [ `nw_stats $tree | grep nodes | cut -f 2` -gt 8 ]; then
	run_treeshrink.py -t $tree -O treeshrink -o . -q $quantile -x $params.OUTGRP
	mv treeshrink.nwk ${name}_shrinked.nwk
	mv treeshrink_summary.txt ${name}_treeshrink.log
	mv treeshrink.txt ${name}_pruned_nodes.log
else
	cat $tree > ${name}_shrinked.nwk
	echo "Shrinjing are useless on such small number of sequences" > ${name}_pruned_nodes.log
fi
"""
}


process extract_terminal_branch_lengths_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "RAxML/$filename"
}

input:
 set val(name), file(tree) from g_317_tree_g_109

output:
 set val(name), file("${name}.branches")  into g_109_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1> ${name}.branches

if [ `grep OUTGRP ${name}.branches | cut -f 2 | python3 -c "import sys; print(float(sys.stdin.readline().strip()) > 0)"` == False ]; then
	cat "${name}.branches"
	echo "Something went wrong: outgroup is not furthest leaf in the tree"
	exit 1
fi

"""
}


process rooting_raxml_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(tree) from g_317_tree_g_301

output:
 set val("${name}_rooted"), file("*.nwk")  into g_301_tree_g_327

"""
nw_reroot -l $tree $params.OUTGRP 1>${name}_rooted.nwk

"""
}


process RAxML_anc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /RAxML_nodeLabelledRootedTree.nwk$/) "RAxML/$filename"
	else if (filename =~ /RAxML_marginalAncestralProbabilities.state$/) "RAxML/$filename"
	else if (filename =~ /RAxML_anc_rec.log$/) "report/$filename"
}

input:
 set val(namet), file(tree) from g_301_tree_g_327
 set val(name), file(mulal) from g_129_phylip_g_327

output:
 set val("raxml"),file("RAxML_nodeLabelledRootedTree.nwk")  into g_327_tree_g_375, g_327_tree_g_381
 set val("raxml"), file("RAxML_marginalAncestralProbabilities.state")  into g_327_state_g_381
 file "RAxML_marginalAncestralStates.fasta"  into g_327_multipleFasta
 file "RAxML_anc_rec.log"  into g_327_logFile

"""
raxmlHPC-PTHREADS-SSE3 -T $THREADS -f A -m $params.RAxML_model -s $mulal -t $tree -n ANCESTORS
mv RAxML_info.ANCESTORS RAxML_anc_rec.log
mv RAxML_marginalAncestralStates.ANCESTORS RAxML_marginalAncestralStates.fasta
# mv RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.txt
# mv RAxML_nodeLabelledRootedTree.ANCESTORS RAxML_nodeLabelledRootedTree.nwk

raxml_states2iqtree_states.py RAxML_marginalAncestralProbabilities.ANCESTORS RAxML_marginalAncestralProbabilities.state
rename_internal_nodes.py $tree RAxML_nodeLabelledRootedTree.ANCESTORS RAxML_nodeLabelledRootedTree.nwk

"""
}

syn4f = params.mutations_raxml.syn4f
use_probabilities = params.mutations_raxml.use_probabilities
mnum192 = params.mutations_raxml.mnum192
proba_min = params.mutations_raxml.proba_min
//* @style @multicolumn:{syn4f, use_probabilities}, {mnum192, proba_min}

syn4f = syn4f == "true" ? "--syn4f" : ""
proba = use_probabilities == "true" ? "--proba" : ""

process mutations_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "mutspec_tables/$filename"
	else if (filename =~ /.*.log$/) "report/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
	else if (filename =~ /ms12syn_${label}.txt$/) "tmp/$filename"
}

input:
 set val(namet), file(tree) from g_327_tree_g_381
 set val(label), file(states1) from g_327_state_g_381
 set val(names2), file(states2) from g_151_state_g_381

output:
 file "*.tsv"  into g_381_outputFileTSV
 file "*.log"  into g_381_logFile
 file "*.pdf"  into g_381_outputFilePdf
 file "ms12syn_${label}.txt"  into g_381_outputFileTxt_g_375

"""
3.collect_mutations.py --tree $tree --states $states1 --states $states2 \
	--gencode $params.gencode --syn $syn4f $proba --no-mutspec --outdir mout

mv mout/* .
mv mutations.tsv observed_mutations_${label}.tsv
mv expected_mutations.tsv expected_mutations.txt
mv run.log ${label}_mut_extraction.log

calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.txt -o . \
	--exclude ${params.OUTGRP},ROOT --mnum192 $mnum192 -l $label $proba --proba_min $proba_min --syn $syn4f --plot -x pdf

cp ms12syn_${label}.tsv ms12syn_${label}.txt

"""
}

when:
params.run_pyvolve == "true"

replics = params.pyvolve_raxml.replics
scale_tree = params.pyvolve_raxml.scale_tree
syn4f = params.pyvolve_raxml.syn4f
mnum192 = params.pyvolve_raxml.mnum192
//* @style @multicolumn:{replics, scale_tree}, {syn4f, mnum192}

syn4f = syn4f == "true" ? "--syn4f" : ""


process pyvolve_raxml {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "simulation/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
	else if (filename =~ /.*.log$/) "simulation/$filename"
}

input:
 file spectra from g_381_outputFileTxt_g_375
 set val(label), file(tree) from g_327_tree_g_375
 set val(name), file(mulal) from g_128_nucl_mulal_g_375

output:
 file "*.tsv"  into g_375_outputFileTSV
 file "*.pdf"  into g_375_outputFilePdf
 file "*.log"  into g_375_logFile

"""
nw_prune $tree $params.OUTGRP | python3 /home/dolphin/dolphin/scripts/resci.py > ${tree}.ingroup
echo "Tree outgroup pruned"
awk '/^>/ {P=index(\$1, "OUTGRP")==0} {if(P) print}' $mulal > ${mulal}.ingroup
echo "Tree outgroup sequence filtered out"

nw_labels -L ${tree}.ingroup | grep -v ROOT | xargs -I{} echo -e "{}\tNode{}" > map.txt
if [ `grep -c NodeNode map.txt` -eq 0 ]; then
	nw_rename ${tree}.ingroup map.txt > ${tree}.tmp
	cat ${tree}.tmp > ${tree}.ingroup
	echo "Internal nodes renamed"
fi

#filter out sequences with ambigous nucleotides
cat ${mulal}.ingroup | perl -e '\$p=\$s="777"; while (<STDIN>) {chomp; if (\$_=~/^>/) {\$h=\$_; if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"} \$p=\$h; \$s="777"} else {\$s=\$_}} if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"}' > ${mulal}.clean
cat ${mulal}.clean 
if [ `grep -c ">" ${mulal}.clean` -lt 1 ]; then
	echo -e "There are no sequences without ambigous nucleotides in the alignment.\nInterrupted" > pyvolve_${label}.log
	exit 0
fi

pyvolve_process.py -a ${mulal}.clean -t ${tree}.ingroup -s $spectra -o seqfile.fasta -r $replics -c $params.gencode -l $scale_tree --write_anc
echo "Mutation samples generated"

for fasta_file in seqfile_sample-*.fasta
do
	echo "Processing \$fasta_file"
	alignment2iqtree_states.py \$fasta_file  \${fasta_file}.state
	3.collect_mutations.py --tree ${tree}.ingroup --states  \${fasta_file}.state --gencode $params.gencode --syn $syn4f --no-mutspec --outdir mout --force
	cat mout/run.log >> pyvolve_${label}.log
	echo -e "\n\n">> pyvolve_${label}.log
	cat mout/mutations.tsv >  \${fasta_file}.mutations
done
echo "Mutations extraction done"

concat_mutations.py seqfile_sample-*.fasta.mutations mutations_${label}_pyvolve.tsv
echo "Mutations concatenation done"
cat mutations_${label}_pyvolve.tsv

calculate_mutspec.py -b mutations_${label}_pyvolve.tsv -e mout/expected_mutations.tsv -o . \
	-l ${label}_pyvolve --exclude ${params.OUTGRP},ROOT --syn $syn4f --mnum192 $mnum192 --plot -x pdf
echo "Mutational spectrum calculated"

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
/opt/dolphin/scripts/mutnumbers.pl $seqs 1>mutnumbers.tsv

"""
}

quantile = params.shrink_iqtree.quantile


process shrink_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}_shrinked.nwk$/) "tmp/$filename"
	else if (filename =~ /.*.log$/) "IQTREE/$filename"
}

input:
 set val(name), file(tree) from g_383_tree_g_315

output:
 set val("${name}_shrinked"), file("${name}_shrinked.nwk")  into g_315_tree_g_302, g_315_tree_g_132
 file "*.log"  into g_315_logFile

"""
if [ `nw_stats $tree | grep nodes | cut -f 2` -gt 8 ]; then
	run_treeshrink.py -t $tree -O treeshrink -o . -q $quantile -x $params.OUTGRP
	mv treeshrink.nwk ${name}_shrinked.nwk
	mv treeshrink_summary.txt ${name}_treeshrink.log
	mv treeshrink.txt ${name}_pruned_nodes.log
else
	cat $tree > ${name}_shrinked.nwk
	echo "Shrinjing are useless on such small number of sequences" > ${name}_pruned_nodes.log
fi
"""
}


process extract_terminal_branch_lengths_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /${name}.branches$/) "IQTREE/$filename"
}

input:
 set val(name), file(tree) from g_315_tree_g_132

output:
 set val(name), file("${name}.branches")  into g_132_branches

"""
nw_distance -m p -s l -n $tree | sort -grk 2 1> ${name}.branches

if [ `grep OUTGRP ${name}.branches | cut -f 2 | python3 -c "import sys; print(float(sys.stdin.readline().strip()) > 0)"` == False ]; then
	cat "${name}.branches"
	echo "Something went wrong: outgroup is not furthest leaf in the tree"
	exit 1
fi

"""
}


process rooting_iqtree_tree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.nwk$/) "tmp/$filename"
}

input:
 set val(name), file(tree) from g_315_tree_g_302

output:
 set val("${name}_rooted"), file("*.nwk")  into g_302_tree_g_326

"""
nw_reroot -l $tree $params.OUTGRP 1>${name}_rooted.nwk

"""
}


process IQTREE_anc {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /iqtree_anc_tree.nwk$/) "IQTREE/$filename"
	else if (filename =~ /iqtree_anc.state$/) "IQTREE/$filename"
	else if (filename =~ /.*.log$/) "report/$filename"
}

input:
 set val(name), file(mulal) from g_129_phylip_g_326
 set val(namet), file(tree) from g_302_tree_g_326

output:
 set val("iqtree"), file("iqtree_anc_tree.nwk")  into g_326_tree_g_374, g_326_tree_g_380
 set val("iqtree"), file("iqtree_anc.state")  into g_326_state_g_380
 file "*.log"  into g_326_logFile

errorStrategy 'retry'
maxRetries 3

script:
"""
iqtree2 -s $mulal -m $params.IQTREE_model -asr -nt $THREADS --prefix anc
mv anc.iqtree iqtree_anc_report.log
# mv anc.state iqtree_anc.state
mv anc.log iqtree_anc.log
nw_reroot anc.treefile OUTGRP | sed 's/;/ROOT;/' > iqtree_anc_tree.nwk

iqtree_states_add_part.py anc.state iqtree_anc.state
"""

}

syn4f = params.mutations_iqtree.syn4f
use_probabilities = params.mutations_iqtree.use_probabilities
mnum192 = params.mutations_iqtree.mnum192
proba_min = params.mutations_iqtree.proba_min
//* @style @multicolumn:{syn4f, use_probabilities}, {mnum192, proba_min}

syn4f = syn4f == "true" ? "--syn4f" : ""
proba = use_probabilities == "true" ? "--proba" : ""

process mutations_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "mutspec_tables/$filename"
	else if (filename =~ /.*.log$/) "report/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
	else if (filename =~ /ms12syn_${label}.txt$/) "tmp/$filename"
}

input:
 set val(namet), file(tree) from g_326_tree_g_380
 set val(label), file(states1) from g_326_state_g_380
 set val(names2), file(states2) from g_151_state_g_380

output:
 file "*.tsv"  into g_380_outputFileTSV
 file "*.log"  into g_380_logFile
 file "*.pdf"  into g_380_outputFilePdf
 file "ms12syn_${label}.txt"  into g_380_outputFileTxt_g_374

"""
3.collect_mutations.py --tree $tree --states $states1 --states $states2 \
	--gencode $params.gencode --syn $syn4f $proba --no-mutspec --outdir mout

mv mout/* .
mv mutations.tsv observed_mutations_${label}.tsv
mv expected_mutations.tsv expected_mutations.txt
mv run.log ${label}_mut_extraction.log

calculate_mutspec.py -b observed_mutations_${label}.tsv -e expected_mutations.txt -o . \
	--exclude ${params.OUTGRP},ROOT --mnum192 $mnum192 -l $label $proba --proba_min $proba_min --syn $syn4f --plot -x pdf

cp ms12syn_${label}.tsv ms12syn_${label}.txt

"""
}

when:
params.run_pyvolve == "true"

replics = params.pyvolve_iqtree.replics
scale_tree = params.pyvolve_iqtree.scale_tree
syn4f = params.pyvolve_iqtree.syn4f
mnum192 = params.pyvolve_iqtree.mnum192
//* @style @multicolumn:{replics, scale_tree}, {syn4f, mnum192}

syn4f = syn4f == "true" ? "--syn4f" : ""


process pyvolve_iqtree {

publishDir params.outdir, overwrite: true, mode: 'copy',
	saveAs: {filename ->
	if (filename =~ /.*.tsv$/) "simulation/$filename"
	else if (filename =~ /.*.pdf$/) "mutspec_images/$filename"
	else if (filename =~ /.*.log$/) "simulation/$filename"
}

input:
 file spectra from g_380_outputFileTxt_g_374
 set val(label), file(tree) from g_326_tree_g_374
 set val(name), file(mulal) from g_128_nucl_mulal_g_374

output:
 file "*.tsv"  into g_374_outputFileTSV
 file "*.pdf"  into g_374_outputFilePdf
 file "*.log"  into g_374_logFile

"""
nw_prune $tree $params.OUTGRP | python3 /home/dolphin/dolphin/scripts/resci.py > ${tree}.ingroup
echo "Tree outgroup pruned"
awk '/^>/ {P=index(\$1, "OUTGRP")==0} {if(P) print}' $mulal > ${mulal}.ingroup
echo "Tree outgroup sequence filtered out"

nw_labels -L ${tree}.ingroup | grep -v ROOT | xargs -I{} echo -e "{}\tNode{}" > map.txt
if [ `grep -c NodeNode map.txt` -eq 0 ]; then
	nw_rename ${tree}.ingroup map.txt > ${tree}.tmp
	cat ${tree}.tmp > ${tree}.ingroup
	echo "Internal nodes renamed"
fi

#filter out sequences with ambigous nucleotides
cat ${mulal}.ingroup | perl -e '\$p=\$s="777"; while (<STDIN>) {chomp; if (\$_=~/^>/) {\$h=\$_; if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"} \$p=\$h; \$s="777"} else {\$s=\$_}} if (\$s!~/[^ACGT]/i) {print "\$p\n\$s\n"}' > ${mulal}.clean
cat ${mulal}.clean 
if [ `grep -c ">" ${mulal}.clean` -lt 1 ]; then
	echo -e "There are no sequences without ambigous nucleotides in the alignment.\nInterrupted" > pyvolve_${label}.log
	exit 0
fi

pyvolve_process.py -a ${mulal}.clean -t ${tree}.ingroup -s $spectra -o seqfile.fasta -r $replics -c $params.gencode -l $scale_tree --write_anc
echo "Mutation samples generated"

for fasta_file in seqfile_sample-*.fasta
do
	echo "Processing \$fasta_file"
	alignment2iqtree_states.py \$fasta_file  \${fasta_file}.state
	3.collect_mutations.py --tree ${tree}.ingroup --states  \${fasta_file}.state --gencode $params.gencode --syn $syn4f --no-mutspec --outdir mout --force
	cat mout/run.log >> pyvolve_${label}.log
	echo -e "\n\n">> pyvolve_${label}.log
	cat mout/mutations.tsv >  \${fasta_file}.mutations
done
echo "Mutations extraction done"

concat_mutations.py seqfile_sample-*.fasta.mutations mutations_${label}_pyvolve.tsv
echo "Mutations concatenation done"
cat mutations_${label}_pyvolve.tsv

calculate_mutspec.py -b mutations_${label}_pyvolve.tsv -e mout/expected_mutations.tsv -o . \
	-l ${label}_pyvolve --exclude ${params.OUTGRP},ROOT --syn $syn4f --mnum192 $mnum192 --plot -x pdf
echo "Mutational spectrum calculated"

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
