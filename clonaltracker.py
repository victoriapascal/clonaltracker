import sys
import subprocess
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from shutil import copyfile
import shutil
import pickle
#import pickle5 as pickle

def check_input_is_fasta(fasta):
	'''Checks that the input is a valid nucleotide fasta sequence'''
	check = ''
	with open(fasta, 'r') as f:
		lines = f.readlines()
		if lines[0].startswith('>') and all(base.upper() in ('A', 'C', 'T', 'G', 'N') for base in lines[1].strip()) ==True:
			check = True
		else:
			check = False
	return check

def run_blastn(db, out_f, fasta):
	'''Run blastn on any given fatsa sequence using the van gene DB'''
	fa = fasta.split('/')[-1]
	print('1. Running blastn for %s...' %(fa))
	output_file = open(out_f, 'w')
	cmd = ['blastn', '-query', str(fasta), '-db', db, '-outfmt', '6', '-evalue', '1e-10']
	torun = subprocess.Popen(cmd, stdout=output_file, stderr=subprocess.PIPE)
	out, err = torun.communicate()
	output_file.close()


def check_blastn_output(out_dir, fasta):
	'''Read the blastn output and choose the van type based on the best candidate gene'''
	output_file = out_dir + os.sep + 'blastn_van_' + str(fasta)
	type_pident = {} #the key is the van gene and the value is a list of region coordinates
	with open(output_file, 'r') as f:
		for l in f:
			line = l.strip().split('\t')
			evalue = float(line[10])
			ident = float(line[2])
			if evalue < 1e-05 and ident > 95:
				record = line[0]
				if line[1] in type_pident.keys():
					entry = [ int(line[6]), int(line[7]), record]
					type_pident[line[1]].append(entry)
				else:
					type_pident[line[1]] = []
					entry = [int(line[6]), int(line[7]), record]
					type_pident[line[1]].append(entry)
	
	return type_pident

def run_poppunk(list_genomes, poppunk_db, out_dir):
	'''Run poppunk usng the vanAB=new_outbreak database and the two input genomes'''
	print('2. Running PopPUNK to cluster inputted genomes into a larger vanA and vanB database')
	#pp_out = out_dir + os.sep + 'poppunk_analysis'
	#os.mkdir(out_dir + os.sep + 'poppunk_analysis')
	cwd = os.getcwd()
	os.chdir(out_dir)
	
	##update DB
	cmd = ['poppunk_assign', '--db', 'vanAB_dataset', '--query', 'list_new_genomes.txt', '--output', 'vanAB_dataset_updated', '--update-db']
	torun = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = torun.communicate()
	
	##fit the model
	cmd2 = ['poppunk', '--fit-model', 'dbscan', '--ref-db', 'vanAB_dataset_updated']
	torun2 = subprocess.Popen(cmd2, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out2, err2= torun2.communicate()

	##visualize
	cmd3 = ['poppunk_visualise', '--ref-db', 'vanAB_dataset_updated', '--output', 'vanAB_dataset_updated/microreact_viz', '--microreact']
	torun3 = subprocess.Popen(cmd3, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out3, err3= torun3.communicate()
	if err:
		print(err)
	if err2:
		print(err2)
	if err3:
		print(err3)
	os.chdir(cwd)

def run_tetyper(fasta, out_dir, type_pident, location_reads):
	'''Running TETyper to assess transposon type'''
	print('3. Running TETyper to assess transposon type')
	sample = fasta.split('/')[-1].split('.')[0] + "_"
	fq1 = [file for file in os.listdir(location_reads) if 'R1.fastq.gz' in file and sample in file][0]
	fq2 = [file for file in os.listdir(location_reads) if 'R2.fastq.gz' in file and sample in file][0]
	#fq1 = [file for file in os.listdir(location_reads) if '_R1_val_1.fq' in file and sample in file][0]
	#fq2 = [file for file in os.listdir(location_reads) if '_R2_val_2.fq' in file and sample in file][0]
	tnp_ref = ''
	script_dir = os.path.abspath( os.path.dirname( __file__ ))
	for a in type_pident.keys():
		if "vanA" in a:
			tnp_ref = script_dir + os.sep + 'van_representatives/tnp_db/vanA_M97297.fa'
		else:
			if "vanB" in a:
				tnp_ref = script_dir + os.sep + 'van_representatives/tnp_db/vanB_AY655721.fa'
	if not os.path.isdir(out_dir + os.sep + 'TETyper_out'):
		os.mkdir(out_dir + os.sep + 'TETyper_out')
	else:
		pass
	tet_out = out_dir + os.sep + 'TETyper_out/' + fasta.split('/')[-1].split('.')[0]
	print(tet_out)
	cmd = ['TETyper.py', '--ref', tnp_ref, '--assembly', fasta, '--outprefix', tet_out, '--flank_len', '5', '--fq1', location_reads + fq1, '--fq2', location_reads + fq2, ] 
	torun = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	out, err = torun.communicate()
	if err:
		print(err)
def parse_tetyper_output(output_dir):
	'''parse TETyper output to check if both genomes have the same transposon type'''
	snps = {}
	names = []
	tet_out = output_dir + os.sep + 'TETyper_out/'
	summary = [files for files in os.listdir(tet_out) if 'summary' in files]
	for f in summary:
		with open(tet_out + f, 'r') as f1:
			next(f1)
			name = f.split('_')[0]
			names.append(name)
			snps[name] = []
			for line in f1:
				lines = line.strip().split('\t')
				snps[name].append(lines[0])
				snps[name].append(lines[2])
#				snps[name].append(lines[3])

	result = ''
	if set(snps[names[0]]) == set(snps[names[1]]):
		result = True
	else:
		result = False
	
	return result

def get_contigs_seqs(output_dir,fna):
	'''
	From the TETyper blast output, get the transposons that align with the transposon
	'''
	#files = output_dir + os.sep + 'TETyper_out/'
	bf = [f for f in os.listdir(output_dir) if 'contig' in f]
	sample_contigs = {}
	for b in bf:
		name = b.replace('blast_contigs_', '').replace('.fasta', '').replace('.fa', '').replace('.fna', '')
		sample_contigs[name] = []
		with open(output_dir + os.sep + b, 'r') as f1:
			for line in f1:
				line = line.strip().split('\t')
				if not line[0] in sample_contigs[name]:
					sample_contigs[name].append(line[0])
	
	for a in sample_contigs.keys():
		print(a)
		sample = [f for f in os.listdir(fna) if a + '.f' in f][0]
		print(sample)
		with open(fna + os.sep + sample, 'r') as f2:
			lines = f2.readlines()
			out = open(output_dir + os.sep + a + "_tnp_contig.fa", 'w')
			header = ''
			for c in sample_contigs[a]:
				header = '>' + c + '\n'
				index = lines.index(header)
				for entry in lines[index:index+2]:
					out.write(entry)
			out.close()
	return sample_contigs

def run_ragtag(output_dir, sample_contigs, van_type):
	'''
	Run ragtag to scaffold contig only when the transposon is found split in more than one contigs	
	'''
	for sample, contigs in sample_contigs.items():
		script_dir = os.path.abspath( os.path.dirname( __file__ ))
		if len(contigs) > 1: ## when there the tnp is split into multiple contigs, run ragtag
			out = output_dir + os.sep + sample + "_ragtag"
			ref_tn = script_dir + os.sep +'van_representatives/tnp_db/' + [f for f in os.listdir(script_dir + os.sep + 'van_representatives/tnp_db/') if van_type in f and f.endswith(".fa")][0]
			isolate = sample + '_tnp_contig.fa'
			query = output_dir + os.sep + isolate
			#query = output_dir + os.sep + [f for f in os.listdir(output_dir) if sample in f and '_tnp_contig.fa' in f][0]
			cmd = ['ragtag.py', 'scaffold','-o', out, ref_tn, query]
			torun = subprocess.Popen(cmd,stdout=subprocess.PIPE, stderr=subprocess.PIPE)			
			out1, err = torun.communicate()
			with open(str(out) + os.sep + 'ragtag.scaffold.fasta', 'r') as f:
				lines = f.readlines()
			try:
				os.remove(output_dir + os.sep + sample + "_tnp_contig.fa") ##remove existing file
				os.remove(output_dir + os.sep + sample + "_tnp_contig.fa.fai")
			except:
				pass
			with open(output_dir + os.sep + sample + "_tnp_contig.fa", 'w') as f2:
				for num, line in zip(range(len(lines)), lines):
					if line.startswith(">"):
						header = ">" + sample + "_contig" + str(num) + "_ragtag" + '\n'
						f2.write(header)
					else:
						f2.write(line.replace("N", ""))
						  
					
				
def run_isescan(out_dir):
	'''
	Run isescan to assess potential tnp IS
	'''
	print("4. Running ISEScan to predict ISs")
	is_out = output_dir + os.sep + "ISEScan_trim"
	if not os.path.isdir(is_out):
		os.mkdir(is_out)
	contigf = [file for file in os.listdir(out_dir) if "_tnp_trimmed.fa" in file and not '.fa.n' in file]
	for cf in contigf:
		od = cf.replace("_tnp_trimmed.fa", '')
		cmd = ['isescan.py', '--seqfile', out_dir + os.sep+ cf, '--output', is_out + os.sep + od]
		torun = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = torun.communicate()

def parse_isescan_output(out_dir):
	'''
	Check ISEScan output folder to see if the two genomes have the same (if any) IS
	'''
	samples_id = {}
	#samples = out_dir.split('_')[2].split('-')
	samples = [f.replace('_tnp_contig.fa', '') for f in os.listdir(out_dir) if 'tnp_contig.fa' in f and '.fai' not in f]
	for sample in samples:
		folders = os.listdir(out_dir + os.sep + "ISEScan_trim/" + sample)
		if out_dir in folders:
			root = out_dir + os.sep + "ISEScan_trim/" + sample 
			table = [file for file in os.listdir(root + os.sep + out_dir) if file.endswith('.tsv')]
			samples_id[sample] = []
			with open(root + os.sep + out_dir + os.sep + table[0], 'r') as f:
				next(f)
				for line in f:
					cluster_id = line.strip().split('\t')[2]
					samples_id[sample].append(cluster_id)
		else:
			samples_id[sample] = ['None']

	match = ''

	if set(sorted(samples_id[samples[0]])) == set(sorted(samples_id[samples[1]])):
		match = True
	else:
		match = False
	
	print("%s has this IS: %s and %s has this IS: %s" %(samples[0], ','.join(samples_id[samples[0]]), samples[1], ','.join(samples_id[samples[1]])))	
	return match

def assess_clonality_with_mash_sketch(out_dir, fasta, fasta2):
	'''
	When the genomes have the same van type and surrounding regions, assess the clonality usng MASH
	'''
	print('5. Running MASH sketch to assess clonality (whole genome similarity)')
	cmd = ['mash', 'sketch', '-k', '32', '-s', '8000', '-o', out_dir + os.sep + 'mash_sketch', fasta, fasta2]
	torun = subprocess.Popen(cmd, stdout=None, stderr=subprocess.PIPE)
	out, err = torun.communicate()

def assess_clonality_with_mash_dist(out_dir):

	print('6. Running MASH dist to get MASH distances')
	output_file = open(out_dir + os.sep + 'mash_dist.txt', 'w')
	cmd = ['mash', 'dist', out_dir + os.sep+ 'mash_sketch.msh', out_dir + os.sep+ 'mash_sketch.msh']
	torun = subprocess.Popen(cmd, stdout=output_file, stderr=subprocess.PIPE)
	out, err = torun.communicate()
	output_file.close()

def evaluate_mash_distances(out_dir):
	sim = set()
	prop = set()
	with open(out_dir + os.sep +'mash_dist.txt', 'r') as out:
		for line in out:
			line = line.strip().split('\t')
			if not line[0] == line[1]:
				sim.add(float(line[2]))
				coef = line[4].split("/")
				prop.add(float(int(coef[0])/int(coef[1])))
	for a in list(prop):
		if a > 0.95:
			print('These two genomes are highly similar and might be clonal. Their MASH distance is %s' %(a))
		else:
			print('These two genomes are unlikely to be clonal as their MASH distance is %s' %(a))

def create_gbk_from_isescan_out(out_dir):
	'''
	From ISEScan ISs, ORFs and gene predictions build a genbank find to use it as input for clinker
	'''
	
	samples = os.listdir(out_dir + os.sep + 'ISEScan_trim')
	if not os.path.isdir(out_dir + os.sep + "tn_trim_gbks"):
		os.mkdir(out_dir + os.sep + "tn_trim_gbks")	
	for sample in samples:
		proteome_dir = out_dir + os.sep + 'ISEScan_trim' + os.sep + sample + os.sep + 'proteome' + os.sep + out_dir
		proteome_files = os.listdir(proteome_dir)
		genes = [f for f in proteome_files if '_tnp_trimmed.fa.ffn' in f][0]
		prots = [f for f in proteome_files if '_tnp_trimmed.fa.faa' in f][0]
		is_dir = out_dir + os.sep + 'ISEScan_trim' + os.sep + sample + os.sep + out_dir
		if os.path.isdir(is_dir):
			is_files = os.listdir(is_dir)
			is_nuc = [f for f in is_files if '_tnp_trimmed.fa.is.fna' in f][0]
			is_prot = [f for f in is_files if '_tnp_trimmed.fa.orf.faa' in f][0] 	
		
		records = {}
		##genes with transposon
		with open(proteome_dir + os.sep + genes, 'r') as f:
			 nuc = f.readlines()

		##proteins with transposon
		with open(proteome_dir + os.sep + prots, 'r') as f2:
 			prot = f2.readlines()
		
		for p in prot:
			 if '>' in p:
  				p_index = prot.index(p)
  				records[p.strip()] = [prot[p_index+1].strip()]

		if os.path.isdir(is_dir):
			##insertion nucleotide sequence
			with open(is_dir + os.sep + is_nuc, 'r') as f4:
				 isn = f4.readlines()

			##insertion protein sequence
			with open(is_dir + os.sep + is_prot, 'r') as f3:
				 isp = f3.readlines()
		
			for i in isp:
				 if ">" in i:
  					isp_index = isp.index(i)
  					records[i.strip()].append(isp[isp_index+1].strip())
		
		contigs = out_dir + os.sep + sample + '_tnp_trimmed.fa'
		header = []
		seq = []
		with open(contigs, 'r') as fasta:
			for line in fasta:
				if not line.startswith(">"):
					seq.append(line.strip())
				else:
   					header.append(line.strip())

		seqs = Seq(''.join(seq))
		#seqs = Seq(''.join(seq), alphabet=IUPAC.unambiguous_dna)
		record = SeqRecord(seqs, annotations={"molecule_type": "DNA"})

		for i in records.keys():
			start = int(i.split(' ')[0].split('_')[-3])
			end = int(i.split(' ')[0].split('_')[-2])
			strand = int(i.split(' ')[0].split('_')[-1].replace('+', '1').replace('-', '-1'))
			my_feature = SeqFeature(FeatureLocation(start, end, strand=strand), type="CDS", id=i)
			my_feature.qualifiers['translation'] = records[i][0]
			record.features.append(my_feature)

		with open(out_dir + os.sep + 'tn_trim_gbks/' + sample + '.gb', 'w') as input_handle:
			SeqIO.write(record, input_handle, "genbank")
def run_clinker(out_dir):
	'''
	From the the tn gbk files, run clincker to visualize the similarities/differences between the isolates transposons
	'''

	output_file1 =  out_dir + os.sep + 'clinker_tn_trim_viz.html'
	output_file2 = out_dir + os.sep + 'clinker_tn_trim_msa.aln'
	gbks = os.listdir(out_dir + os.sep + 'tn_trim_gbks')
	inp = out_dir + '/tn_trim_gbks/*'
	cmd = 'clinker '  + inp + ' -p '+ output_file1 + ' -o ' + output_file2 + ' -i ' + '0.01'
	torun = subprocess.Popen([cmd], stderr=subprocess.PIPE, shell=True)
	err = torun.communicate()

def trim_tnp_region(blast_tn, contig_fasta, out_tn):
	'''
	From the transposon blastn extract the tn coordinats and trim the transposon sequence
	'''
	seq_coords = {}	
	with open(blast_tn, 'r') as f:
		lines = f.readlines()
		isolate = blast_tn.split('/')[-1].replace('blastn_tnp_', '').replace('.fasta', '').replace('.fna', '').replace('.fa', '')
		if len(lines) == 1:
			contig = lines[0].split('\t')[0]
			start = int(lines[0].split('\t')[6])
			end = int(lines[0].split('\t')[7])
			seq_coords[isolate] = [contig, start, end]
		else:
			contig_coords = {}
			for line in lines:
				line = line.split('\t')
				c = line[0]
				if c in contig_coords.keys():
					contig_coords[c].append(int(line[6]))
					contig_coords[c].append(int(line[7]))
				else:	
					contig_coords[c] = []
					contig_coords[c].append(int(line[6]))
					contig_coords[c].append(int(line[7]))
			c_l = 0
			win_c = ''
			for contig, coords in contig_coords.items():
				start = min(coords)
				end = max(coords)
				length = abs(end -start)
				if length > c_l:
					c_l = length
					win_c = contig
			seq_coords[isolate] = [win_c, min(contig_coords[win_c]), max(contig_coords[win_c])]

	##trim the sequence
	iso = contig_fasta.split('/')[-1].replace('_tnp_contig.fa', '')
	with open(contig_fasta, 'r') as f:
		records = f.readlines()
		for r in records:
			if r.startswith(">") and seq_coords[iso][0] in r:
				header = records.index(r)
				region = records[header + 1][int(seq_coords[iso][1]):int(seq_coords[iso][2])]
				out = out_tn
				with open(out, 'w') as outf:
					outf.write(r)
					outf.write(region)


def parse_alignment_scores(out_dir):
	aln_scores = {}
	genes = set()
	selected_scores = []
	selected_genes = set()
	aln = [file for file in os.listdir(out_dir) if file.endswith('.aln')][0]
	with open(out_dir + os.sep + aln, 'r') as f:
		next(f)
		next(f)
		next(f)
		for line in f:
			lines = line.split("  ")
			pair = lines[0] + '*' + lines[1]
			genes.add(lines[0])
			genes.add(lines[1])
			iden = float(lines[2])
			if iden > 0.3:
				aln_scores[pair] = iden
				selected_scores.append(iden)
				selected_genes.add(lines[0])
				selected_genes.add(lines[1])
	
	average = float(sum(selected_scores)/len(selected_scores))
	considered_genes = len(list(selected_genes)) - len(list(genes))
	return average, abs(considered_genes)
	
def compare_tn_sequence_with_DB(out_dir):
	'''
	Compare the trimmed transposon sequence to the vanAB DB ones
	'''
	tn_seqs = [f for f in os.listdir(out_dir) if '_tnp_trimmed' in f and f.endswith('.fa')]
	
	##load tn DB
	script_dir = os.path.abspath( os.path.dirname( __file__ ))
	with open(script_dir + os.sep + 'van_representatives/tnp_db/tnp_db.pickle', 'rb') as handle:
                tnp_db = pickle.load(handle)

	with open(script_dir + os.sep + 'van_representatives/tnp_db/tnp_nomen.pickle', 'rb') as handle2:
                tn_nomen = pickle.load(handle2)
	
	for t in tn_seqs:
		iso = t.replace('_tnp_trimmed.fa', '')
		with open(out_dir + os.sep + t, 'r') as f:
			for line in f:
				if not line.startswith('>'):
					seq = line.strip()
					comp = [a for a in tnp_db.keys() if seq == a]	
					if comp:
						print('Isolate %s has the transposon type %s, also found in isolate %s' %(iso, tn_nomen[seq], ', '.join(tnp_db[seq])))


def makeblastdb(out_dir):
	'''
	Makeblastdb from one of the tnp reference sequences
	'''
	#tnp_seqs = [file for file in os.listdir(out_dir) if '_tnp_trimmed.fa' in file and not '.fa.n' in file]
	ref = out_dir + os.sep + [file.replace('_tnp_contig.fa', '') for file in os.listdir(out_dir) if '_tnp_contig.fa' in file][0] + '_tnp_trimmed.fa'
	#ref = out_dir + os.sep + tnp_seqs[0]
	cmd = ['makeblastdb', '-in', ref, '-out', ref, '-dbtype', 'nucl']
	torun = subprocess.Popen(cmd, stdout=None, stderr=subprocess.PIPE)
	out, err = torun.communicate()

def parse_tnp_synteny(out_dir):
	'''
	Parse the output file blasting one transposon to the other
	'''
	blast_out = [file for file in os.listdir(out_dir) if '_comp_tn_' in file][0]
	lines = []
	with open(out_dir + os.sep + blast_out, 'r') as f:
		lines = f.readlines()
	if len(lines) ==  1:
		coords = set()
		for line in lines:
			line = line.split('\t')
			iden = float(line[2])
			start1 = line[6]
			end1 = line[7]
			start2 = line[8]
			end2 = line[9]
			coords.add(start1)
			coords.add(start2)
			coords.add(end1)
			coords.add(end2)
			if iden == 100.000 and len(list(coords)) == 2:
				print('These transposons are identical')
				return True
			else:
				return False
	else:
		return False
		
def create_output_html(out_dir):
	output_folder = out_dir
	van_blast = [f for f in os.listdir(output_folder) if 'blastn_van_' in f]
	samples = [f.replace('blastn_van_', '').replace('.txt', '') for f in van_blast]
	tnp_blast = [f for f in os.listdir(output_folder) if 'blastn_comp' in f]
	mash = [f for f in os.listdir(output_folder) if 'mash_dist.txt' in f][0]


	with open(output_folder + os.sep + 'results.html', 'w') as out:
		out.write('<!DOCTYPE html>' + '\n' +'<html>' + '\n')
		out.write('<h2 style="font-family:arial">ClonalTracker output results for genomes ' + samples[0] + ' and ' + samples[1] + '</h2>' + '\n')
		out.write('<p style="font-family:arial"><b>1. <i>van</i> typing results using blastN:</b></p>'+ '\n')
		out.write('<p style="font-family:arial"> - ' + van_blast[0] + '\n')
		out.write('    <p><iframe src="' + output_folder + os.sep + van_blast[0] + '" frameborder="0" height="40" width="100%" /></iframe></p>' + '\n')
		out.write('<p style="font-family:arial"> - ' + van_blast[1] + '\n')
		out.write('    <p><iframe src="' + output_folder + os.sep + van_blast[1] + '" frameborder="0" height="40" width="100%" /></iframe></p>' + '\n')
		if 'TETyper_out' in os.listdir(output_folder):
			tetyper = [f for f in os.listdir(output_folder + os.sep + 'TETyper_out') if 'summary' in f]
			out.write('<p style="font-family:arial"><b>2. Transposon typing results</b></p>'+ '\n')
			#out.write('<p style="font-family:arial">2.1 TETyper results:</p>'+ '\n')
			#out.write('<p style="font-family:arial"> - ' + tetyper[0] + '\n')
			#out.write('    <p><iframe src="' + output_folder + os.sep + 'TETyper_out/' + tetyper[0] + '" frameborder="0" height="80" width="100%" /></iframe></p>' + '\n')
			#out.write('<p style="font-family:arial"> - ' + tetyper[1] + '\n')
			#out.write('    <p><iframe src="' + output_folder + os.sep + 'TETyper_out/' + tetyper[1] + '" frameborder="0" height="80" width="100%" /></iframe></p>' + '\n')
			out.write('<p style="font-family:arial">2.1 Clinker results after running RagTag (if needed to scaffold the transposon contigs) & ISEScan:</p>'+ '\n')
			out.write('    <p><iframe src="' + output_folder + os.sep + 'clinker_tn_trim_viz.html' + '" frameborder="0" width=1000 height=800" /></iframe></p>' + '\n')
			out.write('<p style="font-family:arial">2.2 blastN results to check the synteny of both transposons:</p>'+ '\n')
			out.write('    <p><iframe src="' + output_folder + os.sep + tnp_blast[0] + '" frameborder="0" height="40" width="100%" /></iframe></p>' + '\n')
		if mash:
			out.write('<p style="font-family:arial"><b>3. Whole genome comparison using Mash:</b></p>'+ '\n')
			out.write('    <p><iframe src="' + output_folder + os.sep + mash + '" frameborder="0" height="90" width="100%" /></iframe></p>' + '\n')

		out.write('</html>')


if __name__ == '__main__' :
	fa1 = sys.argv[1]
	fa2 = sys.argv[2]
	loc_reads = '/'.join(fa1.split('/')[0:-1])
	fasta1 = fa1.split('/')[-1]
	fasta2 = fa2.split('/')[-1]
	out1 = '.'.join(fasta1.split(".")[0:-1])
	out2 = '.'.join(fasta2.split(".")[0:-1])
	output_dir = "VRE_cwd_" + out1 + "-" + out2
	if not os.path.exists(output_dir):
		os.makedirs(output_dir)
	sys.stdout= open(output_dir + os.sep + 'run.logs', 'w')
	check1 = check_input_is_fasta(fa1)
	check2 = check_input_is_fasta(fa2)
	if check1 == True and check2 == True: ##if both files are fasta
		##copy files for poppunk analysis
		#copyfile(fa1, output_dir + os.sep + fasta1)
		#copyfile(fa2, output_dir + os.sep + fasta2)
		with open(output_dir + os.sep + 'list_new_genomes.txt', 'w') as list1:
			list1.write(out1 + '\t' + fa1 + '\n')
			list1.write(out2 + '\t' + fa2)
		#run blastn for both files
		out_b1 = output_dir + '/blastn_' + 'van_' + str(fasta1)
		script_dir = os.path.abspath( os.path.dirname( __file__ ))
		DB = script_dir + os.sep + 'van_representatives/'
		blast_db = DB + 'van_type_DB/van_nuc_seq_repre.fa'
		run_blastn(blast_db, out_b1, fa1)
		out_b2 = output_dir + '/blastn_' + 'van_' + str(fasta2)
		run_blastn(blast_db, out_b2, fa2)
		#check both blastn outputs
		records1= check_blastn_output(output_dir, fasta1)
		records2 = check_blastn_output(output_dir, fasta2)
		van_type_set = set([a.split('_')[-1] for a in records1.keys()])
		van_type2_set = set([a.split('_')[-1] for a in records2.keys()])
		van_type = ', '.join([a.split('_')[-1] for a in records1.keys()])
		van_type2 = ', '.join([a.split('_')[-1] for a in records2.keys()])
		#run poppunk analysis
		ppdb = script_dir + os.sep + 'vanAB_dataset_poppunk'
		destination = shutil.copytree(ppdb, output_dir + os.sep + 'vanAB_dataset')
		poppunk_folder = 'vanAB_dataset'
		run_poppunk(output_dir + os.sep + 'list_new_genomes.txt', poppunk_folder, output_dir)
		
		if van_type_set == van_type2_set: ##if both genomes have the same van type
			print('The two genomes are %s type' %(van_type))
			tnp_ref = [f for f in os.listdir(DB + 'tnp_db') if van_type in f and f.endswith('.fa')][0]
			tnp_db = DB + 'tnp_db' + os.sep + tnp_ref
			out_blast_tnp1 = output_dir + '/blast_contigs_' + str(fasta1)
			out_blast_tnp2 = output_dir + '/blast_contigs_' + str(fasta2)
			run_blastn(tnp_db, out_blast_tnp1, fa1)
			run_blastn(tnp_db, out_blast_tnp2, fa2)
			##run tetyper
			#run_tetyper(fa1, output_dir, records1, loc_reads)
			#run_tetyper(fa2, output_dir, records2, loc_reads)
			#result  = parse_tetyper_output(output_dir)
			loc_fastas = '/'.join(fa1.split('/')[:-1])
			print(loc_fastas)
			##Run isescan to check IS
			num_contigs = get_contigs_seqs(output_dir, loc_fastas)
			run_ragtag(output_dir, num_contigs, van_type)
			tnp_db = [f for f in os.listdir(DB + 'tnp_db') if van_type in f]
			tn_db = [f for f in os.listdir(DB + 'tnp_db') if van_type in f][0].split('.')[0] + '.fa'
			inf1 = output_dir + os.sep + fasta1.split('.')[0] + '_tnp_contig.fa'
			out_t1 = output_dir + os.sep + 'blastn_tnp_' + fasta1	
			run_blastn(DB + 'tnp_db/' + tn_db, out_t1, inf1)
			inf2 = output_dir + os.sep + fasta2.split('.')[0] + '_tnp_contig.fa'
			out_t2 = output_dir + os.sep + 'blastn_tnp_' + fasta2 
			run_blastn(DB + 'tnp_db/' + tn_db, out_t2, inf2)
			compare_tn_sequence_with_DB(output_dir)
			blast_tn_out = [f for f in os.listdir(output_dir) if 'blastn_tnp_' in f]
			for b in blast_tn_out:
				contig_fa = output_dir + os.sep + b.replace("blastn_tnp_", '').replace('.fna', '').replace('.fasta', '').replace('.fa', '') + '_tnp_contig.fa'
				out_fa = output_dir + os.sep + b.replace('blastn_tnp_', '').replace('.fna', '').replace('.fasta', '').replace('fa', '') + '_tnp_trimmed.fa'
				trim_tnp_region(output_dir + os.sep + b, contig_fa, out_fa)
			makeblastdb(output_dir)
			tnp_ref_db = output_dir + os.sep + [file.replace('_tnp_contig.fa', '') for file in os.listdir(output_dir) if '_tnp_contig.fa' in file][0] + '_tnp_trimmed.fa'
			tnp_query = output_dir + os.sep + [file.replace('_tnp_contig.fa', '') for file in os.listdir(output_dir) if '_tnp_contig.fa' in file and not '.fai' in file][1] +'_tnp_trimmed.fa'
			out_t3 = output_dir + os.sep + 'blastn_comp_tn_' + fasta2
			run_blastn(tnp_ref_db,out_t3, tnp_query)
			run_isescan(output_dir)
			identity = parse_isescan_output(output_dir)
			create_gbk_from_isescan_out(output_dir)
			ref_tn = [f for f in os.listdir(DB + '/tnp_db/') if van_type in f and f.endswith(".gb")]
			##copy tnp reference gbk file
			#copyfile(DB + '/tnp_db/' + ref_tn[0], output_dir + os.sep + 'tn_gbks/' + ref_tn[0])
			run_clinker(output_dir)
			avg, num = parse_alignment_scores(output_dir)
			assess_clonality_with_mash_sketch(output_dir, fa1, fa2)
			assess_clonality_with_mash_dist(output_dir)
			synteny = parse_tnp_synteny(output_dir)
			create_output_html(output_dir)
			if identity == True and avg >= 0.99 and num == 0 and synteny == True: # if they have the same SNPs, deletions and ISs
				#evaluate MASH output to assess clonality
				print("The transposons of these two genomes seem to be identical")
				evaluate_mash_distances(output_dir)
			else:
				print('The transposons of these two genomes are not identical, average gene identity with gene identity > 0.3 is  %s, not considering %s genes for lack of similarity' %(avg, num))

		else:
			create_output_html(output_dir)
			if not van_type == van_type2:
				print('These two genomes are not the same van type, %s is %s and %s is %s' %(fasta1, van_type, fasta2, van_type2))
			elif not van_type or not van_type2:
				print('At least one of these genomes do not have a vancomycin-resistant gene')
	else:
		print("Error: one of the input files seems not to be a nucleotide FASTA sequence. Please, input a correct file")
	sys.stdout.close()

