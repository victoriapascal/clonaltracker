import sys
import subprocess
import os
from Bio import SeqIO, SeqFeature
from Bio.Seq import Seq
from BCBio import GFF
from Bio.SeqRecord import SeqRecord
#from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Align.Applications import ClustalwCommandline
from shutil import copyfile
import shutil
import pickle

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

def run_blastn(db, out_dir, fasta, use):
	'''Run blastn on any given fatsa sequence using the van gene DB'''
	fa = fasta.split('/')[-1]
	print('1. Running blastn for %s...' %(fa))
	output_file = open(out_dir + os.sep + 'blastn_' + str(use) + "_" + str(fa), 'w')
	#db = '/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/van_type_DB/van_nuc_seq_repre.fa'
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
	print(' '.join(cmd2))
	if err:
		print(err)
	if err2:
		print(err2)
	if err3:
		print(err3)
	os.chdir(cwd)

def run_tetyper(fasta, out_dir, type_pident):
	'''Running TETyper to assess transposon type'''
	print('3. Running TETyper to assess transposon type')
	#location_reads = "/".join(fasta.split('/')[:-1]) + '/'
	#location_reads = '/hpc/dla_mm/vpascalandreu/data/vanA_raw_reads/vanA_reads/'
	location_reads = '/hpc/dla_mm/vpascalandreu/data/vanB_raw_reads_renamed2/'
	#location_reads = '/hpc/dla_mm/vpascalandreu/data/new_outbreak_assemblies_vrefidia/'
	#location_reads = '/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/vanB_dataset_wget/'
	sample = fasta.split('/')[-1].split('.')[0] + "_"
	fq1 = [file for file in os.listdir(location_reads) if 'R1.fastq.gz' in file and sample in file][0]
	fq2 = [file for file in os.listdir(location_reads) if 'R2.fastq.gz' in file and sample in file][0]
	#fq1 = [file for file in os.listdir(location_reads) if '_R1_val_1.fq' in file and sample in file][0]
	#fq2 = [file for file in os.listdir(location_reads) if '_R2_val_2.fq' in file and sample in file][0]
	tnp_ref = ''
	for a in type_pident.keys():
		if "vanA" in a:
			tnp_ref = '/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/tnp_db/vanA_M97297.fa'
		else:
			if "vanB" in a:
				tnp_ref = '/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/tnp_db/vanB_AY655721.2.fa'
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
				snps[name].append(lines[3])

	result = ''
	if set(snps[names[0]]) == set(snps[names[1]]):
		result = True
	else:
		result = False
	
	with open('/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/tnp_db/tnp_db.pickle', 'rb') as handle:
		tnp_db = pickle.load(handle)
	with open('/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/tnp_db/tn_nomenclature.pickle', 'rb') as handle2:
		tn_nomen = pickle.load(handle2)
	
	snp1 = '-'.join(snps[names[0]])
	snp2 = '-'.join(snps[names[1]])
	if snp1 in tnp_db.keys():
		print(names[0] + " has the transposon type " + tn_nomen[snp1] + ", the same transposon as " + ','.join(tnp_db[snp1]))

	if snp2 in tnp_db.keys():
                print(names[1] + " has the transposon type " + tn_nomen[snp2] + ", the same transposon as " + ','.join(tnp_db[snp2]))


	return result

def get_contigs_seqs(output_dir):
	'''
	From the TETyper blast output, get the transposons that align with the transposon
	'''
	files = output_dir + os.sep + 'TETyper_out/'
	bf = [f for f in os.listdir(files) if 'blast' in f]
	sample_contigs = {}
	for b in bf:
		name = b.replace('_blast.txt', '')
		sample_contigs[name] = []
		with open(files + os.sep + b, 'r') as f1:
			for line in f1:
				line = line.strip().split('\t')
				if not line[0] in sample_contigs[name]:
					sample_contigs[name].append(line[0])
	
	##remove hard-coded path when the usual one is known
	fna = '/hpc/dla_mm/vpascalandreu/data/vanB_fastas'	
	#fna = '/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/bk_vanB_dataset_wget/202106102320_results/scaffolds/'	
	#fna = '/hpc/dla_mm/vpascalandreu/data/new_outbreak_assemblies_vrefidia/'
	for a in sample_contigs.keys():
		with open(fna + os.sep + a + '.fasta', 'r') as f2:
			lines = f2.readlines()
			out = open(output_dir + os.sep + a + "_tnp_contig.fa", 'w')
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
		if len(contigs) > 1: ## when there the tnp is split into multiple contigs, run ragtag
			out = output_dir + os.sep + sample + "_ragtag"
			ref_tn = '/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/tnp_db/' + [f for f in os.listdir('/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/tnp_db/') if van_type in f and f.endswith(".fa")][0]
			query = output_dir + os.sep + [f for f in os.listdir(output_dir) if sample in f and '_tnp_contig.fa' in f][0]
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
						f2.write(line)
						  
					
				
def run_isescan(out_dir):
	'''
	Run isescan to assess potential tnp IS
	'''
	print("4. Running ISEScan to predict ISs")
	is_out = output_dir + os.sep + "ISEScan_out"
	if not os.path.isdir(is_out):
		os.mkdir(is_out)
	contigf = [file for file in os.listdir(out_dir) if "_tnp_contig.fa" in file]
	for cf in contigf:
		od = cf.replace("_tnp_contig.fa", '')
		cmd = ['isescan.py', '--seqfile', out_dir + os.sep+ cf, '--output', is_out + os.sep + od]
		torun = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		out, err = torun.communicate()

def parse_isescan_output(out_dir):
	'''
	Check ISEScan output folder to see if the two genomes have the same (if any) IS
	'''
	samples_id = {}
	samples = out_dir.split('_')[2].split('-')
	for sample in samples:
		folders = os.listdir(out_dir + os.sep + "ISEScan_out/" + sample)
		if out_dir in folders:
			root = out_dir + os.sep + "ISEScan_out/" + sample 
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
	output_file = open(out_dir + os.sep + 'mash_dist.tsv', 'w')
	cmd = ['mash', 'dist', out_dir + os.sep+ 'mash_sketch.msh', out_dir + os.sep+ 'mash_sketch.msh']
	torun = subprocess.Popen(cmd, stdout=output_file, stderr=subprocess.PIPE)
	out, err = torun.communicate()
	output_file.close()

def evaluate_mash_distances(out_dir):
	sim = set()
	prop = set()
	with open(out_dir + os.sep +'mash_dist.tsv', 'r') as out:
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
	
	samples = os.listdir(out_dir + os.sep + 'ISEScan_out')
	if not os.path.isdir(out_dir + os.sep + "tn_gbks"):
		os.mkdir(out_dir + os.sep + "tn_gbks")	
	for sample in samples:
		proteome_dir = out_dir + os.sep + 'ISEScan_out' + os.sep + sample + os.sep + 'proteome' + os.sep + out_dir
		proteome_files = os.listdir(proteome_dir)
		genes = [f for f in proteome_files if '_tnp_contig.fa.ffn' in f][0]
		prots = [f for f in proteome_files if '_tnp_contig.fa.faa' in f][0]
		is_dir = out_dir + os.sep + 'ISEScan_out' + os.sep + sample + os.sep + out_dir
		if os.path.isdir(is_dir):
			is_files = os.listdir(is_dir)
			is_nuc = [f for f in is_files if '_tnp_contig.fa.is.fna' in f][0]
			is_prot = [f for f in is_files if '_tnp_contig.fa.orf.faa' in f][0] 	
		
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
		
		contigs = out_dir + os.sep + sample + '_tnp_contig.fa'
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

		with open(out_dir + os.sep + 'tn_gbks/' + sample + '.gb', 'w') as input_handle:
			SeqIO.write(record, input_handle, "genbank")
def run_clinker(out_dir):
	'''
	From the the tn gbk files, run clincker to visualize the similarities/differences between the isolates transposons
	'''

	output_file1 =  out_dir + os.sep + 'clinker_tn_viz.html'
	output_file2 = out_dir + os.sep + 'clinker_tn_msa.aln'
	gbks = os.listdir(out_dir + os.sep + 'tn_gbks')
	inp = out_dir + '/tn_gbks/*'
	cmd = 'clinker '  + inp + ' -p '+ output_file1 + ' -o ' + output_file2
	torun = subprocess.Popen([cmd], stderr=subprocess.PIPE, shell=True)
	err = torun.communicate()

if __name__ == '__main__' :
	fa1 = sys.argv[1]
	fa2 = sys.argv[2]
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
		#copy files for poppunk analysis
#		copyfile(fa1, output_dir + os.sep + fasta1)
#		copyfile(fa2, output_dir + os.sep + fasta2)
#		with open(output_dir + os.sep + 'list_new_genomes.txt', 'w') as list1:
#			list1.write(out1 + '\t' + fa1 + '\n')
#			list1.write(out2 + '\t' + fa2)
		#run blastn for both files
		blast_db = '/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/van_type_DB/van_nuc_seq_repre.fa'
		run_blastn(blast_db, output_dir, fa1, 'van')
		run_blastn(blast_db, output_dir, fa2, 'van')
		#check both blastn outputs
		records1= check_blastn_output(output_dir, fasta1)
		records2 = check_blastn_output(output_dir, fasta2)
		van_type_set = set([a.split('_')[-1] for a in records1.keys()])
		van_type2_set = set([a.split('_')[-1] for a in records2.keys()])
		van_type = ', '.join([a.split('_')[-1] for a in records1.keys()])
		van_type2 = ', '.join([a.split('_')[-1] for a in records2.keys()])
		#run poppunk analysis
#		ppdb = '/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/vanAB_dataset_poppunk'
		#destination = shutil.copytree(ppdb, output_dir + os.sep + 'vanAB_dataset')
		#poppunk_folder = 'vanAB_dataset'
		#run_poppunk(output_dir + os.sep + 'list_new_genomes.txt', poppunk_folder, output_dir)
		
		if van_type_set == van_type2_set: ##if both genomes have the same van type
			print('The two genomes are %s type' %(van_type))
			
			##run tetyper
			#run_tetyper(fa1, output_dir, records1)
			#run_tetyper(fa2, output_dir, records2)
			result  = parse_tetyper_output(output_dir)
			##Run isescan to check IS
			num_contigs = get_contigs_seqs(output_dir)
			run_ragtag(output_dir, num_contigs, van_type)	
			run_isescan(output_dir)
			identity = parse_isescan_output(output_dir)
			create_gbk_from_isescan_out(output_dir)
			ref_tn = [f for f in os.listdir('/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/tnp_db/') if van_type in f and f.endswith(".gb")]
			##copy tnp reference gbk file
			copyfile('/hpc/dla_mm/vpascalandreu/VRE_pipeline_validation/pipeline/van_representatives/tnp_db/' + ref_tn[0], output_dir + os.sep + 'tn_gbks/' + ref_tn[0])
			run_clinker(output_dir)
			assess_clonality_with_mash_sketch(output_dir, fa1, fa2)
			assess_clonality_with_mash_dist(output_dir)

			if result == True and identity == True: # if they have the same SNPs, deletions and ISs
				#evaluate MASH output to assess clonality
				evaluate_mash_distances(output_dir)
			else:
				print('The transposons of these two genomes are not identical')
			
		else:
			if not van_type == van_type2:
				print('These two genomes are not the same van type, %s is %s and %s is %s' %(fasta1, van_type, fasta2, van_type2))
			elif not van_type or not van_type2:
				print('At least one of these genomes do not have a vancomycin-resistant gene')
	else:
		print("Error: one of the input files seems not to be a nucleotide FASTA sequence. Please, input a correct file")
	sys.stdout.close()

