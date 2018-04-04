r1='TTGAGAGCCGAGGTTGTATTTTATCGCACTTGCCGCCATACTCACCTGACCTGAACCCCATAGAGAAAGGTTTCAGTGCTTACAAACAGGCACTGCAACTTAACGAGGACCTACTCACAGGAGGTGAAGGTGACTTCTTCATCATTGACGA
   ...: GTTTGTCCGATTGGTTTTCACCTCTGAGCTTACCCAAAAAACTTTTCCGAGGGTCGGGCTATGCTGTGTACTAGGCTTTGGAGGGATTATGTTGCCTCTCACATGTTTCCACCTCACCAACCCAATTTTTTCTACCGCTTGATTCAATGTCATTT
   ...: TTTTTTTCATCTCACCAACTCATGTTTTTCTACCACTCAATCCTATGTGATCTATTTGGTTTATATCAGAAATCTATCATAAATTTTAACCCTTCACATTTTAAGGTGTGTTATCAGATATAACAGGCGGCTCAATTGAGTGAGGTTTTTATAAA
   ...: ACAATCAGACAGATGACAATAGCCTCGAATGAGACAGGTTTTTTATTAACAGTCAGACATGGTTCTTATATAAAAATGAGATAGATAACAGGAGGCTCGATTGAGAGGGATTTTATCAAAAGTCAGACGTATCACAGGCGGCTCCATTGAGCGAG
   ...: GTTTCTTAATAAAAATCAGACAAACAACAGGCAGCTCGATTGAGCGAGGTGTTTTTATAACAATAATTTGAGGGGATGGAAACCAAGAGTCACGAAACAATTGAAAGAATTTTGGTTCAAACCAAAAAACTGGATGCAAACCTATGTATCAATCA
   ...: TGCTGCTCTGCTCACTTGCAAAGGGATTTCAAGTGGCGCTCGTAAGCTTTGACACTTTTGGTTAGAGCACACGCTGTGGTAAGCCAATCAGTGTCACCAGTAACCCAAAGTAGTGTGCAATCACACTAGGATCACAGGCCCCAGGACTCAGAACC
   ...: CAGCAATCACCAAGAATGAAGGTATCACCAACTTTGAAACAGTTTCACAACTGGTATCACAGGTCCCCCCTTCCACATCATCACCACACCCCATTTCCCTGCCTCCTTATTTCTGATACAGTCCATTTTTATGTTTCACTCACCTTTCACACTTG
   ...: TTCCATGTTCCTATCCCATTCCATCTCATTCCACCTACCCTGCATTGACCCCCTTGGAGAATTTTCACCCACATCAACCCCTTCCACATCATGAACCATGTTTCATGCCCTTGTCAACCCACAAGCTGGCCACATTGTTTCCATTCCTCAGTTGA
   ...: AGTTTTTCACTTACCCCTGATAGAATGTTCCATTCCCTGGCACAACCACAAGTTGGCCCCTCTCTTGTTTTACACCTCTTTACTATATAA'
   ...:
   ...:
   ...: r2='GGTTTGCATCCAGTTTTTTGGTTTGAACCAAAATTCTTTCAATTGTTTCGTGACTCTTGGTTTCCATCCCCTCAAATTATTGTTATAAAAACACCTCGCTCAATCGAGCTGCCTGTTGTTTGTCTGATTTTTATTAAGAAACCTCGCTCAA
   ...: TGGAGCCGCCTGTGATACGTCTGACTTTTGATAAAATCCCTCTCAATCGATTGAGAGGGATTTTATCAAAAGTCAGACGTATCACAGGCGGCTCCATTGAGCGAGGTTTCTTAATAAAAATCAGACAAACAACAGGCAGCTCGATTGAGCGAGGT
   ...: GTTTTTATAACAATAATTTGAGGGGATGGAAACCAAGAGTCACGAAACAATTGAAAGAATTTTGGTTCAAACCAAAAAACTGGATGCAAACCTATGTATCAATCATGCTGCTCTGCTCACTTGCAAAGGGATTTCAAGTGGCGCTCGTAAGCTTT
   ...: GACACTTTTGGTTAGAGCACACGCTGTGGTAAGCCAATCAGTGTCACCAGTAACCCAAAGTAGTGTGCAATCACACTAGGATCACAGGCCCCAGGACTCAGAACCCAGCAATCACCAAGAATGAAGGTATCACCAACTTTGAAACAGTTTCACAA
   ...: CTGGTATCACAGGTCCCCCCTTCCACATCATCACCACACCCCATTTCCCTGCCTCCTTATTTCTGATACAGTCCATTTTTATGTTTCACTCACCTTTCACACTTGTTCCATGTTCCTATCCCATTCCATCTCATTCCACCTACCCTGCATCGACC
   ...: CCCTTGGAGAATTTTCACCCACATCAACCCCTTCCACATCATGAACCATGTTTCATGCCCTTGTCAACCCACAAGCTGGCCACATTGTTTCCATTCCTCAGTTGAAGTTTTTCACTTACCCCTGATAGAATGTTCCATTCCCTGGCACAACCACA
   ...: AGTTGGCCCCTCTCTTGTTTTACACCTCTTTACTATATAA'
 chars = ['H', 'I', 'G', 'F','1', '<', '0']

d = ['A', 'C', 'T', 'G']
rep = ''
for i in range(1500):
     rep += random.choice(d)

s1_start = ''
s1_end = ''

s2_start = ''
s2_end = ''
for i in range(2500):
   s1_start += random.choice(d)
   if i < 2478:
     s2_start  += random.choice(d)
     if i < 1437:
       s1_end +=  random.choice(d)
       if i < 1300:
          s2_end+= random.choice(d)
   r1 = s1_start + rep + s1_end
   r2 = s2_start + rep + s2_end



def make_seqs(fname, rep_len=1500, start_len=2000, end_len=1300, start_var=53, end_var=12):
 d = ['A', 'C', 'T', 'G']
 rep = ''
 for i in range(rep_len):
     rep += random.choice(d)
 s1_start = ''
 s1_end = ''
 s2_start = '' 
 s2_end = ''
 for i in range(start_len):
   s1_start += random.choice(d)
   if i < 2478:
     s2_start  += random.choice(d)
     if i < 1437:
       s1_end +=  random.choice(d)
       if i < 1300:
          s2_end+= random.choice(d)
 r1 = s1_start + rep + s1_end
 r2 = s2_start + rep + s2_end
 write_gfa_fasta(fname, rep, s1_start, s2_start, s1_end, s2_end)
 return r1, r2

def write_gfa_fasta(fname, rep, s1_start, s2_start, s1_end, s2_end):
    with open(fname.split('.')[0]+'.gfa', 'w') as f:
    	 f.write('S\t rep\t%s\n' % rep)
    	 f.write('S\t s1_start\t%s\n' % s1_start)
    	 f.write('S\t s2_start\t%s\n' % s2_start)
    	 f.write('S\t s1_end\t%s\n' % s1_end)
    	 f.write('S\t s2_end\t%s\n' % s2_end)
	 f.write('L\t rep \t + \t s1_start \t + \t 199M\n')
	 f.write('L\t rep \t + \t s2_start \t + \t 199M\n')
	 f.write('L\t rep \t - \t s1_end \t + \t 199M\n')
	 f.write('L\t rep \t - \t s2_end \t + \t 199M\n')
    with open(fname.split('.')[0]+'.fasta', 'w') as f:
    	 f.write('>rep\n%s\n' % rep)
    	 f.write('>s1_start\n%s\n' % s1_start)
    	 f.write('>s2_start\n%s\n' % s2_start)
    	 f.write('>s1_end\n%s\n' % s1_end)
    	 f.write('>s2_end\n%s\n' % s2_end)



def write_gfa_fasta_multi_reps(fname, rep, starts, ends):
    with open(fname.split('.')[0]+'.gfa', 'w') as f:
        f.write('S\t rep\t%s\n' % rep)
        for i,  s in enumerate(starts):
            f.write('S\t s%d_start\t%s\n' % (i, s))
            f.write('L\t rep \t + \t s%d_start \t + \t 199M\n' % i)
        for i,  s in enumerate(ends):
            f.write('S\t s%d_end\t%s\n' % (i, s))	 
            f.write('L\t rep \t - \t s%d_end \t + \t 199M\n'% i)
    with open(fname.split('.')[0]+'.fasta', 'w') as f:
    	 f.write('>rep\n%s\n' % rep)
         for i, start in enumerate(starts):
    	   f.write('>s%d_start\n%s\n' % (i, start))
         for i, end in enumerate(ends):
            f.write('>s%d_end\n%s\n' % (i,end))
             	 


def write_fasta(fname, r1, r2):
  with open(fname, 'w') as f:
     for i in range(1000):
     	  r = '@r1_' +fname+'_'+ str(i)
     	  r22 = '@r2_' + str(i)
          b = random.choice([i for i in range(len(r1) - 251)])
          b2 = random.choice([i for i in range(len(r2) - 251)])
          read = r1[b:b+250]
	  read2 = r2[b2:b2+250] 
          qual = ''
   	  for j in range(250):
                 qual += random.choice(chars)
          f.write(r +'\n')
          f.write(read +'\n')
          f.write('+' +'\n')
	  f.write(qual +'\n')
          f.write(r22 +'\n')
          f.write(read2 +'\n')
          f.write('+' +'\n')
	  f.write(qual +'\n')


def write_fasta_reps(fname, read):

  
def write_fasta(fname, r1, r2):
  with open(fname, 'w') as f:
     for i in range(1000):
     	  r = '@r1_' +fname+'_'+ str(i)
     	  r22 = '@r2_' + str(i)
          b = random.choice([i for i in range(len(r1) - 251)])
          b2 = random.choice([i for i in range(len(r2) - 251)])
          read = r1[b:b+250]
          read2 = r2[b2:b2+250] 
          qual = ''
          for j in range(250):
            qual += random.choice(chars)
          f.write(r +'\n')
          f.write(read +'\n')
          f.write('+' +'\n')
          f.write(qual +'\n')
          f.write(r22 +'\n')
          f.write(read2 +'\n')
          f.write('+' +'\n')
          f.write(qual +'\n')):
  

def make_seqs_write_fasta(fname, rep_len=1500, start_len=2000, end_len=1300, start_var=53, end_var=1):
    r1, r2 = make_seqs(fname , rep_len, start_len, end_len, start_var, end_var)
    write_fasta(fname, r1, r2)

make_seqs_write_fasta('test_rep_r1.fastq')

v = no varieies
r = repeats
mask for each variety included
scale al by total bymber of reads:repeated section r*const times

class SampleRef(object):
     def __init__(self, repeat_structure)
    
class Variety(object):
    def __init__(self, repeat_mask,additional_repeats):
        self.repeat_mask = repeat_mask # mask of which repets from each section in the reference are present in this variety, and hence included in reads
        self.additional_repeats =additional_repeats #

-g /tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/bsg_testing/yr/contig_graphs/yr_good_contigs_raw.good.gfa --cidxread1 /tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/strawberry/kmer_compression_index/yr_r1_tail_8m2.fastq --cidxread2 /tgac/workarea/Research-Groups/RG-Bernardo-Clavijo/strawberry/kmer_compression_index/yr_r2_tail_8m2.fastq  --max_mem 140  -o yr
