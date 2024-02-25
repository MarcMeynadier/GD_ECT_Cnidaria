#!/usr/bin/env python

import pysam
import subprocess
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from collections import defaultdict

###########################################

def main():
   hmm_list = ['GSHPx']
   all_dbs = read_dbs('db_list.txt')
   for hmm_name in hmm_list:
      hmm_file = hmm_extract(hmm_name,'/Users/copley/data/pfam/Pfam-A.hmm')
      target_seqs = []
      observed = dict()
      for db in all_dbs:
         print(db[0],db[1],db[2])
         db_name = db[0]
         db_fn = db[1]
         fasta = pysam.Fastafile(db_fn)
         hit_list = do_hmmsearch(hmm_name,hmm_file,db_name,db_fn)
         for hit_id in hit_list:
            if hit_id in observed:
               continue
            seq = fasta.fetch(hit_id)
            seq_obj = Seq(seq)
            seq_rec = SeqRecord(seq_obj,id=hit_id)
            target_seqs.append(seq_rec)
            observed[hit_id] = 1
      target_db_fn = mk_target_db(hmm_name,target_seqs) # targets
      stk_fn = do_hmmalign(hmm_name,hmm_file,target_db_fn) # stk
      alimask_fn = do_esl_alimask(hmm_name,stk_fn) # alimask
      alimanip_fn = do_esl_alimanip(hmm_name,alimask_fn) # alimanip
      mfa_fn = convert_to_mfa(hmm_name,alimanip_fn)
      filtered_fn = filter_mfa(hmm_name,mfa_fn)
#      mk_tree(hmm_name,filtered_fn)

###########################################

def mk_tree(hmm_name,filtered_fn):

   subprocess.call(['iqtree','-s',filtered_fn,
                    '-nt','AUTO','-mset','LG'])

   return

###########################################

def filter_mfa(hmm_name,mfa_fn):

   filtered_out_fn = hmm_name + '_auto_filtered.mfa'

   seq_hash = defaultdict(set)
   reps = []
   for seq_rec in SeqIO.parse(mfa_fn,'fasta'):
      sp_tag = seq_rec.id.split('|')[0]
      seq_str = str(seq_rec.seq)
      if seq_str in seq_hash:
         if sp_tag not in seq_hash[seq_str]:
            seq_hash[seq_str].add(sp_tag)
            reps.append(seq_rec)
      else:
         seq_hash[seq_str].add(sp_tag)
         reps.append(seq_rec)

   with open(filtered_out_fn,'w') as fh:
      for seq_rec in reps:
         print(seq_rec.format('fasta'),end='',file=fh)

   return filtered_out_fn

###########################################

def convert_to_mfa(hmm_name,alimanip_fn):

   mfa_out_fn = hmm_name + '_auto.mfa'
   alignment = AlignIO.read(alimanip_fn,'stockholm')
   AlignIO.write(alignment,mfa_out_fn,'fasta')
   return mfa_out_fn

###########################################

def do_esl_alimask(hmm_name,stk_fn):
   alimask_out_fn = hmm_name + '.alimask'
   subprocess.call(['esl-alimask','-o',alimask_out_fn,
                    '--pavg','0.5','-p',stk_fn])
   return alimask_out_fn

###########################################

def do_esl_alimanip(hmm_name,alimask_fn):
   alimanip_out_fn = hmm_name + '.alimanip'

   subprocess.call(['esl-alimanip','-o',alimanip_out_fn,
                    '--minpp','0.3',alimask_fn])

   # need to do this sequentially - first remove bad residues to gap
   # then filter out sequences that are now too short
   # seem to be able to write to same name as input

   subprocess.call(['esl-alimanip','-o',alimanip_out_fn,
                    '--lnfract','0.7',alimanip_out_fn])

   return alimanip_out_fn

###########################################

def do_hmmalign(hmm_name,hmm_file,target_db_fn):
   stk_out_fn = hmm_name + '.stk'
   subprocess.call(['hmmalign','-o',stk_out_fn,hmm_file,target_db_fn])
   return stk_out_fn

###########################################

def mk_target_db(hmm_name,target_seqs):
   target_db_fn = hmm_name + '_targets.fasta'
   with open(target_db_fn,'w') as fh:
      for seq_rec in target_seqs:
         print(seq_rec.format('fasta'),end='',file=fh)
   return target_db_fn

###########################################

def do_hmmsearch(hmm_name,hmm_fn,db_name,db_fn):
   print(hmm_name,hmm_fn,db_name,db_fn)
   out_fn = hmm_name + '_' + db_name + '.hmmsearch'
   dt_out_fn = hmm_name + '_' + db_name + '_domtbl.out'
   print(out_fn,dt_out_fn)
   subprocess.call(['hmmsearch','--cut_ga','-o',out_fn,
                    '--domtblout',dt_out_fn,hmm_fn,db_fn])
   with open(dt_out_fn) as fh:
      db_ids = set()
      for line in fh:
         if line.startswith('#'): continue
         id = line.split()[0]
         db_ids.add(id)
   return db_ids

###########################################

def hmm_extract(hmm,hmmlib):
   hmm_fn = hmm + '.hmm'
   subprocess.call(['hmmfetch','-o',hmm_fn,hmmlib,hmm])
   return hmm_fn

###########################################

def read_dbs(db_list_file):

   dbs = []
   with open(db_list_file) as fh:
      for line in fh:
         if line.startswith('#'): continue
         line = line.rstrip()
         parts = line.split()
         name = parts[0]
         fasta = parts[1]
         index = parts[2]
         db = (name,fasta,index)
         dbs.append(db)
   return dbs

###########################################

if __name__ == "__main__":

   main()
