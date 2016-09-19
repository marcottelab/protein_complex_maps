#!/usr/bin/env python
import os
import sys

if( len(sys.argv) < 3 ):
    sys.stderr.write("usage: python msblender_out-to-pep_count-mFDRpsm.py <msblender_out> <FDR_cutoff(eg 0.01)> <(optional: eFDRpsm, mFDRpsm)> \n") 
    sys.exit(1)

filename_mb_out = sys.argv[1]
FDR_cutoff = float(sys.argv[2])
FDR_string = sys.argv[2].replace('.','')

error_model = 'mFDRpsm'
error_model_list = ['eFDRpsm','mFDRpsm']
if( len(sys.argv) == 4 and sys.argv[3] in error_model_list ):
    error_model = sys.argv[3]

sys.stderr.write("FDR cutoff: %.3f\nError model:%s\n"%(FDR_cutoff,error_model))

filename_base = filename_mb_out.replace('.msblender_in','').replace('.msblender_out','')

psm_mvScore = dict()
psm_TD = dict()

f_mb_out = open(filename_mb_out,'r')
f_mb_out.readline()
for line in f_mb_out:
    tokens = line.strip().split("\t")
    tmp_mvScore = float(tokens[-1])
    tmp_psm = tokens[0]
    psm_mvScore[tmp_psm] = float(tokens[-1])
    psm_TD[tmp_psm] = tokens[1]
f_mb_out.close()

count_T = 0
count_D = 0

D_pep_count = dict()
pep_count = dict()
sample_list = []
sum_error = 0.0
f_log = open('%s.pep_count_%s%s.log'%(filename_base,error_model,FDR_string),'w')
f_log.write('PSM_id\tFDR\tmvScore\n')
for tmp_psm in sorted(psm_mvScore.keys(),key=psm_mvScore.get,reverse=True):
    if( psm_TD[tmp_psm] == 'D' ):
        count_D += 1
    elif( psm_TD[tmp_psm] == 'F' ):
        count_T += 1
    else:
        sys.stderr.write('No T/D info: %s\n. Exit.'%tmp_psm)
        sys.exit(1)
    
    tmp_FDR = 0.0
    if( error_model == 'eFDRpsm' ):
        if( psm_mvScore[tmp_psm] < 1.0 ):
            tmp_FDR = float(count_D)/(count_T+count_D)
    elif( error_model == 'mFDRpsm' ):
        tmp_score = psm_mvScore[tmp_psm]
        sum_error += (1.0 - psm_mvScore[tmp_psm])
        tmp_FDR = sum_error/count_T
    else:
        sys.stderr.write('No error model is defined. Exit.')
        sys.exit(1)
    
    if( tmp_FDR < FDR_cutoff ):
        tmp_pep = tmp_psm.split('.')[-1]
        if( psm_TD[tmp_psm] == 'F' ):
            f_log.write('%s\t%.3f\t%.2f\n'%(tmp_psm,tmp_FDR,psm_mvScore[tmp_psm]))
            
            sample_name = tmp_psm.split('.')[0]
            sample_list.append(sample_name)
            if( not pep_count.has_key(tmp_pep) ):
                pep_count[tmp_pep] = dict()
            if( not pep_count[tmp_pep].has_key(sample_name) ):
                pep_count[tmp_pep][sample_name] = 0
            pep_count[tmp_pep][sample_name] += 1
        else:
            if( not D_pep_count.has_key(tmp_pep) ):
                D_pep_count[tmp_pep] = 0
            D_pep_count[tmp_pep] += 1
    else:
        break
f_log.close()

sample_list = sorted(list(set(sample_list)))

print "pep_count: %s" % (pep_count)
print "D_pep_count: %s" % (D_pep_count)
sys.stderr.write('Peptide FDR: %.3f\n'%( float(len(D_pep_count))/(len(pep_count)+len(D_pep_count)) ))
f_count = open('%s.pep_count_%s%s'%(filename_base, error_model,FDR_string),'w')
f_count.write('#Peptide FDR: %.3f\n'%( float(len(D_pep_count))/(len(pep_count)+len(D_pep_count)) ))
f_count.write('#PepSeq\tTotalCount\t%s\n'%('\t'.join(sample_list)))
for tmp_pep in sorted(pep_count.keys(),key=pep_count.get):
    out_list = []
    total_count = 0
    for tmp_sample in sample_list:
        if( pep_count[tmp_pep].has_key(tmp_sample) ):
            out_list.append( "%d"%pep_count[tmp_pep][tmp_sample] )
            total_count += pep_count[tmp_pep][tmp_sample]
        else:
            out_list.append('0')
    f_count.write('%s\t%d\t%s\n'%(tmp_pep,total_count,'\t'.join(out_list)))
f_count.close()
