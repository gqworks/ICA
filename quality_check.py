import glob
import re
import pandas as pd
import gzip
from cyvcf2 import VCF

def get_dosage_int(geno):
    x = re.split('\\||/',geno)
    return sum([int(i) for i in x])

nucleotide_map = {
    'A':'T',
    'T':'A',
    'C':'G',
    'G':'C'
}
def flip_dosage(num):
    if num == 0:
        return 2
    elif num == 2:
        return 0
    else:
        return num

def fix_effect_allele(var, effect_allele, dosage_int):
    if effect_allele == var.ALT[0]: #Matched, dont change anything
        return dosage_int
    if effect_allele == var.REF: #Only flip if not palindromic with reference
        if nucleotide_map[effect_allele] == var.ALT[0]:
            #Is palindromic, cannot use this SNP to calculate PRS
            return 0
        else:
            return flip_dosage(dosage_int)

def calc_prs(vcf_path,PRS_table):
    PRS_total = []
    vcf = VCF(vcf_path)
    for i,row in PRS_table.iterrows():
        iterator = vcf(f'{row["hg38_chr"]}:{row["hg38_pos"]}-{row["hg38_pos"]}')
        if all(False for _ in iterator):
            if row['effect_is_ref']:
                PRS_total.append(row['effect_weight']*2)
            continue
        iterator = vcf(f'{row["hg38_chr"]}:{row["hg38_pos"]}-{row["hg38_pos"]}')
        var = next(iterator)
        if var.POS != row["hg38_pos"]: #Sometimes vcf gets neighbour positions
            continue
        hg38_pos_id = f'{var.CHROM}_{var.POS}'
        dosage_int = sum(var.genotypes[0][0:2])
        effect_allele = PRS_table[PRS_table.hg38_pos_id==hg38_pos_id]['effect_allele'].item()
        effect_weight = PRS_table[PRS_table.hg38_pos_id==hg38_pos_id]['effect_weight'].item()
        fixed_dosage = fix_effect_allele(var,effect_allele,dosage_int)
        PRS = effect_weight*fixed_dosage
        PRS_total.append(PRS)
    return(sum(PRS_total))

def calc_mlc(vcf_path,MLC_table):
    vcf = VCF(vcf_path)
    MLC_scores = []
    for v in vcf('chrM'):
        #Skip indels
        if len(v.ALT[0])>1 or len(v.REF)>1:
            continue
        #Check ref
        if v.REF != MLC_table.loc[MLC_table.Position==v.POS,'Reference'].to_list()[0]:
            print(v.POS)
            print("Incorrect reference allele")
            continue
        score = MLC_table[(MLC_table.Position==v.POS)&(MLC_table.Alternate==v.ALT[0])]['MLC_score'].item()
        MLC_scores.append(score)
    return(MLC_scores)
    
# Ploidy
ploidy_files = glob.glob("project/results/germline/*/*/*.ploidy_estimation_metrics.csv")
dat_list = []
for file in ploidy_files:
    sample_name = file.split('/')[4]
    dat = pd.read_csv(file, names = ['1','2','index',sample_name], usecols = ['index',sample_name])
    dat = dat.set_index('index')
    dat_list.append(dat)
ploidy_merged_dat = pd.concat(dat_list, axis=1).T

# Mitochondrial Copy Number 
contig_cov_files = glob.glob("project/results/germline/*/*/*.wgs_contig_mean_cov.csv")
copy_number_dat = []
for file in contig_cov_files:
    sample_name = file.split('/')[4]
    dat = pd.read_csv(file, names = ['contig','read_count','coverage'])
    autosomal_coverage = dat[dat.contig.str.contains("Autosomal regions")]['coverage'].item()
    mitochondrial_coverage = dat[dat.contig=="chrM"]['coverage'].item()
    mito_copy_number = 2*mitochondrial_coverage/autosomal_coverage
    copy_number_dat.append([sample_name,autosomal_coverage,mitochondrial_coverage,mito_copy_number])
copy_number_dat = pd.DataFrame(copy_number_dat, columns = ['index','autosomal_cov','mitochondrial_cov','mitochondrial_copy_number'])
copy_number_dat = copy_number_dat.set_index('index')

# Polygenic Risk Score (PRS)
PRS_table = pd.read_csv('project/aux/PGS000902_hg38.txt', sep='\t', comment = '#')
vcf_files = glob.glob("project/results/germline/*/*/*.hard-filtered.vcf.gz")
prs_dat_list = []
for vcf in vcf_files:
    sample_name = vcf.split('/')[4]
    prs_dat_list.append([sample_name,calc_prs(vcf,PRS_table)])
prs_dat = pd.DataFrame(prs_dat_list, columns = ['index','prs'])
prs_dat = prs_dat.set_index('index')

# Mitochondrial local constraint (MLC)
vcf_files = glob.glob("project/results/germline/*/*/*.hard-filtered.vcf.gz")
MLC_table = pd.read_csv("project/aux/MLC_supplementary_dataset_7.tsv", sep='\t')
mlc_dat_list = []
for vcf in vcf_files:
    sample_name = vcf.split('/')[4]
    MLC_scores = calc_mlc(vcf,MLC_table)
    MSS = sum(MLC_scores)
    mean_MLC = sum(MLC_scores)/len(MLC_scores)
    mlc_dat_list.append([sample_name,MSS,mean_MLC])
mlc_dat = pd.DataFrame(mlc_dat_list, columns = ['index','MSS','Mean_MLC'])
mlc_dat = mlc_dat.set_index('index')

# Tumour Mutational Burden
tmb_files = glob.glob("project/results/somatic/*TMB30*/*/*.tmb.metrics.csv")
tmb_dat_list = []
for file in tmb_files:
    sample_name = file.split('/')[4]
    dat = pd.read_csv(file, names = ['1','2','index',sample_name], usecols = ['index',sample_name])
    dat = dat.set_index('index')
    tmb_dat_list.append(dat)
tmb_merged_dat = pd.concat(tmb_dat_list, axis=1).T

all_merged = ploidy_merged_dat.merge(copy_number_dat, left_index=True, right_index=True)
all_merged = all_merged.merge(prs_dat, left_index=True, right_index=True)
all_merged = all_merged.merge(mlc_dat, left_index=True, right_index=True)
all_merged = all_merged.merge(tmb_merged_dat, how = 'left', left_index=True, right_index=True)
all_merged.to_csv('all_sample_metrics.csv')

    
    
    