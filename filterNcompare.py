import glob, re
import pandas as pd

#Filter SNVs
files = glob.glob("germline_variants/snv/*.csv")
for file in files:
    sample_name = re.sub('_germline.*','',file.split('/')[-1])
    dat = pd.read_csv(file)
    rule1 = (dat.acmg_score > 0)
    rule2 = (dat.SpliceAI_acceptorGainScore > 0.5) & (dat.gnomad_allAf < 0.01)
    rule3 = (dat.SpliceAI_acceptorLossScore > 0.5) & (dat.gnomad_allAf < 0.01)
    rule4 = (dat.SpliceAI_donorGainScore > 0.5) & (dat.gnomad_allAf < 0.01)
    rule5 = (dat.SpliceAI_donorLossScore > 0.5) & (dat.gnomad_allAf < 0.01)
    rule6 = (dat['primateAI-3D_score'] > 0.5) & (dat.gnomad_allAf < 0.01)
    rule7 = (dat['alphamissense_score'] > 0.5) & (dat.gnomad_allAf < 0.01)
    
    filt_dat = dat[rule1|rule2|rule3|rule4|rule5|rule6|rule7]
    filt_dat = filt_dat.sort_values(by=['acmg_score'], ascending=False)
    filt_dat.to_csv(f'germline_variants/snv_filtered/{sample_name}_filtered.csv', index=False)

#Filter CNV
files = glob.glob("germline_variants/cnv/*.csv") 
for file in files:
    sample_name = re.sub('_germline.*','',file.split('/')[-1])
    dat = pd.read_csv(file)
    if(len(dat)==0):
        continue
    rule1 = (dat['clingen'] == 'pathogenic')
    rule2 = (dat['clinvar'] == 'pathogenic')
    
    filt_dat = dat[rule1|rule2]
    filt_dat = filt_dat.sort_values(by=['vid'], ascending=False)
    filt_dat.to_csv(f'germline_variants/cnv_filtered/{sample_name}_filtered.csv', index=False)

#Filter SV
files = glob.glob("germline_variants/sv/*.csv")
for file in files:
    sample_name = re.sub('_germline.*','',file.split('/')[-1])
    dat = pd.read_csv(file)
    if(len(dat)==0):
        continue
    rule1 = (dat['clingen'] == 'pathogenic')
    rule2 = (dat['clinvar'] == 'pathogenic')
    filt_dat = dat[rule1|rule2]
    filt_dat = filt_dat.sort_values(by=['vid'], ascending=False)
    filt_dat.to_csv(f'germline_variants/sv_filtered/{sample_name}_filtered.csv', index=False)


# Get differences betewen parent and progenitor lines
sample_rel = pd.read_csv("sample_relationship.csv")
for i, row in sample_rel.iterrows():
    parent_dat = pd.read_csv(f'germline_variants/snv_filtered/{row["Parent"]}_filtered.csv')
    progenitor_dat = pd.read_csv(f'germline_variants/snv_filtered/{row["Proband"]}_filtered.csv')
    only_parent = parent_dat[~parent_dat.vid.isin(progenitor_dat.vid)]
    only_parent.insert(0, row["Proband"], '(-)')
    only_progenitor = progenitor_dat[~progenitor_dat.vid.isin(parent_dat.vid)]
    only_progenitor.insert(0, row["Proband"], '(+)')
    merge = pd.concat([only_parent,only_progenitor])
    merge.to_csv(f'germline_variants/snv_unique/{row["Proband"]}_unique.csv', index=False)
    