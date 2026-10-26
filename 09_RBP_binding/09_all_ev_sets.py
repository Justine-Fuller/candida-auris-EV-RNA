# Rerun on ALL RNAs found in EVs: detected = reads in >=2 of 3 EV replicates in any strain x treatment group
# comparison = RNAs never detected in EVs but detected in cells (>=2 of 3 Cell sRNA replicates in some group)
import gzip,csv,collections
U='/mnt/user-data/uploads/CandidaEV/v2_downstream_2026-10-12/step6F_whole_RNA_tag/run/'; A='./rbp/all/'
ev=set(); cell=set(); cls={}
for g in ['B8441_untreated','B8441_caspofungin','MMC1_untreated','MMC1_caspofungin']:
    rd=csv.DictReader(gzip.open(U+f'F3_gene_tests/{g}.FX.gene_tests.tsv.gz','rt'),delimiter='\t')
    evc=[c for c in rd.fieldnames if '_EV_' in c and c.split('_')[-1].isdigit()]; cc=[c for c in rd.fieldnames if '_sRNA_' in c]
    for r in rd:
        if sum(float(r[c])>0 for c in evc)>=2: ev.add(r['model_id'])
        if sum(float(r[c])>0 for c in cc)>=2: cell.add(r['model_id'])
seqs={r['model_id']:r for r in csv.DictReader(gzip.open(U+'F1_sequences/gene_sequences.tsv.gz','rt'),delimiter='\t')}
cellonly=cell-ev
with open(A+'ev_all.fa','w') as f:
    for m in sorted(ev): f.write(f'>{m}\n{seqs[m]["seq_whole"]}\n')
with open(A+'cell_only.fa','w') as f:
    for m in sorted(cellonly): f.write(f'>{m}\n{seqs[m]["seq_whole"]}\n')
with open(A+'sets.tsv','w') as f:
    f.write('model_id\tset\trna_class\tlen_whole\n')
    for m in sorted(ev|cellonly): f.write(f"{m}\t{'EV_detected' if m in ev else 'cell_only'}\t{seqs[m]['rna_class']}\t{seqs[m]['len_whole']}\n")
print(len(ev),len(cellonly))
