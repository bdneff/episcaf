#!/usr/bin/env python3
import pandas as pd
df = pd.read_csv('runs/run_test_rfd3_nompmn/04_filter/metrics_full.csv')
ok = df[df['status']=='ok']
print('Total ok rows:', len(ok))
print('With mean_pae:', ok['mean_pae'].notna().sum())
print()
passing = ok[
    (ok['mean_pae'] < 7) & (ok['overall_rmsd'] <= 3) &
    (ok['epitope_chunk_rmsd'] <= 2) & (ok['af3_n_clash_res'] == 1)
]
print('=== PASSING ALL 4 FILTERS:', len(passing), '===')
if len(passing) > 0:
    print(passing[['id','pred','overall_rmsd','epitope_chunk_rmsd','mean_pae','af3_n_clash_res']].sort_values('epitope_chunk_rmsd').to_string(index=False))
