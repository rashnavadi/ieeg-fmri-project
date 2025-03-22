import os
import pandas as pd

# Path where all the subject folders are located
# root_dir = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/Tara/main_IED_electrodes_as_seeds_staticFC_results'
root_dir = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/Tara/all_electrodes_as_seeds_staticFC_results'

# Column name and significance threshold
pval_col = 'p-value (FC)'
alpha = 0.05

# Containers
all_significant_results = []
total_rows_checked = 0
total_significant = 0

# Walk through all files
for dirpath, _, filenames in os.walk(root_dir):
    for file in filenames:
        if file.endswith('.xlsx') and '_FC_analysis_results' in file and not file.startswith('~$') and not file.startswith('._'):
            file_path = os.path.join(dirpath, file)
            try:
                df = pd.read_excel(file_path, engine='openpyxl')
                if pval_col in df.columns:
                    total_rows_checked += len(df)
                    sig_df = df[df[pval_col] < alpha].copy()
                    total_significant += len(sig_df)
                    if not sig_df.empty:
                        sig_df['SourceFile'] = file_path
                        all_significant_results.append(sig_df)
                else:
                    print(f"⚠️ Column '{pval_col}' not found in {file}")
            except Exception as e:
                print(f"❌ Could not process {file_path}: {e}")

# Combine and save
if all_significant_results:
    result_df = pd.concat(all_significant_results, ignore_index=True)
    output_file = 'significant_FC_results_all_subjects.xlsx'
    result_df.to_excel(output_file, index=False)
    print(f"\n✅ Saved significant results to: {output_file}")
else:
    print("⚠️ No significant p-values found in any file.")

# Final summary
print(f"\n🔍 Processed {total_rows_checked} total rows across all Excel files.")
print(f"✅ Found {total_significant} rows with significant p-value (< {alpha}).")
