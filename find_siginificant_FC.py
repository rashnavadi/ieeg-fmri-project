import os
import pandas as pd
from collections import defaultdict

# Path to the subject-level analysis folder
root_dir = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/Tara/all_electrodes_as_seeds_staticFC_results'

# Column name and significance threshold
pval_col = 'p-value (FC)'
alpha = 0.05

# Output directory and text summary file
output_dir = '/Volumes/Rashnavadi/Documents/Data_Analysis/2023/analyses/ICE/Tara/stats_main_channels'
os.makedirs(output_dir, exist_ok=True)
text_summary_path = os.path.join(output_dir, 'stat_results_for_all_EEG_electrodes_as_seeds.txt')

# Containers
all_significant_results = []
total_rows_checked = 0
total_significant = 0
subject_stats = defaultdict(lambda: {'total': 0, 'significant': 0, 'details': set()})
log_lines = []  # Collect messages to write to .txt file

# Walk through all Excel files
for dirpath, _, filenames in os.walk(root_dir):
    for file in filenames:
        if file.endswith('.xlsx') and '_FC_analysis_results' in file and not file.startswith('~$') and not file.startswith('._'):
            file_path = os.path.join(dirpath, file)
            try:
                df = pd.read_excel(file_path, engine='openpyxl')
                if pval_col in df.columns:
                    num_rows = len(df)
                    sig_df = df[df[pval_col] < alpha].copy()
                    num_sig = len(sig_df)

                    total_rows_checked += num_rows
                    total_significant += num_sig

                    subj = df['Subject'].iloc[0] if 'Subject' in df.columns else file.split('_')[0]
                    subject_stats[subj]['total'] += num_rows
                    subject_stats[subj]['significant'] += num_sig

                    details = []
                    if {'Run', 'IED Type', 'Seed'}.issubset(df.columns):
                        unique_info = df[['Run', 'IED Type', 'Seed']].drop_duplicates().astype(str)
                        for _, row in unique_info.iterrows():
                            detail_str = f"Run: {row['Run']} | IED Type: {row['IED Type']} | Seed: {row['Seed']}"
                            subject_stats[subj]['details'].add(detail_str)
                            details.append(detail_str)

                    # Log and print for this subject
                    log_lines.append(f"\n📄 {subj}:")
                    print(f"\n📄 {subj}:")
                    for detail in sorted(details):
                        log_lines.append(f"  {detail}")
                        print(f"  {detail}")
                    summary = f"  ➤ {num_rows} total rows, {num_sig} significant (p < {alpha})"
                    log_lines.append(summary)
                    print(summary)

                    if not sig_df.empty:
                        sig_df['SourceFile'] = file_path
                        all_significant_results.append(sig_df)
                else:
                    msg = f"⚠️ Column '{pval_col}' not found in {file}"
                    print(msg)
                    log_lines.append(msg)
            except Exception as e:
                msg = f"❌ Could not process {file_path}: {e}"
                print(msg)
                log_lines.append(msg)

# Save all significant results
if all_significant_results:
    result_df = pd.concat(all_significant_results, ignore_index=True)
    output_file = os.path.join(output_dir, 'significant_FC_results_all_subjects.xlsx')
    result_df.to_excel(output_file, index=False)
    msg = f"\n✅ Saved significant results to: {output_file}"
    print(msg)
    log_lines.append(msg)
else:
    msg = "⚠️ No significant p-values found in any file."
    print(msg)
    log_lines.append(msg)

# Final summary
summary1 = f"\n🔍 Processed {total_rows_checked} total rows across all Excel files."
summary2 = f"✅ Found {total_significant} rows with significant p-value (< {alpha})."
print(summary1)
print(summary2)
log_lines.append(summary1)
log_lines.append(summary2)

# Ranked per-subject summary
log_lines.append("\n📊 Ranked Subjects by Significant FC Findings:\n")
print("\n📊 Ranked Subjects by Significant FC Findings:\n")
sorted_subjects = sorted(subject_stats.items(), key=lambda x: x[1]['significant'], reverse=True)
for subj, stats in sorted_subjects:
    line = f"⭐ {subj}: {stats['total']} total studies, {stats['significant']} significant FC values (p < {alpha})"
    print(line)
    log_lines.append(line)

# Write all to summary file
with open(text_summary_path, 'w') as f:
    f.write('\n'.join(log_lines))
