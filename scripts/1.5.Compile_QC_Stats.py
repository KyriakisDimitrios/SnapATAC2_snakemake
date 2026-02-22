import sys
import os
import logging
import time
import numpy as np
import pandas as pd

# --- Logger Configuration ---
logging.basicConfig(
    format='%(asctime)s | %(levelname)-7s | %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
    stream=sys.stderr
)

start_time = time.time()
logging.info("Started: Compile Master QC Stats & Outlier Detection")

# --- Argument Parsing ---
try:
    output_master_csv = sys.argv[1]
    input_csvs = sys.argv[2:]
except IndexError:
    logging.error("Not enough arguments provided.")
    sys.exit(1)

logging.info(f"Aggregating {len(input_csvs)} sample QC files...")

# --- Concatenation Logic ---
df_list = []
for file in input_csvs:
    try:
        df = pd.read_csv(file, index_col='sample')
        df_list.append(df)
    except Exception as e:
        logging.error(f"Failed to read {file}: {e}")

if df_list:
    master_df = pd.concat(df_list)
    master_df.sort_index(inplace=True)
    master_df = master_df.round(3)

    master_df.to_csv(output_master_csv)
    logging.info(f"Successfully saved master QC stats: {output_master_csv}")

    # --- Outlier Detection ---
    logging.info("Scanning for statistical outliers using zero-hardcode scaling...")


    def get_iqr_bounds(series):
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        return (q1 - 1.5 * iqr), (q3 + 1.5 * iqr)


    tsse_lower, tsse_upper = get_iqr_bounds(master_df['median_tsse'])
    frag_lower, frag_upper = get_iqr_bounds(master_df['median_nfragments'])
    ncells_lower, ncells_upper = get_iqr_bounds(master_df['ncells'])

    # --- Purely Automated Log10 Ceiling for Doublet Soup ---
    # Transforms to log space to calculate natural variance, then converts back
    log_ncells = np.log10(master_df['ncells'])
    q1_log = log_ncells.quantile(0.25)
    q3_log = log_ncells.quantile(0.75)
    iqr_log = q3_log - q1_log

    # 2.0 multiplier in log-space safely encompasses healthy biological spread
    upper_log = q3_log + (2.0 * iqr_log)
    ncells_extreme_upper = 10 ** upper_log

    logging.info(f"Calculated automated cell ceiling via Log10-IQR: {ncells_extreme_upper:.0f} cells")

    # Create boolean mask
    outlier_mask = (
            (master_df['median_tsse'] < tsse_lower) |
            (master_df['median_tsse'] > tsse_upper) |
            (master_df['median_nfragments'] < frag_lower) |
            (master_df['median_nfragments'] > frag_upper) |
            (master_df['ncells'] < ncells_lower) |
            (master_df['ncells'] > ncells_extreme_upper)
    )

    outliers_df = master_df[outlier_mask].copy()

    if not outliers_df.empty:
        # --- BIOLOGICAL RESCUE ---
        # Rescue if TSSe is great AND cell count is safely below the automated log-ceiling
        rescue_mask = (outliers_df['median_tsse'] >= 8.0) & (outliers_df['ncells'] <= ncells_extreme_upper)

        rescued_count = rescue_mask.sum()
        if rescued_count > 0:
            logging.info(
                f"Biologically rescued {rescued_count} sample(s) that flagged math but show high-quality signal.")

        outliers_df = outliers_df[~rescue_mask]

    if not outliers_df.empty:
        def get_reason(row):
            reasons = []
            if row['median_tsse'] < tsse_lower: reasons.append("Low TSSe")
            if row['median_tsse'] > tsse_upper: reasons.append("High TSSe")
            if row['median_nfragments'] < frag_lower: reasons.append("Low Frags")
            if row['median_nfragments'] > frag_upper: reasons.append("High Frags")
            if row['ncells'] < ncells_lower: reasons.append("Low Cell Count")
            if row['ncells'] > ncells_extreme_upper: reasons.append("Extreme Cell Count (Doublet Soup)")
            return " & ".join(reasons)


        outliers_df['outlier_reason'] = outliers_df.apply(get_reason, axis=1)

        out_dir = os.path.dirname(output_master_csv)
        outlier_csv_path = os.path.join(out_dir, "project_possible_outliers.csv")

        outliers_df.to_csv(outlier_csv_path)
        logging.info(f"Flagged {len(outliers_df)} fatal outliers. Saved to {outlier_csv_path}")
    else:
        logging.info("No fatal outliers detected after biological rescue.")

else:
    logging.error("No valid CSV files were parsed.")
    sys.exit(1)

duration = time.time() - start_time
logging.info(f"Completed: Compile Master QC in {duration:.2f} seconds")