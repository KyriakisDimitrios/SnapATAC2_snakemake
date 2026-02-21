import sys
import os
import logging
import time
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
    input_csvs = sys.argv[2:]  # Captures all 95 file paths seamlessly
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
    # Bind all rows together
    master_df = pd.concat(df_list)

    # Sort alphabetically by index
    master_df.sort_index(inplace=True)

    # Force all numeric columns to 3 decimal places
    master_df = master_df.round(3)

    # Save the master table
    master_df.to_csv(output_master_csv)
    logging.info(f"Successfully saved master QC stats: {output_master_csv}")

    # --- Outlier Detection (IQR Method) ---
    logging.info("Scanning for statistical outliers across the cohort...")


    def get_iqr_bounds(series):
        q1 = series.quantile(0.25)
        q3 = series.quantile(0.75)
        iqr = q3 - q1
        return (q1 - 1.5 * iqr), (q3 + 1.5 * iqr)


    # Calculate bounds for all three critical metrics
    tsse_lower, tsse_upper = get_iqr_bounds(master_df['median_tsse'])
    frag_lower, frag_upper = get_iqr_bounds(master_df['median_nfragments'])
    ncells_lower, ncells_upper = get_iqr_bounds(master_df['ncells'])

    # Create boolean mask for mathematical outliers
    outlier_mask = (
            (master_df['median_tsse'] < tsse_lower) |
            (master_df['median_tsse'] > tsse_upper) |
            (master_df['median_nfragments'] < frag_lower) |
            (master_df['median_nfragments'] > frag_upper) |
            (master_df['ncells'] < ncells_lower) |
            (master_df['ncells'] > ncells_upper)
    )

    outliers_df = master_df[outlier_mask].copy()

    if not outliers_df.empty:
        # --- BIOLOGICAL RESCUE ---
        # Keep samples if TSSe is highly enriched (> 8.0) AND cell count isn't doublet soup (< 60k)
        rescue_mask = (outliers_df['median_tsse'] >= 8.0) & (outliers_df['ncells'] < 60000)

        rescued_count = rescue_mask.sum()
        if rescued_count > 0:
            logging.info(
                f"Biologically rescued {rescued_count} sample(s) that flagged statistical limits but show high-quality signal.")

        # Drop the rescued samples from the outlier dataframe
        outliers_df = outliers_df[~rescue_mask]

    # Re-check if empty after rescue
    if not outliers_df.empty:
        # Tag exactly why the sample was flagged
        def get_reason(row):
            reasons = []
            if row['median_tsse'] < tsse_lower: reasons.append("Low TSSe")
            if row['median_tsse'] > tsse_upper: reasons.append("High TSSe")
            if row['median_nfragments'] < frag_lower: reasons.append("Low Frags")
            if row['median_nfragments'] > frag_upper: reasons.append("High Frags")
            if row['ncells'] < ncells_lower: reasons.append("Low Cell Count")
            if row['ncells'] > ncells_upper: reasons.append("High Cell Count")
            return " & ".join(reasons)


        outliers_df['outlier_reason'] = outliers_df.apply(get_reason, axis=1)

        # Determine output path dynamically
        out_dir = os.path.dirname(output_master_csv)
        outlier_csv_path = os.path.join(out_dir, "project_possible_outliers.csv")

        outliers_df.to_csv(outlier_csv_path)
        logging.info(f"Flagged {len(outliers_df)} fatal outliers. Saved to {outlier_csv_path}")
    else:
        logging.info("No fatal outliers detected in the cohort after biological rescue.")

else:
    logging.error("No valid CSV files were parsed.")
    sys.exit(1)

duration = time.time() - start_time
logging.info(f"Completed: Compile Master QC in {duration:.2f} seconds")