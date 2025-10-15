# create_cellphonedb.py

import os
import sqlite3
from pathlib import Path
import pandas as pd

# --------------------------------------------------
# Configuration
# --------------------------------------------------

# Path to the CellPhoneDB v5.0.0 resource directory (adjust as needed)
CELLPHONEDB_DIR = Path(
    "/home/DX6/carrolllab/cchan4/LabNGS/resources/cellphonedb/v5.0.0"
)

# Path to output SQLite database file in user's home directory
LOCAL_DB_PATH = Path.home() / "databases/cellphonedb.sqlite"

# Mapping of expected table names to CSV filenames
table_files = {
    "multidata_table": "multidata_table.csv",
    "interaction_input": "interaction_input.csv",
    "interaction_table": "interaction_table.csv",
    "gene_input": "gene_input.csv",
    "gene_table": "gene_table.csv",
    "gene_orthologue": "gene_orthologue.csv",
    "gene_synonym_to_gene_name": "gene_synonym_to_gene_name.csv",
    "protein_input": "protein_input.csv",
    "protein_table": "protein_table.csv",
    "complex_input": "complex_input.csv",
    "complex_table": "complex_table.csv",
    "complex_composition_table": "complex_composition_table.csv",
    "receptor_to_transcription_factor": "receptor_to_transcription_factor.csv",
}

# --------------------------------------------------
# Load and sanitize CSVs into DataFrames
# --------------------------------------------------

dataframes = {}

for name, filename in table_files.items():
    try:
        # Load CSV into DataFrame
        df = pd.read_csv(CELLPHONEDB_DIR / filename)

        # Clean column names: strip whitespace, lowercase, replace spaces with underscores
        df.columns = [
            c.strip().replace(" ", "_").lower() for c in df.columns
        ]

        dataframes[name] = df
    except Exception as e:
        print(f"??  Failed to load {filename}: {e}")

# --------------------------------------------------
# Write all tables to SQLite database
# --------------------------------------------------

with sqlite3.connect(str(LOCAL_DB_PATH)) as conn:
    for name, df in dataframes.items():
        try:
            df.to_sql(name, conn, if_exists="replace", index=False)
            print(f"? Loaded: {name}")
        except Exception as e:
            print(f"? Failed to load {name}: {e}")

    # Show all tables in the database
    cursor = conn.cursor()
    cursor.execute("SELECT name FROM sqlite_master WHERE type='table';")
    print("? Tables in DB:", [t[0] for t in cursor.fetchall()])

print(f"SQLite database written to: {LOCAL_DB_PATH}")
