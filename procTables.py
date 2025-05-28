import pandas as pd


# Get original codon table then merge with new tables
codon_tb = pd.read_csv("/Users/la.uqmzardb/projects/misc/hgunter_codon_table/mRNArchitect/data/ref_codon_usage.txt.bak",
                       sep=" ")
# DEBUG
print(codon_tb.head())

xls = pd.read_excel("/Users/la.uqmzardb/projects/misc/hgunter_codon_table/new_codon_models_v3.xlsx", sheet_name=None)
print(xls)
for name, table in xls.items():
    print(table.head())
    print(name)
    _df = table[["codon", "frequency"]].rename(columns={"frequency": name})
    codon_tb = codon_tb.merge(_df, "left", "codon")
    # DEBUG
    print(codon_tb.head())

codon_tb.to_csv("mRNArchitect/data/ref_codon_usage.txt", sep=" ", index=False)