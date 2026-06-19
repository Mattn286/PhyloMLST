import pandas as pd
import networkx as nx
from itertools import combinations
from collections import defaultdict
from Bio import Phylo
import argparse

def main(args):
    assert args.locus_threshold in {1, 2, 3, 4, 5, 6}, \
        "Locus threshold must be an integer between 1 and 6."
    assert isinstance(args.not_ic_threshold, int) and args.not_ic_threshold >= 0, \
        "Not IC threshold must be a non-negative integer."
        
    df = pd.read_csv(args.mlst, sep="\t")
    founder_ST_to_IC_dict = get_founder_STs(args.founder)
    baseline_IC_to_ST_dict = get_ST_IC(args.st_to_ic)
    if args.tree is not None:
        tree = get_tree(args.tree)
        print(f"Parsed tree with {len(tree.get_terminals())} isolates")
    
    df, marker_columns, ST_column = prepare_input(df)
    print(f"Found {len(marker_columns)} marker genes: {', '.join(marker_columns)}")
    
    df = assign_CC(df, marker_columns, args.locus_threshold) 
    print("Assigned isolates to clonal complexes")
    
    df = assign_IC(df, founder_ST_to_IC_dict, baseline_IC_to_ST_dict, ST_column)
    print("Assigned isolates to international clones")
    
    if args.tree is not None:
        df = correct_IC_by_phylogeny_func(df, tree, founder_ST_to_IC_dict, ST_column, baseline_IC_to_ST_dict, IC_column="Assigned IC", not_IC_value="Not IC", not_ic_threshold=args.not_ic_threshold)
        print("Corrected international clone assignment")
        
    if args.mask_non_ST is True:
        df = mask_non_ST(df)
    
    print_output(df)
    df.to_csv(args.output, sep="\t", index=True, index_label='Genome')
    if args.report is not None:
        write_ic_counts_to_tsv(df, args.report)
    
def get_tree(path):
    tree = Phylo.read(path, "newick")
    return tree

def get_founder_STs(path):
    df = pd.read_csv(path, sep="\t", header=None)
    assert df.shape[1] >= 2, "TSV must have at least two columns (ST and IC)."
    
    st_to_ic = {str(row[0]): row[1] for _, row in df.iterrows()}
    return st_to_ic

def get_ST_IC(path):
    if path is None:
        return None

    df = pd.read_csv(path, sep="\t", header=None)

    ic_dict = {}
    for _, row in df.iterrows():
        key = row[0]
        value_str = str(row[1]) if pd.notna(row[1]) else ""
        values = [int(x) for x in value_str.split()] if value_str else []
        ic_dict[key] = values

    return ic_dict

def prepare_input(df):
    df.columns = df.columns.astype(str)
    
    assert "Genome" in df.columns and "ST" in df.columns, \
        "Input DataFrame must contain 'Genome' and 'ST' columns."
    
    df.set_index("Genome", inplace=True)
    ST_column = "ST"
    
    cols_to_remove = {"ST", "Scheme", "clonal_complex", "species"}
    marker_columns = [c for c in df.columns if c not in cols_to_remove]
    
    return df, marker_columns, ST_column
    
def assign_CC(df, marker_columns, locus_threshold=1):
    # Create a dict where the keys are tuples containing unique combinations of marker genes (i.e. a profile) and the values are lists of isolates with this combination 
    profile_groups = defaultdict(list)
    for idx, row in df.iterrows():
        profile_key = tuple(row[col] for col in marker_columns)
        profile_groups[profile_key].append(str(idx))
    
    # Create graph
    G = nx.Graph()
    
    # Add nodes to the network and connect nodes with identical profiles immediately
    for group in profile_groups.values():
        G.add_nodes_from(group)
        for i in range(len(group)):
            for j in range(i + 1, len(group)):
                G.add_edge(group[i], group[j])
    
    # Flatten profiles for SLV comparison
    unique_profiles = list(profile_groups.keys()) #Gives a list of tuples containing all profiles
    
    # Function to calculate locus differences between two profiles
    def profile_diff(p1, p2):
        return sum(a != b for a, b in zip(p1, p2))
    
    # Links profiles which differ by 1 marker gene
    edges_to_add = []
    for i, j in combinations(range(len(unique_profiles)), 2):
        if profile_diff(unique_profiles[i], unique_profiles[j]) <= locus_threshold:
            # Create a list of isolate pairs which differ by 1 marker gene
            edges_to_add.extend([(x, y) for x in profile_groups[unique_profiles[i]] for y in profile_groups[unique_profiles[j]]])
    
    G.add_edges_from(edges_to_add)
    
    # Find connected components (clonal complexes)
    clonal_complexes = list(nx.connected_components(G))
    
    # Sort CCs by size descending
    clonal_complexes.sort(key=lambda x: -len(x))
    
    # Map isolates to their CC ID, largest = 0, next = 1, etc.
    cc_map = {isolate: cc_id for cc_id, cc in enumerate(clonal_complexes) for isolate in cc}
    df["CC"] = df.index.astype(str).map(cc_map).astype(int)
    
    return df

def mask_non_ST(df):
    non_numeric_mask = pd.to_numeric(df["ST"], errors="coerce").isna()
    for col in ["CC", "Assigned IC", "Corrected IC"]:
        if col in df.columns:
            df[col] = df[col].astype(object)
            df.loc[non_numeric_mask, col] = "NA"
    return df

def assign_IC(df, founder_ST_to_IC_dict, baseline_IC_to_ST_dict, ST_column, IC_column="Assigned IC", not_IC_value="Not IC"):
    CC_to_IC_dict = {}
    for _, row in df.iterrows():
        ST = row[ST_column]
        if not (CC_to_IC_dict.get(row["CC"]) in founder_ST_to_IC_dict.values()): # If CC has not already been assigned an IC
            CC_to_IC_dict[row["CC"]] = founder_ST_to_IC_dict.get(ST, not_IC_value) # Assign CC to IC correponding to founder ST
    df[IC_column] = df["CC"].map(CC_to_IC_dict)
    
    for index, row in df.iterrows(): 
        try:
            ST = int(row[ST_column])
        except:
            continue
        for IC, ST_lst in baseline_IC_to_ST_dict.items():
            if ST in ST_lst:
                df.at[index, IC_column] = IC
                break
    
    return df

def correct_IC_by_phylogeny_func(df, tree, founder_ST_to_IC_dict, ST_column,
                                  baseline_IC_to_ST_dict, IC_column="Assigned IC",
                                  not_IC_value="Not IC", not_ic_threshold=0):

    # Build st_to_corrected_ic from baseline, inverting IC->ST_list to ST->IC
    st_to_corrected_ic = {}
    if baseline_IC_to_ST_dict is not None:
        for ic, st_list in baseline_IC_to_ST_dict.items():
            for st in st_list:
                st_to_corrected_ic[st] = ic

    # Apply baseline overrides to df before tree correction
    df["Corrected IC"] = df[ST_column].map(st_to_corrected_ic)
    df[IC_column] = df[IC_column].where(
        ~df[ST_column].isin(st_to_corrected_ic), df["Corrected IC"]
    )

    # Reorder df to match tree leaf order
    tree_leaf_names = [c.name for c in tree.get_terminals()]
    df = df.loc[df.index.intersection(tree_leaf_names)]

    def is_valid_st(st):
        return pd.notna(st) and str(st).isdigit()

    # Isolates explicitly assigned Not IC — these are what make a subclade unclean
    not_ic_isolates = set(df[df[IC_column] == not_IC_value].index.astype(str))

    for ic in df[IC_column].unique():
        if ic == not_IC_value or pd.isna(ic):
            continue

        ic_isolates = set(df[df[IC_column] == ic].index.astype(str))
        if not ic_isolates:
            continue

        # Find MRCA of all isolates assigned to this IC
        mrca = tree.common_ancestor(ic_isolates)
        mrca_leaves = {c.name for c in mrca.get_terminals()}

        # A clade is unclean if it contains more than the threshold Not IC isolates
        non_ic_in_mrca = mrca_leaves & not_ic_isolates
        if len(non_ic_in_mrca) <= not_ic_threshold:
            # Clean MRCA — confirm all valid STs in this IC
            for st in df.loc[df[IC_column] == ic, ST_column].unique():
                if is_valid_st(st):
                    st_to_corrected_ic[st] = ic
            continue

        # Contaminated MRCA — find maximal clean subclades containing a founder ST
        rescued_sts = set()

        def find_clean_subclades(clade):
            leaves = {c.name for c in clade.get_terminals()}
            clade_ic = leaves & ic_isolates
            clade_non_ic = leaves & not_ic_isolates

            if not clade_ic:
                return

            if len(clade_non_ic) <= not_ic_threshold:
                # Clean subclade — rescue valid STs if a founder ST is present
                sts_in_clade = {
                    st for st in df.loc[df.index.isin(clade_ic), ST_column].unique()
                    if is_valid_st(st)
                }
                if any(founder_ST_to_IC_dict.get(st) == ic for st in sts_in_clade):
                    rescued_sts.update(sts_in_clade)
                return

            for child in clade.clades:
                find_clean_subclades(child)

        find_clean_subclades(mrca)

        for st in rescued_sts:
            st_to_corrected_ic[st] = ic

    # Map final ST->IC assignments back to all isolates
    df["Corrected IC"] = df[ST_column].map(st_to_corrected_ic).fillna(not_IC_value)

    return df

def print_output(df):
    assigned_counts = df["Assigned IC"].value_counts(dropna=False).sort_index()
    has_corrected = "Corrected IC" in df.columns
    corrected_counts = df["Corrected IC"].value_counts(dropna=False).sort_index() if has_corrected else None
    
    all_keys = sorted(set(assigned_counts.index).union(corrected_counts.index if has_corrected else []), 
                      key=lambda x: (str(type(x)), x))
    
    header = f"{'IC':<20} {'Assigned':>10}" + (f" {'Corrected':>10}" if has_corrected else "")
    print("\n")
    print(header)
    print("-" * 32)
    
    for ic in all_keys:
        assigned = assigned_counts.get(ic, 0)
        label = str(ic) if pd.notna(ic) else "NaN"
        
        if has_corrected:
            corrected = corrected_counts.get(ic, 0)
            print(f"{label:<20} {assigned:>10} {corrected:>10}")
        else:
            print(f"{label:<20} {assigned:>10}")

def write_ic_counts_to_tsv(df, output_path):
    assigned = df["Assigned IC"].value_counts(dropna=False)
    has_corrected = "Corrected IC" in df.columns
    
    if has_corrected:
        corrected = df["Corrected IC"].value_counts(dropna=False)
        summary_df = pd.DataFrame({"Assigned": assigned, "Corrected": corrected})
    else:
        summary_df = pd.DataFrame({"Assigned": assigned})
    
    summary_df = summary_df.fillna(0).astype(int).sort_index()
    summary_df.index = summary_df.index.map(lambda x: str(x) if pd.notna(x) else "NaN")
    summary_df.index.name = "IC"
    summary_df.to_csv(output_path, sep="\t")

def parse_args():
    parser = argparse.ArgumentParser(
        description="Cluster MLST profiles and assign isolates to ICs using a phylogenetic tree"
    )

    parser.add_argument(
        "-m", "--mlst",
        required=True,
        help="TSV file containing MLST data for all isolates"
    )

    parser.add_argument(
        "-t", "--tree",
        default=None,
        help="Optional: Phylogenetic tree in Newick format (default: None)"
    )

    parser.add_argument(
        "-f", "--founder",
        required=True,
        help="Path to TSV file matching founder ST to IC"
    )

    parser.add_argument(
        "-s", "--st_to_ic",
        default=None,
        help="Optional: path to TSV file matching ICs to a set of STs"
    )
    
    parser.add_argument(
        "-o", "--output",
        default="out.tsv",
        help="Optional: output file name (default: 'out.tsv')"
    )
    
    parser.add_argument(
        "-x", "--mask_non_ST",
        action="store_true",
        help="Mask CC and IC assignments for isolates with non-numeric STs"
    )
        
    parser.add_argument(
        "-r", "--report",
        default=None,
        help="Optional: report file name (default: None)"
    )
    
    parser.add_argument(
        "-l", "--locus_threshold",
        type=int,
        default=1,
        help="Optional: number of locus differences allowed when grouping isolates into CCs (default: 1)"
    )

    parser.add_argument(
        "-n", "--not_ic_threshold",
        type=int,
        default=0,
        help="Optional: number of Not IC isolates permitted in a subtree before it is considered unclean (default: 0)"
    )

    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()
    main(args)
