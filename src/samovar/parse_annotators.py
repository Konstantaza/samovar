"""Module for parsing and processing taxonomic annotation outputs from various tools.

This module provides functionality to read and process taxonomic annotation outputs
from tools like Kraken, Kaiju, and MetaPhlAn. It standardizes the output format
to allow for easy comparison and analysis of results from different tools.
"""

import pandas as pd
import os
import re
from ete3 import NCBITaxa
from typing import Dict, List, Optional
import sqlite3
from samovar.annotators_wrapper import get_annotator_instance

# Initialize NCBI taxonomy database
ncbi = NCBITaxa()

class Annotation:
    """Class for handling taxonomic annotations from multiple tools.
    
    Provides functionality to combine, process, and analyze results from different 
    classification tools at various taxonomic ranks.
    """
    
    def __init__(self, file_path_type: Dict[str, str], get_true_annotation: Optional[str] = None):
        """Initialize the Annotation object and combine input files.
        
        Args:
            file_path_type: Dictionary mapping file paths to their tool types (e.g., {'path/to/file': 'kraken2'}).
            get_true_annotation: Optional regex pattern to extract ground truth taxIDs from sequence headers.
        """
        self.id = 0
        self.DataFrame = pd.DataFrame()
        
        # Iteratively process each annotation file
        for path, tool_type in file_path_type.items():
            # Instantiate the appropriate annotator wrapper
            annotator = get_annotator_instance(tool_type, run_config={}, config={})
            
            # Parse the raw output into a standard DataFrame [seq, taxID]
            df = annotator.parse_output(path)
            
            # Set read IDs as index and ensure taxIDs are strings for consistency
            df = df.set_index("seq").astype({"taxID": "string"})
            
            # Rename columns to include tool type and unique instance ID to avoid collisions
            df.columns = [f"{col}_{tool_type}_{self.id}" for col in df.columns]
            
            # Handle potential duplicate sequence IDs in the tool's output
            if df.index.duplicated().any():
                df = df[~df.index.duplicated(keep="first")]
            
            # Concatenate the processed data into the main DataFrame
            try:
                self.DataFrame = pd.concat([self.DataFrame, df], axis=1)
            except Exception as e:
                raise ValueError(f"Failed to concatenate data from {path}: {e}")
            
            # Increment the ID for the next tool in the dictionary
            self.id += 1

        # Extract ground truth annotations if a regex pattern is provided
        self.true_annotation = []
        if get_true_annotation is not None:
            self.true_annotation = []
            for i in self.DataFrame.index:
                match = re.search(get_true_annotation, str(i))
                if match:
                    self.true_annotation.append(match.group(0))
                else:
                    self.true_annotation.append("")
            print(f"Extracted {len(self.true_annotation)} ground truth annotations.")

        # Build unique lists of annotations and ranks for further analysis
        taxid_columns = self.tr()
        all_found_taxids = []
        for name, column in taxid_columns.items():
            all_found_taxids += list(set(column.dropna()))

        self.annotation_list = self.list2set(all_found_taxids)
        self.true_annotation_list = self.list2set(self.true_annotation)
        self.rank_list = self.list2set([*self.annotation_list, *self.true_annotation_list])

    def true_annotation_unique(self) -> set:
        """Return a set of unique ground truth annotations."""
        return set(self.true_annotation)

    def true_annotation_rename(self, change_dict: Dict[str, str]):
        """Rename ground truth labels using a mapping dictionary."""
        self.true_annotation = [change_dict.get(ta, "") for ta in self.true_annotation]

    def rank_annotation(self, rank: str = "species") -> "RankAnnotation":
        """Map all current annotations to a specific taxonomic rank (e.g., 'genus')."""
        rank_list = [self.rank(j, rank) for j in self.rank_list]
        rank_dict = dict(zip(self.rank_list, rank_list))
        return RankAnnotation(rank).make(self.full(), rank_dict)

    def expand_annotation(self, rank: List[str] = ["species"]) -> "ExpandAnnotation":
        """Generate annotations for multiple taxonomic ranks simultaneously."""
        full_rank_annotation = ExpandAnnotation()
        for r in rank:
            full_rank_annotation.add(self.rank_annotation(r))
        return full_rank_annotation

    def correct_annotations(self, rank: str = "species") -> pd.DataFrame:
        """Calculate counts of correct predictions at a specific rank."""
        return pd.DataFrame(self.rank_annotation(rank).correct_annotation().annotation.value_counts())

    def full(self) -> pd.DataFrame:
        """Return the complete DataFrame including the ground truth column."""
        tmp = self.tr()
        tmp["true"] = self.true_annotation
        return tmp

    def correct_level(self, level: str = "sp") -> None:
        """Adjust all taxIDs in the DataFrame to the specified taxonomic level."""
        taxid_cols = [col for col in self.DataFrame.columns if col.startswith("taxID")]
        unique_taxids = set()
        for col in taxid_cols:
            unique_taxids.update(self.DataFrame[col].dropna().unique())
        
        taxid_map = {}
        for taxid in unique_taxids:
            if taxid == "0" or pd.isna(taxid):
                taxid_map[taxid] = taxid
                continue
            try:
                lineage = ncbi.get_lineage(int(taxid))
                ranks = ncbi.get_rank(lineage)
                for tid, r in ranks.items():
                    if r == level:
                        taxid_map[taxid] = str(tid)
                        break
                if taxid not in taxid_map:
                    taxid_map[taxid] = taxid
            except (ValueError, KeyError):
                taxid_map[taxid] = taxid
        
        for col in taxid_cols:
            self.DataFrame[col] = self.DataFrame[col].map(taxid_map)

    def tr(self) -> pd.DataFrame:
        """Filter the DataFrame to include only taxID columns."""
        return self.DataFrame.copy().filter(regex="taxID.*")

    @staticmethod
    def list2set(a: List) -> List[str]:
        """Convert a list to a unique set of strings."""
        return list(set([str(i) for i in a]))

    @staticmethod
    def list2rank(a: List, at_rank: str) -> List:
        """Convert a list of taxIDs to their corresponding IDs at a target rank."""
        b = Annotation.list2set(a)
        return [Annotation.rank(i, at_rank) for i in b]

    @staticmethod
    def rank(j: str, i: str) -> Optional[str]:
        """Retrieve the taxID at a specific rank for a given input taxID."""
        if j == "0" or pd.isna(j):
            return "0"
        try:
            lineage = ncbi.get_lineage(int(j))
            ranks = ncbi.get_rank(lineage)
            for taxid, r in ranks.items():
                if r == i:
                    return str(taxid)
            return None
        except (ValueError, KeyError):
            return None

    def export(self, file: Optional[str] = None) -> pd.DataFrame:
        """Export the taxID columns and ground truth to a CSV file."""
        df_return = self.DataFrame.loc[:, [col for col in self.DataFrame if col.startswith("taxID")]]
        lencol = [col for col in self.DataFrame if col.startswith("len")]
        if len(lencol) > 0:
            df_return["length"] = self.DataFrame.loc[:, lencol[0]].to_list()
        df_return["true"] = self.true_annotation
        if file is not None:
            df_return.to_csv(file)
        return df_return


class RankAnnotation:
    """Class representing taxonomic annotations at a single fixed rank."""
    
    def __init__(self, rank: str):
        self.rank = rank
        self.annotation = pd.DataFrame()

    def add(self, name: str, annotation: List):
        self.annotation[name] = annotation

    def make(self, annotation: pd.DataFrame, rank_dict: Dict[str, str]) -> "RankAnnotation":
        for name, column in annotation.items():
            self.add(str(name), [rank_dict.get(str(ta)) for ta in column])
        self.reindex(annotation.index)
        return self

    def reindex(self, index):
        self.annotation.index = index

    def y(self) -> pd.Series:
        return self.annotation["true"]

    def x(self) -> pd.DataFrame:
        return self.annotation.copy().drop("true", axis=1)

    def correct_annotation(self) -> "RankAnnotation":
        tmp = RankAnnotation(self.rank)
        for name, column in self.x().items():
            tmp.add(str(name), pd.DataFrame(column == self.y()))
        tmp.reindex(self.annotation.index)
        return tmp


class ExpandAnnotation:
    """Class managing taxonomic annotations across multiple ranks."""
    
    def __init__(self):
        self.rank_annotation = {}

    def add(self, rank_annotation: RankAnnotation):
        self.rank_annotation[rank_annotation.rank] = rank_annotation

    def get(self, rank: str) -> RankAnnotation:
        return self.rank_annotation[rank]


def match_annotation(annotation_name: str) -> Optional[str]:
    """Determine the tool type from a filename (expects .tool_name.out format)."""
    basename = os.path.basename(annotation_name)
    if basename.endswith(".out"):
        parts = basename.split(".")
        if len(parts) < 2:
            return None
        
        tool_tag = parts[-2].lower()
        
        # Mapping filename tags to internal tool types
        mapping = {
            "kraken": "kraken", "kraken1": "kraken",
            "kraken2": "kraken2",
            "krakenuniq": "krakenuniq", "krakenu": "krakenuniq",
            "metaphlan": "metaphlan", "metaphlan4": "metaphlan", "mpa": "metaphlan", "mp4": "metaphlan",
            "kaiju": "kaiju",
            "custom": "custom"
        }
        
        if tool_tag in mapping:
            return mapping[tool_tag]
        else:
            raise ValueError(f"Unknown tool tag '{tool_tag}' in file {basename}.")
    return None