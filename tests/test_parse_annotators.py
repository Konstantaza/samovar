"""Tests for the parse_annotators and annotators_wrapper modules."""

import os
import sys
import pytest
import pandas as pd
import sqlite3

from unittest.mock import MagicMock

# ==============================================================================
# ЗАГЛУШКА
# ==============================================================================
mock_ete3 = MagicMock()
mock_ncbi_instance = MagicMock()

# Создаем функции-обманки, которые понимают, кого у них спрашивают
def fake_get_lineage(taxid):
    if int(taxid) == 9606: 
        return [1, 9605, 9606] # Родословная человека (Homo sapiens)
    if int(taxid) == 511145: 
        return [1, 561, 562, 511145] # Родословная штамма E. coli
    raise ValueError("Taxid not found") # Ошибка для невалидных ID (как 1234567890)

def fake_get_rank(lineages):
    # Словарь рангов для наших фейковых родословных
    ranks = {1: 'no rank', 9605: 'genus', 9606: 'species', 561: 'genus', 562: 'species', 511145: 'strain'}
    return {t: ranks.get(t, 'no rank') for t in lineages}

# Привязываем наши умные функции к моку (используем side_effect вместо return_value)
mock_ncbi_instance.get_lineage.side_effect = fake_get_lineage
mock_ncbi_instance.get_rank.side_effect = fake_get_rank
mock_ete3.NCBITaxa.return_value = mock_ncbi_instance

# Подменяем системный модуль ete3 нашим фейком
sys.modules['ete3'] = mock_ete3

# Import our new OOP architecture
from samovar.annotators_wrapper import get_annotator_instance
from samovar.parse_annotators import (
    Annotation,
    RankAnnotation,
    ExpandAnnotation
)


@pytest.fixture
def test_data_dir():
    """Return path to test data directory."""
    return os.path.join(os.path.dirname(__file__), "data")


@pytest.fixture
def kaiju_file(test_data_dir):
    """Return path to Kaiju test file."""
    return os.path.join(test_data_dir, "kaiju.log")


@pytest.fixture
def kraken1_file(test_data_dir):
    """Return path to Kraken1 test file."""
    return os.path.join(test_data_dir, "kraken.log")


@pytest.fixture
def kraken2_file(test_data_dir):
    """Return path to Kraken2 test file."""
    return os.path.join(test_data_dir, "kraken2.log")


@pytest.fixture
def krakenu_file(test_data_dir):
    """Return path to Kraken Unique test file."""
    return os.path.join(test_data_dir, "krakenunique.log")


@pytest.fixture
def metaphlan_file(test_data_dir):
    """Return path to MetaPhlAn test file."""
    return os.path.join(test_data_dir, "metaphlan4.log")


# Testing individual parsers (OOP approach)
@pytest.mark.parametrize("tool_type, fixture_name", [
    ("kaiju", "kaiju_file"),
    ("kraken1", "kraken1_file"),
    ("kraken2", "kraken2_file"),
    ("krakenuniq", "krakenu_file"),
    ("metaphlan", "metaphlan_file")
])
def test_parse_standard_outputs(tool_type, fixture_name, request):
    """Test reading output files using the new OOP wrappers.
    
    This replaces 5 old separate functions because our new classes
    guarantee a standardized ["seq", "taxID"] output format.
    """
    # Dynamically get the file path from the fixture name
    file_path = request.getfixturevalue(fixture_name)
    
    # Instantiate the correct annotator via our factory
    annotator = get_annotator_instance(tool_type, run_config={}, config={})
    
    # Parse the output
    df = annotator.parse_output(file_path)
    
    # Assertions to ensure our OOP standardization works flawlessly
    assert isinstance(df, pd.DataFrame), f"{tool_type} failed to return a DataFrame"
    
    # Every tool must return exactly these two columns now
    assert list(df.columns) == ["seq", "taxID"], f"{tool_type} returned wrong columns"
    
    # Ensure sequence IDs are unique (no duplicate reads)
    assert len(df.seq.unique()) == len(df.seq), f"Duplicate sequences found in {tool_type}"
    
    # Ensure the dataframe is not empty
    assert len(df) > 0, f"Parsed DataFrame for {tool_type} is empty"


# Testing the main Annotation class
def test_read_annotation(test_data_dir):
    """Test reading multiple annotation files and combining them."""
    # Define the dictionary of files and their corresponding tool types
    file_path_type = {
        os.path.join(test_data_dir, "kaiju.log"): "kaiju",
        os.path.join(test_data_dir, "kraken.log"): "kraken1",
        os.path.join(test_data_dir, "kraken2.log"): "kraken2",
        os.path.join(test_data_dir, "krakenunique.log"): "krakenuniq"
    }
    
    # Initialize the main class. Under the hood, this will use our new OOP 
    # factory to parse each file and concatenate them into a single DataFrame.
    ann = Annotation(file_path_type, get_true_annotation=r".*")
    
    # Ensure export works without crashing
    out_csv = "tests_outs/annotation.csv"
    os.makedirs("tests_outs", exist_ok=True)
    ann.export(out_csv)

    # Assertions
    assert isinstance(ann.DataFrame, pd.DataFrame), "Annotation.DataFrame is not a pandas DataFrame"
    assert len(ann.DataFrame) > 0, "Combined DataFrame is empty"
    
    # Check if columns were properly renamed to avoid collisions (e.g., taxID_kraken2_2)
    assert any(col.startswith("taxID_") for col in ann.DataFrame.columns), "Columns were not properly renamed"


def test_annotation_class(test_data_dir):
    """Test Annotation class rank expansions and subsetting."""
    file_path_type = {
        os.path.join(test_data_dir, "kaiju.log"): "kaiju",
        os.path.join(test_data_dir, "kraken.log"): "kraken1"
    }
    ann = Annotation(file_path_type, get_true_annotation=r".*")
    
    assert isinstance(ann.DataFrame, pd.DataFrame)
    assert len(ann.true_annotation) > 0
    
    rank_ann = ann.rank_annotation("species")
    assert isinstance(rank_ann, RankAnnotation)
    
    # This also hits the real database
    expand_ann = ann.expand_annotation(["species", "genus"])
    assert isinstance(expand_ann, ExpandAnnotation)
    assert len(expand_ann.rank_annotation) == 2


def test_correct_level():
    """Test correct_level method of Annotation class with synthetic data."""
    # We create synthetic data to avoid relying on specific file parsing here.
    # We use known NCBI taxonomy IDs for testing.
    test_data = {
        'seq': ['seq1', 'seq2', 'seq3', 'seq4'],
        'taxID_kraken1_0': ['9606', '511145', '0', '1234567890'],  # Human, E. coli strain, unclassified, invalid
        'taxID_kaiju_1': ['9606', '0', '0', '1234567890']
    }
    # Important: Set 'seq' as index, just like our parser does
    df = pd.DataFrame(test_data).set_index('seq')
    
    # Create an empty Annotation object and inject our fake DataFrame
    ann = Annotation({}, get_true_annotation=None)
    ann.DataFrame = df
    
    # Test correcting to SPECIES level 
    ann.correct_level(level='species')
    
    assert ann.DataFrame['taxID_kraken1_0'].iloc[0] == '9606'  # Human (9606) is already species level
    assert ann.DataFrame['taxID_kraken1_0'].iloc[1] == '562'   # E. coli strain (511145) MUST roll up to species (562)
    assert ann.DataFrame['taxID_kraken1_0'].iloc[2] == '0'     # Unclassified stays 0
    assert ann.DataFrame['taxID_kraken1_0'].iloc[3] == '1234567890' # Invalid ID stays invalid
    
    # Test correcting to GENUS level
    ann.correct_level(level='genus')
    
    assert ann.DataFrame['taxID_kraken1_0'].iloc[0] != '9606'  # Homo (genus) is 9605, so it should change
    assert ann.DataFrame['taxID_kraken1_0'].iloc[1] != '511145' # Escherichia (genus) is 561, so it should change


# Auxiliary Classes and MetaPhlAn DB integration
def test_rank_annotation_class():
    """Test RankAnnotation class initialization."""
    # We initialize the class with a specific rank
    rank_ann = RankAnnotation("species")
    assert rank_ann.rank == "species", "RankAnnotation failed to set rank"
    # By default, its internal annotation should be a DataFrame (even if empty)
    assert hasattr(rank_ann, 'annotation'), "RankAnnotation missing 'annotation' attribute"


def test_expand_annotation_class():
    """Test ExpandAnnotation class initialization."""
    expand_ann = ExpandAnnotation()
    assert isinstance(expand_ann.rank_annotation, dict), "ExpandAnnotation.rank_annotation is not a dict"


def test_metaphlan_with_db(test_data_dir):
    """Test reading MetaPhlAn output with SQLite database mapping using OOP."""
    # Prepare a fake database directory
    db_dir = os.path.join(test_data_dir, "mock_metaphlan_db")
    os.makedirs(db_dir, exist_ok=True)
    db_file = os.path.join(db_dir, "mpa_v30_CHOCOPhlAn_201901_species_map.db")
    
    # Build a tiny SQLite database on the fly
    conn = sqlite3.connect(db_file)
    cursor = conn.cursor()
    cursor.execute("""
        CREATE TABLE mpa_species_map (
            ref_id TEXT,
            tax_id TEXT
        )
    """)
    # Insert some fake mappings (e.g., matching IDs that might be in metaphlan4.log)
    fake_mapping_data = [
        ("M367-c418", "511145"),  # E. coli strain
        ("M1206-c595", "9606"),   # Human
    ]
    cursor.executemany("INSERT INTO mpa_species_map VALUES (?, ?)", fake_mapping_data)
    conn.commit()
    conn.close()
    
    # Test parsing using our OOP architecture
    # We inject the mock DB directory into the run_config
    run_config = {"db_path": db_dir} 
    annotator = get_annotator_instance("metaphlan", run_config=run_config, config={})
    
    # Parse the actual test file
    df = annotator.parse_output(os.path.join(test_data_dir, "metaphlan4.log"))
    
    # Assertions
    assert isinstance(df, pd.DataFrame), "MetaPhlanAnnotator failed to return DataFrame"
    assert "taxID" in df.columns, "taxID column missing"
    
    # Cleanup the fake database to keep the system clean
    os.remove(db_file)
    os.rmdir(db_dir)