# Introdcution

## msgpypulse
* msgpypulse is a python module to analyse MS-GF+ results after PHRP & MASIC merger

* A custom build pipeline and can be run in a single command.

A comprehensive MS-GF+ downstream analysis pipeline for peptide identification and protein quantification.

## Overview

msgpypulse is a Python-based CLI tool that automates the complete MS-GF+ downstream analysis workflow. It processes peptide identification results from MS-GF+ and generates various analytical outputs including peptide cross-tabs, protein quantification matrices, and coverage heatmaps.

## Features

- **Complete Pipeline Automation**: End-to-end processing from SICs to final protein quantification
- **Flexible Configuration**: Customizable filtering thresholds and analysis parameters
- **Multiple Output Formats**: TSV tables, PDF heatmaps, and annotated cross-tabs
- **Robust Error Handling**: Comprehensive logging and error recovery
- **Modular Design**: Each pipeline step can be run independently

## Installation

### Prerequisites

```bash
python >= 3.8
pip
```

### Install from PyPI

```bash
pip install msgpypulse
```

### Install from Source

```bash
git clone https://github.com/thulasis/msgpypulse.git
cd msgpypulse
pip install -e .
```

## Quick Start

### Basic Usage

```bash
msgpypulse_cli.py -i /path/to/sic_input -o /path/to/final_output.tsv --database /path/to/database
```

### With Custom Parameters

```bash
msgpypulse_cli.py \
  -i /path/to/sic_input \
  -o /path/to/final_output.tsv \
  --database /path/to/database \
  --score_field MSMSScore \
  --threshold 10 \
  --num_pep 2 \
  --mode all_matches \
  --rollup sum \
  --cluster_rollup \
  --cluster_output /path/to/cluster_output.tsv
```

## Pipeline Steps

The pipeline executes the following steps in order:

1. **SICs Preparation**: Convert SIC files to standardized format
2. **FDR Estimation**: Estimate false discovery rates
3. **Shared Peptide Merging**: Integrate peptides from PHRPOut results
4. **True Peptide Filtering**: Filter based on configurable score thresholds
5. **Peptide Crosstab Generation**: Create sample-by-peptide matrices
6. **FASTA Clustering**: Cluster database sequences using CD-HIT
7. **Peptide Crosstab Annotation**: Add protein metadata (clusters, functions, genes)
8. **Peptide-Protein Map Building**: Generate coverage matrices for each sample
9. **Protein Rollup**: Aggregate peptide evidence to protein level
10. **Cluster Crosstab (Optional)**: Generate cluster-level summaries

## Configuration Options

| Argument | Description | Default | Required |
|----------|-------------|---------|----------|
| `-i`, `--input` | Input SIC directory | None | Yes |
| `-o`, `--output` | Final output file path | None | Yes |
| `--database` | Database directory path | None | Yes |
| `--score_field` | Score field for filtering | `MSMSScore` | No |
| `--threshold` | Threshold for filtering | `10` | No |
| `--score_field2` | Optional second score field | None | No |
| `--threshold2` | Optional second threshold | None | No |
| `--num_pep` | Minimum peptides per protein | `2` | No |
| `--mode` | Analysis mode | `all_matches` | No |
| `--rollup` | Rollup method | `sum` | No |
| `--outlier_alpha` | Grubbs alpha for outlier detection | None | No |
| `--coverage_tsv` | Path to grouped_coverage.tsv | None | No |
| `--cluster_rollup` | Run cluster crosstab step | False | No |
| `--cluster_output` | Output file for cluster crosstab | None | Required if `--cluster_rollup` |

## Usage Examples

### Example 1: Basic Analysis

```bash
msgpypulse_cli.py \
  -i /path/to/SICs \
  -o ./final_protein_matrix.tsv \
  --database /path/to/database.fasta
```

### Example 2: Advanced Analysis with Clustering

```bash
msgpypulse_cli.py \
  -i ./data/sics \
  -o ./results/final_proteins.tsv \
  --database ./data/database.fasta \
  --score_field MSMSScore \
  --threshold 20 \
  --score_field2 EValue \
  --threshold2 0.01 \
  --num_pep 3 \
  --mode requires_unique \
  --rollup weighted \
  --outlier_alpha 0.05 \
  --cluster_rollup \
  --cluster_output ./results/cluster_proteins.tsv
```

### Example 3: Running Individual Steps

The pipeline can be customized by modifying the `MSGFPDOWNPipeline` class methods. For example, to run only the FDR estimation:

```python
from msgpypulse_cli import MSGFPDOWNPipeline

pipeline = MSGFPDOWNPipeline()
pipeline.fdr_estimator(input_dir="./data/SICs", output_dir="./results/fdr_estd")
```


## Dependencies

- **Python**: >= 3.8
- **Core Libraries**:
  - `pandas` - Data manipulation
  - `matplotlib` - Plotting
  - `seaborn` - Statistical visualization
  - `msgpypulse` - Core analysis functions
- **Optional**:
  - `cd-hit` - Sequence clustering (installed via system package manager)

## Development

### Contributing

1. Fork the repository
2. Create a feature branch
3. Make your changes
4. Test your changes
5. Submit a pull request

## Troubleshooting

### Common Issues

**Issue**: "ModuleNotFoundError: No module named 'msgpypulse'"  
**Solution**: Ensure msgpypulse is installed: `pip install msgpypulse`

**Issue**: "FileNotFoundError: [Errno 2] No such file or directory"  
**Solution**: Check that input paths exist and are accessible

**Issue**: "ImportError: cannot import name 'prepare_SICs'"  
**Solution**: Update msgpypulse to the latest version

**Issue**: Slow performance with large datasets  
**Solution**: 
- Use `--mode unique_only` for faster processing
- Ensure sufficient system memory
- Consider using `--rollup sum` instead of weighted methods

### Getting Help

- **GitHub Issues**: [Report bugs or request features](https://github.com/thulasis/msgpypulse/issues)
- **Discussions**: [Community discussions](https://github.com/thulasis/msgpypulse/discussions)
- **Email**: tulasi_relangi@baylor.edu

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.

## Citation

If you use msgpypulse in your research, please cite:

```
@misc{thulasis2025msgpypulse,
  title={msgpypulse: A comprehensive MS-GF+ downstream analysis pipeline},
  authors={Tulasi Rao Relangi, PhD},{Harrison Hall, Graduate Student}
  group=Wright Lab, Baylor University
  year={2025},
  url={https://github.com/thulasis/msgpypulse}
}
```

## Acknowledgements

- **MS-GF+**: For the initial peptide identification
- **PHRP**: For shared peptide identification
- **CD-HIT**: For sequence clustering
