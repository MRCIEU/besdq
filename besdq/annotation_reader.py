"""Annotation TSV + EBI YAML metadata reader for GWAS-SSF import."""

import datetime
import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

try:
    import yaml
    _YAML_AVAILABLE = True
except ImportError:
    _YAML_AVAILABLE = False


@dataclass
class TraitConfig:
    file_path: str
    trait_id: str
    trait_name: str
    trait_chr: Optional[str]
    trait_bp: Optional[int]
    sample_size: int
    trait_var: float
    gene: Optional[str]
    context: Optional[str]
    study_metadata: dict = field(default_factory=dict)


_YAML_METADATA_KEYS = [
    'gwas_catalog_api',
    'date_metadata_last_modified',
    'genome_assembly',
    'genotyping_technology',
    'imputation_panel',
    'analysis_software',
    'adjusted_covariates',
    'samples',
]


def _to_json_safe(obj):
    """Recursively convert non-JSON-serializable objects to strings."""
    if isinstance(obj, (datetime.date, datetime.datetime)):
        return obj.isoformat()
    if isinstance(obj, dict):
        return {k: _to_json_safe(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [_to_json_safe(i) for i in obj]
    return obj


def _load_yaml_meta(yaml_path: Path) -> dict:
    if not yaml_path.exists():
        return {}
    if not _YAML_AVAILABLE:
        return {}
    with open(yaml_path, 'r') as fh:
        data = yaml.safe_load(fh) or {}
    meta = {}
    for key in _YAML_METADATA_KEYS:
        if key in data:
            meta[key] = _to_json_safe(data[key])
    return meta


def _yaml_sample_size(yaml_path: Path) -> Optional[int]:
    if not yaml_path.exists():
        return None
    if not _YAML_AVAILABLE:
        raise ImportError(
            "PyYAML is required to read sample_size from YAML metadata files. "
            "Install it with: pip install PyYAML"
        )
    with open(yaml_path, 'r') as fh:
        data = yaml.safe_load(fh) or {}
    samples = data.get('samples', [])
    if samples and isinstance(samples, list) and 'sample_size' in samples[0]:
        try:
            return int(samples[0]['sample_size'])
        except (ValueError, TypeError):
            return None
    return None


def read_trait_annotation(tsv_path: str) -> list[TraitConfig]:
    """Parse a trait annotation TSV and return a list of TraitConfig objects.

    Required columns: file_path, trait_id, trait_name
    Optional columns: trait_chr, trait_bp, sample_size, trait_var, gene, context
    """
    path = Path(tsv_path)
    if not path.exists():
        raise FileNotFoundError(f"Annotation file not found: {tsv_path}")

    with open(path, 'r') as fh:
        header = fh.readline().rstrip('\n').split('\t')
        col_idx = {name.strip(): i for i, name in enumerate(header)}

        required = ('file_path', 'trait_id', 'trait_name')
        missing = [c for c in required if c not in col_idx]
        if missing:
            raise ValueError(f"Annotation TSV missing required columns: {missing}")

        traits = []
        for lineno, line in enumerate(fh, start=2):
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')

            def get(col: str) -> Optional[str]:
                i = col_idx.get(col)
                if i is None or i >= len(parts):
                    return None
                v = parts[i].strip()
                return v if v else None

            file_path = get('file_path')
            trait_id = get('trait_id')
            trait_name = get('trait_name')

            if not file_path or not trait_id or not trait_name:
                raise ValueError(f"Line {lineno}: required columns (file_path, trait_id, trait_name) must be non-empty")

            if not Path(file_path).exists():
                raise FileNotFoundError(f"Line {lineno}: file_path does not exist: {file_path}")

            # Optional positional columns
            trait_chr_raw = get('trait_chr')
            trait_bp_raw = get('trait_bp')

            if (trait_chr_raw is None) != (trait_bp_raw is None):
                raise ValueError(
                    f"Line {lineno}: trait_chr and trait_bp must both be present or both absent"
                )

            trait_chr = trait_chr_raw
            trait_bp = int(trait_bp_raw) if trait_bp_raw is not None else None

            # Sample size: TSV takes precedence over YAML
            yaml_path = Path(file_path + '-meta.yaml')
            sample_size_raw = get('sample_size')
            if sample_size_raw is not None:
                try:
                    sample_size = int(sample_size_raw)
                except ValueError:
                    raise ValueError(f"Line {lineno}: sample_size must be an integer, got '{sample_size_raw}'")
            else:
                sample_size = _yaml_sample_size(yaml_path)
                if sample_size is None:
                    raise ValueError(
                        f"Line {lineno}: sample_size not found in TSV or YAML for {file_path}"
                    )

            trait_var_raw = get('trait_var')
            trait_var = float(trait_var_raw) if trait_var_raw is not None else 1.0

            gene = get('gene')
            context = get('context')

            study_metadata = _load_yaml_meta(yaml_path)

            traits.append(TraitConfig(
                file_path=file_path,
                trait_id=trait_id,
                trait_name=trait_name,
                trait_chr=trait_chr,
                trait_bp=trait_bp,
                sample_size=sample_size,
                trait_var=trait_var,
                gene=gene,
                context=context,
                study_metadata=study_metadata,
            ))

    return traits
