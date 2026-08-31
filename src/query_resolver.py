"""
Resolves a comparison request (species pair + chain + gene segment type) into
concrete parameters the comparison engine can use: two IMGT species names,
an IG locus (IGH/IGK/IGL), and a region label to filter sequences on.

Supports both structured input and a loose freeform string like:
    "macaque vs mouse, heavy chain, V genes"
"""
import re
from species_aliases import resolve_species

CHAIN_MAP = {
    'heavy': 'IGH', 'igh': 'IGH', 'h': 'IGH',
    'kappa': 'IGK', 'igk': 'IGK', 'k': 'IGK',
    'lambda': 'IGL', 'igl': 'IGL', 'l': 'IGL',
    'light': 'LIGHT',   # special: means "run both IGK and IGL"
}

GENE_TYPE_MAP = {
    'v': {'V-REGION'}, 'variable': {'V-REGION'}, 'v gene': {'V-REGION'}, 'v genes': {'V-REGION'},
    'd': {'D-REGION'}, 'diversity': {'D-REGION'}, 'd gene': {'D-REGION'}, 'd genes': {'D-REGION'},
    'j': {'J-REGION'}, 'joining': {'J-REGION'}, 'j gene': {'J-REGION'}, 'j genes': {'J-REGION'},
    'c': {'C-REGION', 'CH1'}, 'constant': {'C-REGION', 'CH1'},
}

GENE_TYPE_CANONICAL_LABEL = {
    'v': 'V', 'variable': 'V', 'v gene': 'V', 'v genes': 'V',
    'd': 'D', 'diversity': 'D', 'd gene': 'D', 'd genes': 'D',
    'j': 'J', 'joining': 'J', 'j gene': 'J', 'j genes': 'J',
    'c': 'C', 'constant': 'C',
}

class QueryError(ValueError):
    pass

def resolve_chain(chain_str):
    key = chain_str.strip().lower()
    if key not in CHAIN_MAP:
        raise QueryError(
            f"Could not interpret chain '{chain_str}'. "
            f"Use one of: heavy, kappa, lambda, light."
        )
    return CHAIN_MAP[key]

def resolve_gene_type(gene_str):
    key = gene_str.strip().lower()
    if key not in GENE_TYPE_MAP:
        raise QueryError(
            f"Could not interpret gene type '{gene_str}'. "
            f"Use one of: V, D, J, C."
        )
    return GENE_TYPE_MAP[key]

def resolve_query(species_a_raw, species_b_raw, chain_raw, gene_type_raw, available_species):
    species_a, match_a = resolve_species(species_a_raw, available_species)
    species_b, match_b = resolve_species(species_b_raw, available_species)
    if species_a is None:
        raise QueryError(f"Could not resolve species '{species_a_raw}' to any species in the loaded data.")
    if species_b is None:
        raise QueryError(f"Could not resolve species '{species_b_raw}' to any species in the loaded data.")

    loci = resolve_chain(chain_raw)
    loci = ['IGK', 'IGL'] if loci == 'LIGHT' else [loci]

    region_set = resolve_gene_type(gene_type_raw)

    return dict(
        species_a=species_a, species_a_match=match_a,
        species_b=species_b, species_b_match=match_b,
        loci=loci, region_set=region_set,
        gene_type_label=GENE_TYPE_CANONICAL_LABEL[gene_type_raw.strip().lower()],
    )

_VS_SPLIT = re.compile(r'\s+vs\.?\s+|\s+versus\s+', re.IGNORECASE)

def parse_freeform_query(text):
    """
    Loosely parse a string like: "macaque vs mouse, heavy chain, V genes"
    Returns (species_a_raw, species_b_raw, chain_raw, gene_type_raw).
    Any part that can't be found defaults to: chain='heavy', gene_type='V'.
    """
    parts = [p.strip() for p in text.split(',')]
    if not parts or _VS_SPLIT.search(parts[0]) is None:
        raise QueryError(
            "Freeform query must start with '<species A> vs <species B>', "
            "e.g. \"macaque vs mouse, heavy chain, V genes\""
        )
    sp = _VS_SPLIT.split(parts[0], maxsplit=1)
    species_a_raw, species_b_raw = sp[0].strip(), sp[1].strip()

    chain_raw = 'heavy'
    gene_type_raw = 'V'
    rest = ' '.join(parts[1:]).lower()
    for key in CHAIN_MAP:
        if re.search(rf'\b{re.escape(key)}\b', rest):
            chain_raw = key
            break
    for key in ['v genes', 'd genes', 'j genes', 'v gene', 'd gene', 'j gene',
                'variable', 'diversity', 'joining', 'constant', ' v ', ' d ', ' j ', ' c ']:
        if key.strip() and re.search(rf'\b{re.escape(key.strip())}\b', rest):
            gene_type_raw = key.strip()
            break

    return species_a_raw, species_b_raw, chain_raw, gene_type_raw
