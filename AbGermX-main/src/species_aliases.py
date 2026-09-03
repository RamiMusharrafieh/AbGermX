"""
Covers the species present in the IMGT/GENE-DB bulk reference file.
"""
import difflib

ALIASES = {
    'human': 'Homo sapiens',
    'mouse': 'Mus musculus',
    'house mouse': 'Mus musculus',
    'rat': 'Rattus norvegicus',
    'norway rat': 'Rattus norvegicus',
    'black rat': 'Rattus rattus',
    'rabbit': 'Oryctolagus cuniculus',
    'chicken': 'Gallus gallus',
    'cow': 'Bos taurus',
    'cattle': 'Bos taurus',
    'sheep': 'Ovis aries',
    'goat': 'Capra hircus',
    'pig': 'Sus scrofa',
    'swine': 'Sus scrofa',
    'horse': 'Equus caballus',
    'dog': 'Canis lupus familiaris',
    'cat': 'Felis catus',
    'ferret': 'Mustela putorius furo',
    'mink': 'Neogale vison',
    'alpaca': 'Vicugna pacos',
    'camel': 'Camelus dromedarius',
    'dromedary': 'Camelus dromedarius',
    'gorilla': 'Gorilla gorilla gorilla',
    'chimp': 'Pan troglodytes',
    'chimpanzee': 'Pan troglodytes',
    'orangutan': 'Pongo abelii',
    'sumatran orangutan': 'Pongo abelii',
    'bornean orangutan': 'Pongo pygmaeus',
    'lemur': 'Lemur catta',
    'ring-tailed lemur': 'Lemur catta',
    'naked mole rat': 'Heterocephalus glaber',
    'naked mole-rat': 'Heterocephalus glaber',
    'platypus': 'Ornithorhynchus anatinus',
    'zebrafish': 'Danio rerio',
    'salmon': 'Salmo salar',
    'atlantic salmon': 'Salmo salar',
    'trout': 'Oncorhynchus mykiss',
    'rainbow trout': 'Oncorhynchus mykiss',
    'grass carp': 'Ctenopharyngodon idella',
    'dolphin': 'Tursiops truncatus',
    'bottlenose dolphin': 'Tursiops truncatus',
    'hamster': 'Mesocricetus auratus',
    'golden hamster': 'Mesocricetus auratus',
    'baboon': 'Papio anubis anubis',
    'olive baboon': 'Papio anubis anubis',
    'sooty mangabey': 'Cercocebus atys',
    # macaques: several species share the common name "macaque"
    'macaque': 'Macaca mulatta',              # default: rhesus
    'rhesus macaque': 'Macaca mulatta',
    'rhesus': 'Macaca mulatta',
    'cynomolgus': 'Macaca fascicularis',
    'cynomolgus macaque': 'Macaca fascicularis',
    'crab-eating macaque': 'Macaca fascicularis',
    'pig-tailed macaque': 'Macaca nemestrina',
    'pigtail macaque': 'Macaca nemestrina',
}

def resolve_species(query, available_species):
    """
    Resolve a free-text species query to an exact IMGT species string
    present in `available_species` (an iterable of species names actually
    found in the loaded data).
    Resolution order: exact match -> alias table -> fuzzy match.
    Returns (resolved_name, match_type) or (None, 'unresolved').
    """
    q = query.strip()
    q_lower = q.lower()
    available = list(available_species)
    available_lower = {s.lower(): s for s in available}

    # 1. exact match against real data
    if q_lower in available_lower:
        return available_lower[q_lower], 'exact'

    # 2. alias table
    if q_lower in ALIASES:
        alias_target = ALIASES[q_lower]
        if alias_target in available:
            return alias_target, 'alias'

    # 3. match against both real species names and alias keys
    candidates = available + list(ALIASES.keys())
    close = difflib.get_close_matches(q_lower, [c.lower() for c in candidates], n=1, cutoff=0.6)
    if close:
        match = close[0]
        if match in available_lower:
            return available_lower[match], 'fuzzy'
        if match in ALIASES and ALIASES[match] in available:
            return ALIASES[match], 'fuzzy-alias'

    return None, 'unresolved'
