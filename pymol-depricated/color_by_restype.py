import colorsys,sys
from pymol import cmd

aa_1_3 = {
  'A': 'ALA', 
  'C': 'CYS', 
  'D': 'ASP', 
  'E': 'GLU', 
  'F': 'PHE', 
  'G': 'GLY', 
  'H': 'HIS', 
  'I': 'ILE', 
  'K': 'LYS', 
  'L': 'LEU', 
  'M': 'MET', 
  'N': 'ASN', 
  'P': 'PRO', 
  'Q': 'GLN', 
  'R': 'ARG', 
  'S': 'SER', 
  'T': 'THR', 
  'V': 'VAL', 
  'W': 'TRP', 
  'Y': 'TYR', 
}

aa_3_1 = {
  'ALA' : 'A', 
  'CYS' : 'C', 
  'ASP' : 'D', 
  'GLU' : 'E', 
  'PHE' : 'F', 
  'GLY' : 'G', 
  'HIS' : 'H', 
  'ILE' : 'I', 
  'LYS' : 'K', 
  'LEU' : 'L', 
  'MET' : 'M', 
  'ASN' : 'N', 
  'PRO' : 'P', 
  'GLN' : 'Q', 
  'ARG' : 'R', 
  'SER' : 'S', 
  'THR' : 'T', 
  'VAL' : 'V', 
  'TRP' : 'W', 
  'TYR' : 'Y', 
}

aa_types = {
  'A': 'hydrophobic',
  'C': 'cysteine',
  'D': 'negative',
  'E': 'negative',
  'F': 'aromatic',
  'G': 'hydrophobic',
  'H': 'polar',
  'I': 'hydrophobic',
  'K': 'positive',
  'L': 'hydrophobic',
  'M': 'hydrophobic',
  'N': 'polar',
  'P': 'proline',
  'Q': 'polar',
  'R': 'positive',
  'S': 'polar',
  'T': 'polar',
  'V': 'hydrophobic',
  'W': 'aromatic',
  'Y': 'aromatic',
}


def color_by_restype(selection="all",
        hydrophobic='white',
        aromatic='magenta',
        polar='cyan',
        positive='blue',
        negative='red',
        cysteine='yellow',
        proline='green',
        ):

  """
  usage: color_by_restype <selection>, <optional overrides of default colors>

  e.g. color_by_restype protein and chain A, hydrophobic=wheat

  Residue groups:
    hydrophobic: AGILMV
    aromatic: FWY
    polar: HNQST
    positive: KR
    negative: DE
    cysteine: C
    proline: P
  """
  colors = {
    'hydrophobic': hydrophobic,
    'aromatic': aromatic,
    'polar': polar,
    'positive': positive,
    'negative': negative,
    'cysteine': cysteine,
    'proline': proline,
  }

  for aa in aa_types:
    sel = selection + " and r. %s" % aa_1_3[aa]
#    print sel,"-->", colors[aa_types[aa]]
    cmd.color(colors[aa_types[aa]],sel)

cmd.extend("color_by_restype",color_by_restype)
