import colorsys,sys
from pymol import cmd


cmd.set_color('lightblue', [0.6, 0.8, 1.0])
def color_by_attype(selection="all",
        hydrophobic='white',
        positive='blue',
        negative='red',
        polar_n='lightblue',
        polar_o='salmon',
        sulfur='yellow'):

  """

  usage: color_by_attype <selection>, <optional overrides of default colors>

  e.g. color_by_attype protein and chain A, hydrophobic=white

  Atom types:            Default colours:
    hydrophobic            white
    polar_n                lightblue
    polar_o                salmon
    positive               blue
    negative               red
    sulfur                 yellow
  """
# used a list to maintain order!
  at_list = ['positive','negative','polar_o','polar_n','sulfur']

  at_types = {
    'positive': '(r. arg+lys and e. n)',
    'negative': '(n. oxt+ot1+ot2 or (r. asp+glu and e. o))',
    'polar_o': '(n. o or (r. ser+thr+asn+gln+tyr+tip+sol+wat+hoh and e. o))',
    'polar_n': '(n. n or (r. asn+gln+his+trp and e. n))',
    'sulfur': '(e. s)',
  }

  colors = {
    'hydrophobic': hydrophobic,
    'positive': positive,
    'negative': negative,
    'polar_n': polar_n,
    'polar_o': polar_o,
    'sulfur': sulfur,
  }

  cmd.color(hydrophobic,selection)
  for at in at_list:
    sel = selection + ' and ' + at_types[at]
#    print sel,"-->", colors[at]
    cmd.color(colors[at],sel)

cmd.extend("color_by_attype",color_by_attype)
