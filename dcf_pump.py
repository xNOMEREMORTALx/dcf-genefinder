"""dcf_pump (L3449): P-type ion PUMPS sorted by the metal they move, label-free. These are the 'battery chargers':
ATP energy is stored as an ion gradient that channels later release (nerve impulse, signalling).
Author: Gary NOMEREMORTAL Mason (DCF). Not for medical use.
Rule: P-type ATPase phosphorylation signature D-K-T-G-T-[LIVM]-[TI] (the aspartate that takes the phosphate), then
  TM4 metal-binding loop [NY]-V-P-E-G-L      -> NaK  (Na/K pump, incl. the gastric H/K pump)
  TM4 metal-binding loop A-[IV]-P-E-G-L      -> Ca   (SERCA/PMCA/SPCA calcium pumps)
  (development note: the residue BEFORE P-E-G-L decides; fish Na/K paralogs atp1a1a.2/.3/.5 end P-E-G-L-P)
  C-P-C in TM6 and >= 2 copper-binding GM[TSHN]C..C repeats -> Cu (copper pumps ATP7A/B)
  otherwise                                    -> other (lipid flippases ATP8-11, P5 ATP13A)
Functions: pump_class(protein) -> class or None."""
import re
def pump_class(p):
    if not re.search(r'DKTGT[LIVM][TI]', p): return None
    if re.search(r'[NY]VPEGL', p): return 'NaK'
    if re.search(r'A[IV]PEGL', p): return 'Ca'
    if 'CPC' in p and len(re.findall(r'GM[TSHN]C..C', p)) >= 2: return 'Cu'
    return 'other'
