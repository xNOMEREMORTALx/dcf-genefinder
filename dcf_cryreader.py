"""dcf_cryreader (L3801; vectorised traceback L3813; family rule L3829): read a protein's cryptochrome/photolyase family and its four-tryptophan electron chain (TrpA-D).
Gary NOMEREMORTAL Mason, Digital Cosmos Framework. Not for medical use.
Method (frozen in L3777/L3781, blind-tested in L3785 and L3797):
  - the protein (400-900 aa) is scored by local alignment (dcf_localalign, BLOSUM62, gap 8) against six references;
  - best score >= 350 = family member; family = best reference;
  - the tetrad is read by traceback alignment to robin CRY4 at TrpA 395, TrpB 372, TrpC 318, TrpD 369;
  - marker: a confident (score >= 900) cryptochrome/6-4 protein carries W, Y or F at TrpD (blind: 31 of 31 confident hits,
    L3785 + L3797); weak hits (350-899) broke this 3 times in 17 (L3785); no DASH/photolyase has been aromatic at TrpD."""
import os, runpy
import numpy as np
_L = runpy.run_path(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'dcf_localalign.py'))
score, M, enc = _L['score'], _L['M'], _L['enc']
GAP = 8
POS = {'A': 395, 'B': 372, 'C': 318, 'D': 369}
CRY = ('CRY4(robin)', 'phr6-4(fly)')
REFS = {
    'CRY4(robin)': 'MLHRTIHLFRKELRLHDNPVLLAALQSSEALYPVYILDRAFLTSSMHIGALRWHFLLQSLEDLHKNLCQLGSCLLVIQGEYETVLRDHIQKWSITQVTLDAEMEPFYKEMEANIQCLGAELGFEVLSLGSHSLYDTQRILDINGGSPPLTYKRFLHILSLLGDPEVPVRNLTAEDFQRCSAPDPDLAECYRVPLPVDLKISPENLSPWRGGETEGLQRLEQHLTDQGWVASFTKPRTIPNSLLPSTTGLSPYFSMGCLSVRTFFYRLSNIYAQAKHHSLPPVSLQGQLLWREFFYTVASATPNFTQMAGNPICLQICWYKDAERLHKWKMAQTGFPWIDAIMTQLRQEGWIHHLARHAVACFLTRGDLWISWEEGMKVFEELLLDADYSINAGNWMWLSASAFFHQYTRIFCPVRFGKRTDPQGNYIRKYLPILKNFPSKYIYEPWTASEEEQKQAGCIIGRDYPFPMVNHKEASDHNLQLMRQVREEQHRTAQLTRDDADDPMEIKVKRDHTEENISKGKVARTTE',
    'phr6-4(fly)': 'MDSQRSTLVHWFRKGLRLHDNPALSHIFTAANAAPGKYFVRPIFILDPGILDWMQVGANRWRFLQQTLEDLDNQLRKLNSRLFVVRGKPAEVFPRIFKSWRVEMLTFETDIEPYSVTRDAAVQKLAKAEGVRVETHCSHTIYNPELVIAKNLGKAPITYQKFLGIVEQLKVPKVLGVPEKLKNMPTPPKDEVEQKDSAAYDCPTMKQLVKRPEELGPNKFPGGETEALRRMEESLKDEIWVARFEKPNTAPNSLEPSTTVLSPYLKFGCLSARLFNQKLKEIIKRQPKHSQPPVSLIGQLMWREFYYTVAAAEPNFDRMLGNVYCMQIPWQEHPDHLEAWTHGRTGYPFIDAIMRQLRQEGWIHHLARHAVACFLTRGDLWISWEEGQRVFEQLLLDQDWALNAGNWMWLSASAFFHQYFRVYSPVAFGKKTDPQGHYIRKYVPELSKYPAGCIYEPWKASLVDQRAYGCVLGTDYPHRIVKHEVVHKENIKRMGAAYKVNREVRTGKEEESSFEEKSETSTSGKRKVRRATGSAPKRKR',
    'phr-classII(fly)': 'MFTLASYWRESFKIVLPLQAMKRTKAQKAGPSKKAAKNEKASSEPKSDQESSDEEASTSKALLVSKPDYQNFEQFLTHLEHQRVCTAANIQEFSFRKKRVRVLSKTEDVKESSLGGVVYWMSRDGRVQDNWALLFAQRLALKLELPLTVVFCLVPKFLNATIRHYKFMMGGLQEVEQQCRALDIPFHLLMGSAVEKLPQFVKSKDIGAVVCDFAPLRLPRQWVEDVGKALPKSVPLVQVDAHNVVPLWVASDKQEYAARTIRNKINSKLGEYLSVFPPVVRHPHGTGCKNVNTVDWSAAYASLQCDMEVDEVQWAKPGYKAACQQLYEFCSRRLRHFNDKRNDPTADALSGLSPWLHFGHISAQRCALEVQRFRGQHKASADAFCEEAIVRRELADNFCFYNEHYDSLKGLSSWAYQTLDAHRKDKRDPCYSLEELEKSLTYDDLWNSAQLQLVREGKMHGFLRMYWAKKILEWTATPEHALEYAILLNDKYSLDGRDPNGYVGCMWSIGGVHDMGWKERAIFGKVRYMNYQGCRRKFDVNAFVMRYGGKVHKKK',
    'cry-dash(zf)': 'MMSASRTVICLLRNDLRLHDNEVFHWAQRNAEHIIPLYCFDPRHYQGTYHYNFPKTGPFRLRFLLDSVKDLRALLKKHGSTLLVRQGKPEDVVCELIKQLGSVSTVAFHEEVASEEKSVEEKLKEICCQNKVRVQTFWGSTLYHRDDLPFSHIGGLPDVYTQFRKAVEAQGRVRPVLSTPEQVKSPPSGLEEGPIPTFDSLGQTEPLDDCRSAFPCRGGETEALARLKHYFWDTNAVATYKETRNGMIGVDFSTKFSPWLALGCISPRYIYEQIKKYEVERTANQSTYWVIFELLWRDYFKFVALKYGNRIFYMNGLQDKHVPWNTDMKMFDAWKEGRTGVPFVDANMRELALTGFMSNRGRQNVASFLTKDLGLDWRLGAEWFEYLLVDHDVCSNYGNWLYSAGIGNDPRENRKFNMIKQGLDYDNNGDYVRQWVPELRGIKGGDVHTPWTLSNSALSHAQVSLNQTYPCPIITAPEWSRHVNNKSSGSSSSKGRKGSSYTARQHKDRGIDFYFSKNKHF',
    'CPD-photolyase(Neurospora)': 'MAPSKRKASAPPQTSHVNGNPSADKKRKTTTDAPPTNPNTSSDPLRAPHPFYKDSETHGIVLRKFYPHEMSNARAQAYNDNELPRPIETLSAALAETAALRKSLPVRQAVVHWFKMDLRLHDNRSLWLASQKAKEAGVPLICLYVLSPEDLEAHLRAPIRVDFMLRTLEVLKTDLEDLGIPLWVETVEKRKEVPTKIKELMKSWGASHLFCAMEYEVDELRREAKLVKLLAEGEKGEKMAADVVHDTCVVMPGALQSGSGGQYAVYSPWFRAWIKHIEENPECLEIYEKPGPNPPGTKEKHENLFACSIPEAPEGKRLRDDEKARYHSLWPAGEHEALKRLEKFCDEAIGKYAERRNIPAMQGTSNLSVHFASGTLSARTAIRTARDRNNTKKLNGGNEGIQRWISEVAWRDFYKHVLVHWPYVCMNKPFKPTYSNIEWSYNVDHFHAWTQGRTGFPIIDAAMRQVLSTGYMHNRLRMIVASFLAKDLLVDWRMGERYFMEHLIDGDFASNNGGWGFAASVGVDPQPYFRVFNPLLQSEKFDPDGDYIRKWVEELRDLPELKGGKGGEIHDPYGRGSEKVKKKLEEKGYPRPIVEHSGARDRALDAYKRGLARDL',
    'DASH(diatom)': 'MSSSRSKQSFWFLIFVHTVLILTFAAGAATAMRGASVSTVTKATSRASSQVVLHWFRHGDLRLLDNPALIHSSKTAESCVPVFCFDDSVYGNDNRTPDTRAPHSNDRGQLKCGPRRAQFVLDSVQDLRRSLQSRGSALYVAHGKPAQVFQRLVDAWPAVPAADTAAPNGSLLTIVCQREVVREENDAVRAVQSVLRRRFPQAKVQQIWGSTMYELDDLPFATDLANMPDTFTPFRNKVEKNCQIGTPLPVPKQLSLPENFPSALKQGLEYLPTLKELGYTDAQIQQVETHDERGVLHFMVAKQPDSHDCLKDYFETRNGMLGPNYSTKFSPWLAHGNVSPRYVAAQCRKYEEERVENKSTYWVVFELLWRDFCKFFATKHGDAIFYPYGTTERTDRHKPWSTFGRNLQAWQEGRTGYPLVDANMRELVATGFMSNRGRQNVASFLAINLNHDWRCGGDFFESHLLDYDVYSNWVNWCAAAGMTGGRLNRFNISKQSKDYDQHGDYVRHWLPELAKVPNEFVHEPWKMTSFQQMEYECKLGVDYPNPIVPPSRPNPHTDRNNRGRGGHQSKGNNRHGPPDKSRKANASGSNRHQKYEMKSLQPGSFRVKES'}
def sw_map(ref, q):
    """Smith-Waterman with traceback (gap 8); returns {reference position (1-based): query residue}, best score.
    Row-vectorised (L3813): within a row, horizontal gaps are resolved exactly by a running maximum,
    H[j] = max_k<=j (V[k] - GAP*(j-k)) with V = max(0, diagonal, up). Same scores and tie rules as the L3753 loop version."""
    a = np.array(enc(ref)); b = np.array(enc(q)); n, m = len(a), len(b)
    H = np.zeros((n + 1, m + 1), np.int32); T = np.zeros((n + 1, m + 1), np.int8); idx = GAP * np.arange(m + 1, dtype=np.int64)
    for i in range(1, n + 1):
        d = np.empty(m + 1, np.int64); d[0] = -1 << 40; d[1:] = H[i - 1, :-1] + M[a[i - 1], b]
        u = H[i - 1].astype(np.int64) - GAP
        V = np.maximum(0, np.maximum(d, u)); V[0] = 0
        h = np.maximum.accumulate(V + idx) - idx; h[0] = 0; H[i] = h
        l = np.empty(m + 1, np.int64); l[0] = -1 << 40; l[1:] = h[:-1] - GAP
        T[i] = np.where(h == 0, 0, np.where(h == d, 1, np.where(h == u, 2, 3)))
    i, j = np.unravel_index(np.argmax(H), H.shape); mp = {}
    while i > 0 and j > 0 and T[i, j] != 0:
        t = T[i, j]
        if t == 1: mp[i] = q[j - 1]; i -= 1; j -= 1
        elif t == 2: i -= 1
        else: j -= 1
    return mp, int(H.max())
def read(seq):
    """Return dict: length, scores per reference, best, family, tetrad (4 letters, '-' = not aligned), call (plain words)."""
    seq = ''.join(c for c in seq.upper() if c.isalpha())
    r = {'length': len(seq)}
    if not 400 <= len(seq) <= 900:
        r['call'] = 'outside the 400-900 residue window the method was tested on (fragment or fusion?)'; return r
    sc = {k: score(seq, v) for k, v in REFS.items()}; v, best = max((s, k) for k, s in sc.items())  # same tie rule as the scan (L3777)
    r.update(scores=sc, best=best, score=v)
    if v < 350:
        r['call'] = 'not detected as a cryptochrome/photolyase (best score %d < 350)' % v; return r
    mp, _ = sw_map(REFS['CRY4(robin)'], seq); t = ''.join(mp.get(POS[k], '-') for k in 'ABCD'); r['tetrad'] = t
    fam = 'cryptochrome / 6-4 photolyase' if best in CRY else 'DASH / DNA photolyase'; r['family'] = fam
    conf = v >= 900; arom = t[3] in 'WYF'
    # L3829 reader update, family rule from L3809: TrpD aromatic -> cryptochrome-type chain; non-aromatic -> photolyase-type.
    # Record: 261 of 264 insect cryptochromes (PANTHER SF17) aromatic; 2,266 of 2,300 class I photolyases (SF18) non-aromatic;
    # 40 of 40 confident cryptochromes in 11 blind animals. Score only decides superfamily membership (>= 350).
    r['chain_type'] = 'cryptochrome-type' if arom else 'photolyase-type'
    rule = ('cryptochrome-type chain (aromatic TrpD; rule right for 98.6% of 2,564 PANTHER-classified proteins)' if arom
            else 'photolyase-type chain (non-aromatic TrpD, as in DNA-repair photolyases and CRY-DASH; rule right for 98.6% of 2,564 PANTHER-classified proteins)')
    agree = (best in CRY) == arom
    r['call'] = 'nearest %s (score %d, %s); TrpD = %s -> %s%s' % (best, v, 'strong' if conf else 'weak', t[3], rule,
                 '' if agree else '; NOTE: nearest reference and TrpD disagree - check by tree or PANTHER')
    return r
